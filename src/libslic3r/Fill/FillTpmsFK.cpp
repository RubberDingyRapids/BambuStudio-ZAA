#include "../ClipperUtils.hpp"
#include "../MarchingSquares.hpp"
#include "FillTpmsFK.hpp"
#include <cmath>
#include <algorithm>

namespace marchsq {
using namespace Slic3r;

using coordr_t = long; // length type for (r, c) raster coordinates.
using Pointf   = Vec2d; // (x, y) field point in coordf_t.

// Scalar field for the Fischer-Koch S triply-periodic minimal surface,
// sampled and contoured with marching squares (see MarchingSquares.hpp).
struct ScalarField
{
    static constexpr float gsizef = 0.40f;                        // grid cell size in mm (roughly line segment length).
    static constexpr float rsizef = 0.004f;                       // raster pixel size in mm (roughly point accuracy).
    const coord_t          rsize  = scaled(rsizef);              // raster pixel size in coord_t.
    const coordr_t         gsize  = std::round(gsizef / rsizef); // grid cell size in coordr_t.
    Point                  size;                                 // field size in coord_t.
    Point                  offs;                                 // field offset in coord_t.
    coordf_t               z;                                    // z offset as a float.
    float                  freq;                                 // field frequency in cycles per mm.
    float                  isoval = 0.0f;                        // iso value threshold to use.

    explicit ScalarField(const BoundingBox bb, const coordf_t z = 0.0, const float period = 10.0f)
        : size{bb.size()}, offs{bb.min}, z{z}, freq{float(2 * PI) / period}
    {}

    // Fischer - Koch S equation:
    // cos(2x)sin(y)cos(z) + cos(2y)sin(z)cos(x) + cos(2z)sin(x)cos(y) = 0
    float get_scalar(coordf_t x, coordf_t y, coordf_t z_arg) const
    {
        const float fx = freq * float(x);
        const float fy = freq * float(y);
        const float fz = freq * float(z_arg);
        return std::cos(2 * fx) * std::sin(fy) * std::cos(fz) + std::cos(2 * fy) * std::sin(fz) * std::cos(fx) +
               std::cos(2 * fz) * std::sin(fx) * std::cos(fy);
    }

    float get_scalar(Coord p) const
    {
        Pointf pf = to_Pointf(p);
        return get_scalar(pf.x(), pf.y(), z);
    }

    inline coord_t  to_coord (const coordr_t& x) const { return x * rsize; }
    inline coordr_t to_coordr(const coord_t& x)  const { return x / rsize; }
    inline Point  to_Point (const Coord& p) const { return Point(to_coord(p.c) + offs.x(), to_coord(p.r) + offs.y()); }
    inline Coord  to_Coord (const Point& p) const { return Coord(to_coordr(p.y() - offs.y()), to_coordr(p.x() - offs.x())); }
    inline Pointf to_Pointf(const Point& p) const { return Pointf(unscaled(p.x()), unscaled(p.y())); }
    inline Pointf to_Pointf(const Coord& p) const { return to_Pointf(to_Point(p)); }
};

// Register ScalarField as a raster type for MarchingSquares.
template<> struct _RasterTraits<ScalarField>
{
    using ValueType = float;
    static float  get (const ScalarField& sf, size_t row, size_t col) { return sf.get_scalar(Coord(row, col)); }
    static size_t rows(const ScalarField& sf) { return sf.to_coordr(sf.size.y()); }
    static size_t cols(const ScalarField& sf) { return sf.to_coordr(sf.size.x()); }
};

// Extract the polylines for the scalar field's iso-zero contour.
inline Polylines get_polylines(const ScalarField& sf, const double tolerance = SCALED_EPSILON)
{
    std::vector<Ring> rings = execute(sf, sf.isoval, {sf.gsize, sf.gsize});

    Polylines polys;
    polys.reserve(rings.size());
    for (const Ring& ring : rings) {
        Polyline poly;
        Points&  pts = poly.points;
        pts.reserve(ring.size() + 1);
        for (const Coord& crd : ring)
            pts.emplace_back(sf.to_Point(crd));
        // Marching squares rings are closed loops; duplicate the first point to
        // turn the ring into a closed Polyline.
        pts.push_back(pts.front());
        if (tolerance >= 0.0)
            poly.simplify(tolerance);
        polys.emplace_back(poly);
    }
    return polys;
}

} // namespace marchsq

namespace Slic3r {

void FillTpmsFK::_fill_surface_single(const FillParams&              params,
                                      unsigned int                   thickness_layers,
                                      const std::pair<float, Point>& direction,
                                      ExPolygon                      expolygon,
                                      Polylines&                     polylines_out)
{
    auto infill_angle = float(this->angle + (CorrectionAngle * 2 * M_PI) / 360.);
    if (std::abs(infill_angle) >= EPSILON)
        expolygon.rotate(-infill_angle);

    float density_factor = std::min(0.9f, params.density);
    // Density (field period) adjusted to have a good %of weight.
    const float vari_T = 4.18f * spacing * params.multiline / density_factor;

    BoundingBox bbox = expolygon.contour.bounding_box();
    // Enlarge the bounding box by the multi-line width to avoid artifacts at the edges.
    bbox.offset(scale_((params.multiline + 1) * spacing));
    marchsq::ScalarField sf = marchsq::ScalarField(bbox, this->z, vari_T);
    // Get simplified lines using a coarse tolerance (this is infill).
    Polylines polylines = marchsq::get_polylines(sf, SCALED_RESOLUTION);

    // Apply multiline offset if needed
    multiline_fill(polylines, params, spacing);

    // Prune the lines within the expolygon.
    polylines = intersection_pl(std::move(polylines), expolygon);

    if (!polylines.empty()) {
        // Remove very small bits, but be careful to not remove infill lines connecting thin walls!
        // The infill perimeter lines should be separated by around a single infill line width.
        const double minlength = scale_(0.8 * this->spacing);
        polylines.erase(std::remove_if(polylines.begin(), polylines.end(),
                                       [minlength](const Polyline& pl) { return pl.length() < minlength; }),
                        polylines.end());
    }

    if (!polylines.empty()) {
        // connect lines
        size_t polylines_out_first_idx = polylines_out.size();

        // chain_infill not suitable for this pattern due to internal "islands"; this also affects performance a lot.
        this->connect_infill(std::move(polylines), expolygon, polylines_out, this->spacing, params);

        // new paths must be rotated back
        if (std::abs(infill_angle) >= EPSILON) {
            for (auto it = polylines_out.begin() + polylines_out_first_idx; it != polylines_out.end(); ++it)
                it->rotate(infill_angle);
        }
    }
}

} // namespace Slic3r
