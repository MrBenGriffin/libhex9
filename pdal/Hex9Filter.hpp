// Hex9Filter.hpp — PDAL filter: WGS84 lon/lat -> Hex9 warped octahedral
// cartesian (w_oct), the seamless 3D storage CRS.
//
// Expects points already in a WGS84 geographic CRS (EPSG:4979 lon/lat/h, or
// 4326 lon/lat) — reproject upstream with filters.reprojection. X is longitude,
// Y is latitude (degrees); Z is the orthogonal WGS84 (ellipsoidal) altitude and
// is NEVER folded into the octahedral coordinate.
//
//   mode=augment (default): add WoctX/WoctY/WoctZ dimensions, keep X/Y/Z.
//   mode=replace          : set X/Y/Z = w_oct, move elevation to Altitude — the
//                           output cloud is then IN w_oct, so a COPC/EPT octree
//                           over it is a Hex9 spatial index (level 1 = octant).
#pragma once

#include <pdal/Filter.hpp>

namespace pdal {

class PDAL_EXPORT Hex9Filter : public Filter
{
public:
    Hex9Filter() = default;
    std::string getName() const override;

private:
    bool m_replace = false;                 // "replace" arg (else augment)
    Dimension::Id m_wx, m_wy, m_wz, m_alt;

    void addArgs(ProgramArgs& args) override;
    void addDimensions(PointLayoutPtr layout) override;
    void ready(PointTableRef table) override;
    void filter(PointView& view) override;

    Hex9Filter& operator=(const Hex9Filter&) = delete;
    Hex9Filter(const Hex9Filter&) = delete;
};

} // namespace pdal
