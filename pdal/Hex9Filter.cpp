// Hex9Filter.cpp — filters.hex9 : WGS84 lon/lat -> Hex9 w_oct (see header).
#include "Hex9Filter.hpp"

#include <pdal/PluginHelper.hpp>
#include <pdal/util/ProgramArgs.hpp>

#include "hex9_c.h"   // libhex9 C ABI (extern "C" guarded)

namespace pdal {

static PluginInfo const s_info {
    "filters.hex9",
    "Project WGS84 lon/lat (EPSG:4979/4326) to Hex9 warped octahedral "
    "cartesian (w_oct) — the seamless 3D storage CRS. Altitude (Z) is kept "
    "orthogonal (WGS84 ellipsoidal height).",
    "https://github.com/MrBenGriffin/hex9"
};

CREATE_SHARED_STAGE(Hex9Filter, s_info)

std::string Hex9Filter::getName() const { return s_info.name; }

void Hex9Filter::addArgs(ProgramArgs& args)
{
    args.add("replace", "Replace X/Y/Z with w_oct and move elevation to the "
             "Altitude dimension (default: add WoctX/WoctY/WoctZ, keep X/Y/Z).",
             m_replace, false);
}

void Hex9Filter::addDimensions(PointLayoutPtr layout)
{
    if (m_replace)
        m_alt = layout->registerOrAssignDim("Altitude", Dimension::Type::Double);
    else
    {
        m_wx = layout->registerOrAssignDim("WoctX", Dimension::Type::Double);
        m_wy = layout->registerOrAssignDim("WoctY", Dimension::Type::Double);
        m_wz = layout->registerOrAssignDim("WoctZ", Dimension::Type::Double);
    }
}

void Hex9Filter::ready(PointTableRef)
{
    char err[256] = {0};
    if (hex9_warp_init(err, sizeof err))
        log()->get(LogLevel::Warning)
            << "hex9_warp_init: " << err << " (falling back to identity warp)\n";

    const SpatialReference& srs = getSpatialReference();
    if (!srs.empty() && !srs.isGeographic())
        log()->get(LogLevel::Warning)
            << "filters.hex9 expects a geographic CRS (lon/lat); got a "
            << "projected one. Reproject to EPSG:4979 upstream "
            << "(filters.reprojection).\n";
}

void Hex9Filter::filter(PointView& view)
{
    const point_count_t n = view.size();
    for (PointId i = 0; i < n; ++i)
    {
        const double lon = view.getFieldAs<double>(Dimension::Id::X, i);
        const double lat = view.getFieldAs<double>(Dimension::Id::Y, i);
        double xyz[3];
        hex9_to_woct(lon, lat, xyz);

        if (m_replace)
        {
            const double h = view.getFieldAs<double>(Dimension::Id::Z, i);
            view.setField(m_alt, i, h);
            view.setField(Dimension::Id::X, i, xyz[0]);
            view.setField(Dimension::Id::Y, i, xyz[1]);
            view.setField(Dimension::Id::Z, i, xyz[2]);
        }
        else
        {
            view.setField(m_wx, i, xyz[0]);
            view.setField(m_wy, i, xyz[1]);
            view.setField(m_wz, i, xyz[2]);
        }
    }
}

} // namespace pdal
