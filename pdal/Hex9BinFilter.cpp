// Hex9BinFilter.cpp — filters.hex9bin (see header).
#include "Hex9BinFilter.hpp"

#include <pdal/PluginHelper.hpp>
#include <pdal/util/ProgramArgs.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <sstream>
#include <unordered_map>
#include <vector>

#include "hex9_c.h"

namespace pdal {

static PluginInfo const s_info {
    "filters.hex9bin",
    "Aggregate a point cloud into Hex9 cells (population-ceiling adaptive "
    "digest); one point per cell at its centroid with count/value/layer/bin "
    "and mean/min/max/mode rollups of named dimensions.",
    "https://github.com/MrBenGriffin/hex9"
};

CREATE_SHARED_STAGE(Hex9BinFilter, s_info)

std::string Hex9BinFilter::getName() const { return s_info.name; }

void Hex9BinFilter::addArgs(ProgramArgs& args)
{
    args.add("min_layer", "Finest digest layer (points binned here first).",
             m_minLayer, 12);
    args.add("max_layer", "Coarsest layer; <0 = fixed single layer = min_layer.",
             m_maxLayer, -1);
    args.add("ceiling", "Population ceiling per cell (weight units).",
             m_ceiling, 256.0);
    args.add("floor", "Minimum accumulated weight for a cell to emit.",
             m_floor, 1.0);
    args.add("weight", "Per-point weight dimension (default: count, 1 each).",
             m_weightDim, "");
    args.add("aggregates", "Comma list 'Dim:func' (func=mean|min|max|mode).",
             m_aggSpec, "");
}

static Hex9BinFilter::Agg parseFunc(const std::string& f)
{
    if (f == "mean") return Hex9BinFilter::Agg::Mean;
    if (f == "min")  return Hex9BinFilter::Agg::Min;
    if (f == "max")  return Hex9BinFilter::Agg::Max;
    if (f == "mode") return Hex9BinFilter::Agg::Mode;
    throw pdal_error("filters.hex9bin: unknown aggregate func '" + f + "'");
}

void Hex9BinFilter::addDimensions(PointLayoutPtr layout)
{
    m_dLayer = layout->registerOrAssignDim("H9Layer", Dimension::Type::Unsigned8);
    m_dCount = layout->registerOrAssignDim("H9Count", Dimension::Type::Unsigned64);
    m_dValue = layout->registerOrAssignDim("H9Value", Dimension::Type::Double);
    m_dBinHi = layout->registerOrAssignDim("H9BinHi", Dimension::Type::Unsigned64);
    m_dBinLo = layout->registerOrAssignDim("H9BinLo", Dimension::Type::Unsigned64);

    // parse "Dim:func,Dim:func,..."
    std::stringstream ss(m_aggSpec);
    std::string tok;
    while (std::getline(ss, tok, ','))
    {
        if (tok.empty()) continue;
        auto colon = tok.find(':');
        if (colon == std::string::npos)
            throw pdal_error("filters.hex9bin: bad aggregate '" + tok +
                             "' (want Dim:func)");
        Spec sp;
        sp.srcName = tok.substr(0, colon);
        sp.func = parseFunc(tok.substr(colon + 1));
        sp.src = layout->findDim(sp.srcName);
        if (sp.src == Dimension::Id::Unknown)
            throw pdal_error("filters.hex9bin: no such dimension '" +
                             sp.srcName + "'");
        switch (sp.func)
        {
            case Agg::Mean: case Agg::Mode: sp.outName = sp.srcName; break;
            case Agg::Min:  sp.outName = sp.srcName + "Min"; break;
            case Agg::Max:  sp.outName = sp.srcName + "Max"; break;
        }
        sp.out = layout->registerOrAssignDim(sp.outName, Dimension::Type::Double);
        m_specs.push_back(sp);
    }
}

void Hex9BinFilter::ready(PointTableRef table)
{
    char err[256] = {0};
    if (hex9_warp_init(err, sizeof err))
        log()->get(LogLevel::Warning) << "hex9_warp_init: " << err << "\n";
    if (!m_weightDim.empty())
        m_weight = table.layout()->findDim(m_weightDim);
}

PointViewSet Hex9BinFilter::run(PointViewPtr view)
{
    PointViewSet out;
    PointViewPtr ov = view->makeNew();
    out.insert(ov);
    const point_count_t n = view->size();
    if (n == 0) return out;

    // 1. lon/lat + weights -> full L29 uuids
    std::vector<double> lon(n), lat(n);
    for (PointId i = 0; i < n; ++i) {
        lon[i] = view->getFieldAs<double>(Dimension::Id::X, i);
        lat[i] = view->getFieldAs<double>(Dimension::Id::Y, i);
    }
    std::vector<uint8_t> uu(n * 16);
    hex9_encode_many(lon.data(), lat.data(), n, uu.data());

    std::vector<double> w;
    double* wp = nullptr;
    if (m_weight != Dimension::Id::Unknown) {
        w.resize(n);
        for (PointId i = 0; i < n; ++i) w[i] = view->getFieldAs<double>(m_weight, i);
        wp = w.data();
    }

    // 2. adaptive digest
    const int maxL = (m_maxLayer < 0) ? m_minLayer : m_maxLayer;
    char err[256] = {0};
    hex9_adaptive* a = hex9_adaptive_create(uu.data(), wp, n, m_minLayer, maxL,
                                            m_ceiling, m_floor, err, sizeof err);
    if (!a) throwError(err);
    const int nc = hex9_adaptive_count(a);
    std::vector<int64_t> assign(n);
    hex9_adaptive_assign(a, assign.data());

    // 3. attribute rollup (second pass over the assignment)
    struct Acc { double sum = 0, mn = std::numeric_limits<double>::infinity(),
                 mx = -std::numeric_limits<double>::infinity(); long cnt = 0;
                 std::unordered_map<long, long> hist; };
    std::vector<std::vector<Acc>> acc(m_specs.size(), std::vector<Acc>(nc));
    for (PointId p = 0; p < n; ++p) {
        const int64_t c = assign[p];
        if (c < 0 || c >= nc) continue;
        for (size_t s = 0; s < m_specs.size(); ++s) {
            const double v = view->getFieldAs<double>(m_specs[s].src, p);
            Acc& ac = acc[s][c];
            ac.sum += v; ac.cnt++;
            ac.mn = std::min(ac.mn, v); ac.mx = std::max(ac.mx, v);
            if (m_specs[s].func == Agg::Mode) ac.hist[std::lround(v)]++;
        }
    }

    // 4. emit one point per cell
    for (int c = 0; c < nc; ++c) {
        uint8_t binu[16], full[16];
        int layer; double value; int64_t npts;
        hex9_adaptive_cell(a, c, binu, &layer, &value, &npts);
        hex9_adaptive_cell_full(a, c, full);
        double clon, clat;
        hex9_decode(full, &clon, &clat);

        const PointId id = ov->size();
        ov->setField(Dimension::Id::X, id, clon);
        ov->setField(Dimension::Id::Y, id, clat);
        ov->setField(Dimension::Id::Z, id, 0.0);
        ov->setField(m_dLayer, id, (uint8_t)layer);
        ov->setField(m_dCount, id, (uint64_t)npts);
        ov->setField(m_dValue, id, value);

        uint64_t hi = 0, lo = 0;              // big-endian split of the bin uuid
        for (int k = 0; k < 8; ++k) hi = (hi << 8) | binu[k];
        for (int k = 8; k < 16; ++k) lo = (lo << 8) | binu[k];
        ov->setField(m_dBinHi, id, hi);
        ov->setField(m_dBinLo, id, lo);

        for (size_t s = 0; s < m_specs.size(); ++s) {
            const Acc& ac = acc[s][c];
            double r = 0.0;
            switch (m_specs[s].func) {
                case Agg::Mean: r = ac.cnt ? ac.sum / ac.cnt : 0.0; break;
                case Agg::Min:  r = ac.cnt ? ac.mn : 0.0; break;
                case Agg::Max:  r = ac.cnt ? ac.mx : 0.0; break;
                case Agg::Mode: {
                    long best = 0, bestN = -1;
                    for (auto& kv : ac.hist)
                        if (kv.second > bestN) { bestN = kv.second; best = kv.first; }
                    r = (double)best; break;
                }
            }
            ov->setField(m_specs[s].out, id, r);
        }
    }
    hex9_adaptive_destroy(a);
    return out;
}

} // namespace pdal
