// Hex9BinFilter.hpp — filters.hex9bin : aggregate a point cloud into Hex9 cells.
//
// Encodes each point (X=lon, Y=lat, WGS84) to its Hex9 address, runs the
// population-ceiling adaptive digest (hex9_adaptive) to a mixed-layer cell set
// (dense areas resolve fine, sparse coarse), and emits ONE point per cell at
// the cell centroid, carrying count/value/layer/bin plus mean/min/max/mode
// rollups of named input dimensions (colour, DEM/Z, intensity, class...).
//
// Attribute aggregation is a second pass over the per-point cell assignment
// (hex9_adaptive_assign), so the split policy (weight/ceiling) and the carried
// data are independent — count is just weight=1.
//
//   min_layer / max_layer : digest layer range (max_layer<0 => fixed=min_layer).
//   ceiling / floor       : population-ceiling params (weight units).
//   weight                : per-point weight dimension (default: count=1 each).
//   aggregates            : comma list "Dim:func[,...]", func = mean|min|max|mode.
//                           Output dim = Dim (mean/mode), DimMin, DimMax.
#pragma once

#include <pdal/Filter.hpp>
#include <string>
#include <vector>

namespace pdal {

class PDAL_EXPORT Hex9BinFilter : public Filter
{
public:
    Hex9BinFilter() = default;
    std::string getName() const override;
    enum class Agg { Mean, Min, Max, Mode };

private:
    struct Spec { std::string srcName, outName; Agg func;
                  Dimension::Id src{Dimension::Id::Unknown},
                                out{Dimension::Id::Unknown}; };

    int m_minLayer = 12;
    int m_maxLayer = -1;                 // <0 => fixed layer = m_minLayer
    double m_ceiling = 256.0;
    double m_floor = 1.0;
    std::string m_weightDim;
    std::string m_aggSpec;
    std::vector<Spec> m_specs;
    Dimension::Id m_weight{Dimension::Id::Unknown};
    Dimension::Id m_dLayer, m_dCount, m_dValue, m_dBinHi, m_dBinLo;

    void addArgs(ProgramArgs& args) override;
    void addDimensions(PointLayoutPtr layout) override;
    void ready(PointTableRef table) override;
    PointViewSet run(PointViewPtr view) override;

    Hex9BinFilter& operator=(const Hex9BinFilter&) = delete;
    Hex9BinFilter(const Hex9BinFilter&) = delete;
};

} // namespace pdal
