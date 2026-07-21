/* h9_table.cpp — set-returning surface: grid enumeration and the two exact
 * hierarchy enumerators (curve descendants / owned sub-zones), as DuckDB
 * table functions.
 *
 * h9_grid enumerates by bounding box (cells whose geographic CENTROID lies
 * within it — the core's semantic; PostGIS additionally filters to a real
 * bounds geometry in its glue, which a WHERE over the centroid replicates
 * here). Emits both the full identity uuid (h9_id, reversible, safe to
 * re-bin at any layer <= the grid layer) and the layer-scoped bin key
 * (h9_bin, for joining), plus WKB polygon/centroid — the F4 doctrine: render
 * from geom/centroid, use h9_id as key.
 *
 * max_cells defaults to 708588 everywhere, matching the postgis_hex9
 * hex9.grid_max_cells GUC; override per-call with the named parameter.
 */
#include "duckdb.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

#include "h9_uuid.hpp"
#include "h9_wkb.hpp"
#include "hex9_c.h"

#include <memory>
#include <vector>

namespace duckdb {

static constexpr int64_t H9_DEFAULT_MAX_CELLS = 708588;

static int64_t H9MaxCellsParam(const named_parameter_map_t &params, const char *fname) {
	auto it = params.find("max_cells");
	if (it == params.end()) {
		return H9_DEFAULT_MAX_CELLS;
	}
	auto v = it->second.GetValue<int64_t>();
	if (v <= 0) {
		throw BinderException("%s: max_cells must be > 0, got %lld", fname, (long long)v);
	}
	return v;
}

/* ── h9_grid(lon_min, lat_min, lon_max, lat_max, layer) ─────────────────── */

struct H9GridBindData : public FunctionData {
	double lon_min, lat_min, lon_max, lat_max;
	int32_t layer;
	int32_t densify = 0;
	int64_t max_cells = H9_DEFAULT_MAX_CELLS;

	unique_ptr<FunctionData> Copy() const override {
		return make_uniq<H9GridBindData>(*this);
	}
	bool Equals(const FunctionData &other_p) const override {
		auto &o = other_p.Cast<H9GridBindData>();
		return lon_min == o.lon_min && lat_min == o.lat_min && lon_max == o.lon_max && lat_max == o.lat_max &&
		       layer == o.layer && densify == o.densify && max_cells == o.max_cells;
	}
};

struct H9GridGlobalState : public GlobalTableFunctionState {
	hex9_grid *grid = nullptr;
	int count = 0;
	int next = 0;
	int n_ring = 0;
	std::vector<double> ring;

	~H9GridGlobalState() override {
		if (grid) {
			hex9_grid_destroy(grid);
		}
	}
	idx_t MaxThreads() const override {
		return 1;
	}
};

static unique_ptr<FunctionData> H9GridBind(ClientContext &, TableFunctionBindInput &input,
                                           vector<LogicalType> &return_types, vector<string> &names) {
	auto result = make_uniq<H9GridBindData>();
	for (auto &v : input.inputs) {
		if (v.IsNull()) {
			throw BinderException("h9_grid: arguments must not be NULL");
		}
	}
	result->lon_min = input.inputs[0].GetValue<double>();
	result->lat_min = input.inputs[1].GetValue<double>();
	result->lon_max = input.inputs[2].GetValue<double>();
	result->lat_max = input.inputs[3].GetValue<double>();
	result->layer = input.inputs[4].GetValue<int32_t>();
	if (result->layer < 1 || result->layer > hex9_lmax()) {
		throw BinderException("h9_grid: layer must be 1..%d, got %d", hex9_lmax(), result->layer);
	}
	auto densify_it = input.named_parameters.find("densify");
	if (densify_it != input.named_parameters.end()) {
		result->densify = densify_it->second.GetValue<int32_t>();
		if (result->densify < 0 || result->densify > 9 || result->layer + result->densify > hex9_lmax()) {
			throw BinderException("h9_grid: densify must be 0..9 with layer + densify <= %d", hex9_lmax());
		}
	}
	result->max_cells = H9MaxCellsParam(input.named_parameters, "h9_grid");

	return_types = {LogicalType::UUID, LogicalType::UUID, LogicalType::BLOB, LogicalType::BLOB};
	names = {"h9_id", "h9_bin", "geom", "centroid"};
	return std::move(result);
}

static unique_ptr<GlobalTableFunctionState> H9GridInit(ClientContext &, TableFunctionInitInput &input) {
	auto &bind_data = input.bind_data->Cast<H9GridBindData>();
	auto state = make_uniq<H9GridGlobalState>();
	char errbuf[256] = {0};
	state->grid = hex9_grid_create(bind_data.lon_min, bind_data.lat_min, bind_data.lon_max, bind_data.lat_max,
	                               bind_data.layer, bind_data.densify, bind_data.max_cells, errbuf, sizeof(errbuf));
	if (!state->grid) {
		throw InvalidInputException("h9_grid: %s", errbuf);
	}
	state->count = hex9_grid_count(state->grid);
	state->n_ring = hex9_ring_npoints(bind_data.densify);
	state->ring.resize(size_t(state->n_ring) * 2);
	return std::move(state);
}

static void H9GridFunction(ClientContext &, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<H9GridBindData>();
	auto &state = data_p.global_state->Cast<H9GridGlobalState>();
	auto id_data = FlatVector::GetData<hugeint_t>(output.data[0]);
	auto bin_data = FlatVector::GetData<hugeint_t>(output.data[1]);
	idx_t out = 0;
	while (out < STANDARD_VECTOR_SIZE && state.next < state.count) {
		const int i = state.next++;
		uint8_t u[16];
		hex9_grid_cell_id(state.grid, i, u);
		id_data[out] = H9BytesToUuid(u);
		hex9_grid_cell_uuid(state.grid, i, u);
		bin_data[out] = H9BytesToUuid(u);
		int n = hex9_grid_cell_ring(state.grid, i, bind_data.densify, state.ring.data(), state.n_ring);
		if (n != state.n_ring) {
			throw InternalException("h9_grid: ring build failed for cell %d", i);
		}
		FlatVector::GetData<string_t>(output.data[2])[out] = H9WkbPolygon(output.data[2], state.ring.data(), n);
		double lon, lat;
		hex9_grid_cell_centroid(state.grid, i, &lon, &lat);
		FlatVector::GetData<string_t>(output.data[3])[out] = H9WkbPoint(output.data[3], lon, lat);
		out++;
	}
	output.SetCardinality(out);
}

/* ── h9_curve_cells / h9_owned_cells (uuid, layer) ──────────────────────── */

using H9EnumFn = int64_t (*)(const uint8_t[16], int, uint8_t *, uint8_t *, int64_t);

struct H9EnumBindData : public FunctionData {
	uint8_t uuid[16];
	int32_t layer;
	int64_t max_cells;
	const char *fname;
	H9EnumFn fn;

	unique_ptr<FunctionData> Copy() const override {
		return make_uniq<H9EnumBindData>(*this);
	}
	bool Equals(const FunctionData &other_p) const override {
		auto &o = other_p.Cast<H9EnumBindData>();
		return std::memcmp(uuid, o.uuid, 16) == 0 && layer == o.layer && max_cells == o.max_cells && fn == o.fn;
	}
};

struct H9EnumGlobalState : public GlobalTableFunctionState {
	std::vector<uint8_t> bins;
	std::vector<uint8_t> curves;
	int64_t count = 0;
	int64_t next = 0;

	idx_t MaxThreads() const override {
		return 1;
	}
};

static unique_ptr<FunctionData> H9EnumBindCommon(TableFunctionBindInput &input, vector<LogicalType> &return_types,
                                                 vector<string> &names, const char *fname, H9EnumFn fn) {
	auto result = make_uniq<H9EnumBindData>();
	if (input.inputs[0].IsNull() || input.inputs[1].IsNull()) {
		throw BinderException("%s: arguments must not be NULL", fname);
	}
	H9UuidToBytes(input.inputs[0].GetValueUnsafe<hugeint_t>(), result->uuid);
	result->layer = input.inputs[1].GetValue<int32_t>();
	if (result->layer < 0 || result->layer > hex9_lmax()) {
		throw BinderException("%s: layer must be 0..%d, got %d", fname, hex9_lmax(), result->layer);
	}
	result->max_cells = H9MaxCellsParam(input.named_parameters, fname);
	result->fname = fname;
	result->fn = fn;
	return_types = {LogicalType::UUID, LogicalType::UUID};
	names = {"h9_bin", "h9_curve"};
	return std::move(result);
}

static unique_ptr<FunctionData> H9CurveCellsBind(ClientContext &, TableFunctionBindInput &input,
                                                 vector<LogicalType> &return_types, vector<string> &names) {
	return H9EnumBindCommon(input, return_types, names, "h9_curve_cells", hex9_curve_cells);
}

static unique_ptr<FunctionData> H9OwnedCellsBind(ClientContext &, TableFunctionBindInput &input,
                                                 vector<LogicalType> &return_types, vector<string> &names) {
	return H9EnumBindCommon(input, return_types, names, "h9_owned_cells", hex9_owned_cells);
}

static unique_ptr<GlobalTableFunctionState> H9EnumInit(ClientContext &, TableFunctionInitInput &input) {
	auto &bind_data = input.bind_data->Cast<H9EnumBindData>();
	auto state = make_uniq<H9EnumGlobalState>();
	/* Size the arrays from the exact descendant count (cheap arithmetic): the
	 * cell's layer is its curve address's layer, and 9^(to - from) follows. */
	uint8_t curve[16];
	if (hex9_curve(bind_data.uuid, curve) != 0) {
		throw InvalidInputException("%s: not a canonical H9 or curve UUID", bind_data.fname);
	}
	const int from_layer = hex9_curve_layer(curve);
	const int64_t expected = hex9_curve_ncells(from_layer, bind_data.layer);
	if (expected < 0) {
		throw InvalidInputException("%s: layer %d is above the cell's own layer %d, or the descendant count "
		                            "overflows",
		                            bind_data.fname, bind_data.layer, from_layer);
	}
	if (expected > bind_data.max_cells) {
		throw InvalidInputException("%s: %lld cells exceeds max_cells=%lld; use a coarser layer or raise "
		                            "max_cells",
		                            bind_data.fname, (long long)expected, (long long)bind_data.max_cells);
	}
	state->bins.resize(size_t(expected) * 16);
	state->curves.resize(size_t(expected) * 16);
	state->count = bind_data.fn(bind_data.uuid, bind_data.layer, state->bins.data(), state->curves.data(), expected);
	if (state->count < 0) {
		throw InvalidInputException("%s: enumeration failed (malformed UUID?)", bind_data.fname);
	}
	return std::move(state);
}

static void H9EnumFunction(ClientContext &, TableFunctionInput &data_p, DataChunk &output) {
	auto &state = data_p.global_state->Cast<H9EnumGlobalState>();
	auto bin_data = FlatVector::GetData<hugeint_t>(output.data[0]);
	auto curve_data = FlatVector::GetData<hugeint_t>(output.data[1]);
	idx_t out = 0;
	while (out < STANDARD_VECTOR_SIZE && state.next < state.count) {
		const auto i = size_t(state.next++);
		bin_data[out] = H9BytesToUuid(state.bins.data() + i * 16);
		curve_data[out] = H9BytesToUuid(state.curves.data() + i * 16);
		out++;
	}
	output.SetCardinality(out);
}

/* ── registration ───────────────────────────────────────────────────────── */

void RegisterH9Table(ExtensionLoader &loader) {
	const auto DBL = LogicalType::DOUBLE;
	const auto I32 = LogicalType::INTEGER;

	TableFunction grid("h9_grid", {DBL, DBL, DBL, DBL, I32}, H9GridFunction, H9GridBind, H9GridInit);
	grid.named_parameters["densify"] = I32;
	grid.named_parameters["max_cells"] = LogicalType::BIGINT;
	loader.RegisterFunction(grid);

	TableFunction curve_cells("h9_curve_cells", {LogicalType::UUID, I32}, H9EnumFunction, H9CurveCellsBind,
	                          H9EnumInit);
	curve_cells.named_parameters["max_cells"] = LogicalType::BIGINT;
	loader.RegisterFunction(curve_cells);

	TableFunction owned_cells("h9_owned_cells", {LogicalType::UUID, I32}, H9EnumFunction, H9OwnedCellsBind,
	                          H9EnumInit);
	owned_cells.named_parameters["max_cells"] = LogicalType::BIGINT;
	loader.RegisterFunction(owned_cells);
}

} // namespace duckdb
