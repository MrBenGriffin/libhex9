/* h9_scalar.cpp — core addressing surface: encode/decode/bin, cell geometry,
 * labels, lineage, adjacency, and the UUID<->UHUGEINT reinterpretation pair.
 * Semantics (bounds, messages, doctrine guards) mirror extension/postgis_hex9.
 */
#include "duckdb.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

#include "h9_uuid.hpp"
#include "h9_wkb.hpp"
#include "hex9_c.h"

#include <memory>
#include <vector>

namespace duckdb {

/* ── version / lmax ─────────────────────────────────────────────────────── */

static void H9VersionFn(DataChunk &args, ExpressionState &, Vector &result) {
	result.SetVectorType(VectorType::CONSTANT_VECTOR);
	auto data = ConstantVector::GetData<string_t>(result);
	data[0] = StringVector::AddString(result, hex9_version());
}

static void H9LmaxFn(DataChunk &args, ExpressionState &, Vector &result) {
	result.SetVectorType(VectorType::CONSTANT_VECTOR);
	ConstantVector::GetData<int32_t>(result)[0] = hex9_lmax();
}

/* ── encode / decode / bin ──────────────────────────────────────────────── */

static void H9EncodeFn(DataChunk &args, ExpressionState &, Vector &result) {
	const auto count = args.size();
	auto &lon = args.data[0];
	auto &lat = args.data[1];
	/* Flat fast path: one hex9_encode_many call per chunk. */
	if (lon.GetVectorType() == VectorType::FLAT_VECTOR && lat.GetVectorType() == VectorType::FLAT_VECTOR &&
	    FlatVector::Validity(lon).AllValid() && FlatVector::Validity(lat).AllValid()) {
		auto lon_d = FlatVector::GetData<double>(lon);
		auto lat_d = FlatVector::GetData<double>(lat);
		std::unique_ptr<uint8_t[]> buf(new uint8_t[count * 16]);
		if (hex9_encode_many(lon_d, lat_d, count, buf.get()) == 0) {
			auto out = FlatVector::GetData<hugeint_t>(result);
			for (idx_t i = 0; i < count; i++) {
				out[i] = H9BytesToUuid(buf.get() + i * 16);
			}
			result.SetVectorType(VectorType::FLAT_VECTOR);
			return;
		}
		/* Batch failed: fall through per-row for a precise error. */
	}
	BinaryExecutor::Execute<double, double, hugeint_t>(lon, lat, result, count, [](double lo, double la) {
		uint8_t u[16];
		if (hex9_encode(lo, la, u) != 0) {
			throw InvalidInputException("h9_encode: invalid coordinate (lon %f, lat %f)", lo, la);
		}
		return H9BytesToUuid(u);
	});
}

static void H9DecodeFn(DataChunk &args, ExpressionState &, Vector &result) {
	const auto count = args.size();
	UnifiedVectorFormat uf;
	args.data[0].ToUnifiedFormat(count, uf);
	auto in = UnifiedVectorFormat::GetData<hugeint_t>(uf);
	auto &entries = StructVector::GetEntries(result);
	auto lon_d = FlatVector::GetData<double>(*entries[0]);
	auto lat_d = FlatVector::GetData<double>(*entries[1]);
	for (idx_t i = 0; i < count; i++) {
		auto idx = uf.sel->get_index(i);
		if (!uf.validity.RowIsValid(idx)) {
			FlatVector::SetNull(result, i, true);
			continue;
		}
		uint8_t u[16];
		H9UuidToBytes(in[idx], u);
		double lo, la;
		if (hex9_decode(u, &lo, &la) != 0) {
			throw InvalidInputException("h9_decode: not a valid H9 UUID");
		}
		lon_d[i] = lo;
		lat_d[i] = la;
	}
}

static void H9DecodeWkbFn(DataChunk &args, ExpressionState &, Vector &result) {
	UnaryExecutor::Execute<hugeint_t, string_t>(args.data[0], result, args.size(), [&](hugeint_t v) {
		uint8_t u[16];
		H9UuidToBytes(v, u);
		double lo, la;
		if (hex9_decode(u, &lo, &la) != 0) {
			throw InvalidInputException("h9_decode_wkb: not a valid H9 UUID");
		}
		return H9WkbPoint(result, lo, la);
	});
}

static void H9BinFn(DataChunk &args, ExpressionState &, Vector &result) {
	BinaryExecutor::Execute<hugeint_t, int32_t, hugeint_t>(
	    args.data[0], args.data[1], result, args.size(), [](hugeint_t v, int32_t layer) {
		    if (layer < 0 || layer > hex9_lmax()) {
			    throw InvalidInputException("h9_bin: layer must be 0..%d, got %d", hex9_lmax(), layer);
		    }
		    uint8_t u[16], out[16];
		    H9UuidToBytes(v, u);
		    if (hex9_bin(u, layer, out) != 0) {
			    throw InvalidInputException("h9_bin: not a valid H9 UUID at layer %d", layer);
		    }
		    return H9BytesToUuid(out);
	    });
}

/* ── cell polygon ───────────────────────────────────────────────────────── */

static string_t H9CellRing(Vector &result, hugeint_t v, int32_t layer, int32_t densify) {
	if (layer < 1 || layer > hex9_lmax()) {
		throw InvalidInputException("h9_cell: layer must be 1..%d, got %d", hex9_lmax(), layer);
	}
	if (densify < 0) {
		throw InvalidInputException("h9_cell: densify must be >= 0, got %d", densify);
	}
	if (layer + densify > hex9_lmax()) {
		throw InvalidInputException("h9_cell: layer + densify must be <= %d (layer=%d, densify=%d -> %d)",
		                            hex9_lmax(), layer, densify, layer + densify);
	}
	/* Hard cap to keep ring sizes sane (6*3^densify+1 = 118099 at densify=9). */
	if (densify > 9) {
		throw InvalidInputException("h9_cell: densify must be <= 9, got %d", densify);
	}
	const int n_ring = hex9_ring_npoints(densify);
	std::vector<double> lonlat(size_t(n_ring) * 2);
	uint8_t u[16];
	H9UuidToBytes(v, u);
	if (hex9_cell_ring(u, layer, densify, lonlat.data(), n_ring) != n_ring) {
		throw InvalidInputException("h9_cell: not a valid H9 UUID at layer %d", layer);
	}
	return H9WkbPolygon(result, lonlat.data(), n_ring);
}

static void H9CellFn(DataChunk &args, ExpressionState &, Vector &result) {
	BinaryExecutor::Execute<hugeint_t, int32_t, string_t>(
	    args.data[0], args.data[1], result, args.size(),
	    [&](hugeint_t v, int32_t layer) { return H9CellRing(result, v, layer, 0); });
}

static void H9CellDensifyFn(DataChunk &args, ExpressionState &, Vector &result) {
	TernaryExecutor::Execute<hugeint_t, int32_t, int32_t, string_t>(
	    args.data[0], args.data[1], args.data[2], result, args.size(),
	    [&](hugeint_t v, int32_t layer, int32_t densify) { return H9CellRing(result, v, layer, densify); });
}

/* ── labels ─────────────────────────────────────────────────────────────── */

using H9LabelCFn = int (*)(const uint8_t[16], int, char *, size_t);

static void H9LabelExec(DataChunk &args, Vector &result, H9LabelCFn fn, const char *fname) {
	BinaryExecutor::Execute<hugeint_t, int32_t, string_t>(
	    args.data[0], args.data[1], result, args.size(), [&](hugeint_t v, int32_t layer) {
		    if (layer < 0 || layer > hex9_lmax()) {
			    throw InvalidInputException("%s: layer must be 0..%d, got %d", fname, hex9_lmax(), layer);
		    }
		    uint8_t u[16];
		    H9UuidToBytes(v, u);
		    char buf[48];
		    int n = fn(u, layer, buf, sizeof(buf));
		    if (n < 0) {
			    throw InvalidInputException("%s: not a valid H9 UUID at layer %d", fname, layer);
		    }
		    return StringVector::AddString(result, buf, idx_t(n));
	    });
}

static void H9LabelFn(DataChunk &args, ExpressionState &, Vector &result) {
	H9LabelExec(args, result, hex9_label, "h9_label");
}

static void H9LabelKeyFn(DataChunk &args, ExpressionState &, Vector &result) {
	H9LabelExec(args, result, hex9_label_key, "h9_label_key");
}

static void H9ParseLabelFn(DataChunk &args, ExpressionState &, Vector &result) {
	UnaryExecutor::Execute<string_t, hugeint_t>(args.data[0], result, args.size(), [](string_t s) {
		uint8_t out[16];
		if (hex9_parse_label(s.GetString().c_str(), out) < 0) {
			throw InvalidInputException("h9_parse_label: not a valid H9 label: '%s'", s.GetString());
		}
		return H9BytesToUuid(out);
	});
}

static void H9LabelCentroidFn(DataChunk &args, ExpressionState &, Vector &result) {
	UnaryExecutor::Execute<string_t, string_t>(args.data[0], result, args.size(), [&](string_t s) {
		double lo, la;
		if (hex9_label_centroid(s.GetString().c_str(), &lo, &la) != 0) {
			throw InvalidInputException("h9_label_centroid: not a valid H9 label: '%s'", s.GetString());
		}
		return H9WkbPoint(result, lo, la);
	});
}

/* ── lineage ────────────────────────────────────────────────────────────── */

static void H9CellParentFn(DataChunk &args, ExpressionState &, Vector &result) {
	UnaryExecutor::Execute<hugeint_t, hugeint_t>(args.data[0], result, args.size(), [](hugeint_t v) {
		uint8_t u[16], out[16];
		H9UuidToBytes(v, u);
		if (hex9_cell_parent(u, out) != 0) {
			throw InvalidInputException("h9_cell_parent: no parent (L0 cell) or malformed H9 UUID");
		}
		return H9BytesToUuid(out);
	});
}

static void H9CellAncestorFn(DataChunk &args, ExpressionState &, Vector &result) {
	BinaryExecutor::Execute<hugeint_t, int32_t, hugeint_t>(
	    args.data[0], args.data[1], result, args.size(), [](hugeint_t v, int32_t layer) {
		    if (layer < 0 || layer > hex9_lmax()) {
			    throw InvalidInputException("h9_cell_ancestor: layer must be 0..%d, got %d", hex9_lmax(), layer);
		    }
		    uint8_t u[16], out[16];
		    H9UuidToBytes(v, u);
		    if (hex9_cell_ancestor(u, layer, out) != 0) {
			    throw InvalidInputException(
			        "h9_cell_ancestor: input coarser than layer %d, or malformed H9 UUID", layer);
		    }
		    return H9BytesToUuid(out);
	    });
}

static void H9CellChildrenFn(DataChunk &args, ExpressionState &, Vector &result) {
	const auto count = args.size();
	UnifiedVectorFormat uf;
	args.data[0].ToUnifiedFormat(count, uf);
	auto in = UnifiedVectorFormat::GetData<hugeint_t>(uf);
	auto list_entries = FlatVector::GetData<list_entry_t>(result);
	ListVector::Reserve(result, count * 9);
	auto child_data = FlatVector::GetData<hugeint_t>(ListVector::GetEntry(result));
	idx_t offset = 0;
	for (idx_t i = 0; i < count; i++) {
		auto idx = uf.sel->get_index(i);
		if (!uf.validity.RowIsValid(idx)) {
			FlatVector::SetNull(result, i, true);
			list_entries[i] = list_entry_t {offset, 0};
			continue;
		}
		uint8_t u[16], out[9 * 16];
		H9UuidToBytes(in[idx], u);
		if (hex9_cell_children(u, out) != 0) {
			throw InvalidInputException("h9_cell_children: cell at lmax (no children) or malformed H9 UUID");
		}
		for (idx_t j = 0; j < 9; j++) {
			child_data[offset + j] = H9BytesToUuid(out + j * 16);
		}
		list_entries[i] = list_entry_t {offset, 9};
		offset += 9;
	}
	ListVector::SetListSize(result, offset);
}

static void H9CommonAncestorFn(DataChunk &args, ExpressionState &, Vector &result) {
	const auto count = args.size();
	auto &list = args.data[0];
	UnifiedVectorFormat list_uf;
	list.ToUnifiedFormat(count, list_uf);
	auto list_data = UnifiedVectorFormat::GetData<list_entry_t>(list_uf);
	auto &elem = ListVector::GetEntry(list);
	const auto elem_count = ListVector::GetListSize(list);
	UnifiedVectorFormat elem_uf;
	elem.ToUnifiedFormat(elem_count, elem_uf);
	auto elem_data = UnifiedVectorFormat::GetData<hugeint_t>(elem_uf);

	UnifiedVectorFormat layer_uf;
	args.data[1].ToUnifiedFormat(count, layer_uf);
	auto layer_data = UnifiedVectorFormat::GetData<int32_t>(layer_uf);

	auto &entries = StructVector::GetEntries(result);
	auto label_vec = entries[0].get();
	auto bin_data = FlatVector::GetData<hugeint_t>(*entries[1]);
	auto out_layer_data = FlatVector::GetData<int32_t>(*entries[2]);
	auto label_data = FlatVector::GetData<string_t>(*label_vec);

	std::vector<uint8_t> staging;
	for (idx_t i = 0; i < count; i++) {
		auto lidx = list_uf.sel->get_index(i);
		auto layer_idx = layer_uf.sel->get_index(i);
		if (!list_uf.validity.RowIsValid(lidx) || !layer_uf.validity.RowIsValid(layer_idx)) {
			FlatVector::SetNull(result, i, true);
			continue;
		}
		const int32_t layer = layer_data[layer_idx];
		if (layer < 0 || layer > hex9_lmax()) {
			throw InvalidInputException("h9_common_ancestor: layer must be 0..%d, got %d", hex9_lmax(), layer);
		}
		const auto entry = list_data[lidx];
		if (entry.length == 0) {
			FlatVector::SetNull(result, i, true);
			continue;
		}
		staging.resize(entry.length * 16);
		for (idx_t j = 0; j < entry.length; j++) {
			auto eidx = elem_uf.sel->get_index(entry.offset + j);
			if (!elem_uf.validity.RowIsValid(eidx)) {
				throw InvalidInputException("h9_common_ancestor: NULL element in UUID list");
			}
			H9UuidToBytes(elem_data[eidx], staging.data() + j * 16);
		}
		char buf[48];
		uint8_t out_uuid[16];
		int anc = hex9_common_ancestor(staging.data(), entry.length, layer, buf, sizeof(buf), out_uuid);
		if (anc < 0) {
			/* No common ancestor (cells span L0 hexes): NULL, not an error. */
			FlatVector::SetNull(result, i, true);
			continue;
		}
		label_data[i] = StringVector::AddString(*label_vec, buf);
		bin_data[i] = H9BytesToUuid(out_uuid);
		out_layer_data[i] = anc;
	}
}

/* ── adjacency (full-UUID input only; bins rejected — doctrine) ─────────── */

static void H9NeighborsFn(DataChunk &args, ExpressionState &, Vector &result) {
	const auto count = args.size();
	UnifiedVectorFormat uuid_uf, layer_uf;
	args.data[0].ToUnifiedFormat(count, uuid_uf);
	args.data[1].ToUnifiedFormat(count, layer_uf);
	auto in = UnifiedVectorFormat::GetData<hugeint_t>(uuid_uf);
	auto layer_data = UnifiedVectorFormat::GetData<int32_t>(layer_uf);
	auto list_entries = FlatVector::GetData<list_entry_t>(result);
	ListVector::Reserve(result, count * 6);
	auto child_data = FlatVector::GetData<hugeint_t>(ListVector::GetEntry(result));
	idx_t offset = 0;
	for (idx_t i = 0; i < count; i++) {
		auto uidx = uuid_uf.sel->get_index(i);
		auto lidx = layer_uf.sel->get_index(i);
		if (!uuid_uf.validity.RowIsValid(uidx) || !layer_uf.validity.RowIsValid(lidx)) {
			FlatVector::SetNull(result, i, true);
			list_entries[i] = list_entry_t {offset, 0};
			continue;
		}
		const int32_t layer = layer_data[lidx];
		if (layer < 1 || layer > hex9_lmax()) {
			throw InvalidInputException("h9_neighbors: layer must be 1..%d, got %d", hex9_lmax(), layer);
		}
		uint8_t u[16], out[6 * 16];
		H9UuidToBytes(in[uidx], u);
		H9RejectBinInput(u, "h9_neighbors");
		int n = hex9_neighbors(u, layer, out);
		if (n < 0) {
			throw InvalidInputException("h9_neighbors: not a valid H9 UUID at layer %d", layer);
		}
		for (int j = 0; j < n; j++) {
			child_data[offset + idx_t(j)] = H9BytesToUuid(out + j * 16);
		}
		list_entries[i] = list_entry_t {offset, idx_t(n)};
		offset += idx_t(n);
	}
	ListVector::SetListSize(result, offset);
}

using H9KFillFn = int64_t (*)(const uint8_t[16], int, int, uint8_t *, int64_t);

static void H9KExec(DataChunk &args, Vector &result, const char *fname, H9KFillFn fn) {
	const auto count = args.size();
	UnifiedVectorFormat uuid_uf, layer_uf, k_uf;
	args.data[0].ToUnifiedFormat(count, uuid_uf);
	args.data[1].ToUnifiedFormat(count, layer_uf);
	args.data[2].ToUnifiedFormat(count, k_uf);
	auto in = UnifiedVectorFormat::GetData<hugeint_t>(uuid_uf);
	auto layer_data = UnifiedVectorFormat::GetData<int32_t>(layer_uf);
	auto k_data = UnifiedVectorFormat::GetData<int32_t>(k_uf);
	auto list_entries = FlatVector::GetData<list_entry_t>(result);
	std::vector<uint8_t> buf;
	idx_t offset = 0;
	for (idx_t i = 0; i < count; i++) {
		auto uidx = uuid_uf.sel->get_index(i);
		auto lidx = layer_uf.sel->get_index(i);
		auto kidx = k_uf.sel->get_index(i);
		if (!uuid_uf.validity.RowIsValid(uidx) || !layer_uf.validity.RowIsValid(lidx) ||
		    !k_uf.validity.RowIsValid(kidx)) {
			FlatVector::SetNull(result, i, true);
			list_entries[i] = list_entry_t {offset, 0};
			continue;
		}
		const int32_t layer = layer_data[lidx];
		const int32_t k = k_data[kidx];
		if (layer < 1 || layer > hex9_lmax()) {
			throw InvalidInputException("%s: layer must be 1..%d, got %d", fname, hex9_lmax(), layer);
		}
		if (k < 0) {
			throw InvalidInputException("%s: k must be >= 0, got %d", fname, k);
		}
		uint8_t u[16];
		H9UuidToBytes(in[uidx], u);
		H9RejectBinInput(u, fname);
		/* Nominal disk size 1 + 3k(k+1) is an upper bound for both ring and disk. */
		const int64_t ncells = hex9_disk_ncells(k);
		if (ncells > 60000000) {
			throw OutOfRangeException("%s: k=%d implies up to %lld cells; use a smaller k", fname, k,
			                          (long long)ncells);
		}
		buf.resize(size_t(ncells) * 16);
		int64_t n = fn(u, layer, k, buf.data(), ncells);
		if (n < 0) {
			throw InvalidInputException("%s: not a valid H9 UUID at layer %d", fname, layer);
		}
		ListVector::Reserve(result, offset + idx_t(n));
		auto child_data = FlatVector::GetData<hugeint_t>(ListVector::GetEntry(result));
		for (int64_t j = 0; j < n; j++) {
			child_data[offset + idx_t(j)] = H9BytesToUuid(buf.data() + j * 16);
		}
		list_entries[i] = list_entry_t {offset, idx_t(n)};
		offset += idx_t(n);
	}
	ListVector::SetListSize(result, offset);
}

static void H9KringFn(DataChunk &args, ExpressionState &, Vector &result) {
	H9KExec(args, result, "h9_kring", hex9_k_ring);
}

static void H9KdiskFn(DataChunk &args, ExpressionState &, Vector &result) {
	H9KExec(args, result, "h9_kdisk", hex9_k_disk);
}

/* ── UUID <-> UHUGEINT reinterpretation ─────────────────────────────────── */

static void H9IdIntFn(DataChunk &args, ExpressionState &, Vector &result) {
	UnaryExecutor::Execute<hugeint_t, uhugeint_t>(args.data[0], result, args.size(), [](hugeint_t v) {
		uint8_t u[16];
		H9UuidToBytes(v, u);
		return H9BytesToUhuge(u);
	});
}

static void H9IdFromIntFn(DataChunk &args, ExpressionState &, Vector &result) {
	UnaryExecutor::Execute<uhugeint_t, hugeint_t>(args.data[0], result, args.size(), [](uhugeint_t v) {
		uint8_t u[16];
		H9UhugeToBytes(v, u);
		return H9BytesToUuid(u);
	});
}

/* ── registration ───────────────────────────────────────────────────────── */

void RegisterH9Scalar(ExtensionLoader &loader) {
	const auto UUID = LogicalType::UUID;
	const auto I32 = LogicalType::INTEGER;
	const auto DBL = LogicalType::DOUBLE;
	const auto STR = LogicalType::VARCHAR;
	const auto BLOB = LogicalType::BLOB;
	const auto lonlat_struct = LogicalType::STRUCT({{"lon", DBL}, {"lat", DBL}});
	const auto ancestor_struct = LogicalType::STRUCT({{"label", STR}, {"h9_bin", UUID}, {"layer", I32}});

	loader.RegisterFunction(ScalarFunction("h9_version", {}, STR, H9VersionFn));
	loader.RegisterFunction(ScalarFunction("h9_lmax", {}, I32, H9LmaxFn));
	loader.RegisterFunction(ScalarFunction("h9_encode", {DBL, DBL}, UUID, H9EncodeFn));
	loader.RegisterFunction(ScalarFunction("h9_decode", {UUID}, lonlat_struct, H9DecodeFn));
	loader.RegisterFunction(ScalarFunction("h9_decode_wkb", {UUID}, BLOB, H9DecodeWkbFn));
	loader.RegisterFunction(ScalarFunction("h9_bin", {UUID, I32}, UUID, H9BinFn));

	ScalarFunctionSet cell("h9_cell");
	cell.AddFunction(ScalarFunction({UUID, I32}, BLOB, H9CellFn));
	cell.AddFunction(ScalarFunction({UUID, I32, I32}, BLOB, H9CellDensifyFn));
	loader.RegisterFunction(cell);

	loader.RegisterFunction(ScalarFunction("h9_label", {UUID, I32}, STR, H9LabelFn));
	loader.RegisterFunction(ScalarFunction("h9_label_key", {UUID, I32}, STR, H9LabelKeyFn));
	loader.RegisterFunction(ScalarFunction("h9_parse_label", {STR}, UUID, H9ParseLabelFn));
	loader.RegisterFunction(ScalarFunction("h9_label_centroid", {STR}, BLOB, H9LabelCentroidFn));

	loader.RegisterFunction(ScalarFunction("h9_cell_parent", {UUID}, UUID, H9CellParentFn));
	loader.RegisterFunction(ScalarFunction("h9_cell_ancestor", {UUID, I32}, UUID, H9CellAncestorFn));
	loader.RegisterFunction(ScalarFunction("h9_cell_children", {UUID}, LogicalType::LIST(UUID), H9CellChildrenFn));
	loader.RegisterFunction(
	    ScalarFunction("h9_common_ancestor", {LogicalType::LIST(UUID), I32}, ancestor_struct, H9CommonAncestorFn));

	loader.RegisterFunction(ScalarFunction("h9_neighbors", {UUID, I32}, LogicalType::LIST(UUID), H9NeighborsFn));
	loader.RegisterFunction(ScalarFunction("h9_kring", {UUID, I32, I32}, LogicalType::LIST(UUID), H9KringFn));
	loader.RegisterFunction(ScalarFunction("h9_kdisk", {UUID, I32, I32}, LogicalType::LIST(UUID), H9KdiskFn));

	loader.RegisterFunction(ScalarFunction("h9_id_int", {UUID}, LogicalType::UHUGEINT, H9IdIntFn));
	loader.RegisterFunction(ScalarFunction("h9_id_from_int", {LogicalType::UHUGEINT}, UUID, H9IdFromIntFn));
}

} // namespace duckdb
