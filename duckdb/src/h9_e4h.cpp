/* h9_e4h.cpp — the E4H aperture-4 structural tail and the grid verbs
 * (aim / walk_to / vision_cone). Semantics mirror extension/postgis_hex9:
 * E4H uuids are ADDRESSES under the universality contract (libhex9
 * docs/universality.md — the exact classifier mints bit-identical uuids on
 * every platform at every depth), the 0xE marker is decisive, and the h9_*
 * machinery rejects E4H input. The verbs are cell-first: bin KEYS in and
 * out (full-uuid re-derivation internal), bearings advisory FP in degrees
 * clockwise from north.
 */
#include "duckdb.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

#include "h9_uuid.hpp"
#include "h9_wkb.hpp"
#include "hex9_c.h"

#include <vector>

namespace duckdb {

/* ── E4H scalars ────────────────────────────────────────────────────────── */

static void H9eEncodeFn(DataChunk &args, ExpressionState &, Vector &result) {
	const auto count = args.size();
	UnifiedVectorFormat lon_uf, lat_uf, layer_uf, depth_uf;
	args.data[0].ToUnifiedFormat(count, lon_uf);
	args.data[1].ToUnifiedFormat(count, lat_uf);
	args.data[2].ToUnifiedFormat(count, layer_uf);
	args.data[3].ToUnifiedFormat(count, depth_uf);
	auto lon_d = UnifiedVectorFormat::GetData<double>(lon_uf);
	auto lat_d = UnifiedVectorFormat::GetData<double>(lat_uf);
	auto layer_d = UnifiedVectorFormat::GetData<int32_t>(layer_uf);
	auto depth_d = UnifiedVectorFormat::GetData<int32_t>(depth_uf);
	auto out = FlatVector::GetData<hugeint_t>(result);
	for (idx_t i = 0; i < count; i++) {
		auto oidx = lon_uf.sel->get_index(i);
		auto aidx = lat_uf.sel->get_index(i);
		auto lidx = layer_uf.sel->get_index(i);
		auto didx = depth_uf.sel->get_index(i);
		if (!lon_uf.validity.RowIsValid(oidx) || !lat_uf.validity.RowIsValid(aidx) ||
		    !layer_uf.validity.RowIsValid(lidx) || !depth_uf.validity.RowIsValid(didx)) {
			FlatVector::SetNull(result, i, true);
			continue;
		}
		const int32_t layer = layer_d[lidx], depth = depth_d[didx];
		uint8_t u[16];
		const int rc = hex9_e4h_encode(lon_d[oidx], lat_d[aidx], layer, depth, u);
		if (rc == 3) {
			throw InvalidInputException("h9e_encode: point two seams from its host");
		}
		if (rc != 0) {
			throw InvalidInputException(
			    "h9e_encode: need layer >= 0, depth >= 0 and layer + 2 + depth <= %d; got layer %d depth %d",
			    hex9_lmax(), layer, depth);
		}
		out[i] = H9BytesToUuid(u);
	}
	result.SetVectorType(VectorType::FLAT_VECTOR);
}

using H9ePointCFn = int (*)(const uint8_t[16], double *, double *);

static void H9ePointExec(DataChunk &args, Vector &result, H9ePointCFn fn, const char *fname) {
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
		if (fn(u, &lo, &la) != 0) {
			throw InvalidInputException("%s: not an E4H uuid (or bad tail grammar)", fname);
		}
		lon_d[i] = lo;
		lat_d[i] = la;
	}
}

static void H9eDecodeFn(DataChunk &args, ExpressionState &, Vector &result) {
	H9ePointExec(args, result, hex9_e4h_decode, "h9e_decode");
}

static void H9ePartnerFn(DataChunk &args, ExpressionState &, Vector &result) {
	H9ePointExec(args, result, hex9_e4h_partner, "h9e_partner");
}

static void H9eDecodeWkbFn(DataChunk &args, ExpressionState &, Vector &result) {
	UnaryExecutor::Execute<hugeint_t, string_t>(args.data[0], result, args.size(), [&](hugeint_t v) {
		uint8_t u[16];
		H9UuidToBytes(v, u);
		double lo, la;
		if (hex9_e4h_decode(u, &lo, &la) != 0) {
			throw InvalidInputException("h9e_decode_wkb: not an E4H uuid");
		}
		return H9WkbPoint(result, lo, la);
	});
}

static void H9eBinFn(DataChunk &args, ExpressionState &, Vector &result) {
	BinaryExecutor::Execute<hugeint_t, int32_t, hugeint_t>(
	    args.data[0], args.data[1], result, args.size(), [](hugeint_t v, int32_t depth) {
		    uint8_t u[16], out[16];
		    H9UuidToBytes(v, u);
		    if (hex9_e4h_bin(u, depth, out) != 0) {
			    throw InvalidInputException(
			        "h9e_bin: not an E4H uuid, or depth %d exceeds the tail (truncation only)", depth);
		    }
		    return H9BytesToUuid(out);
	    });
}

static void H9eDepthFn(DataChunk &args, ExpressionState &, Vector &result) {
	UnaryExecutor::Execute<hugeint_t, int32_t>(args.data[0], result, args.size(), [](hugeint_t v) {
		uint8_t u[16];
		H9UuidToBytes(v, u);
		return int32_t(hex9_e4h_depth(u));
	});
}

static void H9IsE4hFn(DataChunk &args, ExpressionState &, Vector &result) {
	UnaryExecutor::Execute<hugeint_t, bool>(args.data[0], result, args.size(), [](hugeint_t v) {
		uint8_t u[16];
		H9UuidToBytes(v, u);
		return hex9_is_e4h(u) != 0;
	});
}

static void H9eHostFn(DataChunk &args, ExpressionState &, Vector &result) {
	UnaryExecutor::Execute<hugeint_t, hugeint_t>(args.data[0], result, args.size(), [](hugeint_t v) {
		uint8_t u[16], host[16], digits[28];
		int half, nd;
		H9UuidToBytes(v, u);
		if (hex9_e4h_split(u, host, &half, digits, &nd) != 0) {
			throw InvalidInputException("h9e_host: not an E4H uuid");
		}
		return H9BytesToUuid(host);
	});
}

static void H9eLabelFn(DataChunk &args, ExpressionState &, Vector &result) {
	UnaryExecutor::Execute<hugeint_t, string_t>(args.data[0], result, args.size(), [&](hugeint_t v) {
		uint8_t u[16];
		H9UuidToBytes(v, u);
		char buf[64];
		const int n = hex9_e4h_label(u, buf, sizeof(buf));
		if (n < 0) {
			throw InvalidInputException("h9e_label: not an E4H uuid");
		}
		return StringVector::AddString(result, buf, idx_t(n));
	});
}

static void H9eParseLabelFn(DataChunk &args, ExpressionState &, Vector &result) {
	UnaryExecutor::Execute<string_t, hugeint_t>(args.data[0], result, args.size(), [](string_t s) {
		uint8_t out[16];
		if (hex9_e4h_parse_label(s.GetString().c_str(), out) != 0) {
			throw InvalidInputException("h9e_parse_label: not a valid E4H label: '%s'", s.GetString());
		}
		return H9BytesToUuid(out);
	});
}

/* hexagon binning: canonical key = the matched pair's MODE-0 member
 * (ruling 2026-08-05; mode = descent rotation parity, pure address
 * arithmetic). Both halves map to one key; every host owns exactly
 * 4^depth hexagons. */
static void H9eHexFn(DataChunk &args, ExpressionState &, Vector &result) {
	UnaryExecutor::Execute<hugeint_t, hugeint_t>(args.data[0], result, args.size(), [](hugeint_t v) {
		uint8_t u[16], out[16];
		H9UuidToBytes(v, u);
		if (hex9_e4h_hex(u, out) != 0) {
			throw InvalidInputException("h9e_hex: not an E4H uuid");
		}
		return H9BytesToUuid(out);
	});
}

static void H9eModeFn(DataChunk &args, ExpressionState &, Vector &result) {
	UnaryExecutor::Execute<hugeint_t, int32_t>(args.data[0], result, args.size(), [](hugeint_t v) {
		uint8_t u[16];
		H9UuidToBytes(v, u);
		return int32_t(hex9_e4h_mode(u));
	});
}

/* ── grid verbs ─────────────────────────────────────────────────────────── */

/* TRUE (geodesic, WGS84) azimuth <-> the verbs' bearing convention (the
 * flattening term at the latitude where the heading holds; identity at
 * poles/cardinals, ~0.19 deg max at equatorial diagonals). */
static void H9BearingFromTrueFn(DataChunk &args, ExpressionState &, Vector &result) {
	BinaryExecutor::Execute<double, double, double>(
	    args.data[0], args.data[1], result, args.size(),
	    [](double b, double lat) { return hex9_bearing_from_true(b, lat); });
}

static void H9BearingToTrueFn(DataChunk &args, ExpressionState &, Vector &result) {
	BinaryExecutor::Execute<double, double, double>(
	    args.data[0], args.data[1], result, args.size(),
	    [](double b, double lat) { return hex9_bearing_to_true(b, lat); });
}

static void H9AimFn(DataChunk &args, ExpressionState &, Vector &result) {
	TernaryExecutor::Execute<hugeint_t, int32_t, double, hugeint_t>(
	    args.data[0], args.data[1], args.data[2], result, args.size(),
	    [](hugeint_t v, int32_t layer, double bearing) {
		    uint8_t u[16], out[16];
		    H9UuidToBytes(v, u);
		    if (hex9_aim(u, layer, bearing, out) != 0) {
			    throw InvalidInputException("h9_aim: not a valid H9 key at layer %d (or non-finite bearing)",
			                                layer);
		    }
		    return H9BytesToUuid(out);
	    });
}

/* gather an optional LIST(UUID) argument's row entry into packed bytes */
static size_t H9GatherList(DataChunk &args, idx_t arg_no, idx_t row, std::vector<uint8_t> &bytes) {
	bytes.clear();
	if (args.ColumnCount() <= arg_no) {
		return 0;
	}
	auto &list = args.data[arg_no];
	UnifiedVectorFormat list_uf;
	list.ToUnifiedFormat(args.size(), list_uf);
	auto lidx = list_uf.sel->get_index(row);
	if (!list_uf.validity.RowIsValid(lidx)) {
		return 0;
	}
	auto entry = UnifiedVectorFormat::GetData<list_entry_t>(list_uf)[lidx];
	auto &elem = ListVector::GetEntry(list);
	const auto elem_count = ListVector::GetListSize(list);
	UnifiedVectorFormat elem_uf;
	elem.ToUnifiedFormat(elem_count, elem_uf);
	auto elem_data = UnifiedVectorFormat::GetData<hugeint_t>(elem_uf);
	bytes.resize(entry.length * 16);
	size_t m = 0;
	for (idx_t j = 0; j < entry.length; j++) {
		auto eidx = elem_uf.sel->get_index(entry.offset + j);
		if (!elem_uf.validity.RowIsValid(eidx)) {
			continue;
		}
		H9UuidToBytes(elem_data[eidx], bytes.data() + m * 16);
		m++;
	}
	return m;
}

static void H9WalkToFn(DataChunk &args, ExpressionState &, Vector &result) {
	const auto count = args.size();
	UnifiedVectorFormat src_uf, dst_uf, layer_uf;
	args.data[0].ToUnifiedFormat(count, src_uf);
	args.data[1].ToUnifiedFormat(count, dst_uf);
	args.data[2].ToUnifiedFormat(count, layer_uf);
	auto src_d = UnifiedVectorFormat::GetData<hugeint_t>(src_uf);
	auto dst_d = UnifiedVectorFormat::GetData<hugeint_t>(dst_uf);
	auto layer_d = UnifiedVectorFormat::GetData<int32_t>(layer_uf);
	auto list_entries = FlatVector::GetData<list_entry_t>(result);
	std::vector<uint8_t> obs, path;
	idx_t offset = 0;
	for (idx_t i = 0; i < count; i++) {
		auto sidx = src_uf.sel->get_index(i);
		auto didx = dst_uf.sel->get_index(i);
		auto lidx = layer_uf.sel->get_index(i);
		if (!src_uf.validity.RowIsValid(sidx) || !dst_uf.validity.RowIsValid(didx) ||
		    !layer_uf.validity.RowIsValid(lidx)) {
			FlatVector::SetNull(result, i, true);
			list_entries[i] = list_entry_t {offset, 0};
			continue;
		}
		uint8_t src[16], dst[16];
		H9UuidToBytes(src_d[sidx], src);
		H9UuidToBytes(dst_d[didx], dst);
		const size_t nobs = H9GatherList(args, 3, i, obs);
		const int64_t max_expand = 5000;
		path.resize(size_t(max_expand + 2) * 16);
		const int64_t n = hex9_walk_to(src, dst, layer_d[lidx], obs.data(), nobs, max_expand, path.data(),
		                               max_expand + 2);
		if (n < 0) {
			throw InvalidInputException("h9_walk_to: invalid keys or layer %d", layer_d[lidx]);
		}
		if (n == 0) { /* no path within budget: NULL, not an error */
			FlatVector::SetNull(result, i, true);
			list_entries[i] = list_entry_t {offset, 0};
			continue;
		}
		ListVector::Reserve(result, offset + idx_t(n));
		auto out_data = FlatVector::GetData<hugeint_t>(ListVector::GetEntry(result));
		for (int64_t j = 0; j < n; j++) {
			out_data[offset + idx_t(j)] = H9BytesToUuid(path.data() + j * 16);
		}
		list_entries[i] = list_entry_t {offset, idx_t(n)};
		offset += idx_t(n);
	}
	ListVector::SetListSize(result, offset);
}

static void H9VisionConeFn(DataChunk &args, ExpressionState &, Vector &result) {
	const auto count = args.size();
	UnifiedVectorFormat u_uf, layer_uf, bear_uf, half_uf, k_uf;
	args.data[0].ToUnifiedFormat(count, u_uf);
	args.data[1].ToUnifiedFormat(count, layer_uf);
	args.data[2].ToUnifiedFormat(count, bear_uf);
	args.data[3].ToUnifiedFormat(count, half_uf);
	args.data[4].ToUnifiedFormat(count, k_uf);
	auto u_d = UnifiedVectorFormat::GetData<hugeint_t>(u_uf);
	auto layer_d = UnifiedVectorFormat::GetData<int32_t>(layer_uf);
	auto bear_d = UnifiedVectorFormat::GetData<double>(bear_uf);
	auto half_d = UnifiedVectorFormat::GetData<double>(half_uf);
	auto k_d = UnifiedVectorFormat::GetData<int32_t>(k_uf);
	auto list_entries = FlatVector::GetData<list_entry_t>(result);
	std::vector<uint8_t> obs, buf;
	idx_t offset = 0;
	for (idx_t i = 0; i < count; i++) {
		auto uidx = u_uf.sel->get_index(i);
		auto lidx = layer_uf.sel->get_index(i);
		auto bidx = bear_uf.sel->get_index(i);
		auto hidx = half_uf.sel->get_index(i);
		auto kidx = k_uf.sel->get_index(i);
		if (!u_uf.validity.RowIsValid(uidx) || !layer_uf.validity.RowIsValid(lidx) ||
		    !bear_uf.validity.RowIsValid(bidx) || !half_uf.validity.RowIsValid(hidx) ||
		    !k_uf.validity.RowIsValid(kidx)) {
			FlatVector::SetNull(result, i, true);
			list_entries[i] = list_entry_t {offset, 0};
			continue;
		}
		const int32_t k = k_d[kidx];
		if (k < 0 || hex9_disk_ncells(k) > 60000000) {
			throw OutOfRangeException("h9_vision_cone: k=%d out of range", k);
		}
		uint8_t u[16];
		H9UuidToBytes(u_d[uidx], u);
		const size_t nobs = H9GatherList(args, 5, i, obs);
		const int64_t cap = hex9_disk_ncells(k);
		buf.resize(size_t(cap) * 16);
		const int64_t n = hex9_vision_cone(u, layer_d[lidx], bear_d[bidx], half_d[hidx], k, obs.data(), nobs,
		                                   buf.data(), cap);
		if (n < 0) {
			throw InvalidInputException("h9_vision_cone: invalid key, layer %d, or non-finite angles",
			                            layer_d[lidx]);
		}
		ListVector::Reserve(result, offset + idx_t(n));
		auto out_data = FlatVector::GetData<hugeint_t>(ListVector::GetEntry(result));
		for (int64_t j = 0; j < n; j++) {
			out_data[offset + idx_t(j)] = H9BytesToUuid(buf.data() + j * 16);
		}
		list_entries[i] = list_entry_t {offset, idx_t(n)};
		offset += idx_t(n);
	}
	ListVector::SetListSize(result, offset);
}

/* ── registration ───────────────────────────────────────────────────────── */

void RegisterH9E4h(ExtensionLoader &loader) {
	const auto UUID = LogicalType::UUID;
	const auto I32 = LogicalType::INTEGER;
	const auto DBL = LogicalType::DOUBLE;
	const auto STR = LogicalType::VARCHAR;
	const auto BLOB = LogicalType::BLOB;
	const auto lonlat_struct = LogicalType::STRUCT({{"lon", DBL}, {"lat", DBL}});

	loader.RegisterFunction(ScalarFunction("h9e_encode", {DBL, DBL, I32, I32}, UUID, H9eEncodeFn));
	loader.RegisterFunction(ScalarFunction("h9e_decode", {UUID}, lonlat_struct, H9eDecodeFn));
	loader.RegisterFunction(ScalarFunction("h9e_decode_wkb", {UUID}, BLOB, H9eDecodeWkbFn));
	loader.RegisterFunction(ScalarFunction("h9e_partner", {UUID}, lonlat_struct, H9ePartnerFn));
	loader.RegisterFunction(ScalarFunction("h9e_bin", {UUID, I32}, UUID, H9eBinFn));
	loader.RegisterFunction(ScalarFunction("h9e_depth", {UUID}, I32, H9eDepthFn));
	loader.RegisterFunction(ScalarFunction("h9_is_e4h", {UUID}, LogicalType::BOOLEAN, H9IsE4hFn));
	loader.RegisterFunction(ScalarFunction("h9e_host", {UUID}, UUID, H9eHostFn));
	loader.RegisterFunction(ScalarFunction("h9e_label", {UUID}, STR, H9eLabelFn));
	loader.RegisterFunction(ScalarFunction("h9e_parse_label", {STR}, UUID, H9eParseLabelFn));
	loader.RegisterFunction(ScalarFunction("h9e_hex", {UUID}, UUID, H9eHexFn));
	loader.RegisterFunction(ScalarFunction("h9e_mode", {UUID}, I32, H9eModeFn));

	loader.RegisterFunction(ScalarFunction("h9_bearing_from_true", {DBL, DBL}, DBL, H9BearingFromTrueFn));
	loader.RegisterFunction(ScalarFunction("h9_bearing_to_true", {DBL, DBL}, DBL, H9BearingToTrueFn));
	loader.RegisterFunction(ScalarFunction("h9_aim", {UUID, I32, DBL}, UUID, H9AimFn));

	ScalarFunctionSet walk("h9_walk_to");
	walk.AddFunction(ScalarFunction({UUID, UUID, I32}, LogicalType::LIST(UUID), H9WalkToFn));
	walk.AddFunction(ScalarFunction({UUID, UUID, I32, LogicalType::LIST(UUID)}, LogicalType::LIST(UUID), H9WalkToFn));
	loader.RegisterFunction(walk);

	ScalarFunctionSet cone("h9_vision_cone");
	cone.AddFunction(ScalarFunction({UUID, I32, DBL, DBL, I32}, LogicalType::LIST(UUID), H9VisionConeFn));
	cone.AddFunction(
	    ScalarFunction({UUID, I32, DBL, DBL, I32, LogicalType::LIST(UUID)}, LogicalType::LIST(UUID), H9VisionConeFn));
	loader.RegisterFunction(cone);
}

} // namespace duckdb
