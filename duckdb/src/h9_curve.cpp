/* h9_curve.cpp — Hamiltonian curve surface. Curve indexes are HUGEINT
 * (int128): they reach 12*9^L - 1, overflowing int64 above L18 but fitting
 * int128 with room (the PostGIS shim needs numeric for the same values).
 */
#include "duckdb.hpp"
#include "duckdb/common/operator/cast_operators.hpp"
#include "duckdb/common/types/hugeint.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

#include "h9_uuid.hpp"
#include "hex9_c.h"

namespace duckdb {

static void H9CurveFn(DataChunk &args, ExpressionState &, Vector &result) {
	UnaryExecutor::Execute<hugeint_t, hugeint_t>(args.data[0], result, args.size(), [](hugeint_t v) {
		uint8_t u[16], out[16];
		H9UuidToBytes(v, u);
		if (hex9_curve(u, out) != 0) {
			throw InvalidInputException("h9_curve: not a canonical H9 UUID");
		}
		return H9BytesToUuid(out);
	});
}

static void H9CurveDecodeFn(DataChunk &args, ExpressionState &, Vector &result) {
	UnaryExecutor::Execute<hugeint_t, hugeint_t>(args.data[0], result, args.size(), [](hugeint_t v) {
		uint8_t u[16], out[16];
		H9UuidToBytes(v, u);
		if (hex9_curve_decode(u, out) != 0) {
			throw InvalidInputException("h9_curve_decode: not a valid curve UUID");
		}
		return H9BytesToUuid(out);
	});
}

static void H9IsCurveFn(DataChunk &args, ExpressionState &, Vector &result) {
	UnaryExecutor::Execute<hugeint_t, bool>(args.data[0], result, args.size(), [](hugeint_t v) {
		uint8_t u[16];
		H9UuidToBytes(v, u);
		return hex9_is_curve(u) != 0;
	});
}

static void H9CurveLayerFn(DataChunk &args, ExpressionState &, Vector &result) {
	UnaryExecutor::Execute<hugeint_t, int32_t>(args.data[0], result, args.size(), [](hugeint_t v) {
		uint8_t u[16];
		H9UuidToBytes(v, u);
		int layer = hex9_curve_layer(u);
		if (layer < 0) {
			throw InvalidInputException("h9_curve_layer: not a curve UUID");
		}
		return layer;
	});
}

static void H9CurveBinFn(DataChunk &args, ExpressionState &, Vector &result) {
	BinaryExecutor::Execute<hugeint_t, int32_t, hugeint_t>(
	    args.data[0], args.data[1], result, args.size(), [](hugeint_t v, int32_t layer) {
		    if (layer < 0 || layer > hex9_lmax()) {
			    throw InvalidInputException("h9_curve_bin: layer must be 0..%d, got %d", hex9_lmax(), layer);
		    }
		    uint8_t u[16], out[16];
		    H9UuidToBytes(v, u);
		    if (hex9_curve_bin(u, layer, out) != 0) {
			    throw InvalidInputException(
			        "h9_curve_bin: input must be a curve UUID at layer >= %d", layer);
		    }
		    return H9BytesToUuid(out);
	    });
}

static void H9CurveIndexFn(DataChunk &args, ExpressionState &, Vector &result) {
	UnaryExecutor::Execute<hugeint_t, hugeint_t>(args.data[0], result, args.size(), [](hugeint_t v) {
		uint8_t u[16];
		H9UuidToBytes(v, u);
		char buf[48];
		int n = hex9_curve_index(u, buf, sizeof(buf));
		if (n < 0) {
			throw InvalidInputException("h9_curve_index: not a valid H9 or curve UUID");
		}
		hugeint_t out;
		if (!TryCast::Operation<string_t, hugeint_t>(string_t(buf, idx_t(n)), out, true)) {
			throw InternalException("h9_curve_index: index numeral '%s' did not fit HUGEINT", buf);
		}
		return out;
	});
}

static void H9CurvePackFn(DataChunk &args, ExpressionState &, Vector &result) {
	BinaryExecutor::Execute<hugeint_t, int32_t, hugeint_t>(
	    args.data[0], args.data[1], result, args.size(), [](hugeint_t index, int32_t layer) {
		    if (layer < 0 || layer > hex9_lmax()) {
			    throw InvalidInputException("h9_curve_pack: layer must be 0..%d, got %d", hex9_lmax(), layer);
		    }
		    const string dec = Hugeint::ToString(index);
		    uint8_t out[16];
		    if (hex9_curve_pack(dec.c_str(), layer, out) != 0) {
			    throw InvalidInputException("h9_curve_pack: index %s out of range for layer %d", dec, layer);
		    }
		    return H9BytesToUuid(out);
	    });
}

static void H9CurveLabelFn(DataChunk &args, ExpressionState &, Vector &result) {
	UnaryExecutor::Execute<hugeint_t, string_t>(args.data[0], result, args.size(), [&](hugeint_t v) {
		uint8_t u[16];
		H9UuidToBytes(v, u);
		char buf[48];
		int n = hex9_curve_label(u, buf, sizeof(buf));
		if (n < 0) {
			throw InvalidInputException("h9_curve_label: not a valid H9 or curve UUID");
		}
		return StringVector::AddString(result, buf, idx_t(n));
	});
}

static void H9CurveFromLabelFn(DataChunk &args, ExpressionState &, Vector &result) {
	UnaryExecutor::Execute<string_t, hugeint_t>(args.data[0], result, args.size(), [](string_t s) {
		uint8_t out[16];
		if (hex9_curve_parse_label(s.GetString().c_str(), out) != 0) {
			throw InvalidInputException("h9_curve_from_label: not a valid curve label: '%s'", s.GetString());
		}
		return H9BytesToUuid(out);
	});
}

void RegisterH9Curve(ExtensionLoader &loader) {
	const auto UUID = LogicalType::UUID;
	const auto I32 = LogicalType::INTEGER;

	loader.RegisterFunction(ScalarFunction("h9_curve", {UUID}, UUID, H9CurveFn));
	loader.RegisterFunction(ScalarFunction("h9_curve_decode", {UUID}, UUID, H9CurveDecodeFn));
	loader.RegisterFunction(ScalarFunction("h9_is_curve", {UUID}, LogicalType::BOOLEAN, H9IsCurveFn));
	loader.RegisterFunction(ScalarFunction("h9_curve_layer", {UUID}, I32, H9CurveLayerFn));
	loader.RegisterFunction(ScalarFunction("h9_curve_bin", {UUID, I32}, UUID, H9CurveBinFn));
	loader.RegisterFunction(ScalarFunction("h9_curve_index", {UUID}, LogicalType::HUGEINT, H9CurveIndexFn));
	loader.RegisterFunction(ScalarFunction("h9_curve_pack", {LogicalType::HUGEINT, I32}, UUID, H9CurvePackFn));
	loader.RegisterFunction(ScalarFunction("h9_curve_label", {UUID}, LogicalType::VARCHAR, H9CurveLabelFn));
	loader.RegisterFunction(ScalarFunction("h9_curve_from_label", {LogicalType::VARCHAR}, UUID, H9CurveFromLabelFn));
}

} // namespace duckdb
