/* h9_wkb.hpp — minimal little-endian WKB writer (POINT, single-ring POLYGON).
 *
 * WKB BLOBs are the cross-extension geometry interchange: the spatial
 * extension ingests them via ST_GeomFromWKB, and they are valid standalone
 * payloads (GDAL, GeoPandas, ...) when spatial is absent. Coordinates are
 * (lon, lat) degrees, WGS84, matching the core's conventions.
 */
#pragma once

#include "duckdb.hpp"

#include <cstring>

namespace duckdb {

inline void H9WkbPutU32(uint8_t *p, uint32_t v) {
	std::memcpy(p, &v, 4); /* little-endian hosts only (x86, arm64, wasm) */
}

inline void H9WkbPutD(uint8_t *p, double v) {
	std::memcpy(p, &v, 8);
}

/* WKB POINT: 1 (byte order) + 4 (type) + 16 (x, y) = 21 bytes. */
inline string_t H9WkbPoint(Vector &result, double lon, double lat) {
	uint8_t buf[21];
	buf[0] = 1;
	H9WkbPutU32(buf + 1, 1);
	H9WkbPutD(buf + 5, lon);
	H9WkbPutD(buf + 13, lat);
	return StringVector::AddStringOrBlob(result, const_char_ptr_cast(buf), sizeof(buf));
}

/* WKB POLYGON with one ring of `npts` closed (lon, lat) points:
 * 1 + 4 (type) + 4 (nrings) + 4 (npoints) + 16*npts bytes. */
inline string_t H9WkbPolygon(Vector &result, const double *lonlat, int npts) {
	const idx_t size = 13 + idx_t(npts) * 16;
	auto blob = StringVector::EmptyString(result, size);
	auto p = reinterpret_cast<uint8_t *>(blob.GetDataWriteable());
	p[0] = 1;
	H9WkbPutU32(p + 1, 3);
	H9WkbPutU32(p + 5, 1);
	H9WkbPutU32(p + 9, uint32_t(npts));
	std::memcpy(p + 13, lonlat, size_t(npts) * 16);
	blob.Finalize();
	return blob;
}

} // namespace duckdb
