/* h9_uuid.hpp — hex9 uuid <-> DuckDB UUID/UHUGEINT conversions + doctrine guard.
 *
 * hex9 ids are raw 16-byte big-endian buffers. DuckDB's UUID is a hugeint_t
 * whose upper word has bit 63 flipped so that signed hugeint order equals
 * bytewise uuid order (see duckdb src/common/types/uuid.cpp). UHUGEINT is the
 * unflipped raw 128-bit value, so uhugeint order == uuid order too.
 */
#pragma once

#include "duckdb.hpp"
#include "hex9_c.h"

namespace duckdb {

inline uint64_t H9Be64Get(const uint8_t *p) {
	uint64_t v = 0;
	for (int i = 0; i < 8; i++) {
		v = (v << 8) | p[i];
	}
	return v;
}

inline void H9Be64Put(uint64_t v, uint8_t *p) {
	for (int i = 7; i >= 0; i--) {
		p[i] = uint8_t(v & 0xFF);
		v >>= 8;
	}
}

inline hugeint_t H9BytesToUuid(const uint8_t u[16]) {
	hugeint_t r;
	r.upper = int64_t(H9Be64Get(u) ^ (uint64_t(1) << 63));
	r.lower = H9Be64Get(u + 8);
	return r;
}

inline void H9UuidToBytes(hugeint_t v, uint8_t out[16]) {
	H9Be64Put(uint64_t(v.upper) ^ (uint64_t(1) << 63), out);
	H9Be64Put(v.lower, out + 8);
}

inline uhugeint_t H9BytesToUhuge(const uint8_t u[16]) {
	uhugeint_t r;
	r.upper = H9Be64Get(u);
	r.lower = H9Be64Get(u + 8);
	return r;
}

inline void H9UhugeToBytes(uhugeint_t v, uint8_t out[16]) {
	H9Be64Put(v.upper, out);
	H9Be64Put(v.lower, out + 8);
}

/* Adjacency input must be a full UUID from h9_encode — bins are layer-scoped
 * keys, not addresses (their key tail cannot carry the meta the resolution
 * needs); give the doctrine error rather than a generic one. Mirrors
 * postgis_hex9's h9_reject_bin_input. */
inline bool H9IsBin(const uint8_t u[16]) {
	return (u[15] >> 4) == 0x0Fu;
}

inline void H9RejectBinInput(const uint8_t u[16], const char *fname) {
	if (H9IsBin(u)) {
		throw InvalidInputException(
		    "%s: bin UUIDs are layer keys, not addresses — pass the full UUID from h9_encode. "
		    "Re-derive adjacency from the original point; output bins are for joining, "
		    "not for further traversal.",
		    fname);
	}
}

} // namespace duckdb
