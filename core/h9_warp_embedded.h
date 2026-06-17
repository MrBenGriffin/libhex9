/* h9_warp_embedded.h — declarations for the .incbin-embedded warp blob.
 *
 * The bytes are defined in h9_warp_embedded.cpp via a compile-time .incbin of
 * the blob named by H9_WARP_BLOB (since F6: WGS84_l5_warp_f6.full.f64g.h9warp,
 * 12.8 MB — v3 format with per-vertex deltas AND hhg9's final CT gradients).
 * Baked into the library's data section — no runtime file load, fully
 * self-contained. (Replaces the former ~13 MB auto-generated ASCII header.)
 */
#pragma once
#include <cstddef>

namespace h9 {
extern const unsigned char *const EMBEDDED_WARP_DATA;  /* the warp blob */
extern const std::size_t          EMBEDDED_WARP_SIZE;  /* byte count */
}
