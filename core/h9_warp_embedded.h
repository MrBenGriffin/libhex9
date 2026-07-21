/* h9_warp_embedded.h — declarations for the .incbin-embedded warp blob.
 *
 * The bytes are defined in h9_warp_embedded.cpp via a compile-time .incbin of
 * the blob named by H9_WARP_BLOB: since 2.0.0 that is the Sphere-L6
 * fundamental-domain v4 wedge field (Sphere_l6_fund.f64g.h9warp, ~19 MB —
 * per-vertex deltas AND shipped gradients, so no gradient-estimation pass at
 * init). Baked into the library's data section — no runtime file load, fully
 * self-contained.
 */
#pragma once
#include <cstddef>

namespace h9 {
extern const unsigned char *const EMBEDDED_WARP_DATA;  /* the warp blob */
extern const std::size_t          EMBEDDED_WARP_SIZE;  /* byte count */
}
