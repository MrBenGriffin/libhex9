#pragma once
/*
 * h9_layout.h — addressing-layout dimensions, gated by HEX9_USE_L29.
 *
 * The 16-byte UUID is 32 nibbles. The legacy layout (HEX9_USE_L29, the default,
 * set by the CMake option of the same name) spends nibble 30 on h_term — the
 * L30 look-ahead RID — and addresses cells to L29. Undefining HEX9_USE_L29
 * reclaims nibble 30 as body[30] — proven cell-identity-redundant by
 * test/canonical_invariants.py (bins omit h_term yet stay globally unique and
 * traversable) — extending addressing to L30 for the same 16 bytes. key_tail
 * (nibble 31) is fixed in both layouts.
 *
 * SCOPE: this header + the H9_LMAX/H9_DESC_MAX cap refactor parameterise the
 * dimensions only, and are behaviour-neutral on the legacy layout (a
 * USE_L29=ON build must still reproduce test/golden_l29.py). The per-nibble
 * encode/decode/bin logic for the L30 branch (h_term -> body[30], the
 * canonical-terminal decode seed, the key_tail shift) is a separate step and
 * is NOT implemented here.
 */
#if defined(HEX9_USE_L29)
#  define H9_LMAX 29              /* deepest addressable layer */
#else
#  define H9_LMAX 30
#endif

#define H9_NLEVELS     (H9_LMAX + 1)   /* body levels including L0 */
#define H9_NIB_TAIL    31              /* key_tail nibble — fixed in both layouts */
#define H9_NIB_BODYTOP (H9_LMAX)       /* deepest body nibble */
#define H9_NIB_HTERM   (H9_LMAX + 1)   /* h_term nibble; == NIB_TAIL (absent) when L30 */
#define H9_HAS_HTERM   (H9_NIB_HTERM < H9_NIB_TAIL)
#define H9_DESC_MAX    (H9_LMAX + 7)   /* deepest cid/rid index: body + 6-level look-ahead */
