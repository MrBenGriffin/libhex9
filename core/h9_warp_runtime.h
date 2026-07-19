/* h9_warp_runtime.h — process-global WarpState + thin h9_warp_fwd/inv
 * delegators that match the existing legacy API.
 *
 * Lifecycle:
 *   - The state is a single h9::WarpState owned by one TU (see
 *     h9_warp_runtime.cpp). Declared `extern` here so any header that
 *     wants to use it via the thin h9_warp_fwd/inv wrappers can.
 *   - Caller responsibility: call h9_warp_init_embedded() exactly once,
 *     early, before any h9_warp_fwd/inv invocation. For the PostGIS
 *     extension that means _PG_init(); for the CLI it means main().
 *   - h9_warp_init_embedded() reads the .h9warp blob baked into
 *     h9_warp_embedded.h and runs the full WarpState build (~0.9 s).
 *     (Was ~13 s until the 729² grid fill switched from a per-cell
 *     bin-spiral to a multi-source BFS flood — see h9_ct.h build_ct_state.)
 *
 * The h9_warp_fwd / h9_warp_inv shims preserve the legacy call-site
 * signature so h9_addressing.h, main.cpp, lwgeom_hex9.cpp need no
 * source edits. They are only compiled when H9_WARP_ENABLE is set;
 * otherwise the warp is a no-op identity.
 */
#ifndef H9_WARP_RUNTIME_H
#define H9_WARP_RUNTIME_H

#include "h9_warp.h"

#include <string>

namespace h9 {

/* Global, lazily-but-explicitly initialised. Owned by the TU that
 * defines it in h9_warp_runtime.cpp. */
extern WarpState g_warp_state;
extern bool      g_warp_state_ready;

/* Via-sphere counterpart: the Sphere-L6 wedge-fold field (v4 fund blob),
 * selected by the shims when the via-sphere mode is on. Built lazily by
 * hex9_set_via_sphere / h9_warp_init_embedded_sph. */
extern WarpState g_warp_state_sph;
extern bool      g_warp_state_sph_ready;

/* Runtime apply toggle. true (default) ⇒ h9_warp_fwd/inv apply the
 * authalic correction when state is ready; false ⇒ identity. Exposed
 * in PostGIS as the `hex9.use_warp` GUC; CLI / standalone builds get
 * the default value and need not touch it. Decoupled from state-ready
 * so the GUC can pre-disable warp without requiring the 13 s init to
 * be skipped (or vice-versa). */
extern bool      g_warp_use;

/* Build g_warp_state from the embedded sidecar blob baked into
 * h9_warp_embedded.h. Idempotent (no-op after first success). Returns
 * true on success; on failure `err` is filled and g_warp_state_ready
 * stays false. */
bool h9_warp_init_embedded(std::string* err = nullptr,
                           int    grad_maxiter = 2000,
                           double grad_tol     = 1e-12);

/* Build g_warp_state_sph from the embedded Sphere-L6 fund blob.
 * Idempotent; ~1 s (no gradient estimation — the blob ships them). */
bool h9_warp_init_embedded_sph(std::string* err = nullptr);

/* Optional alternate init from a runtime sidecar file (e.g. for tests
 * or operator-supplied warpfiles). */
bool h9_warp_init_from_path(const std::string& path,
                            std::string* err = nullptr,
                            int    grad_maxiter = 2000,
                            double grad_tol     = 1e-12);

/* Lateral-edge treatment (see WarpState in h9_warp.h):
 *   Feather — product default: delta ramps smoothly to exact identity
 *             on the seam (F6 fix).
 *   Bypass  — legacy hard identity band (hhg9 parity runs).
 *   Raw     — no edge treatment at all: the trained field speaks for
 *             itself. REQUIRED when validating new/retrained .h9warp
 *             data (Ben's carve-out, 2026-06-11) — feather or bypass
 *             would mask the field's own edge behaviour and give
 *             misleading results.
 * Both init functions honour the H9_WARP_EDGE environment variable
 * ("feather" | "bypass" | "raw") so tools and test runs can select the
 * mode without recompiling. */
enum class WarpEdgeMode { Feather, Bypass, Raw };
void         h9_warp_set_edge_mode(WarpEdgeMode m);
WarpEdgeMode h9_warp_edge_mode();

} /* namespace h9 */

/* ── Legacy-API delegators, header-only -─────────────────────────────
 * Replace the static-data versions formerly in h9_math.h:281-413. The
 * shims pass through to identity when (a) the state isn't built (CLI
 * forgot to call h9_warp_init_embedded), or (b) the operator has
 * disabled warp at runtime via `SET hex9.use_warp = off`. */
#if H9_WARP_ENABLE
/* Via-sphere: when included via h9_math.h (which defines the per-TU mode
 * flag) the shims route to the Sphere-L6 wedge-fold state; standalone
 * includers (h9_warp_runtime.cpp) compile the classic-only fallback. */
#ifdef H9_VIA_SPHERE_FLAG
static inline const h9::WarpState* h9_warp_pick(bool* ready)
{
    if (h9_g_via_sphere) {
        *ready = h9::g_warp_state_sph_ready;
        return &h9::g_warp_state_sph;
    }
    *ready = h9::g_warp_state_ready;
    return &h9::g_warp_state;
}
#else
static inline const h9::WarpState* h9_warp_pick(bool* ready)
{
    *ready = h9::g_warp_state_ready;
    return &h9::g_warp_state;
}
#endif

static inline void h9_warp_fwd(double rx, double ry, int oct_mode,
                               double *wx, double *wy)
{
    bool ready;
    const h9::WarpState* ws = h9_warp_pick(&ready);
    if (!h9::g_warp_use || !ready) { *wx = rx; *wy = ry; return; }
    h9::warp_do(*ws, rx, ry, oct_mode, wx, wy);
}

/* Hint overload retained for API compatibility — hints are unused now
 * that the inverse seed is identity-based. */
static inline void h9_warp_fwd(double rx, double ry, int oct_mode,
                               double *wx, double *wy,
                               double /*hint_x*/, double /*hint_y*/)
{
    h9_warp_fwd(rx, ry, oct_mode, wx, wy);
}

static inline void h9_warp_inv(double wx, double wy, int oct_mode,
                               double *rx, double *ry)
{
    bool ready;
    const h9::WarpState* ws = h9_warp_pick(&ready);
    if (!h9::g_warp_use || !ready) { *rx = wx; *ry = wy; return; }
    h9::warp_undo(*ws, wx, wy, oct_mode, rx, ry);
}
#endif /* H9_WARP_ENABLE */

#endif /* H9_WARP_RUNTIME_H */
