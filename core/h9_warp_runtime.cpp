/* h9_warp_runtime.cpp — single TU that owns the global WarpState and
 * the 2 MB embedded sidecar blob. Including h9_warp_embedded.h here
 * (and only here) keeps the 13 MB ASCII expansion out of other TUs'
 * preprocessor work.
 */
#include "h9_warp_runtime.h"
#include "h9_warp_embedded.h"

#include <cstdlib>
#include <cstring>

namespace h9 {

WarpState g_warp_state;
bool      g_warp_state_ready = false;
bool      g_warp_use         = true;  /* default ON; PostGIS GUC may override */

void h9_warp_set_edge_mode(WarpEdgeMode m)
{
    g_warp_state.edge_feather = (m == WarpEdgeMode::Feather);
    g_warp_state.edge_bypass  = (m == WarpEdgeMode::Bypass);
}

WarpEdgeMode h9_warp_edge_mode()
{
    if (g_warp_state.edge_feather) return WarpEdgeMode::Feather;
    if (g_warp_state.edge_bypass)  return WarpEdgeMode::Bypass;
    return WarpEdgeMode::Raw;
}

/* Default edge mode follows the field, then the env may override, then
 * the combination is validated: the F6 edge-tangent field (v3, gradients
 * shipped) is boundary-preserving by construction and MUST run raw —
 * feather would compress its legitimate km-scale tangential seam slide
 * into the 350 m ramp, bypass would reintroduce the sliver bug. Legacy
 * (v1/v2) fields default to feather as before. */
static bool apply_edge_mode_env(std::string* err)
{
    h9_warp_set_edge_mode(g_warp_state.field_has_grads ? WarpEdgeMode::Raw
                                                       : WarpEdgeMode::Feather);
    const char* e = std::getenv("H9_WARP_EDGE");
    if (e) {
        if      (!std::strcmp(e, "feather")) h9_warp_set_edge_mode(WarpEdgeMode::Feather);
        else if (!std::strcmp(e, "bypass"))  h9_warp_set_edge_mode(WarpEdgeMode::Bypass);
        else if (!std::strcmp(e, "raw"))     h9_warp_set_edge_mode(WarpEdgeMode::Raw);
        else {
            if (err) *err = std::string("H9_WARP_EDGE must be feather|bypass|raw, got ") + e;
            return false;
        }
    }
    if (g_warp_state.field_has_grads &&
        h9_warp_edge_mode() != WarpEdgeMode::Raw) {
        if (err) *err = "H9_WARP_EDGE=feather/bypass is invalid for the F6 "
                        "edge-tangent warp field (v3, gradients shipped); "
                        "this field must run raw";
        return false;
    }
    return true;
}

bool h9_warp_init_embedded(std::string* err,
                           int    grad_maxiter,
                           double grad_tol)
{
    if (g_warp_state_ready) return true;
    H9WarpData data;
    if (!load_h9warp(EMBEDDED_WARP_DATA, EMBEDDED_WARP_SIZE, data, err))
        return false;
    g_warp_state.data = std::move(data);
    if (!finish_warp_state(g_warp_state, grad_maxiter, grad_tol)) {
        if (err) *err = "h9warp: gradient count mismatch vs mesh";
        return false;
    }
    if (!apply_edge_mode_env(err)) return false;
    g_warp_state_ready = true;
    return true;
}

bool h9_warp_init_from_path(const std::string& path,
                            std::string*       err,
                            int                grad_maxiter,
                            double             grad_tol)
{
    if (g_warp_state_ready) return true;
    if (!build_warp_state(path, g_warp_state, err, grad_maxiter, grad_tol))
        return false;
    if (!apply_edge_mode_env(err)) return false;
    g_warp_state_ready = true;
    return true;
}

} /* namespace h9 */
