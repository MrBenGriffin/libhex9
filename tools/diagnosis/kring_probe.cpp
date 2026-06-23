/* kring_probe.cpp — empirical derivation/validation of the symbolic k-ring
 * neighbour algebra, against the geometric truth of h9grid::enumerate_global.
 *
 * Validates, per layer:
 *   1. SOUNDNESS  — every enumerated cell's origin descent re-converges
 *                   (tree membership == exact integer-descent convergence).
 *   2. ADJACENCY  — geometric edge-adjacency (2 shared ring vertices) gives
 *                   degree 6 everywhere (or reports the defect histogram).
 *   3. STEPS      — same-frame neighbour centroid deltas are exactly
 *                   {(0,±2),(±3,±1)} in integer UV.
 *   4. SEAMS      — every cross-frame adjacent pair is explained by either
 *                   the v-mirror map (u,v)→(u,−v) (edge neighbours) or the
 *                   identity map (corner/vertex neighbours), and the mirror
 *                   target oid agrees with H9_OID_NB.
 *   5. UNIQUENESS — brute-force resolution (all 8 oids × {identity, mirror})
 *                   of every neighbour centroid converges in EXACTLY one
 *                   frame, and the resolved identity matches the geometric
 *                   neighbour.
 *
 * Build:  c++ -O2 -std=c++17 -I core -o build/kring_probe tools/kring_probe.cpp
 * Run:    build/kring_probe [max_layer]
 */
#define H9_WARP_ENABLE 0
#include "h9_math.h"
#include "h9_addressing.h"
#include "h9_uv_lattice.h"
#include "h9_grid.h"
#include "h9_kring.h"

#include <cstdio>
#include <cstdint>
#include <map>
#include <set>
#include <array>
#include <vector>
#include <algorithm>

/* c2 hex centroid offset from supercell origin (template vertex 6). */
static const int64_t C2_DU[3] = { 1, 1, -2 };
static const int64_t C2_DV[3] = { 1, -1, 0 };

/* The 6 hex neighbour steps on the centroid lattice (integer UV). */
static const int64_t STEP[6][2] = {
    {0, 2}, {3, 1}, {3, -1}, {0, -2}, {-3, -1}, {-3, 1}
};

/* ── Derived seam algebra (fitted + verified by this probe at L1..L3) ───────
 * Crossing edge e (0..2) of octant `oid` into H9_OID_NB[oid][e]:
 *   p' = R[e]·p + s·T[mo][e]      (s = 3^layer, mo = H9_OID_MO[oid])
 * R[0] = identity; R[1], R[2] = rotations by ∓120° (numerators over 2).
 * All maps are det=+1 lattice rotations; T depends only on (mo, e). */
static const int SEAM_R[3][4] = {           /* {n00,n01,n10,n11} over /2 */
    { 2,  0, 0,  2 },
    {-1,  3, -1, -1 },
    {-1, -3,  1, -1 },
};
static const int SEAM_T[2][3][2] = {
    { {0, -2}, {3, 1},  {-3, 1} },          /* source octant mode 0 */
    { {0,  2}, {-3, -1}, {3, -1} },         /* source octant mode 1 */
};

static void seam_apply(int e, int mo, int64_t s, int64_t u, int64_t v,
                       int64_t *ou, int64_t *ov) {
    const int64_t nu = SEAM_R[e][0]*u + SEAM_R[e][1]*v;
    const int64_t nv = SEAM_R[e][2]*u + SEAM_R[e][3]*v;
    *ou = nu/2 + s*SEAM_T[mo][e][0];
    *ov = nv/2 + s*SEAM_T[mo][e][1];
}

/* Exact integer descent from the octant root towards supercell origin
 * (ia, ib) at `layer`.  Mirrors cids_from_iauv's walk (classification +
 * NN fallback) but only tracks convergence.  Returns true iff the walk
 * lands exactly on (ia, ib) at a mode-0 leaf — i.e. (oid, ia, ib) is a
 * real supercell of oid's BFS tree. */
static bool descend_converges(int oid, int64_t ia, int64_t ib, int layer) {
    int64_t scale = 1;
    for (int i = 0; i < layer; i++) scale *= 3;
    int     p_mo = (int)H9_OID_MO[oid];
    int64_t ou = 0, ov = 0;
    for (int k = 0; k < layer; k++) {
        const int (*ofs)[2] = (p_mo == 1) ? H9_MODE1_OFS        : H9_MODE0_OFS;
        const int *cm       = (p_mo == 1) ? H9_MODE1_CHILD_MODE : H9_MODE0_CHILD_MODE;
        const uint8_t *ctab = (p_mo == 1) ? H9_UP_CIDS          : H9_DN_CIDS;
        const int64_t resu = ia - ou, resv = ib - ov;
        const uint8_t cid  = h9grid::classify_band_int(resu, resv, scale, p_mo);
        const int64_t cs   = scale / 3;
        int j = -1;
        for (int kk = 0; kk < 9; kk++) if (ctab[kk] == cid) { j = kk; break; }
        if (j < 0) {
            int64_t best = INT64_MAX; int bj = 0;
            for (int kk = 0; kk < 9; kk++) {
                const int64_t du = resu - (int64_t)ofs[kk][0] * cs;
                const int64_t dv = resv - (int64_t)ofs[kk][1] * cs;
                const int64_t d2 = du*du + 3*dv*dv;
                if (d2 < best) { best = d2; bj = kk; }
            }
            j = bj;
        }
        ou += (int64_t)ofs[j][0] * cs;
        ov += (int64_t)ofs[j][1] * cs;
        p_mo = cm[j];
        scale = cs;
    }
    return ou == ia && ov == ib && p_mo == 0;
}

/* Resolve a centroid (cu, cv) in `oid`'s frame to a cell identity, or fail. */
static bool resolve_centroid(int oid, int64_t cu, int64_t cv, int layer,
                             int *c2_out, int64_t *ia_out, int64_t *ib_out) {
    int hits = 0;
    for (int c2 = 0; c2 < 3; c2++) {
        const int64_t ia = cu - C2_DU[c2], ib = cv - C2_DV[c2];
        if (descend_converges(oid, ia, ib, layer)) {
            *c2_out = c2; *ia_out = ia; *ib_out = ib;
            hits++;
        }
    }
    if (hits > 1) std::printf("  !! resolve ambiguity: oid=%d cu=%lld cv=%lld (%d hits)\n",
                              oid, (long long)cu, (long long)cv, hits);
    return hits == 1;
}

struct VKey { int64_t x, y, z; bool operator<(const VKey &o) const {
    if (x != o.x) return x < o.x;
    if (y != o.y) return y < o.y;
    return z < o.z; } };

static VKey vkey(double lon_deg, double lat_deg, double shift) {
    double x, y, z;
    h9_rad_lonlat_to_ecef(lon_deg * M_PI / 180.0, lat_deg * M_PI / 180.0, &x, &y, &z);
    const double r = std::sqrt(x*x + y*y + z*z);   /* metres → unit sphere */
    const double Q = 1e9;                          /* tolerance 1e-9 rel */
    return { (int64_t)std::llround(x / r * Q + shift),
             (int64_t)std::llround(y / r * Q + shift),
             (int64_t)std::llround(z / r * Q + shift) };
}

int main(int argc, char **argv) {
    const int max_layer = (argc > 1) ? std::atoi(argv[1]) : 3;

    for (int L = 1; L <= max_layer; L++) {
        std::vector<H9GridCell> cells;
        h9grid::enumerate_global(L, cells);
        const int n = (int)cells.size();
        std::printf("=== layer %d: %d cells (expect %d) ===\n", L, n, 12 * (int)std::pow(9.0, L));

        /* 1. soundness: every cell's own origin re-converges in its frame */
        int bad_conv = 0;
        for (int i = 0; i < n; i++) {
            if (!descend_converges(cells[i].oid, cells[i].ia, cells[i].ib, L)) {
                if (bad_conv < 5)
                    std::printf("  !! no convergence: oid=%d c2=%d ia=%lld ib=%lld ext=%d\n",
                                cells[i].oid, cells[i].c2,
                                (long long)cells[i].ia, (long long)cells[i].ib, cells[i].ext);
                bad_conv++;
            }
        }
        std::printf("  soundness: %d/%d origins re-converge\n", n - bad_conv, n);

        /* 2. geometric adjacency via shared ring vertices (two offset grids) */
        std::map<std::pair<int,int>, int> shared;
        for (int g = 0; g < 2; g++) {
            std::map<VKey, std::vector<int>> vmap;
            for (int i = 0; i < n; i++)
                for (int v = 0; v < 6; v++)
                    vmap[vkey(cells[i].vlon[v], cells[i].vlat[v], g * 0.5)].push_back(i);
            std::map<std::pair<int,int>, int> sh_g;
            for (auto &kv : vmap) {
                auto &lst = kv.second;
                std::sort(lst.begin(), lst.end());
                lst.erase(std::unique(lst.begin(), lst.end()), lst.end());
                for (size_t a = 0; a < lst.size(); a++)
                    for (size_t b = a + 1; b < lst.size(); b++)
                        sh_g[{lst[a], lst[b]}]++;
            }
            for (auto &kv : sh_g) {
                auto it = shared.find(kv.first);
                if (it == shared.end()) shared[kv.first] = kv.second;
                else it->second = std::max(it->second, kv.second);
            }
        }
        std::vector<std::vector<int>> adj(n);
        for (auto &kv : shared)
            if (kv.second >= 2) {
                adj[kv.first.first].push_back(kv.first.second);
                adj[kv.first.second].push_back(kv.first.first);
            }
        std::map<int,int> deg_hist;
        for (int i = 0; i < n; i++) deg_hist[(int)adj[i].size()]++;
        std::printf("  degree histogram:");
        for (auto &kv : deg_hist) std::printf("  %d:%d", kv.first, kv.second);
        std::printf("\n");

        /* 2b. fit the exact affine seam map per ordered (oidA → oidB) pair from
         * shared-vertex correspondences.  Model: vB = M·vA + s·t  with integer
         * M (2×2) and integer t (×s).  Verified exactly on every sample. */
        {
            const int64_t s_int = (int64_t)std::llround(std::pow(3.0, (double)L));
            /* gather correspondences keyed by geographic vertex */
            std::map<VKey, std::vector<std::array<int64_t,3>>> vown;  /* (oid,u,v) */
            for (int i = 0; i < n; i++)
                for (int v = 0; v < 6; v++)
                    vown[vkey(cells[i].vlon[v], cells[i].vlat[v], 0.0)]
                        .push_back({(int64_t)cells[i].poid[v], cells[i].pu[v], cells[i].pv[v]});
            std::map<std::pair<int,int>, std::vector<std::array<int64_t,4>>> corr;
            for (auto &kv : vown) {
                auto &lst = kv.second;
                std::sort(lst.begin(), lst.end());
                lst.erase(std::unique(lst.begin(), lst.end()), lst.end());
                for (auto &a : lst)
                    for (auto &b : lst)
                        if (a[0] != b[0])
                            corr[{(int)a[0], (int)b[0]}].push_back({a[1], a[2], b[1], b[2]});
            }
            /* The 12 lattice isometries: R60^k (k=0..5) optionally composed with
             * the v-mirror.  In integer (u,v) coords (cartesian (u·U1, v·√3·U1)):
             * R60 = [(u−3v)/2, (u+v)/2] — integer on the u+v-even sublattice.
             * Stored as 2×M numerators over denominator 2. */
            struct Iso { int n00, n01, n10, n11; };   /* (n·p)/2 */
            std::vector<Iso> isos;
            {
                Iso r = {2, 0, 0, 2};                /* identity ×2 */
                for (int k = 0; k < 6; k++) {
                    isos.push_back(r);
                    isos.push_back({r.n00, -r.n01, r.n10, -r.n11});   /* ∘ v-mirror */
                    /* r ← R60·r  (R60 numerators over 2: [[1,-3],[1,1]]) */
                    Iso nr = { (1*r.n00 + -3*r.n10) / 2, (1*r.n01 + -3*r.n11) / 2,
                               (1*r.n00 +  1*r.n10) / 2, (1*r.n01 +  1*r.n11) / 2 };
                    r = nr;
                }
            }
            /* collect adjacent cross pairs per ordered oid pair (for disambiguation) */
            std::map<std::pair<int,int>, std::vector<std::pair<int,int>>> xadj;
            for (int i = 0; i < n; i++)
                for (int j : adj[i])
                    if (cells[i].oid != cells[j].oid)
                        xadj[{cells[i].oid, cells[j].oid}].push_back({i, j});

            std::printf("  seam maps (vB = (N.vA)/2 + s*t):\n");
            for (auto &kv : corr) {
                auto &ps = kv.second;
                auto &aps = xadj[kv.first];
                int nfound = 0;
                for (auto &iso : isos) {
                    /* fit t from first corr; require integer t (×s) */
                    const int64_t fu = iso.n00*ps[0][0] + iso.n01*ps[0][1];
                    const int64_t fv = iso.n10*ps[0][0] + iso.n11*ps[0][1];
                    if (fu % 2 || fv % 2) continue;
                    const int64_t du = ps[0][2] - fu/2, dv = ps[0][3] - fv/2;
                    if (du % s_int || dv % s_int) continue;
                    const int64_t t0 = du / s_int, t1 = dv / s_int;
                    /* verify on all vertex correspondences */
                    bool ok = true;
                    for (auto &p : ps) {
                        const int64_t gu = iso.n00*p[0] + iso.n01*p[1];
                        const int64_t gv = iso.n10*p[0] + iso.n11*p[1];
                        if (gu % 2 || gv % 2 ||
                            gu/2 + t0*s_int != p[2] || gv/2 + t1*s_int != p[3]) { ok = false; break; }
                    }
                    if (!ok) continue;
                    /* disambiguate: for every adjacent pair, mapping the stepped
                     * target centroid (centA + step, or centA itself for half-hex
                     * partners) must land EXACTLY on centB for some step. */
                    int badp = 0;
                    for (auto &ap : aps) {
                        const H9GridCell &A = cells[ap.first], &B = cells[ap.second];
                        const int64_t cuA = A.ia + C2_DU[A.c2], cvA = A.ib + C2_DV[A.c2];
                        const int64_t cuB = B.ia + C2_DU[B.c2], cvB = B.ib + C2_DV[B.c2];
                        bool step_ok = false;
                        const int s_first = (A.ext && B.ext) ? -1 : 0;  /* s=0 only for half-hex partners */
                        for (int s2 = s_first; s2 < 6 && !step_ok; s2++) {
                            const int64_t tu = cuA + (s2 < 0 ? 0 : STEP[s2][0]);
                            const int64_t tv = cvA + (s2 < 0 ? 0 : STEP[s2][1]);
                            const int64_t gu = iso.n00*tu + iso.n01*tv;
                            const int64_t gv = iso.n10*tu + iso.n11*tv;
                            if (gu % 2 || gv % 2) continue;
                            if (gu/2 + t0*s_int == cuB && gv/2 + t1*s_int == cvB)
                                step_ok = true;
                        }
                        if (!step_ok) badp++;
                    }
                    if (badp == 0) {
                        std::printf("    o%d->o%d  N=[%d %d; %d %d] t=[%lld %lld]  "
                                    "(%zu corr, %zu pairs)\n",
                                    kv.first.first, kv.first.second,
                                    iso.n00, iso.n01, iso.n10, iso.n11,
                                    (long long)t0, (long long)t1, ps.size(), aps.size());
                        nfound++;
                    }
                }
                if (nfound != 1)
                    std::printf("    o%d->o%d  !! %d candidate maps (corr=%zu pairs=%zu)\n",
                                kv.first.first, kv.first.second, nfound, ps.size(), aps.size());
            }
        }

        /* 3+4. classify each directed adjacent pair */
        std::set<std::pair<int64_t,int64_t>> same_frame_deltas;
        std::map<std::array<int,3>, int> cross_rel;   /* {oidA, oidB, type} type:0=mirror,1=ident */
        int unexplained = 0, nb_table_viol = 0;
        for (int i = 0; i < n; i++) {
            const int64_t cuA = cells[i].ia + C2_DU[cells[i].c2];
            const int64_t cvA = cells[i].ib + C2_DV[cells[i].c2];
            for (int j : adj[i]) {
                const int64_t cuB = cells[j].ia + C2_DU[cells[j].c2];
                const int64_t cvB = cells[j].ib + C2_DV[cells[j].c2];
                if (cells[i].oid == cells[j].oid) {
                    same_frame_deltas.insert({cuB - cuA, cvB - cvA});
                    continue;
                }
                bool expl = false;
                for (int s = 0; s < 6 && !expl; s++) {
                    const int64_t tu = cuA + STEP[s][0], tv = cvA + STEP[s][1];
                    if (tu == cuB && -tv == cvB) {            /* v-mirror */
                        cross_rel[{cells[i].oid, cells[j].oid, 0}]++;
                        bool in_nb = false;
                        for (int e = 0; e < 3; e++)
                            if (H9_OID_NB[cells[i].oid][e] == cells[j].oid) in_nb = true;
                        if (!in_nb) nb_table_viol++;
                        expl = true;
                    } else if (tu == cuB && tv == cvB) {      /* identity */
                        cross_rel[{cells[i].oid, cells[j].oid, 1}]++;
                        expl = true;
                    }
                }
                if (!expl) {
                    if (unexplained < 12)
                        std::printf("  ?? unexplained: A(oid=%d c2=%d ia=%lld ib=%lld ext=%d) "
                                    "B(oid=%d c2=%d ia=%lld ib=%lld ext=%d)\n",
                                    cells[i].oid, cells[i].c2, (long long)cells[i].ia,
                                    (long long)cells[i].ib, cells[i].ext,
                                    cells[j].oid, cells[j].c2, (long long)cells[j].ia,
                                    (long long)cells[j].ib, cells[j].ext);
                    unexplained++;
                }
            }
        }
        std::printf("  same-frame deltas (%zu):", same_frame_deltas.size());
        for (auto &d : same_frame_deltas)
            std::printf(" (%lld,%lld)", (long long)d.first, (long long)d.second);
        std::printf("\n  cross-frame relations:\n");
        for (auto &kv : cross_rel)
            std::printf("    oid %d -> oid %d  %s  x%d\n", kv.first[0], kv.first[1],
                        kv.first[2] ? "IDENT " : "MIRROR", kv.second);
        std::printf("  unexplained cross pairs: %d   NB-table violations: %d\n",
                    unexplained, nb_table_viol);

        /* 5. brute-force resolution of every neighbour step; compare to truth */
        std::map<std::array<int64_t,4>, int> by_id;   /* (oid,c2,ia,ib) -> index */
        for (int i = 0; i < n; i++)
            by_id[{(int64_t)cells[i].oid, (int64_t)cells[i].c2, cells[i].ia, cells[i].ib}] = i;

        long steps_total = 0, steps_multi = 0, steps_none = 0, steps_wrong = 0;
        std::map<int,int> resolved_deg_hist;
        for (int i = 0; i < n; i++) {
            const int64_t cuA = cells[i].ia + C2_DU[cells[i].c2];
            const int64_t cvA = cells[i].ib + C2_DV[cells[i].c2];
            std::set<int> resolved_nbrs;
            const int64_t s_int2 = (int64_t)std::llround(std::pow(3.0, (double)L));
            for (int s = -1; s < 6; s++) {
                /* s == -1 probes the zero step (half-hex partner search). */
                const int64_t tu = cuA + (s < 0 ? 0 : STEP[s][0]);
                const int64_t tv = cvA + (s < 0 ? 0 : STEP[s][1]);
                steps_total++;
                /* candidate frames: own; 3 edge maps; 6 corner compositions */
                struct Att { int oid; int64_t u, v; };
                Att att[10];
                int natt = 0;
                const int oidA = cells[i].oid, moA = (int)H9_OID_MO[oidA];
                if (s >= 0) att[natt++] = { oidA, tu, tv };
                for (int e = 0; e < 3; e++) {
                    const int nb = H9_OID_NB[oidA][e];
                    const int moB = (int)H9_OID_MO[nb];
                    int64_t u1, v1;
                    seam_apply(e, moA, s_int2, tu, tv, &u1, &v1);
                    att[natt++] = { nb, u1, v1 };
                    for (int e2 = 0; e2 < 3; e2++) {
                        const int nc = H9_OID_NB[nb][e2];
                        if (nc == oidA) continue;
                        int64_t u2, v2;
                        seam_apply(e2, moB, s_int2, u1, v1, &u2, &v2);
                        att[natt++] = { nc, u2, v2 };
                    }
                }
                int nhits = 0, hit_idx = -1;
                std::set<int> hit_set;
                for (int a = 0; a < natt; a++) {
                    int c2r; int64_t iar, ibr;
                    if (resolve_centroid(att[a].oid, att[a].u, att[a].v, L, &c2r, &iar, &ibr)) {
                        auto it = by_id.find({(int64_t)att[a].oid, (int64_t)c2r, iar, ibr});
                        const int idx = (it == by_id.end()) ? -2 : it->second;
                        if (idx == i) continue;        /* self — not a neighbour */
                        if (hit_set.insert(idx).second) nhits++;
                        hit_idx = idx;
                    }
                }
                if (s < 0 && nhits == 0) { steps_total--; continue; }  /* no partner: fine for non-ext */
                if (nhits == 0) steps_none++;
                else if (nhits > 1) {
                    steps_multi++;
                    static int shown_multi = 0;
                    if (shown_multi++ < 8) {
                        std::printf("  ** ambiguous step: cell(oid=%d c2=%d ia=%lld ib=%lld) "
                                    "step=%d hits:", cells[i].oid, cells[i].c2,
                                    (long long)cells[i].ia, (long long)cells[i].ib, s);
                        for (int h : hit_set) std::printf(" %d", h);
                        std::printf("\n");
                    }
                }
                else if (hit_idx == -2) steps_wrong++;       /* resolved id not in mesh */
                else resolved_nbrs.insert(hit_idx);
            }
            resolved_deg_hist[(int)resolved_nbrs.size()]++;
            /* compare to geometric adjacency */
            std::set<int> geo(adj[i].begin(), adj[i].end());
            if (geo != resolved_nbrs && steps_wrong < 12 && i < 100000) {
                static int shown = 0;
                if (shown++ < 8) {
                    std::printf("  ## mismatch at cell %d (oid=%d c2=%d ia=%lld ib=%lld ext=%d): "
                                "geo=%zu resolved=%zu\n",
                                i, cells[i].oid, cells[i].c2, (long long)cells[i].ia,
                                (long long)cells[i].ib, cells[i].ext,
                                geo.size(), resolved_nbrs.size());
                }
            }
        }
        std::printf("  resolution: %ld steps  none=%ld multi=%ld not-in-mesh=%ld\n",
                    steps_total, steps_none, steps_multi, steps_wrong);
        std::printf("  resolved-degree histogram:");
        for (auto &kv : resolved_deg_hist) std::printf("  %d:%d", kv.first, kv.second);
        std::printf("\n");

        /* 6. validate the production header (core/h9_kring.h) end-to-end:
         * uuid -> identity -> neighbors == geometric adjacency, every cell. */
        {
            long id_bad = 0, nb_bad = 0;
            for (int i = 0; i < n; i++) {
                h9kring::H9CellId id;
                if (!h9kring::identity_from_uuid(cells[i].uuid, L, &id) ||
                    id.oid != cells[i].oid || id.c2 != cells[i].c2 ||
                    id.ia != cells[i].ia || id.ib != cells[i].ib ||
                    id.ext != (bool)cells[i].ext) {
                    if (id_bad < 5)
                        std::printf("  !! identity_from_uuid mismatch at cell %d "
                                    "(oid=%d c2=%d ia=%lld ib=%lld ext=%d)\n",
                                    i, cells[i].oid, cells[i].c2, (long long)cells[i].ia,
                                    (long long)cells[i].ib, cells[i].ext);
                    id_bad++;
                    continue;
                }
                h9kring::H9CellId nbs[6];
                const int nn = h9kring::neighbors(id, L, nbs);
                std::set<int> got;
                for (int v = 0; v < nn; v++) {
                    auto it = by_id.find({(int64_t)nbs[v].oid, (int64_t)nbs[v].c2,
                                          nbs[v].ia, nbs[v].ib});
                    got.insert(it == by_id.end() ? -2 : it->second);
                }
                std::set<int> geo(adj[i].begin(), adj[i].end());
                if (got != geo) {
                    if (nb_bad < 5)
                        std::printf("  !! h9_kring neighbors mismatch at cell %d: "
                                    "got %zu geo %zu\n", i, got.size(), geo.size());
                    nb_bad++;
                }
            }
            std::printf("  h9_kring.h: identity %ld bad, neighbors %ld bad (of %d)\n",
                        id_bad, nb_bad, n);
            /* k-disk count spot check on an interior cell */
            for (int i = 0; i < n; i++) {
                if (cells[i].oid != 0 || cells[i].ext) continue;
                h9kring::H9CellId id{cells[i].oid, cells[i].c2, cells[i].ia,
                                     cells[i].ib, (bool)cells[i].ext};
                std::vector<std::pair<h9kring::H9CellId,int>> disk;
                for (int k = 1; k <= 3 && k < (L > 1 ? 4 : 2); k++) {
                    h9kring::k_disk(id, L, k, disk);
                    std::printf("  k_disk(k=%d) from cell %d: %zu cells (interior nominal %d)\n",
                                k, i, disk.size(), 1 + 3*k*(k+1));
                }
                break;
            }
        }
        std::printf("\n");
    }
    return 0;
}
