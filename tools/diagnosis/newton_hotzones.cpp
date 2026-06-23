/* newton_hotzones.cpp — stress the PROTOTYPE 2D-Newton inversion at the octant
 * HOT ZONES, where the AK-inversion Jacobian degenerates (design memo:
 * "Octahedron vertices / poles — most singular, Jacobian degenerates").
 *
 * A uniform equal-area sample (newton_probe) essentially never lands near a
 * vertex or seam, so its 0% fallback only characterises the easy interior.
 * The octahedron geometry that stresses the map:
 *   - 6 VERTICES = 2 poles (±Z) + 4 equatorial-axis points (±X,±Y at lon
 *     0/90/180/-90, lat 0). One barycentric coord → 1, tan → ∞, ∂lat/∂fy → 0,
 *     det → 0 (the beam itself guards this with skip_nr at barycentric > 0.85).
 *   - 12 EDGES = the seams: the EQUATOR (the 4 ±X/±Y edges) and the
 *     lon ∈ {0,±90,180} MERIDIANS (the 8 pole↔equator edges) — the same seam
 *     meridians characterised by the polar/Jacobian warp work.
 * Vertices themselves are invariant (the pole/axis shortcuts nail them); the
 * hazard is the APPROACH, so we sweep a geodesic distance ladder toward each.
 *
 * Per (zone × distance band) we report, Newton-vs-beam: fallback rate, Newton
 * iters, max |Δboct|, roundtrip geodesic (both), and full-L29 / bin-L16
 * address disagreement. Warp is ON throughout (addressing applies warp_inv).
 *
 * Build:  cmake --build build --target newton_hotzones
 *         ./build/newton_hotzones [K_per_cell]
 */
#include "hex9_c.h"
#include "h9_grid.h"
#include "newton_invert.h"
#include "newton_invert_aj.h"

#include <cstdio>
#include <cstring>
#include <cstdint>
#include <cmath>
#include <cstdlib>
#include <algorithm>

static uint64_t g_s = 0x243f6a8885a308d3ULL;
static double urand() {
    g_s += 0x9e3779b97f4a7c15ULL;
    uint64_t z = g_s;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    z ^= z >> 31;
    return (double)(z >> 11) / 9007199254740992.0;
}
static double haversine_m(double lon1, double lat1, double lon2, double lat2) {
    const double R = 6371008.8, d2r = M_PI/180.0;
    const double dlat = (lat2-lat1)*d2r, dlon = (lon2-lon1)*d2r;
    const double a = std::sin(dlat/2)*std::sin(dlat/2) +
                     std::cos(lat1*d2r)*std::cos(lat2*d2r)*std::sin(dlon/2)*std::sin(dlon/2);
    return 2.0*R*std::asin(std::min(1.0, std::sqrt(a)));
}

/* Geodesic destination on a sphere: from (lon1,lat1) deg, bearing brg rad,
 * angular distance ad rad. Used to place a point a fixed GROUND distance from a
 * vertex/edge with random bearing (covers every approach direction). */
static void dest(double lon1, double lat1, double brg, double ad,
                 double *lon2, double *lat2) {
    const double d2r = M_PI/180.0, r2d = 180.0/M_PI;
    const double la1 = lat1*d2r, lo1 = lon1*d2r;
    const double la2 = std::asin(std::sin(la1)*std::cos(ad) +
                                 std::cos(la1)*std::sin(ad)*std::cos(brg));
    const double lo2 = lo1 + std::atan2(std::sin(brg)*std::sin(ad)*std::cos(la1),
                                        std::cos(ad) - std::sin(la1)*std::sin(la2));
    *lat2 = la2*r2d;
    *lon2 = lo2*r2d;
    while (*lon2 >  180.0) *lon2 -= 360.0;
    while (*lon2 < -180.0) *lon2 += 360.0;
}

enum ZoneKind { Z_VERTEX, Z_MERIDIAN, Z_EQUATOR };
struct Zone { const char *name; ZoneKind kind; double lon0, lat0; };

/* point a ground distance `dist_m` into the given zone, random approach */
static void gen_point(const Zone &z, double dist_m, double *lon, double *lat) {
    const double R = 6371008.8;
    const double ad = dist_m / R;                 /* angular distance, rad */
    if (z.kind == Z_VERTEX) {
        dest(z.lon0, z.lat0, urand()*2.0*M_PI, ad, lon, lat);
    } else if (z.kind == Z_MERIDIAN) {
        /* random lat on the seam meridian, away from the poles/equator
         * vertices; offset perpendicular (±) by dist_m */
        const double la = (15.0 + urand()*60.0) * (urand() < 0.5 ? 1.0 : -1.0);
        const double brg = (urand() < 0.5 ? M_PI*0.5 : -M_PI*0.5);   /* ±east = across meridian */
        dest(z.lon0, la, brg, ad, lon, lat);
    } else { /* Z_EQUATOR: random lon away from the 4 vertices, offset ±N/S */
        double lo;
        do { lo = urand()*360.0 - 180.0; }
        while (std::fabs(std::fmod(std::fabs(lo), 90.0) - 45.0) > 44.0);  /* keep ≥1° from a vertex */
        const double brg = (urand() < 0.5 ? 0.0 : M_PI);             /* ±north = across equator */
        dest(lo, 0.0, brg, ad, lon, lat);
    }
}

int main(int argc, char **argv) {
    char err[256] = {0};
    if (hex9_warp_init(err, sizeof err)) { printf("warp init FAILED: %s\n", err); return 1; }
    const long K = (argc > 1) ? atol(argv[1]) : 1000;
    /* solver under test: arg2 "aj" → analytic-Jacobian variant, else FD */
    const bool use_aj = (argc > 2) && std::strcmp(argv[2], "aj") == 0;
    auto solve = [&](double lr, double ar, NStat *st) {
        return use_aj ? lonlat_to_boct_newton_aj(lr, ar, st)
                      : lonlat_to_boct_newton(lr, ar, st);
    };
    printf("solver: %s\n", use_aj ? "AJ-Newton (analytic Jacobian)" : "FD-Newton");

    const Zone zones[] = {
        {"vtx +Z (N pole)", Z_VERTEX,    0.0,  90.0},
        {"vtx -Z (S pole)", Z_VERTEX,    0.0, -90.0},
        {"vtx +X",          Z_VERTEX,    0.0,   0.0},
        {"vtx +Y",          Z_VERTEX,   90.0,   0.0},
        {"vtx -X",          Z_VERTEX,  180.0,   0.0},
        {"vtx -Y",          Z_VERTEX,  -90.0,   0.0},
        {"seam lon=0",      Z_MERIDIAN,  0.0,   0.0},
        {"seam lon=90",     Z_MERIDIAN, 90.0,   0.0},
        {"seam lon=180",    Z_MERIDIAN,180.0,   0.0},
        {"seam lon=-90",    Z_MERIDIAN,-90.0,   0.0},
        {"equator edge",    Z_EQUATOR,   0.0,   0.0},
    };
    const int NZ = (int)(sizeof(zones)/sizeof(zones[0]));
    const double dist_m[] = { 1e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0, 1e3, 1e4, 1e5 };
    const char  *dist_l[] = { "1mm", "1cm", "10cm","1m", "10m", "100m","1km","10km","100km"};
    const int ND = (int)(sizeof(dist_m)/sizeof(dist_m[0]));

    printf("K=%ld points per (zone × band), warp ON\n", K);
    printf("%-16s %-6s  %5s %5s  %8s  %7s %7s  %5s %5s\n",
           "zone","dist","fb","itMx","maxDuvw","rtN_nm","rtB_nm","L29","L16");
    printf("%.*s\n", 86, "----------------------------------------------------------------------------------------");

    long tot_fb = 0, tot_l29 = 0, tot_l16 = 0, tot_n = 0;

    for (int zi = 0; zi < NZ; ++zi) {
        for (int di = 0; di < ND; ++di) {
            long fb = 0, l29 = 0, l16 = 0; int itmax = 0;
            double maxduvw = 0.0, rtN = 0.0, rtB = 0.0;
            for (long k = 0; k < K; ++k) {
                double lon_d, lat_d;
                gen_point(zones[zi], dist_m[di], &lon_d, &lat_d);
                const double lon_r = lon_d*M_PI/180.0, lat_r = lat_d*M_PI/180.0;

                NStat st;
                const H9BOct nb = solve(lon_r, lat_r, &st);
                const H9BOct bb = h9_lonlat_to_boct_beam(lon_r, lat_r);

                const double dm = std::max(std::fabs(nb.u-bb.u),
                                  std::max(std::fabs(nb.v-bb.v), std::fabs(nb.w-bb.w)));
                if (dm > maxduvw) maxduvw = dm;
                if (st.fellback) fb++;
                if (st.iters > itmax) itmax = st.iters;

                uint8_t un[16], ub[16];
                h9_boct_to_uuid_beam(nb, 6, 6, un);
                h9_boct_to_uuid_beam(bb, 6, 6, ub);
                if (std::memcmp(un, ub, 16) != 0) l29++;
                uint8_t bn[16], bbn[16];
                hex9_bin(un, 16, bn); hex9_bin(ub, 16, bbn);
                if (std::memcmp(bn, bbn, 16) != 0) l16++;

                double dlon, dlat;
                hex9_decode(un, &dlon, &dlat);
                const double rn = haversine_m(lon_d, lat_d, dlon, dlat);
                hex9_decode(ub, &dlon, &dlat);
                const double rb = haversine_m(lon_d, lat_d, dlon, dlat);
                if (rn > rtN) rtN = rn;
                if (rb > rtB) rtB = rb;
            }
            printf("%-16s %-6s  %5ld %5d  %8.1e  %7.1f %7.1f  %5ld %5ld\n",
                   zones[zi].name, dist_l[di], fb, itmax, maxduvw,
                   rtN*1e9, rtB*1e9, l29, l16);
            tot_fb += fb; tot_l29 += l29; tot_l16 += l16; tot_n += K;
        }
    }
    printf("%.*s\n", 86, "----------------------------------------------------------------------------------------");
    printf("TOTAL  n=%ld   fallback=%ld (%.4f%%)   L29 diff=%ld (%.3f%%)   L16 diff=%ld (%.4f%%)\n",
           tot_n, tot_fb, 100.0*tot_fb/tot_n, tot_l29, 100.0*tot_l29/tot_n,
           tot_l16, 100.0*tot_l16/tot_n);
    return 0;
}
