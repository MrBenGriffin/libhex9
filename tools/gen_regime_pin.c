/* gen_regime_pin.c — produce the regime pin goldens.
 *
 * The pin freezes what the addressing chain (authalic series + unit-sphere
 * core + Sphere-L6 wedge-fold warp) PRODUCES, so that any future change to
 * the projection or the warp field fails loudly and deliberately instead of
 * slipping through.
 *
 * Why this exists: the 2.0.0 regime change (WGS84-trained field → via-sphere)
 * moved addresses from around layer 7 downward, and the entire test suite
 * passed through it unchanged — the Edinburgh anchor happens to agree to
 * nibble 9, and nothing else pinned the projection. Westminster's layer-8 bin
 * moved (435878503 → 435878530) with no test to say so.
 *
 * The goldens are only as good as the chain that made them. Generate ONLY
 * from a build where test/via_sphere passes — that is the independent oracle
 * (4,010 rows against the frozen hhg9 Python reference); this tool is a
 * recorder, not a validator.
 *
 * Regenerating is a deliberate act. If a diff appears here and you did not
 * intend to change the projection, the change is a bug.
 *
 *   cmake --build build --target gen_regime_pin
 *   ./build/gen_regime_pin            # rewrites test_data/regime_pin*.tsv
 */
#include "hex9_c.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef H9_REPO_DIR
#define H9_REPO_DIR "."
#endif

/* Named points: landmarks, the anchors the suite already used, and the
 * geometric edge cases (poles, antimeridian, octant corners/seams) where a
 * regime change is most likely to show up first. */
static const struct { double lon, lat; const char *name; } NAMED[] = {
    { -3.190000,  55.950000, "edinburgh"        },
    { -0.127600,  51.507200, "westminster"      },
    {  0.000000,   0.000000, "null-island"      },
    {  139.691700, 35.689500, "tokyo"           },
    { -74.006000,  40.712800, "new-york"        },
    {  151.209300,-33.868800, "sydney"          },
    {  18.423300, -33.918900, "cape-town"       },
    { -43.172900, -22.906800, "rio"             },
    {  37.617800,  55.755800, "moscow"          },
    {  77.209000,  28.613900, "delhi"           },
    {  116.407400, 39.904200, "beijing"         },
    { -99.133200,  19.432600, "mexico-city"     },
    {  31.235700,  30.044400, "cairo"           },
    {  103.819800,  1.352100, "singapore"       },
    { -58.381600, -34.603700, "buenos-aires"    },
    {  0.000000,  90.000000, "north-pole"       },
    {  0.000000, -90.000000, "south-pole"       },
    {  0.000000,  89.999900, "near-north-pole"  },
    {  0.000000, -89.999900, "near-south-pole"  },
    {  180.000000, 0.000000, "antimeridian-eq"  },
    { -180.000000, 0.000000, "antimeridian-neg" },
    {  179.999900, 0.000000, "near-antimeridian"},
    {  90.000000,  0.000000, "equator-90e"      },
    { -90.000000,  0.000000, "equator-90w"      },
    {  45.000000,  35.264390, "octant-corner"    },
    { -45.000000,  35.264390, "octant-corner-nw" },
    {  135.000000,-35.264390, "octant-corner-se" },
    {  0.000000,  45.000000, "prime-45n"        },
    {  0.000000, -45.000000, "prime-45s"        },
    {  90.000000, 45.000000, "seam-90e-45n"     },
};
#define N_NAMED (sizeof NAMED / sizeof NAMED[0])

/* Deterministic equal-area sphere sample. A plain LCG, written out in full so
 * the sequence is identical on every platform and compiler — the goldens must
 * be regenerable anywhere, and rand() is not portable. */
static unsigned long long lcg_state = 0x2545F4914F6CDD1DULL;
static double lcg_next(void)
{
    lcg_state = lcg_state * 6364136223846793005ULL + 1442695040888963407ULL;
    return (double)((lcg_state >> 11) & ((1ULL << 53) - 1)) / (double)(1ULL << 53);
}

#define N_RANDOM 200

/* Points whose full ownership ladder (layer 0..29) is pinned. Ownership is
 * hex9_cell_ancestor — the direct deep fold, NOT transitive, the relation
 * that partitions each layer at exactly 9^d per zone. It is the one most
 * sensitive to a fold-table or projection change. */
static const int OWNER_PTS[] = { 0, 1, 3, 15, 19, 24 };
#define N_OWNER (sizeof OWNER_PTS / sizeof OWNER_PTS[0])
#define OWNER_MAX_LAYER 29

static void hexify(const unsigned char *u, char *out)
{
    static const char *H = "0123456789abcdef";
    for (int i = 0; i < 16; ++i) { out[2*i] = H[u[i] >> 4]; out[2*i+1] = H[u[i] & 15]; }
    out[32] = '\0';
}

int main(void)
{
    char err[256] = {0};
    if (hex9_warp_init(err, sizeof err)) {
        fprintf(stderr, "warp init failed: %s\n", err);
        return 1;
    }

    double lon[N_NAMED + N_RANDOM], lat[N_NAMED + N_RANDOM];
    const char *name[N_NAMED + N_RANDOM];
    static char genname[N_RANDOM][16];

    for (size_t i = 0; i < N_NAMED; ++i) {
        lon[i] = NAMED[i].lon; lat[i] = NAMED[i].lat; name[i] = NAMED[i].name;
    }
    for (int i = 0; i < N_RANDOM; ++i) {
        const size_t k = N_NAMED + (size_t)i;
        lon[k] = lcg_next() * 360.0 - 180.0;
        /* asin of a uniform gives equal-area latitude coverage — a uniform
         * latitude would oversample the poles. */
        lat[k] = asin(lcg_next() * 2.0 - 1.0) * (180.0 / M_PI);
        snprintf(genname[i], sizeof genname[i], "s%03d", i);
        name[k] = genname[i];
    }
    const size_t n = N_NAMED + N_RANDOM;

    /* ── address + curve pin ──────────────────────────────────────────── */
    FILE *fh = fopen(H9_REPO_DIR "/test_data/regime_pin.tsv", "w");
    if (!fh) { fprintf(stderr, "cannot write regime_pin.tsv\n"); return 1; }
    fprintf(fh, "# regime_pin.tsv — GENERATED by tools/gen_regime_pin.c. Do not hand-edit.\n");
    fprintf(fh, "# Freezes the addressing chain: lon/lat -> full uuid -> curve uuid.\n");
    fprintf(fh, "# A diff here means the projection changed. That is either deliberate\n");
    fprintf(fh, "# (regenerate, and bump the major version — every stored address moves)\n");
    fprintf(fh, "# or a bug. There is no third case.\n");
    fprintf(fh, "#\n# name\tlon\tlat\tuuid\tcurve\n");
    for (size_t i = 0; i < n; ++i) {
        unsigned char u[16], c[16];
        char uh[33], ch[33];
        if (hex9_encode(lon[i], lat[i], u)) {
            fprintf(stderr, "encode failed at %s\n", name[i]); fclose(fh); return 1;
        }
        if (hex9_curve(u, c)) {
            fprintf(stderr, "curve failed at %s\n", name[i]); fclose(fh); return 1;
        }
        hexify(u, uh); hexify(c, ch);
        /* %.17g — full round-trip precision. NOT cosmetic: an L23 cell is
         * ~70 µm across and 1e-9 deg is ~110 µm, so a prettier %.9f makes the
         * re-parsed point land in a different cell and the pin fails against
         * itself in the deep tail. The coordinate must reload bit-identical. */
        fprintf(fh, "%s\t%.17g\t%.17g\t%s\t%s\n", name[i], lon[i], lat[i], uh, ch);
    }
    fclose(fh);

    /* ── ownership ladder pin ─────────────────────────────────────────── */
    fh = fopen(H9_REPO_DIR "/test_data/regime_pin_owners.tsv", "w");
    if (!fh) { fprintf(stderr, "cannot write regime_pin_owners.tsv\n"); return 1; }
    fprintf(fh, "# regime_pin_owners.tsv — GENERATED by tools/gen_regime_pin.c. Do not hand-edit.\n");
    fprintf(fh, "# The OWNERSHIP ladder (hex9_cell_ancestor) at every layer 0..%d for a\n", OWNER_MAX_LAYER);
    fprintf(fh, "# few points. Ownership is the non-transitive deep fold that partitions\n");
    fprintf(fh, "# each layer at exactly 9^d per zone — the relation most sensitive to a\n");
    fprintf(fh, "# fold-table or projection change.\n");
    fprintf(fh, "#\n# name\tlayer\towner_uuid\n");
    for (size_t j = 0; j < N_OWNER; ++j) {
        const size_t i = (size_t)OWNER_PTS[j];
        unsigned char u[16];
        if (hex9_encode(lon[i], lat[i], u)) {
            fprintf(stderr, "encode failed at %s\n", name[i]); fclose(fh); return 1;
        }
        for (int L = 0; L <= OWNER_MAX_LAYER; ++L) {
            unsigned char a[16];
            char ah[33];
            if (hex9_cell_ancestor(u, L, a)) {
                fprintf(stderr, "ancestor failed at %s layer %d\n", name[i], L);
                fclose(fh); return 1;
            }
            hexify(a, ah);
            fprintf(fh, "%s\t%d\t%s\n", name[i], L, ah);
        }
    }
    fclose(fh);

    printf("wrote test_data/regime_pin.tsv (%zu points)\n", n);
    printf("wrote test_data/regime_pin_owners.tsv (%zu points x %d layers)\n",
           N_OWNER, OWNER_MAX_LAYER + 1);
    return 0;
}
