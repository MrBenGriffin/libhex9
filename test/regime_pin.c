/* regime_pin.c — assert the addressing chain still produces the frozen
 * addresses in test_data/regime_pin*.tsv.
 *
 * This is the guard that was missing at 2.0.0: the WGS84 → via-sphere regime
 * change moved addresses from around layer 7 downward and the whole suite
 * passed through it, because the one anchor it pinned (Edinburgh, layer 8)
 * happened to agree across both regimes while Westminster's did not.
 *
 * A failure here is NOT a flaky test. It means one of:
 *   - the projection or warp field changed (deliberate ⇒ regenerate the pin
 *     with tools/gen_regime_pin and bump the MAJOR version, because every
 *     stored address in every downstream database has just moved), or
 *   - something changed it by accident, which is a bug.
 *
 * Correctness of the pinned values rests on test/via_sphere (4,010 rows vs
 * the frozen hhg9 Python reference). This test pins; that test validates.
 * If both fail, believe via_sphere.
 */
#include "hex9_c.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef H9_REPO_DIR
#define H9_REPO_DIR "."
#endif

#define MAX_REPORT 8

static void hexify(const unsigned char *u, char *out)
{
    static const char *H = "0123456789abcdef";
    for (int i = 0; i < 16; ++i) { out[2*i] = H[u[i] >> 4]; out[2*i+1] = H[u[i] & 15]; }
    out[32] = '\0';
}

/* First differing nibble, or -1 if equal — reported so a failure says how
 * DEEP the divergence is, which is what tells you whether it is a projection
 * change (shallow, layer 5-10) or float noise in the deep tail. */
static int first_diff(const char *a, const char *b)
{
    for (int i = 0; i < 32; ++i) if (a[i] != b[i]) return i;
    return -1;
}

static int check_addresses(void)
{
    FILE *fh = fopen(H9_REPO_DIR "/test_data/regime_pin.tsv", "r");
    if (!fh) { printf("FAIL: cannot open regime_pin.tsv\n"); return 1; }

    char line[512];
    int fails = 0, rows = 0, reported = 0;
    while (fgets(line, sizeof line, fh)) {
        if (line[0] == '#' || line[0] == '\n') continue;
        char name[64], want_u[64], want_c[64];
        double lon, lat;
        if (sscanf(line, "%63s %lf %lf %63s %63s",
                   name, &lon, &lat, want_u, want_c) != 5) {
            printf("FAIL: malformed row: %s", line);
            ++fails; continue;
        }
        ++rows;

        unsigned char u[16], c[16];
        char got_u[33], got_c[33];
        if (hex9_encode(lon, lat, u)) {
            printf("FAIL: %s: encode error\n", name); ++fails; continue;
        }
        if (hex9_curve(u, c)) {
            printf("FAIL: %s: curve error\n", name); ++fails; continue;
        }
        hexify(u, got_u); hexify(c, got_c);

        if (strcmp(got_u, want_u)) {
            if (reported++ < MAX_REPORT)
                printf("FAIL uuid  %-18s (%.6f, %.6f)\n  want %s\n  got  %s\n  diverges at nibble %d\n",
                       name, lon, lat, want_u, got_u, first_diff(want_u, got_u));
            ++fails;
        }
        if (strcmp(got_c, want_c)) {
            if (reported++ < MAX_REPORT)
                printf("FAIL curve %-18s (%.6f, %.6f)\n  want %s\n  got  %s\n  diverges at nibble %d\n",
                       name, lon, lat, want_c, got_c, first_diff(want_c, got_c));
            ++fails;
        }
    }
    fclose(fh);
    if (rows == 0) { printf("FAIL: regime_pin.tsv had no rows\n"); return 1; }
    printf("addresses+curves: %d rows, %d mismatches %s\n",
           rows, fails, fails ? "FAIL" : "PASS");
    return fails;
}

static int check_owners(void)
{
    FILE *fh = fopen(H9_REPO_DIR "/test_data/regime_pin_owners.tsv", "r");
    if (!fh) { printf("FAIL: cannot open regime_pin_owners.tsv\n"); return 1; }

    /* The ladder file gives name+layer+owner; re-derive the point's uuid from
     * the address pin so the two files are cross-checked against each other. */
    char line[512];
    int fails = 0, rows = 0, reported = 0;
    char cur_name[64] = {0};
    unsigned char cur_uuid[16];
    int have_uuid = 0;

    while (fgets(line, sizeof line, fh)) {
        if (line[0] == '#' || line[0] == '\n') continue;
        char name[64], want[64];
        int layer;
        if (sscanf(line, "%63s %d %63s", name, &layer, want) != 3) {
            printf("FAIL: malformed owner row: %s", line);
            ++fails; continue;
        }
        ++rows;

        if (strcmp(name, cur_name)) {
            /* new point — look its lon/lat up in the address pin */
            FILE *pf = fopen(H9_REPO_DIR "/test_data/regime_pin.tsv", "r");
            if (!pf) { printf("FAIL: cannot reopen regime_pin.tsv\n"); fclose(fh); return 1; }
            char pline[512];
            have_uuid = 0;
            while (fgets(pline, sizeof pline, pf)) {
                if (pline[0] == '#' || pline[0] == '\n') continue;
                char pname[64], pu[64], pc[64];
                double plon, plat;
                if (sscanf(pline, "%63s %lf %lf %63s %63s",
                           pname, &plon, &plat, pu, pc) != 5) continue;
                if (!strcmp(pname, name)) {
                    if (!hex9_encode(plon, plat, cur_uuid)) have_uuid = 1;
                    break;
                }
            }
            fclose(pf);
            snprintf(cur_name, sizeof cur_name, "%s", name);
            if (!have_uuid) {
                printf("FAIL: owner point %s not found in regime_pin.tsv\n", name);
                ++fails; continue;
            }
        }
        if (!have_uuid) { ++fails; continue; }

        unsigned char a[16];
        char got[33];
        if (hex9_cell_ancestor(cur_uuid, layer, a)) {
            printf("FAIL: %s layer %d: ancestor error\n", name, layer);
            ++fails; continue;
        }
        hexify(a, got);
        if (strcmp(got, want)) {
            if (reported++ < MAX_REPORT)
                printf("FAIL owner %-18s layer %2d\n  want %s\n  got  %s\n",
                       name, layer, want, got);
            ++fails;
        }
    }
    fclose(fh);
    if (rows == 0) { printf("FAIL: regime_pin_owners.tsv had no rows\n"); return 1; }
    printf("ownership ladder: %d rows, %d mismatches %s\n",
           rows, fails, fails ? "FAIL" : "PASS");
    return fails;
}

int main(void)
{
    char err[256] = {0};
    if (hex9_warp_init(err, sizeof err)) {
        printf("FAIL: warp init: %s\n", err);
        return 1;
    }

    int fails = check_addresses() + check_owners();

    printf(fails ? "=== regime_pin: FAIL (%d) ===\n"
                 : "=== regime_pin: PASS ===\n", fails);
    return fails ? 1 : 0;
}
