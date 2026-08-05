/* e4h_parity.c — byte-exact E4H parity vs the hhg9 reference corpus.
 *
 * test_data/e4h_pin.tsv rows were minted by hhg9/h9/e4h.py (the normative
 * reference; generator tools/gen_e4h_pin.py, knife-edge points filtered).
 * The C exact classifier must reproduce every uuid BYTE-IDENTICALLY, spell
 * the same label, and place decode/partner representatives within 1e-8
 * degrees. A uuid mismatch is a regime violation or a broken port — never a
 * tolerance issue; the corpus margin filter (1e-5 in the residual frame)
 * already excludes every point where the frozen snap could legitimately
 * disagree with the reference's double descent.
 */
#include "hex9_c.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef H9_REPO_DIR
#define H9_REPO_DIR "."
#endif

static void hexify(const unsigned char *u, char *out)
{
    static const char *H = "0123456789abcdef";
    for (int i = 0; i < 16; ++i) { out[2*i] = H[u[i] >> 4]; out[2*i+1] = H[u[i] & 15]; }
    out[32] = '\0';
}

static double wrap180(double a)
{
    a = fmod(a + 180.0, 360.0);
    if (a < 0) a += 360.0;
    return a - 180.0;
}

/* angular closeness in degrees: lat direct, lon scaled by cos(lat) */
static int close_deg(double lon_a, double lat_a, double lon_b, double lat_b,
                     double tol)
{
    const double c = cos(lat_a * (3.14159265358979323846 / 180.0));
    return fabs(lat_a - lat_b) <= tol &&
           fabs(wrap180(lon_a - lon_b)) * (fabs(c) > 1e-6 ? fabs(c) : 1.0)
               <= tol;
}

int main(void)
{
    char err[256];
    if (hex9_init(err, sizeof err)) { printf("FAIL: init: %s\n", err); return 1; }

    FILE *fh = fopen(H9_REPO_DIR "/test_data/e4h_pin.tsv", "r");
    if (!fh) { printf("FAIL: cannot open e4h_pin.tsv\n"); return 1; }

    char line[1024];
    int rows = 0, fails = 0, parsed = 0;
    while (fgets(line, sizeof line, fh)) {
        if (line[0] == '#' || line[0] == '\n') continue;
        char name[64], want_uuid[64], want_label[64];
        double lon, lat, dlon, dlat, plon, plat;
        int layer, depth;
        if (sscanf(line, "%63s %lf %lf %d %d %63s %63s %lf %lf %lf %lf",
                   name, &lon, &lat, &layer, &depth, want_uuid, want_label,
                   &dlon, &dlat, &plon, &plat) != 11) {
            printf("FAIL: bad row: %s", line);
            return 1;
        }
        rows++;

        uint8_t u[16];
        char got[40];
        if (hex9_e4h_encode(lon, lat, layer, depth, u)) {
            printf("FAIL %s: encode rc\n", name);
            fails++;
            continue;
        }
        hexify(u, got);
        if (strcmp(got, want_uuid) != 0) {
            printf("FAIL %s: uuid %s != %s\n", name, got, want_uuid);
            fails++;
            continue;
        }
        char lab[64];
        if (hex9_e4h_label(u, lab, sizeof lab) < 0 ||
            strcmp(lab, want_label) != 0) {
            printf("FAIL %s: label '%s' != '%s'\n", name, lab, want_label);
            fails++;
            continue;
        }
        /* label parse: bare host labels carry the F1 split-hex ambiguity
         * (doctrine); when the parse resolves, it must resolve HERE. */
        uint8_t v[16];
        if (hex9_e4h_parse_label(lab, v) == 0) {
            parsed++;
            if (memcmp(u, v, 16) != 0) {
                printf("FAIL %s: parse_label roundtrip\n", name);
                fails++;
                continue;
            }
        }
        double glon, glat;
        if (hex9_e4h_decode(u, &glon, &glat) ||
            !close_deg(glon, glat, dlon, dlat, 1e-8)) {
            printf("FAIL %s: decode (%.12f,%.12f) != (%.12f,%.12f)\n",
                   name, glon, glat, dlon, dlat);
            fails++;
            continue;
        }
        if (hex9_e4h_partner(u, &glon, &glat) ||
            !close_deg(glon, glat, plon, plat, 1e-8)) {
            printf("FAIL %s: partner (%.12f,%.12f) != (%.12f,%.12f)\n",
                   name, glon, glat, plon, plat);
            fails++;
            continue;
        }
    }
    fclose(fh);

    if (rows < 100) { printf("FAIL: only %d corpus rows\n", rows); return 1; }
    if (parsed < rows / 2) {
        printf("FAIL: bare-label parse resolved only %d/%d\n", parsed, rows);
        return 1;
    }
    if (fails) { printf("FAIL: %d/%d rows\n", fails, rows); return 1; }
    printf("OK: e4h parity — %d rows byte-exact (labels %d parsed), "
           "decode/partner <= 1e-8 deg\n", rows, parsed);
    return 0;
}
