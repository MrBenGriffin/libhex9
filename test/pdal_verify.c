/* pdal_verify.c — validate writers.text output of filters.hex9 against the C
 * ABI: WoctX/Y/Z == hex9_to_woct(X,Y), |L1|==1, and from_woct round-trips X,Y.
 * argv[1] = the CSV written by the pipeline (X,Y,Z,WoctX,WoctY,WoctZ header). */
#include "hex9_c.h"
#include <stdio.h>
#include <math.h>

int main(int argc, char **argv) {
    if (argc < 2) { printf("usage: pdal_verify out.csv\n"); return 2; }
    char err[256] = {0};
    hex9_warp_init(err, sizeof err);
    FILE *f = fopen(argv[1], "r");
    if (!f) { printf("cannot open %s\n", argv[1]); return 1; }
    char line[512];
    if (!fgets(line, sizeof line, f)) { printf("empty file\n"); return 1; }  /* header */
    int rows = 0, bad = 0;
    double max_dw = 0, max_l1 = 0, max_rt = 0;
    double lon, lat, z, wx, wy, wz;
    while (fgets(line, sizeof line, f)) {
        if (sscanf(line, "%lf,%lf,%lf,%lf,%lf,%lf",
                   &lon, &lat, &z, &wx, &wy, &wz) != 6) continue;
        rows++;
        double exp[3];
        hex9_to_woct(lon, lat, exp);
        double dw = fabs(exp[0]-wx) + fabs(exp[1]-wy) + fabs(exp[2]-wz);
        double l1 = fabs(wx) + fabs(wy) + fabs(wz);
        double rlon, rlat, xyz[3] = {wx, wy, wz};
        hex9_from_woct(xyz, &rlon, &rlat);
        double dlo = fabs(rlon - lon); if (dlo > 180) dlo = 360 - dlo;
        double rt = (fabs(fabs(lat)-90) < 1e-6) ? fabs(rlat-lat)
                    : hypot(dlo*cos(lat*M_PI/180), fabs(rlat-lat));
        if (dw > 1e-9 || fabs(l1-1) > 1e-9 || rt > 1e-4) bad++;
        if (dw > max_dw) max_dw = dw;
        if (fabs(l1-1) > max_l1) max_l1 = fabs(l1-1);
        if (rt > max_rt) max_rt = rt;
    }
    fclose(f);
    printf("rows=%d  max|plugin-abi|=%.2e  max|L1-1|=%.2e  max_rt=%.2e deg  bad=%d\n",
           rows, max_dw, max_l1, max_rt, bad);
    printf(bad == 0 && rows > 0 ? "PDAL PLUGIN OK\n" : "PDAL PLUGIN FAILED\n");
    return (bad == 0 && rows > 0) ? 0 : 1;
}
