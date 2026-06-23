/* newton_aj_probe.cpp — verify + benchmark the ANALYTIC-JACOBIAN inversion
 * (newton_invert_aj.h) against the FD-Newton (newton_invert.h) and the beam.
 *
 *   1. Jacobian verification: analytic J3 vs central differences (interior).
 *   2. Output-preserving: AJ boct vs beam (max |Δuvw|), bins L8/12/16/20, full L29.
 *   3. Roundtrip fidelity (encode→decode geodesic, warp ON), AJ vs beam.
 *   4. Timing: beam vs FD-Newton vs AJ-Newton (inversion only), mean iters.
 *
 * Build:  cmake --build build --target newton_aj_probe && ./build/newton_aj_probe [N]
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
#include <chrono>
#include <vector>
#include <algorithm>

static uint64_t g_s = 0x9e3779b97f4a7c15ULL;
static double urand() {
    g_s += 0x9e3779b97f4a7c15ULL; uint64_t z = g_s;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL; z ^= z >> 31;
    return (double)(z >> 11) / 9007199254740992.0;
}
static double haversine_m(double lo1,double la1,double lo2,double la2){
    const double R=6371008.8,d2r=M_PI/180.0; const double dla=(la2-la1)*d2r,dlo=(lo2-lo1)*d2r;
    const double a=std::sin(dla/2)*std::sin(dla/2)+std::cos(la1*d2r)*std::cos(la2*d2r)*std::sin(dlo/2)*std::sin(dlo/2);
    return 2.0*R*std::asin(std::min(1.0,std::sqrt(a)));
}

int main(int argc, char **argv) {
    char err[256] = {0};
    if (hex9_warp_init(err, sizeof err)) { printf("warp init FAILED: %s\n", err); return 1; }
    const long N = (argc > 1) ? atol(argv[1]) : 100000;

    /* 1. Jacobian verification */
    const double jerr = aj_jacobian_max_relerr(50000, 0xcafef00dULL);
    printf("── analytic Jacobian vs central differences ──\n");
    printf("max relative error (50k interior pts): %.3e   %s\n\n",
           jerr, jerr < 1e-6 ? "OK" : "*** SUSPECT ***");

    long boct_bad=0, full_diff=0; double max_duvw=0.0;
    const int bl[]={8,12,16,20}; long bin_diff[4]={0,0,0,0};
    long fell=0; double sum_it=0.0;
    double rtA_max=0, rtA_sum=0, rtB_max=0, rtB_sum=0;

    for (long i = 0; i < N; ++i) {
        const double lon_d=(urand()*2-1)*180.0, lat_d=std::asin(urand()*2-1)*180.0/M_PI;
        const double lon_r=lon_d*M_PI/180.0, lat_r=lat_d*M_PI/180.0;
        NStat st;
        const H9BOct ab = lonlat_to_boct_newton_aj(lon_r, lat_r, &st);
        const H9BOct bb = h9_lonlat_to_boct_beam(lon_r, lat_r);
        const double dm = std::max(std::fabs(ab.u-bb.u), std::max(std::fabs(ab.v-bb.v), std::fabs(ab.w-bb.w)));
        if (dm > max_duvw) max_duvw = dm;
        if (dm > 1e-12) boct_bad++;
        if (st.fellback) fell++;
        sum_it += st.iters;

        uint8_t ua[16], ub[16];
        h9_boct_to_uuid_beam(ab, 6, 6, ua);
        h9_boct_to_uuid_beam(bb, 6, 6, ub);
        if (std::memcmp(ua, ub, 16) != 0) full_diff++;
        for (int b=0;b<4;++b){ uint8_t x[16],y[16]; hex9_bin(ua,bl[b],x); hex9_bin(ub,bl[b],y);
            if (std::memcmp(x,y,16)!=0) bin_diff[b]++; }
        double dlon,dlat;
        hex9_decode(ua,&dlon,&dlat); const double ra=haversine_m(lon_d,lat_d,dlon,dlat);
        hex9_decode(ub,&dlon,&dlat); const double rb=haversine_m(lon_d,lat_d,dlon,dlat);
        rtA_sum+=ra; if(ra>rtA_max)rtA_max=ra; rtB_sum+=rb; if(rb>rtB_max)rtB_max=rb;
    }
    printf("── AJ-Newton vs beam (N=%ld equal-area) ──\n", N);
    printf("max |Δ(u,v,w)|        : %.3e\n", max_duvw);
    printf("boct |Δuvw|>1e-12     : %ld\n", boct_bad);
    printf("full L29 disagree     : %ld  (%.3f%%)\n", full_diff, 100.0*full_diff/N);
    for (int b=0;b<4;++b) printf("bin L%-2d disagree      : %ld\n", bl[b], bin_diff[b]);
    printf("would fall back       : %ld  (%.4f%%)\n", fell, 100.0*fell/N);
    printf("mean iters            : %.3f\n", sum_it/N);
    printf("roundtrip AJ   max %.2f nm  mean %.2f nm\n", rtA_max*1e9, rtA_sum/N*1e9);
    printf("roundtrip beam max %.2f nm  mean %.2f nm\n", rtB_max*1e9, rtB_sum/N*1e9);

    /* 4. timing */
    const long T=(N<40000)?N:40000;
    std::vector<double> lr(T),ar(T);
    for(long i=0;i<T;++i){ lr[i]=((urand()*2-1)*180.0)*M_PI/180.0; ar[i]=(std::asin(urand()*2-1)*180.0/M_PI)*M_PI/180.0; }
    using clk=std::chrono::steady_clock; volatile double sink=0; NStat st;
    auto bench=[&](const char*nm, H9BOct(*fn)(double,double,NStat*)){
        fn(lr[0],ar[0],&st);
        auto t0=clk::now(); for(long i=0;i<T;++i){H9BOct b=fn(lr[i],ar[i],&st); sink+=b.u;} auto t1=clk::now();
        return std::chrono::duration<double,std::micro>(t1-t0).count()/T;
    };
    h9_lonlat_to_boct_beam(lr[0],ar[0]);
    auto t0=clk::now(); for(long i=0;i<T;++i){H9BOct b=h9_lonlat_to_boct_beam(lr[i],ar[i]); sink+=b.u;} auto t1=clk::now();
    const double beam_us=std::chrono::duration<double,std::micro>(t1-t0).count()/T;
    const double fd_us = bench("fd",  lonlat_to_boct_newton);
    const double aj_us = bench("aj",  lonlat_to_boct_newton_aj);
    printf("\n── timing (inversion only) ──\n");
    printf("beam      : %8.3f µs/pt  (%.0f pts/s)\n", beam_us, 1e6/beam_us);
    printf("FD-Newton : %8.3f µs/pt  (%.0f pts/s)  %.1f×\n", fd_us, 1e6/fd_us, beam_us/fd_us);
    printf("AJ-Newton : %8.3f µs/pt  (%.0f pts/s)  %.1f×   (sink %.1f)\n", aj_us, 1e6/aj_us, beam_us/aj_us, (double)sink);
    return 0;
}
