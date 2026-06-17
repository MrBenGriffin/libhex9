/* hex9_ext.cpp — nanobind bindings for libhex9.
 *
 * Thin numpy wrapper over the C ABI batch functions. Arrays cross the boundary
 * once; the GIL is released around the OpenMP loop in libhex9, so this is both
 * parallel and non-blocking to other Python threads.
 */
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/optional.h>

#include <optional>

#include "hex9_c.h"

#include <cstdint>
#include <stdexcept>
#include <string>

namespace nb = nanobind;

using f64_1d  = nb::ndarray<const double, nb::ndim<1>, nb::c_contig, nb::device::cpu>;
using u8_2d_in = nb::ndarray<const uint8_t, nb::ndim<2>, nb::c_contig, nb::device::cpu>;

/* Initialise the warp once; raises on failure. */
static void warp_init() {
    char err[256] = {0};
    if (hex9_warp_init(err, sizeof err))
        throw std::runtime_error(std::string("hex9 warp init failed: ") + err);
}

/* encode(lon[n], lat[n]) -> uint8[n,16] */
static nb::ndarray<nb::numpy, uint8_t, nb::ndim<2>>
encode(f64_1d lon, f64_1d lat) {
    const size_t n = lon.shape(0);
    if (lat.shape(0) != n) throw std::runtime_error("lon and lat must be the same length");
    uint8_t *out = new uint8_t[n * 16];
    {
        nb::gil_scoped_release release;
        hex9_encode_many(lon.data(), lat.data(), n, out);
    }
    nb::capsule owner(out, [](void *p) noexcept { delete[] static_cast<uint8_t *>(p); });
    return nb::ndarray<nb::numpy, uint8_t, nb::ndim<2>>(out, {n, 16}, owner);
}

/* decode(uuid[n,16]) -> (lon[n], lat[n]) */
static nb::tuple decode(u8_2d_in uuid) {
    const size_t n = uuid.shape(0);
    if (uuid.shape(1) != 16) throw std::runtime_error("uuid array must be (n, 16)");
    double *lon = new double[n];
    double *lat = new double[n];
    {
        nb::gil_scoped_release release;
        hex9_decode_many(uuid.data(), n, lon, lat);
    }
    nb::capsule lon_o(lon, [](void *p) noexcept { delete[] static_cast<double *>(p); });
    nb::capsule lat_o(lat, [](void *p) noexcept { delete[] static_cast<double *>(p); });
    return nb::make_tuple(
        nb::ndarray<nb::numpy, double, nb::ndim<1>>(lon, {n}, lon_o),
        nb::ndarray<nb::numpy, double, nb::ndim<1>>(lat, {n}, lat_o));
}

/* bin(uuid[n,16], layer) -> uint8[n,16] */
static nb::ndarray<nb::numpy, uint8_t, nb::ndim<2>>
bin(u8_2d_in uuid, int layer) {
    const size_t n = uuid.shape(0);
    if (uuid.shape(1) != 16) throw std::runtime_error("uuid array must be (n, 16)");
    uint8_t *out = new uint8_t[n * 16];
    int rc;
    {
        nb::gil_scoped_release release;
        rc = hex9_bin_many(uuid.data(), layer, n, out);
    }
    if (rc) { delete[] out; throw std::runtime_error("hex9.bin: layer out of range (0..29)"); }
    nb::capsule owner(out, [](void *p) noexcept { delete[] static_cast<uint8_t *>(p); });
    return nb::ndarray<nb::numpy, uint8_t, nb::ndim<2>>(out, {n, 16}, owner);
}

using u8_1d = nb::ndarray<const uint8_t, nb::ndim<1>, nb::c_contig, nb::device::cpu>;

/* helper: wrap a heap double[] as a numpy ndarray that owns + frees it */
template <int Ndim>
static nb::ndarray<nb::numpy, double, nb::ndim<Ndim>>
own_f64(double *p, std::initializer_list<size_t> shape) {
    nb::capsule owner(p, [](void *q) noexcept { delete[] static_cast<double *>(q); });
    return nb::ndarray<nb::numpy, double, nb::ndim<Ndim>>(p, shape, owner);
}

/* cell(uuid[16], layer, densify=0) -> (npoints, 2) closed lon/lat ring */
static nb::ndarray<nb::numpy, double, nb::ndim<2>>
cell(u8_1d uuid, int layer, int densify) {
    if (uuid.shape(0) != 16) throw std::runtime_error("uuid must be a (16,) uint8 array");
    int P = hex9_ring_npoints(densify);
    if (P < 0) throw std::runtime_error("hex9.cell: densify out of range (0..9)");
    double *r = new double[(size_t)P * 2];
    int got = hex9_cell_ring(uuid.data(), layer, densify, r, P);
    if (got < 0) { delete[] r; throw std::runtime_error("hex9.cell: invalid layer/densify"); }
    return own_f64<2>(r, {(size_t)P, 2});
}

/* grid(lon_min,lat_min,lon_max,lat_max, layer, densify=0, max_cells=0)
 *   -> (uuids[n,16], centroids[n,2], rings)
 * rings: densify==0 -> flat (n,6,2) of the 6 hex corners;
 *        densify>0  -> list of n (npoints,2) closed-ring arrays. */
static nb::tuple grid(double lon_min, double lat_min, double lon_max, double lat_max,
                      int layer, int densify, int64_t max_cells) {
    char err[256] = {0};
    hex9_grid *g = nullptr;
    {
        nb::gil_scoped_release release;
        g = hex9_grid_create(lon_min, lat_min, lon_max, lat_max,
                             layer, densify, max_cells, err, sizeof err);
    }
    if (!g) throw std::runtime_error(std::string("hex9.grid: ") + err);
    const int n = hex9_grid_count(g);

    uint8_t *uu  = new uint8_t[(size_t)n * 16];
    double  *cen = new double[(size_t)n * 2];
    for (int i = 0; i < n; ++i) {
        hex9_grid_cell_uuid(g, i, uu + (size_t)i * 16);
        hex9_grid_cell_centroid(g, i, &cen[2*i], &cen[2*i + 1]);
    }
    nb::capsule uu_o(uu, [](void *q) noexcept { delete[] static_cast<uint8_t *>(q); });
    auto uuids = nb::ndarray<nb::numpy, uint8_t, nb::ndim<2>>(uu, {(size_t)n, 16}, uu_o);
    auto cents = own_f64<2>(cen, {(size_t)n, 2});

    nb::object rings;
    const int P = hex9_ring_npoints(densify);
    if (densify == 0) {
        /* flat (n,6,2): first 6 of the 7-point closed ring (drop the dup close) */
        double *r = new double[(size_t)n * 6 * 2];
        double tmp[7 * 2];
        for (int i = 0; i < n; ++i) {
            hex9_grid_cell_ring(g, i, 0, tmp, 7);
            std::memcpy(r + (size_t)i * 12, tmp, 12 * sizeof(double));
        }
        rings = nb::cast(own_f64<3>(r, {(size_t)n, 6, 2}));
    } else {
        nb::list lst;
        for (int i = 0; i < n; ++i) {
            double *r = new double[(size_t)P * 2];
            hex9_grid_cell_ring(g, i, densify, r, P);
            lst.append(own_f64<2>(r, {(size_t)P, 2}));
        }
        rings = lst;
    }
    hex9_grid_destroy(g);
    return nb::make_tuple(uuids, cents, rings);
}

/* neighbors(uuid[n,16], layer) -> (uint8[n,6,16], int32[n] counts).
 * Unused neighbour slots are zero-filled; counts[i] is 5 or 6 (or -1 on a
 * per-row error). OpenMP across rows — the kring core is read-only. */
static nb::tuple neighbors(u8_2d_in uuid, int layer) {
    const size_t n = uuid.shape(0);
    if (uuid.shape(1) != 16) throw std::runtime_error("uuid array must be (n, 16)");
    uint8_t *out = new uint8_t[n * 6 * 16]();
    int32_t *cnt = new int32_t[n];
    {
        nb::gil_scoped_release release;
        const ptrdiff_t N = (ptrdiff_t)n;
        #pragma omp parallel for schedule(static)
        for (ptrdiff_t i = 0; i < N; ++i)
            cnt[i] = hex9_neighbors(uuid.data() + (size_t)i * 16, layer,
                                    out + (size_t)i * 6 * 16);
    }
    nb::capsule out_o(out, [](void *p) noexcept { delete[] static_cast<uint8_t *>(p); });
    nb::capsule cnt_o(cnt, [](void *p) noexcept { delete[] static_cast<int32_t *>(p); });
    return nb::make_tuple(
        nb::ndarray<nb::numpy, uint8_t, nb::ndim<3>>(out, {n, 6, 16}, out_o),
        nb::ndarray<nb::numpy, int32_t, nb::ndim<1>>(cnt, {n}, cnt_o));
}

/* k_ring/k_disk(uuid[16], layer, k) -> uint8[m,16] */
static nb::ndarray<nb::numpy, uint8_t, nb::ndim<2>>
kring_disk(u8_1d uuid, int layer, int k, bool ring_only) {
    if (uuid.shape(0) != 16) throw std::runtime_error("uuid must be a (16,) uint8 array");
    const int64_t cap = hex9_disk_ncells(k);
    if (cap < 0) throw std::runtime_error("k must be >= 0");
    uint8_t *buf = new uint8_t[(size_t)cap * 16];
    int64_t m;
    {
        nb::gil_scoped_release release;
        m = ring_only ? hex9_k_ring(uuid.data(), layer, k, buf, cap)
                      : hex9_k_disk(uuid.data(), layer, k, buf, cap);
    }
    if (m < 0) { delete[] buf; throw std::runtime_error("hex9 k_ring/k_disk: invalid uuid/layer/k"); }
    nb::capsule owner(buf, [](void *p) noexcept { delete[] static_cast<uint8_t *>(p); });
    return nb::ndarray<nb::numpy, uint8_t, nb::ndim<2>>(buf, {(size_t)m, 16}, owner);
}

/* label(uuid[16], layer, key=False) -> str */
static std::string label(u8_1d uuid, int layer, bool key) {
    if (uuid.shape(0) != 16) throw std::runtime_error("uuid must be a (16,) uint8 array");
    char buf[40];
    const int len = key ? hex9_label_key(uuid.data(), layer, buf, sizeof buf)
                        : hex9_label(uuid.data(), layer, buf, sizeof buf);
    if (len < 0) throw std::runtime_error("hex9.label: invalid layer");
    return std::string(buf, (size_t)len);
}

/* parse_label(str) -> (uuid[16], layer) */
static nb::tuple parse_label(const std::string &lbl) {
    uint8_t *u = new uint8_t[16];
    const int layer = hex9_parse_label(lbl.c_str(), u);
    if (layer < 0) { delete[] u; throw std::runtime_error("hex9.parse_label: invalid label"); }
    nb::capsule owner(u, [](void *p) noexcept { delete[] static_cast<uint8_t *>(p); });
    return nb::make_tuple(
        nb::ndarray<nb::numpy, uint8_t, nb::ndim<1>>(u, {16}, owner), layer);
}

/* label_centroid(str) -> (lon, lat) */
static nb::tuple label_centroid(const std::string &lbl) {
    double lon, lat;
    if (hex9_label_centroid(lbl.c_str(), &lon, &lat))
        throw std::runtime_error("hex9.label_centroid: invalid label");
    return nb::make_tuple(lon, lat);
}

/* common_ancestor(uuids[n,16], layer) -> (label, layer, uuid[16]) | None */
static nb::object common_ancestor(u8_2d_in uuids, int layer) {
    if (uuids.shape(1) != 16) throw std::runtime_error("uuids array must be (n, 16)");
    char buf[40];
    uint8_t *u = new uint8_t[16];
    const int al = hex9_common_ancestor(uuids.data(), uuids.shape(0), layer,
                                        buf, sizeof buf, u);
    if (al < 0) { delete[] u; return nb::none(); }
    nb::capsule owner(u, [](void *p) noexcept { delete[] static_cast<uint8_t *>(p); });
    return nb::make_tuple(std::string(buf), al,
                          nb::ndarray<nb::numpy, uint8_t, nb::ndim<1>>(u, {16}, owner));
}

/* adaptive(uuids[n,16], min_layer, max_layer, ceiling, floor, weight=None)
 *   -> (uuids[m,16], layers[m] i32, values[m] f64, npoints[m] i64, assign[n] i64)
 * Input is FULL uuids (from encode) — bin input is rejected: the digest
 * re-bins across layers, guaranteed only from the full uuid. */
static nb::tuple adaptive(u8_2d_in uuids, int min_layer, int max_layer,
                          double ceiling, double floor_,
                          std::optional<f64_1d> weight) {
    if (uuids.shape(1) != 16) throw std::runtime_error("uuids must be (n, 16)");
    const size_t n = uuids.shape(0);
    const double *w = nullptr;
    if (weight) {
        if (weight->shape(0) != n) throw std::runtime_error("weight must match uuids length");
        w = weight->data();
    }
    char err[256] = {0};
    hex9_adaptive *a = nullptr;
    {
        nb::gil_scoped_release release;
        a = hex9_adaptive_create(uuids.data(), w, n,
                                 min_layer, max_layer, ceiling, floor_, err, sizeof err);
    }
    if (!a) throw std::runtime_error(std::string("hex9.adaptive: ") + err);
    const size_t m = (size_t)hex9_adaptive_count(a);
    uint8_t *uu  = new uint8_t[m * 16];
    int32_t *ly  = new int32_t[m];
    double  *val = new double[m];
    int64_t *np  = new int64_t[m];
    int64_t *as  = new int64_t[n];
    for (size_t i = 0; i < m; i++) {
        int layer_i;
        hex9_adaptive_cell(a, (int)i, uu + i * 16, &layer_i, &val[i], &np[i]);
        ly[i] = (int32_t)layer_i;
    }
    hex9_adaptive_assign(a, as);
    hex9_adaptive_destroy(a);
    nb::capsule uu_o(uu, [](void *p) noexcept { delete[] static_cast<uint8_t *>(p); });
    nb::capsule ly_o(ly, [](void *p) noexcept { delete[] static_cast<int32_t *>(p); });
    nb::capsule va_o(val, [](void *p) noexcept { delete[] static_cast<double *>(p); });
    nb::capsule np_o(np, [](void *p) noexcept { delete[] static_cast<int64_t *>(p); });
    nb::capsule as_o(as, [](void *p) noexcept { delete[] static_cast<int64_t *>(p); });
    return nb::make_tuple(
        nb::ndarray<nb::numpy, uint8_t, nb::ndim<2>>(uu, {m, 16}, uu_o),
        nb::ndarray<nb::numpy, int32_t, nb::ndim<1>>(ly, {m}, ly_o),
        nb::ndarray<nb::numpy, double,  nb::ndim<1>>(val, {m}, va_o),
        nb::ndarray<nb::numpy, int64_t, nb::ndim<1>>(np, {m}, np_o),
        nb::ndarray<nb::numpy, int64_t, nb::ndim<1>>(as, {n}, as_o));
}

NB_MODULE(hex9_ext, m) {
    m.doc() = "libhex9 — Hex9 DGGS fast backend (nanobind + OpenMP).";
    m.def("version", &hex9_version);
    m.def("warp_init", &warp_init,
          "Rebuild the authalic-warp state. Done automatically at import; only "
          "needed again after a failed init.");
    m.def("set_use_warp", [](bool on) { hex9_set_use_warp(on ? 1 : 0); }, nb::arg("on"),
          "Toggle the authalic warp. ON by default at import; set False only to "
          "inspect the raw (non-equal-area) lattice.");
    m.def("set_encoder",  [](int mode) { hex9_set_encoder(mode); }, nb::arg("mode"));
    m.def("encode", &encode, nb::arg("lon"), nb::arg("lat"),
          "Encode lon/lat arrays to an (n,16) uint8 array of cell UUIDs.");
    m.def("decode", &decode, nb::arg("uuid"),
          "Decode an (n,16) uint8 UUID array to (lon, lat) arrays.");
    m.def("bin", &bin, nb::arg("uuid"), nb::arg("layer"),
          "Bin an (n,16) UUID array to cell keys at the given layer.");
    m.def("cell", &cell, nb::arg("uuid"), nb::arg("layer"), nb::arg("densify") = 0,
          "Hexagon ring (npoints,2) lon/lat for a UUID at layer/densify.");
    m.def("grid", &grid,
          nb::arg("lon_min"), nb::arg("lat_min"), nb::arg("lon_max"), nb::arg("lat_max"),
          nb::arg("layer"), nb::arg("densify") = 0, nb::arg("max_cells") = 0,
          "Enumerate cells in a lon/lat bbox -> (uuids[n,16], centroids[n,2], rings). "
          "rings is (n,6,2) when densify=0, else a list of (npoints,2) arrays.");
    m.def("neighbors", &neighbors, nb::arg("uuid"), nb::arg("layer"),
          "Edge-adjacent cells of each (n,16) FULL UUID at layer -> "
          "(uint8[n,6,16], int32[n] counts). Counts are 6, or 5 on octahedron-"
          "vertex half-hexes; unused slots are zero. Bin input is rejected "
          "(count -1): bins are layer keys, not addresses.");
    m.def("k_ring",
          [](u8_1d u, int layer, int k) { return kring_disk(u, layer, k, true); },
          nb::arg("uuid"), nb::arg("layer"), nb::arg("k"),
          "Cells at graph distance exactly k from the FULL uuid -> (m,16), "
          "uuid-sorted. Bin input raises: bins are layer keys, not addresses.");
    m.def("k_disk",
          [](u8_1d u, int layer, int k) { return kring_disk(u, layer, k, false); },
          nb::arg("uuid"), nb::arg("layer"), nb::arg("k"),
          "Cells within graph distance k (centre included) of the FULL uuid -> "
          "(m,16), uuid-sorted. Bin input raises: bins are layer keys, not "
          "addresses.");
    m.def("disk_ncells", &hex9_disk_ncells, nb::arg("k"),
          "Nominal k-disk cell count 1+3k(k+1) (upper bound near octahedron vertices).");
    m.def("label", &label, nb::arg("uuid"), nb::arg("layer"), nb::arg("key") = false,
          "Cell label at layer; key=True appends '.<key_tail>'.");
    m.def("parse_label", &parse_label, nb::arg("label"),
          "Label (bare canonical, or keyed '.k' any flavour) -> (canonical bin uuid[16], layer).");
    m.def("label_centroid", &label_centroid, nb::arg("label"),
          "Geographic centroid (lon, lat) of the labelled cell (grid convention).");
    m.def("common_ancestor", &common_ancestor, nb::arg("uuids"), nb::arg("layer"),
          "Deepest common address-tree ancestor of (n,16) cells at layer -> "
          "(label, layer, uuid[16]), or None when cells span L0 hexes.");
    m.def("adaptive", &adaptive,
          nb::arg("uuids"), nb::arg("min_layer"), nb::arg("max_layer"),
          nb::arg("ceiling"), nb::arg("floor"), nb::arg("weight") = nb::none(),
          "Population-digest multi-layer grid over FULL uuids (n,16) from "
          "encode() — bin input rejected (re-binning across layers is only "
          "guaranteed from the full uuid). Bin weighted addresses at max_layer; "
          "cells reaching `floor` emit (first-fit to `ceiling`), excess re-bins to "
          "the parent layer; min_layer absorbs the rest (partial fills allowed). -> "
          "(uuids[m,16], layers[m], values[m], npoints[m], assign[n]).");

    // Enable the authalic warp by default. Without it, encode/decode/cell return
    // the RAW un-warped lattice, which is NOT equal-area (~7% cell-area spread) —
    // the warp is hex9's equal-area mechanism. Raising here (rather than silently
    // serving raw cells) surfaces a missing/corrupt warp blob at import time.
    // Callers wanting the raw grid can still hex9_ext.set_use_warp(False).
    warp_init();
    hex9_set_use_warp(1);
}
