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
#include <cstring>   /* std::memcpy — manylinux gcc does not provide it transitively */
#include <stdexcept>
#include <string>

namespace nb = nanobind;

using f64_1d  = nb::ndarray<const double, nb::ndim<1>, nb::c_contig, nb::device::cpu>;
using i32_1d  = nb::ndarray<const int32_t, nb::ndim<1>, nb::c_contig, nb::device::cpu>;
using u8_2d_in = nb::ndarray<const uint8_t, nb::ndim<2>, nb::c_contig, nb::device::cpu>;

/* Initialise the warp once; raises on failure. */
static void warp_init() {   /* C-side canonical name is hex9_init since 2.1.0 */
    char err[256] = {0};
    if (hex9_init(err, sizeof err))
        throw std::runtime_error(std::string("hex9 warp init failed: ") + err);
}

/* encode(lon[n], lat[n]) -> uint8[n,16] */
static nb::ndarray<nb::numpy, uint8_t, nb::ndim<2>>
encode(f64_1d lon, f64_1d lat, bool sphere) {
    const size_t n = lon.shape(0);
    if (lat.shape(0) != n) throw std::runtime_error("lon and lat must be the same length");
    uint8_t *out = new uint8_t[n * 16];
    {
        nb::gil_scoped_release release;
        if (sphere) hex9_encode_many_sphere(lon.data(), lat.data(), n, out);
        else        hex9_encode_many(lon.data(), lat.data(), n, out);
    }
    nb::capsule owner(out, [](void *p) noexcept { delete[] static_cast<uint8_t *>(p); });
    return nb::ndarray<nb::numpy, uint8_t, nb::ndim<2>>(out, {n, 16}, owner);
}

/* decode(uuid[n,16]) -> (lon[n], lat[n]) */
static nb::tuple decode(u8_2d_in uuid, bool sphere) {
    const size_t n = uuid.shape(0);
    if (uuid.shape(1) != 16) throw std::runtime_error("uuid array must be (n, 16)");
    double *lon = new double[n];
    double *lat = new double[n];
    int rc;
    {
        nb::gil_scoped_release release;
        if (sphere) rc = hex9_decode_many_sphere(uuid.data(), n, lon, lat);
        else        rc = hex9_decode_many(uuid.data(), n, lon, lat);
    }
    if (rc) { delete[] lon; delete[] lat;
              throw std::runtime_error("hex9.decode: curve or E4H uuid input — "
                                       "use curve_decode / e4h_decode"); }
    nb::capsule lon_o(lon, [](void *p) noexcept { delete[] static_cast<double *>(p); });
    nb::capsule lat_o(lat, [](void *p) noexcept { delete[] static_cast<double *>(p); });
    return nb::make_tuple(
        nb::ndarray<nb::numpy, double, nb::ndim<1>>(lon, {n}, lon_o),
        nb::ndarray<nb::numpy, double, nb::ndim<1>>(lat, {n}, lat_o));
}

/* project(lon[n], lat[n]) -> (cx[n], cy[n], oid[n]) — the b_oct backend
 * surface (continuous, no descent; see hex9_c.h hex9_project). */
static nb::tuple project(f64_1d lon, f64_1d lat, bool sphere) {
    const size_t n = lon.shape(0);
    if (lat.shape(0) != n) throw std::runtime_error("lon and lat must be the same length");
    double *cx = new double[n];
    double *cy = new double[n];
    int32_t *oid = new int32_t[n];
    {
        nb::gil_scoped_release release;
        static_assert(sizeof(int) == sizeof(int32_t), "oid width");
        if (sphere) hex9_project_many_sphere(lon.data(), lat.data(), n,
                                             cx, cy, (int *)oid);
        else        hex9_project_many(lon.data(), lat.data(), n,
                                      cx, cy, (int *)oid);
    }
    nb::capsule cx_o(cx, [](void *p) noexcept { delete[] static_cast<double *>(p); });
    nb::capsule cy_o(cy, [](void *p) noexcept { delete[] static_cast<double *>(p); });
    nb::capsule oid_o(oid, [](void *p) noexcept { delete[] static_cast<int32_t *>(p); });
    return nb::make_tuple(
        nb::ndarray<nb::numpy, double, nb::ndim<1>>(cx, {n}, cx_o),
        nb::ndarray<nb::numpy, double, nb::ndim<1>>(cy, {n}, cy_o),
        nb::ndarray<nb::numpy, int32_t, nb::ndim<1>>(oid, {n}, oid_o));
}

/* unproject(cx[n], cy[n], oid[n]) -> (lon[n], lat[n]) — exact inverse. */
static nb::tuple unproject(f64_1d cx, f64_1d cy, i32_1d oid, bool sphere) {
    const size_t n = cx.shape(0);
    if (cy.shape(0) != n || oid.shape(0) != n)
        throw std::runtime_error("cx, cy and oid must be the same length");
    double *lon = new double[n];
    double *lat = new double[n];
    {
        nb::gil_scoped_release release;
        if (sphere) hex9_unproject_many_sphere(cx.data(), cy.data(),
                                               (const int *)oid.data(), n, lon, lat);
        else        hex9_unproject_many(cx.data(), cy.data(),
                                        (const int *)oid.data(), n, lon, lat);
    }
    nb::capsule lon_o(lon, [](void *p) noexcept { delete[] static_cast<double *>(p); });
    nb::capsule lat_o(lat, [](void *p) noexcept { delete[] static_cast<double *>(p); });
    return nb::make_tuple(
        nb::ndarray<nb::numpy, double, nb::ndim<1>>(lon, {n}, lon_o),
        nb::ndarray<nb::numpy, double, nb::ndim<1>>(lat, {n}, lat_o));
}

/* encode_boct(cx[n], cy[n], oid[n]) -> uint8[n,16] — mint addresses from
 * face coordinates directly (bring-your-own-projection; see hex9_c.h
 * doctrine block: the resulting addresses are NON-CANONICAL unless the
 * inputs came from hex9's own projection). */
static nb::ndarray<nb::numpy, uint8_t, nb::ndim<2>>
encode_boct(f64_1d cx, f64_1d cy, i32_1d oid) {
    const size_t n = cx.shape(0);
    if (cy.shape(0) != n || oid.shape(0) != n)
        throw std::runtime_error("cx, cy and oid must be the same length");
    uint8_t *out = new uint8_t[n * 16];
    int rc;
    {
        nb::gil_scoped_release release;
        rc = hex9_encode_boct_many(cx.data(), cy.data(),
                                   (const int *)oid.data(), n, out);
    }
    if (rc) { delete[] out; throw std::runtime_error("encode_boct: invalid input (oid range / non-finite)"); }
    nb::capsule owner(out, [](void *p) noexcept { delete[] static_cast<uint8_t *>(p); });
    return nb::ndarray<nb::numpy, uint8_t, nb::ndim<2>>(out, {n, 16}, owner);
}

/* decode_boct(uuid[n,16]) -> (cx[n], cy[n], oid[n]) — cell centroids in
 * face coordinates (projection-free). */
static nb::tuple decode_boct(u8_2d_in uuid) {
    const size_t n = uuid.shape(0);
    if (uuid.shape(1) != 16) throw std::runtime_error("uuid array must be (n, 16)");
    double *cx = new double[n];
    double *cy = new double[n];
    int32_t *oid = new int32_t[n];
    int rc;
    {
        nb::gil_scoped_release release;
        rc = hex9_decode_boct_many(uuid.data(), n, cx, cy, (int *)oid);
    }
    if (rc) { delete[] cx; delete[] cy; delete[] oid;
              throw std::runtime_error("decode_boct: malformed uuid"); }
    nb::capsule cx_o(cx, [](void *p) noexcept { delete[] static_cast<double *>(p); });
    nb::capsule cy_o(cy, [](void *p) noexcept { delete[] static_cast<double *>(p); });
    nb::capsule oid_o(oid, [](void *p) noexcept { delete[] static_cast<int32_t *>(p); });
    return nb::make_tuple(
        nb::ndarray<nb::numpy, double, nb::ndim<1>>(cx, {n}, cx_o),
        nb::ndarray<nb::numpy, double, nb::ndim<1>>(cy, {n}, cy_o),
        nb::ndarray<nb::numpy, int32_t, nb::ndim<1>>(oid, {n}, oid_o));
}

/* cell_ring_boct(uuid[16], layer, densify) -> (cx[k], cy[k], oid[k]) —
 * one cell's closed boundary ring in face coordinates, per-vertex frames. */
static nb::tuple cell_ring_boct(u8_2d_in uuid, int layer, int densify) {
    if (uuid.shape(0) != 1 || uuid.shape(1) != 16)
        throw std::runtime_error("uuid array must be (1, 16) — one cell per call");
    const int k = hex9_ring_npoints(densify);
    if (k < 0) throw std::runtime_error("densify out of range");
    double *cx = new double[(size_t)k];
    double *cy = new double[(size_t)k];
    int32_t *oid = new int32_t[(size_t)k];
    int got;
    {
        nb::gil_scoped_release release;
        got = hex9_cell_ring_boct(uuid.data(), layer, densify,
                                  cx, cy, (int *)oid, k);
    }
    if (got != k) { delete[] cx; delete[] cy; delete[] oid;
                    throw std::runtime_error("cell_ring_boct: bad layer/uuid"); }
    nb::capsule cx_o(cx, [](void *p) noexcept { delete[] static_cast<double *>(p); });
    nb::capsule cy_o(cy, [](void *p) noexcept { delete[] static_cast<double *>(p); });
    nb::capsule oid_o(oid, [](void *p) noexcept { delete[] static_cast<int32_t *>(p); });
    return nb::make_tuple(
        nb::ndarray<nb::numpy, double, nb::ndim<1>>(cx, {(size_t)k}, cx_o),
        nb::ndarray<nb::numpy, double, nb::ndim<1>>(cy, {(size_t)k}, cy_o),
        nb::ndarray<nb::numpy, int32_t, nb::ndim<1>>(oid, {(size_t)k}, oid_o));
}

/* cell_uv(uuid[n,16], layers[n]) -> (c_ia, c_ib, c_oid, v_ia[n,6],
 * v_ib[n,6], v_oid[n,6], ext) — the integer lattice identity surface.
 * Datum-free; vertex keys are canonical pool keys (see hex9_c.h). */
static nb::tuple cell_uv(u8_2d_in uuid, i32_1d layers) {
    const size_t n = uuid.shape(0);
    if (uuid.shape(1) != 16) throw std::runtime_error("uuid array must be (n, 16)");
    if (layers.shape(0) != n) throw std::runtime_error("layers must match uuid rows");
    int64_t *cia = new int64_t[n];
    int64_t *cib = new int64_t[n];
    int32_t *coid = new int32_t[n];
    int64_t *via = new int64_t[n * 6];
    int64_t *vib = new int64_t[n * 6];
    int32_t *void_ = new int32_t[n * 6];
    int32_t *ext = new int32_t[n];
    int rc;
    {
        nb::gil_scoped_release release;
        rc = hex9_cell_uv_many(uuid.data(), layers.data(), n,
                               cia, cib, coid, via, vib, void_, ext);
    }
    if (rc) {
        delete[] cia; delete[] cib; delete[] coid;
        delete[] via; delete[] vib; delete[] void_; delete[] ext;
        throw std::runtime_error("hex9.cell_uv: malformed uuid/layer in batch");
    }
    auto own64 = [](int64_t *p) {
        return nb::capsule(p, [](void *q) noexcept { delete[] static_cast<int64_t *>(q); });
    };
    auto own32 = [](int32_t *p) {
        return nb::capsule(p, [](void *q) noexcept { delete[] static_cast<int32_t *>(q); });
    };
    return nb::make_tuple(
        nb::ndarray<nb::numpy, int64_t, nb::ndim<1>>(cia, {n}, own64(cia)),
        nb::ndarray<nb::numpy, int64_t, nb::ndim<1>>(cib, {n}, own64(cib)),
        nb::ndarray<nb::numpy, int32_t, nb::ndim<1>>(coid, {n}, own32(coid)),
        nb::ndarray<nb::numpy, int64_t, nb::ndim<2>>(via, {n, 6}, own64(via)),
        nb::ndarray<nb::numpy, int64_t, nb::ndim<2>>(vib, {n, 6}, own64(vib)),
        nb::ndarray<nb::numpy, int32_t, nb::ndim<2>>(void_, {n, 6}, own32(void_)),
        nb::ndarray<nb::numpy, int32_t, nb::ndim<1>>(ext, {n}, own32(ext)));
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

/* cell_parent(uuid[n,16]) -> uint8[n,16] — canonical one-level cell parent */
static nb::ndarray<nb::numpy, uint8_t, nb::ndim<2>>
cell_parent(u8_2d_in uuid) {
    const size_t n = uuid.shape(0);
    if (uuid.shape(1) != 16) throw std::runtime_error("uuid array must be (n, 16)");
    uint8_t *out = new uint8_t[n * 16];
    int rc;
    {
        nb::gil_scoped_release release;
        rc = hex9_cell_parent_many(uuid.data(), n, out);
    }
    if (rc) { delete[] out; throw std::runtime_error("hex9.cell_parent: L0 cells have no parent (or malformed uuid)"); }
    nb::capsule owner(out, [](void *p) noexcept { delete[] static_cast<uint8_t *>(p); });
    return nb::ndarray<nb::numpy, uint8_t, nb::ndim<2>>(out, {n, 16}, owner);
}

/* cell_ancestor(uuid[n,16], layer) -> uint8[n,16] — the OWNERSHIP relation:
 * direct deep fold, NOT iterated cell_parent (they differ on hexagon-band
 * cells beyond one generation). */
static nb::ndarray<nb::numpy, uint8_t, nb::ndim<2>>
cell_ancestor(u8_2d_in uuid, int layer) {
    const size_t n = uuid.shape(0);
    if (uuid.shape(1) != 16) throw std::runtime_error("uuid array must be (n, 16)");
    uint8_t *out = new uint8_t[n * 16];
    int rc;
    {
        nb::gil_scoped_release release;
        rc = hex9_cell_ancestor_many(uuid.data(), layer, n, out);
    }
    if (rc) { delete[] out; throw std::runtime_error("hex9.cell_ancestor: input coarser than target layer (or layer/uuid invalid)"); }
    nb::capsule owner(out, [](void *p) noexcept { delete[] static_cast<uint8_t *>(p); });
    return nb::ndarray<nb::numpy, uint8_t, nb::ndim<2>>(out, {n, 16}, owner);
}

/* curve(uuid[n,16]) -> uint8[n,16] — packed curve-uuids (0xC marker). */
static nb::ndarray<nb::numpy, uint8_t, nb::ndim<2>>
curve(u8_2d_in uuid) {
    const size_t n = uuid.shape(0);
    if (uuid.shape(1) != 16) throw std::runtime_error("uuid array must be (n, 16)");
    uint8_t *out = new uint8_t[n * 16];
    int rc;
    {
        nb::gil_scoped_release release;
        rc = hex9_curve_many(uuid.data(), n, out);
    }
    if (rc) { delete[] out; throw std::runtime_error("hex9.curve: malformed or non-canonical uuid in batch"); }
    nb::capsule owner(out, [](void *p) noexcept { delete[] static_cast<uint8_t *>(p); });
    return nb::ndarray<nb::numpy, uint8_t, nb::ndim<2>>(out, {n, 16}, owner);
}

/* curve_decode(curve[n,16]) -> uint8[n,16] — canonical bins (h9-uuids pass). */
static nb::ndarray<nb::numpy, uint8_t, nb::ndim<2>>
curve_decode(u8_2d_in cu) {
    const size_t n = cu.shape(0);
    if (cu.shape(1) != 16) throw std::runtime_error("uuid array must be (n, 16)");
    uint8_t *out = new uint8_t[n * 16];
    int rc;
    {
        nb::gil_scoped_release release;
        rc = hex9_curve_decode_many(cu.data(), n, out);
    }
    if (rc) { delete[] out; throw std::runtime_error("hex9.curve_decode: malformed curve uuid in batch"); }
    nb::capsule owner(out, [](void *p) noexcept { delete[] static_cast<uint8_t *>(p); });
    return nb::ndarray<nb::numpy, uint8_t, nb::ndim<2>>(out, {n, 16}, owner);
}

/* cell_children(uuid[n,16]) -> uint8[n,9,16] — curve-rank order. */
static nb::ndarray<nb::numpy, uint8_t, nb::ndim<3>>
cell_children(u8_2d_in uuid) {
    const size_t n = uuid.shape(0);
    if (uuid.shape(1) != 16) throw std::runtime_error("uuid array must be (n, 16)");
    uint8_t *out = new uint8_t[n * 9 * 16];
    int rc = 0;
    {
        nb::gil_scoped_release release;
        for (size_t i = 0; i < n; ++i)
            rc |= hex9_cell_children(uuid.data() + i * 16, out + i * 9 * 16);
    }
    if (rc) { delete[] out; throw std::runtime_error("hex9.cell_children: malformed uuid or cell at lmax"); }
    nb::capsule owner(out, [](void *p) noexcept { delete[] static_cast<uint8_t *>(p); });
    return nb::ndarray<nb::numpy, uint8_t, nb::ndim<3>>(out, {n, 9, 16}, owner);
}

/* Shared body for the two descendant enumerators: one zone row in,
 * (bins[m,16], curves[m,16]) out. */
static nb::tuple cells_of(const uint8_t *zone, int layer,
                          int64_t (*fn)(const uint8_t *, int, uint8_t *,
                                        uint8_t *, int64_t),
                          const char *name) {
    uint8_t cu[16];
    if (hex9_curve(zone, cu) != 0)
        throw std::runtime_error(std::string(name) + ": not a valid H9 or curve uuid");
    const int from = hex9_curve_layer(cu);
    const int64_t want = hex9_curve_ncells(from, layer);
    if (want < 0)
        throw std::runtime_error(std::string(name) + ": bad layer (own layer "
                                 + std::to_string(from) + ")");
    uint8_t *bins = new uint8_t[(size_t)want * 16];
    uint8_t *curves = new uint8_t[(size_t)want * 16];
    int64_t m;
    {
        nb::gil_scoped_release release;
        m = fn(zone, layer, bins, curves, want);
    }
    if (m != want) {
        delete[] bins; delete[] curves;
        throw std::runtime_error(std::string(name) + ": enumeration failed");
    }
    nb::capsule bo(bins, [](void *p) noexcept { delete[] static_cast<uint8_t *>(p); });
    nb::capsule co(curves, [](void *p) noexcept { delete[] static_cast<uint8_t *>(p); });
    return nb::make_tuple(
        nb::ndarray<nb::numpy, uint8_t, nb::ndim<2>>(bins, {(size_t)m, 16}, bo),
        nb::ndarray<nb::numpy, uint8_t, nb::ndim<2>>(curves, {(size_t)m, 16}, co));
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
cell(u8_1d uuid, int layer, int densify, bool sphere) {
    if (uuid.shape(0) != 16) throw std::runtime_error("uuid must be a (16,) uint8 array");
    int P = hex9_ring_npoints(densify);
    if (P < 0) throw std::runtime_error("hex9.cell: densify out of range (0..9)");
    double *r = new double[(size_t)P * 2];
    int got = sphere ? hex9_cell_ring_sphere(uuid.data(), layer, densify, r, P)
                     : hex9_cell_ring(uuid.data(), layer, densify, r, P);
    if (got < 0) { delete[] r; throw std::runtime_error("hex9.cell: invalid layer/densify"); }
    return own_f64<2>(r, {(size_t)P, 2});
}

/* grid(lon_min,lat_min,lon_max,lat_max, layer, densify=0, max_cells=0)
 *   -> (uuids[n,16], centroids[n,2], rings)
 * rings: densify==0 -> flat (n,6,2) of the 6 hex corners;
 *        densify>0  -> list of n (npoints,2) closed-ring arrays. */
static nb::tuple grid(double lon_min, double lat_min, double lon_max, double lat_max,
                      int layer, int densify, int64_t max_cells, bool sphere) {
    char err[256] = {0};
    hex9_grid *g = nullptr;
    {
        nb::gil_scoped_release release;
        g = sphere ? hex9_grid_create_sphere(lon_min, lat_min, lon_max, lat_max,
                                             layer, densify, max_cells, err, sizeof err)
                   : hex9_grid_create(lon_min, lat_min, lon_max, lat_max,
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
    if ((uuid.data()[0] >> 4) == 0xC)
        throw std::runtime_error("hex9.label: curve-uuid input (nibble 0 = 0xC) — "
                                 "cell labels of curve addresses are meaningless; "
                                 "use curve_label, or curve_decode first");
    char buf[40];
    const int len = key ? hex9_label_key(uuid.data(), layer, buf, sizeof buf)
                        : hex9_label(uuid.data(), layer, buf, sizeof buf);
    if (len < 0) throw std::runtime_error("hex9.label: invalid layer");
    return std::string(buf, (size_t)len);
}

/* label(uuids[n,16], layer, key=False) -> list[str] — batch overload */
static nb::list label_many(u8_2d_in uuids, int layer, bool key) {
    if (uuids.shape(1) != 16) throw std::runtime_error("uuids array must be (n, 16)");
    const size_t n = uuids.shape(0);
    nb::list out;
    char buf[40];
    for (size_t i = 0; i < n; ++i) {
        const uint8_t *u = uuids.data() + i * 16;
        if ((u[0] >> 4) == 0xC)
            throw std::runtime_error("hex9.label: curve-uuid input at row " +
                                     std::to_string(i) + " (nibble 0 = 0xC) — "
                                     "use curve_label, or curve_decode first");
        const int len = key ? hex9_label_key(u, layer, buf, sizeof buf)
                            : hex9_label(u, layer, buf, sizeof buf);
        if (len < 0) throw std::runtime_error("hex9.label: invalid layer");
        out.append(nb::str(buf, (size_t)len));
    }
    return out;
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

/* curve_label(uuid[16]) -> str */
static std::string curve_label(u8_1d uuid) {
    if (uuid.shape(0) != 16) throw std::runtime_error("uuid must be a (16,) uint8 array");
    char buf[40];
    const int len = hex9_curve_label(uuid.data(), buf, sizeof buf);
    if (len < 0) throw std::runtime_error("hex9.curve_label: invalid uuid");
    return std::string(buf, (size_t)len);
}

/* curve_label(uuids[n,16]) -> list[str] — batch overload */
static nb::list curve_label_many(u8_2d_in uuids) {
    if (uuids.shape(1) != 16) throw std::runtime_error("uuids array must be (n, 16)");
    const size_t n = uuids.shape(0);
    nb::list out;
    char buf[40];
    for (size_t i = 0; i < n; ++i) {
        const int len = hex9_curve_label(uuids.data() + i * 16, buf, sizeof buf);
        if (len < 0) throw std::runtime_error("hex9.curve_label: invalid uuid at row " +
                                              std::to_string(i));
        out.append(nb::str(buf, (size_t)len));
    }
    return out;
}

/* curve_parse_label(str) -> curve uuid[16] */
static nb::ndarray<nb::numpy, uint8_t, nb::ndim<1>>
curve_parse_label(const std::string &lbl) {
    uint8_t *u = new uint8_t[16];
    if (hex9_curve_parse_label(lbl.c_str(), u)) {
        delete[] u;
        throw std::runtime_error("hex9.curve_parse_label: invalid curve label");
    }
    nb::capsule owner(u, [](void *p) noexcept { delete[] static_cast<uint8_t *>(p); });
    return nb::ndarray<nb::numpy, uint8_t, nb::ndim<1>>(u, {16}, owner);
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

/* ── E4H: the aperture-4 structural tail (see hex9_c.h doctrine block) ── */

static nb::ndarray<nb::numpy, uint8_t, nb::ndim<2>>
e4h_encode(f64_1d lon, f64_1d lat, int layer, int depth, bool sphere) {
    const size_t n = lon.shape(0);
    if (lat.shape(0) != n) throw std::runtime_error("lon and lat must be the same length");
    uint8_t *out = new uint8_t[n * 16];
    int rc;
    {
        nb::gil_scoped_release release;
        rc = sphere
            ? hex9_e4h_encode_many_sphere(lon.data(), lat.data(), n, layer, depth, out)
            : hex9_e4h_encode_many(lon.data(), lat.data(), n, layer, depth, out);
    }
    if (rc) { delete[] out;
              throw std::runtime_error("e4h_encode: invalid layer/depth "
                                       "(need layer + 2 + depth <= lmax())"); }
    nb::capsule owner(out, [](void *p) noexcept { delete[] static_cast<uint8_t *>(p); });
    return nb::ndarray<nb::numpy, uint8_t, nb::ndim<2>>(out, {n, 16}, owner);
}

static nb::tuple e4h_decode_impl(u8_2d_in uuid, bool sphere, bool partner) {
    const size_t n = uuid.shape(0);
    if (uuid.shape(1) != 16) throw std::runtime_error("uuid array must be (n, 16)");
    double *lon = new double[n];
    double *lat = new double[n];
    int rc = 0;
    {
        nb::gil_scoped_release release;
        if (partner) {
            for (size_t i = 0; i < n; ++i)
                rc |= sphere
                    ? hex9_e4h_partner_sphere(uuid.data() + i * 16, &lon[i], &lat[i])
                    : hex9_e4h_partner(uuid.data() + i * 16, &lon[i], &lat[i]);
        } else {
            rc = sphere ? hex9_e4h_decode_many_sphere(uuid.data(), n, lon, lat)
                        : hex9_e4h_decode_many(uuid.data(), n, lon, lat);
        }
    }
    if (rc) { delete[] lon; delete[] lat;
              throw std::runtime_error("e4h_decode: not an E4H uuid / bad grammar"); }
    nb::capsule lon_o(lon, [](void *p) noexcept { delete[] static_cast<double *>(p); });
    nb::capsule lat_o(lat, [](void *p) noexcept { delete[] static_cast<double *>(p); });
    return nb::make_tuple(
        nb::ndarray<nb::numpy, double, nb::ndim<1>>(lon, {n}, lon_o),
        nb::ndarray<nb::numpy, double, nb::ndim<1>>(lat, {n}, lat_o));
}

/* e4h_split(uuid[n,16]) -> (host[n,16], half[n], digits[n,28] int8, -1 pad) */
static nb::tuple e4h_split(u8_2d_in uuid) {
    const size_t n = uuid.shape(0);
    if (uuid.shape(1) != 16) throw std::runtime_error("uuid array must be (n, 16)");
    uint8_t *host = new uint8_t[n * 16];
    int32_t *half = new int32_t[n];
    int8_t  *dig  = new int8_t[n * 28];
    for (size_t i = 0; i < n; ++i) {
        uint8_t d[28];
        int hf, nd;
        if (hex9_e4h_split(uuid.data() + i * 16, host + i * 16, &hf, d, &nd)) {
            delete[] host; delete[] half; delete[] dig;
            throw std::runtime_error("e4h_split: not an E4H uuid at row "
                                     + std::to_string(i));
        }
        half[i] = hf;
        for (int j = 0; j < 28; ++j)
            dig[i * 28 + j] = j < nd ? (int8_t)d[j] : (int8_t)-1;
    }
    nb::capsule h_o(host, [](void *p) noexcept { delete[] static_cast<uint8_t *>(p); });
    nb::capsule f_o(half, [](void *p) noexcept { delete[] static_cast<int32_t *>(p); });
    nb::capsule d_o(dig,  [](void *p) noexcept { delete[] static_cast<int8_t *>(p); });
    return nb::make_tuple(
        nb::ndarray<nb::numpy, uint8_t, nb::ndim<2>>(host, {n, 16}, h_o),
        nb::ndarray<nb::numpy, int32_t, nb::ndim<1>>(half, {n}, f_o),
        nb::ndarray<nb::numpy, int8_t,  nb::ndim<2>>(dig, {n, 28}, d_o));
}

static nb::ndarray<nb::numpy, uint8_t, nb::ndim<2>>
e4h_bin(u8_2d_in uuid, int depth) {
    const size_t n = uuid.shape(0);
    if (uuid.shape(1) != 16) throw std::runtime_error("uuid array must be (n, 16)");
    uint8_t *out = new uint8_t[n * 16];
    for (size_t i = 0; i < n; ++i)
        if (hex9_e4h_bin(uuid.data() + i * 16, depth, out + i * 16)) {
            delete[] out;
            throw std::runtime_error("e4h_bin: bad depth or not E4H at row "
                                     + std::to_string(i));
        }
    nb::capsule owner(out, [](void *p) noexcept { delete[] static_cast<uint8_t *>(p); });
    return nb::ndarray<nb::numpy, uint8_t, nb::ndim<2>>(out, {n, 16}, owner);
}

static nb::ndarray<nb::numpy, uint8_t, nb::ndim<2>> e4h_hex(u8_2d_in uuid) {
    const size_t n = uuid.shape(0);
    if (uuid.shape(1) != 16) throw std::runtime_error("uuid array must be (n, 16)");
    uint8_t *out = new uint8_t[n * 16];
    int rc;
    {
        nb::gil_scoped_release release;
        rc = hex9_e4h_hex_many(uuid.data(), n, out);
    }
    if (rc) { delete[] out;
              throw std::runtime_error("e4h_hex: not an E4H uuid / bad grammar"); }
    nb::capsule owner(out, [](void *p) noexcept { delete[] static_cast<uint8_t *>(p); });
    return nb::ndarray<nb::numpy, uint8_t, nb::ndim<2>>(out, {n, 16}, owner);
}

static nb::ndarray<nb::numpy, int32_t, nb::ndim<1>> e4h_mode(u8_2d_in uuid) {
    const size_t n = uuid.shape(0);
    if (uuid.shape(1) != 16) throw std::runtime_error("uuid array must be (n, 16)");
    int32_t *out = new int32_t[n];
    for (size_t i = 0; i < n; ++i)
        out[i] = hex9_e4h_mode(uuid.data() + i * 16);
    nb::capsule owner(out, [](void *p) noexcept { delete[] static_cast<int32_t *>(p); });
    return nb::ndarray<nb::numpy, int32_t, nb::ndim<1>>(out, {n}, owner);
}

static nb::ndarray<nb::numpy, int32_t, nb::ndim<1>> e4h_depth(u8_2d_in uuid) {
    const size_t n = uuid.shape(0);
    if (uuid.shape(1) != 16) throw std::runtime_error("uuid array must be (n, 16)");
    int32_t *out = new int32_t[n];
    for (size_t i = 0; i < n; ++i)
        out[i] = hex9_e4h_depth(uuid.data() + i * 16);
    nb::capsule owner(out, [](void *p) noexcept { delete[] static_cast<int32_t *>(p); });
    return nb::ndarray<nb::numpy, int32_t, nb::ndim<1>>(out, {n}, owner);
}

static std::string e4h_label_one(u8_1d uuid) {
    if (uuid.shape(0) != 16) throw std::runtime_error("uuid must be a (16,) uint8 array");
    char buf[64];
    const int len = hex9_e4h_label(uuid.data(), buf, sizeof buf);
    if (len < 0) throw std::runtime_error("e4h_label: not an E4H uuid");
    return std::string(buf, (size_t)len);
}

static nb::list e4h_label_many(u8_2d_in uuids) {
    if (uuids.shape(1) != 16) throw std::runtime_error("uuids array must be (n, 16)");
    nb::list out;
    char buf[64];
    for (size_t i = 0; i < uuids.shape(0); ++i) {
        const int len = hex9_e4h_label(uuids.data() + i * 16, buf, sizeof buf);
        if (len < 0) throw std::runtime_error("e4h_label: not an E4H uuid at row "
                                              + std::to_string(i));
        out.append(nb::str(buf, (size_t)len));
    }
    return out;
}

static nb::ndarray<nb::numpy, uint8_t, nb::ndim<1>>
e4h_parse_label(const std::string &label) {
    uint8_t *out = new uint8_t[16];
    if (hex9_e4h_parse_label(label.c_str(), out)) {
        delete[] out;
        throw std::runtime_error("e4h_parse_label: bad label '" + label + "'");
    }
    nb::capsule owner(out, [](void *p) noexcept { delete[] static_cast<uint8_t *>(p); });
    return nb::ndarray<nb::numpy, uint8_t, nb::ndim<1>>(out, {16}, owner);
}

static bool is_e4h_one(u8_1d uuid) {
    if (uuid.shape(0) != 16) throw std::runtime_error("uuid must be a (16,) uint8 array");
    return hex9_is_e4h(uuid.data()) != 0;
}

static nb::ndarray<nb::numpy, bool, nb::ndim<1>> is_e4h_many(u8_2d_in uuid) {
    const size_t n = uuid.shape(0);
    if (uuid.shape(1) != 16) throw std::runtime_error("uuid array must be (n, 16)");
    bool *out = new bool[n];
    for (size_t i = 0; i < n; ++i)
        out[i] = hex9_is_e4h(uuid.data() + i * 16) != 0;
    nb::capsule owner(out, [](void *p) noexcept { delete[] static_cast<bool *>(p); });
    return nb::ndarray<nb::numpy, bool, nb::ndim<1>>(out, {n}, owner);
}

/* Module name is compile-time selectable: the in-tree module stays hex9_ext
 * (tests, embedders); the PyPI wheel builds this same TU as hex9._core
 * (-DHEX9_PY_MODULE=_core). Double-expansion so the macro argument resolves
 * before NB_MODULE's token pasting. */
#ifndef HEX9_PY_MODULE
#define HEX9_PY_MODULE hex9_ext
#endif
#define H9_NB_MODULE_(name, variable) NB_MODULE(name, variable)
H9_NB_MODULE_(HEX9_PY_MODULE, m) {
    m.doc() = "libhex9 — Hex9 DGGS fast backend (nanobind + OpenMP).";
    m.def("version", &hex9_version);
    m.def("lmax", &hex9_lmax, "Deepest addressable layer (29 legacy / 30 reclaimed).");
    m.def("init", &warp_init,
          "Build the addressing chain (authalic series + warp field). Done "
          "automatically at import; only needed again after a failed init. "
          "Canonical name since 2.1.0.");
    m.def("warp_init", &warp_init,
          "Historic alias of init() (pre-2.1.0 name).");
    m.def("set_encoder",  [](int mode) { hex9_set_encoder(mode); }, nb::arg("mode"));
    m.def("encode", &encode, nb::arg("lon"), nb::arg("lat"),
          nb::arg("sphere") = false,
          "Encode lon/lat arrays to an (n,16) uint8 array of cell UUIDs. "
          "sphere=True: lon/lat are already-spherical degrees (no WGS84 "
          "authalic reduction) — the datum is dataset metadata; never mix "
          "within one dataset.");
    m.def("decode", &decode, nb::arg("uuid"), nb::arg("sphere") = false,
          "Decode an (n,16) uint8 UUID array to (lon, lat) arrays. "
          "sphere=True emits spherical degrees (see encode).");
    m.def("project", &project, nb::arg("lon"), nb::arg("lat"),
          nb::arg("sphere") = false,
          "Continuous projection: lon/lat arrays -> (cx, cy, oid) in the "
          "b_oct frame (no descent — the backend surface for embedders and "
          "renderers). sphere=True: inputs are already-spherical degrees.");
    m.def("unproject", &unproject, nb::arg("cx"), nb::arg("cy"), nb::arg("oid"),
          nb::arg("sphere") = false,
          "Exact inverse of project: (cx, cy, oid) arrays -> (lon, lat). "
          "sphere=True emits spherical degrees.");
    m.def("encode_boct", &encode_boct, nb::arg("cx"), nb::arg("cy"), nb::arg("oid"),
          "Mint addresses from face coordinates (cx, cy, oid) directly — the "
          "bring-your-own-projection entry (the (face)(cell) layers without "
          "AKW). Addresses minted from a non-hex9 projection are NON-CANONICAL: "
          "same uuid space, different meaning. Projection identity is dataset "
          "metadata; never mix projections within one dataset.");
    m.def("decode_boct", &decode_boct, nb::arg("uuid"),
          "Cell centroids in face coordinates (cx, cy, oid) — the same "
          "centroid decode() unprojects, emitted before unprojection. "
          "Projection-free.");
    m.def("cell_ring_boct", &cell_ring_boct, nb::arg("uuid"),
          nb::arg("layer"), nb::arg("densify") = 0,
          "One cell's closed boundary ring in face coordinates: (cx, cy, oid) "
          "arrays, per-vertex frames (each vertex in the chart of ITS oid). "
          "Projection-free — render through your own inverse.");
    m.def("cell_uv", &cell_uv, nb::arg("uuid"), nb::arg("layers"),
          "Integer lattice identity of each (n,16) bin at layers[n]: "
          "(c_ia, c_ib, c_oid, v_ia[n,6], v_ib[n,6], v_oid[n,6], ext). "
          "Exact integers at scale 3^layer (b_oct via uv_units); vertex keys "
          "are canonical pool keys for shared-vertex meshes. Datum-free — "
          "upstream of b_oct and lon/lat.");
    m.def("uv_units", []() {
        double u1, v3;
        hex9_uv_units(&u1, &v3);
        return nb::make_tuple(u1, v3);
    }, "The lattice unit lengths (u1, v3): cx = ia*u1/3^L, cy = ib*v3/3^L.");
    m.def("bin", &bin, nb::arg("uuid"), nb::arg("layer"),
          "Bin an (n,16) UUID array to cell keys at the given layer.");
    m.def("cell_parent", &cell_parent, nb::arg("uuid"),
          "Canonical one-level CELL parent of each (n,16) layer-L bin (L>=1): "
          "the layer-(L-1) cell containing the cell's mode-0 d_cell. Every "
          "parent has exactly 9 canonical children. Distinct from bin(), which "
          "answers the point question from a FULL uuid.");
    m.def("cell_ancestor", &cell_ancestor, nb::arg("uuid"), nb::arg("layer"),
          "Canonical layer-`layer` CELL ancestor of each (n,16) bin — the "
          "direct deep re-bin of the cell's mode-0 d_cell interior (leaf-only "
          "reification; NOT iterated cell_parent, which decoheres at nested "
          "splits). Cells already at `layer` pass through unchanged.");
    m.def("cell", &cell, nb::arg("uuid"), nb::arg("layer"), nb::arg("densify") = 0,
          nb::arg("sphere") = false,
          "Hexagon ring (npoints,2) lon/lat for a UUID at layer/densify. "
          "sphere=True emits spherical degrees.");
    m.def("grid", &grid,
          nb::arg("lon_min"), nb::arg("lat_min"), nb::arg("lon_max"), nb::arg("lat_max"),
          nb::arg("layer"), nb::arg("densify") = 0, nb::arg("max_cells") = 0,
          nb::arg("sphere") = false,
          "Enumerate cells in a lon/lat bbox -> (uuids[n,16], centroids[n,2], rings). "
          "rings is (n,6,2) when densify=0, else a list of (npoints,2) arrays. "
          "sphere=True: the bbox AND all emitted lon/lat are spherical degrees.");
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
    m.def("label", &label_many, nb::arg("uuids"), nb::arg("layer"), nb::arg("key") = false,
          "Batch form: (n,16) uuids -> list of n labels.");
    m.def("parse_label", &parse_label, nb::arg("label"),
          "Label (bare canonical, or keyed '.k' any flavour) -> (canonical bin uuid[16], layer).");
    m.def("label_centroid", &label_centroid, nb::arg("label"),
          "Geographic centroid (lon, lat) of the labelled cell (grid convention).");
    m.def("common_ancestor", &common_ancestor, nb::arg("uuids"), nb::arg("layer"),
          "Deepest common address-tree ancestor of (n,16) cells at layer -> "
          "(label, layer, uuid[16]), or None when cells span L0 hexes.");
    m.def("curve", &curve, nb::arg("uuid"),
          "Hamiltonian curve address of each (n,16) bin/full uuid -> packed "
          "curve-uuids (n,16), nibble 0 = 0xC. Byte order at one layer IS "
          "curve order; a curve-uuid truncates EXACTLY to the lineage "
          "ancestor's address. Curve-uuid input passes through.");
    m.def("curve_decode", &curve_decode, nb::arg("uuid"),
          "Constructive inverse: (n,16) curve-uuids -> canonical bin uuids at "
          "each curve address's layer. h9-uuid input passes through.");
    m.def("curve_label", &curve_label, nb::arg("uuid"),
          "Human-readable curve label 'c<slot><base-9 ranks>' (e.g. 'c112504') "
          "of one (16,) uuid — purely positional, length carries the layer. "
          "Accepts a cell uuid (curved first) or a curve-uuid.");
    m.def("curve_label", &curve_label_many, nb::arg("uuids"),
          "Batch form: (n,16) uuids -> list of n curve labels.");
    m.def("e4h_encode", &e4h_encode, nb::arg("lon"), nb::arg("lat"),
          nb::arg("layer") = 6, nb::arg("depth") = 2, nb::arg("sphere") = false,
          "E4H addresses: h9 host bin at `layer` + 0xE marker + half + `depth` "
          "aperture-4 tail digits (exact classifier — bit-identical on every "
          "platform; see docs/universality.md). (n,16) uint8.");
    m.def("e4h_decode",
          [](u8_2d_in u, bool sphere) { return e4h_decode_impl(u, sphere, false); },
          nb::arg("uuid"), nb::arg("sphere") = false,
          "Representative (lon, lat) of each E4H leaf trapezoid centroid.");
    m.def("e4h_partner",
          [](u8_2d_in u, bool sphere) { return e4h_decode_impl(u, sphere, true); },
          nb::arg("uuid"), nb::arg("sphere") = false,
          "Representative (lon, lat) of each leaf's PARTNER half — the other "
          "half of the same fine hexagon (matched pairs share the final digit).");
    m.def("e4h_split", &e4h_split, nb::arg("uuid"),
          "(host[n,16], half[n], digits[n,28] int8, -1-padded) of E4H uuids.");
    m.def("e4h_bin", &e4h_bin, nb::arg("uuid"), nb::arg("depth"),
          "Truncate the aperture-4 tail to `depth` digits — suffix-local "
          "binning (exact, unlike the a9 body digits).");
    m.def("e4h_hex", &e4h_hex, nb::arg("uuid"),
          "Canonical HEXAGON key of each E4H half address — the matched "
          "pair's mode-0 member (both halves map to the same key; "
          "idempotent; every host owns exactly 4^depth hexagons). The GIS "
          "binning surface: E4H addresses are half-hex trapezoids, the "
          "hexagons are the pairs.");
    m.def("e4h_mode", &e4h_mode, nb::arg("uuid"),
          "Transport mode (0/1) per E4H uuid — parity of the descent's "
          "rotation accumulator; -1 where not E4H. Pair members always "
          "carry opposite modes.");
    m.def("e4h_depth", &e4h_depth, nb::arg("uuid"),
          "Tail digit count per uuid; -1 where not E4H.");
    m.def("e4h_label", &e4h_label_one, nb::arg("uuid"),
          "Label '<h9-label>E<half><digits>' for a (16,) E4H uuid.");
    m.def("e4h_label", &e4h_label_many, nb::arg("uuids"),
          "Batch overload: (n,16) -> list of labels.");
    m.def("e4h_parse_label", &e4h_parse_label, nb::arg("label"),
          "Parse an E4H label back to its (16,) uuid.");
    m.def("is_e4h", &is_e4h_one, nb::arg("uuid"),
          "True iff the (16,) uuid is E4H (any nibble == 0xE — decisive: 0xE "
          "cannot occur in h9 or curve uuids).");
    m.def("is_e4h", &is_e4h_many, nb::arg("uuids"),
          "Batch overload: (n,16) -> bool[n].");
    m.def("curve_parse_label", &curve_parse_label, nb::arg("label"),
          "Curve label -> curve-uuid (16,) (inverse of curve_label).");
    m.def("cell_children", &cell_children, nb::arg("uuid"),
          "The 9 canonical children of each (n,16) cell -> (n,9,16), in "
          "CURVE-RANK order (one generation: lineage == ownership).");
    m.def("curve_cells",
          [](u8_1d u, int layer) {
              if (u.shape(0) != 16) throw std::runtime_error("uuid must be (16,)");
              return cells_of(u.data(), layer, hex9_curve_cells, "hex9.curve_cells");
          },
          nb::arg("uuid"), nb::arg("layer"),
          "Every layer-`layer` LINEAGE descendant of one cell (16,), in curve "
          "order -> (bins[m,16], curves[m,16]), m = 9^depth. Transitive, "
          "curve-prefix-exact; 1/6 of descendant area displaced — use "
          "owned_cells for the geometrically-bounded partition.");
    m.def("owned_cells",
          [](u8_1d u, int layer) {
              if (u.shape(0) != 16) throw std::runtime_error("uuid must be (16,)");
              return cells_of(u.data(), layer, hex9_owned_cells, "hex9.owned_cells");
          },
          nb::arg("uuid"), nb::arg("layer"),
          "Every layer-`layer` OWNED sub-zone of one cell (16,) — cells whose "
          "cell_ancestor IS the cell — curve-sorted -> (bins[m,16], "
          "curves[m,16]), m = exactly 9^depth. Owned sets partition each "
          "layer (the OGC API-DGGS #108 aggregation relation); bounded "
          "excursion, not transitive.");
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
    warp_init();
}
