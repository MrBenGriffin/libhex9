/**********************************************************************
 * hex9_c.h — C ABI for libhex9, the Hex9 (H9) DGGS core.
 *
 * libhex9 is a bundled C++ dependency (compiled in-tree like deps/wagyu and
 * deps/flatgeobuf). This header is the ONLY surface the pure-C extension
 * (extensions/postgis_hex9) sees; the C++ implementation lives behind it.
 *
 * DESIGN BOUNDARY
 *   - The core is geometry-library-agnostic: it speaks lon/lat degrees
 *     (WGS84) and 16-byte cell UUIDs. It does NOT use liblwgeom or any
 *     PostgreSQL type. (This differs deliberately from wagyu, which passes
 *     LWGEOM; keeping hex9 geometry-free lets the standalone CLI link the
 *     same core, and keeps all GSERIALIZED handling in the C glue.)
 *   - The C extension glue is responsible for: extracting lon/lat from a
 *     POINT geometry, building GSERIALIZED polygons from the ring coordinates
 *     this API returns, and the precise centroid-in-bounds containment test
 *     (lwpoly_contains_point). The core enumerates a grid by bounding box;
 *     the glue filters that to the actual bounds geometry.
 *
 * CONVENTIONS
 *   - All coordinates are degrees, SRID 4326 (lon, lat) order.
 *   - UUIDs are raw 16-byte big-endian buffers (uint8_t[16]); the glue wraps
 *     them in pg_uuid_t.
 *   - Integer-returning functions return 0 on success, non-zero on error
 *     (error text, when provided, is written to a caller-supplied buffer).
 *   - Layer is 0..29 (h9_bin), 1..29 (cell/grid). The core validates and
 *     returns an error rather than asserting.
 *
 * Part of the Hex9 (H9) Project. Licensed under the GNU GPL v2 or later
 * (to match PostGIS; supersedes the development-era Apache-2.0 header).
 **********************************************************************/

#ifndef HEX9_C_H
#define HEX9_C_H

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ── Lifecycle / configuration ─────────────────────────────────────────────
 *
 * hex9_warp_init() builds the authalic-warp interpolation state once. It is
 * idempotent and should be called from the extension's _PG_init(). Returns 0
 * on success; on failure writes a message into errbuf (if errlen > 0) and
 * returns non-zero (the caller may fall back to the identity warp).
 */
const char *hex9_version(void);                       /* libhex9 build/version string */
int         hex9_lmax(void);                          /* deepest addressable layer (29 legacy / 30) */
int         hex9_warp_init(char *errbuf, size_t errlen);

/* Runtime toggles, mirroring the extension GUCs. Safe to call any time. */
void        hex9_set_use_warp(int on);                /* hex9.use_warp: 0/1 */
void        hex9_set_encoder(int mode);               /* hex9.encoder: 0=nn, 1=containment */

/* ── Point addressing ──────────────────────────────────────────────────────*/

/* Encode (lon, lat) to its layer-29 cell UUID. Returns 0 on success. */
int  hex9_encode(double lon, double lat, uint8_t out_uuid[16]);

/* Decode a UUID to the representative (lon, lat). Returns 0 on success.
 * GUARANTEED for full UUIDs only. Bins are layer-scoped keys, not
 * addresses: decoding one resolves it through the identity machinery,
 * which mis-locates meta-bearing (split-hex / seam-flavour) cells —
 * docs/addressing-doctrine.md, F2. */
int  hex9_decode(const uint8_t uuid[16], double *lon, double *lat);

/* Bin a UUID to its CANONICAL cell key at `layer` (0..29). Returns 0 on
 * success. The bin is the identity coarsening h9_grid enumerates with, so the
 * cell geometrically CONTAINS the point and hex9_bin == hex9_k_disk(.,0) ==
 * h9_grid (the bin-keyed-geometry-by-JOIN guarantee). At a split-hex (6/7/8)
 * leaf this resolves to the canonical MODE-0 parent — NOT the address's
 * reversible mode-1 nibble ancestry, which names a sibling cell that does not
 * contain the point (~20% of split-hex cells). The full UUID stays the
 * reversible load-bearer; the bin is derived from it.
 * Always re-bin from the FULL UUID. Bin input: same layer = identity; coarser
 * layers are a FOSSIL — bins are layer-scoped keys, not nested addresses
 * (docs/addressing-doctrine.md, F3); see test/bin_prefix_guard. */
int  hex9_bin(const uint8_t uuid[16], int layer, uint8_t out_uuid[16]);

/* ── Batch addressing (the fast path) ──────────────────────────────────────
 *
 * Process `n` items in one call; parallelised with OpenMP when libhex9 was
 * built with it (serial otherwise). The point work is independent and the warp
 * state is read-only after hex9_warp_init(), so this scales near-linearly.
 * Python bindings call these with whole numpy arrays (GIL released).
 *
 *   lon/lat    : n doubles each.
 *   out_uuid   : n*16 bytes (row-major: item i at out_uuid + i*16).
 * Return 0 on success.
 */
int  hex9_encode_many(const double *lon, const double *lat, size_t n,
                      uint8_t *out_uuid);
int  hex9_decode_many(const uint8_t *uuid, size_t n,
                      double *lon, double *lat);
int  hex9_bin_many(const uint8_t *uuid, int layer, size_t n,
                   uint8_t *out_uuid);

/* ── Continuous projection (b_oct backend) ─────────────────────────────────
 *
 * Forward map ONLY — no descent, no layer, no L29 cap. Writes the continuous
 * WARPED octahedral coordinate b_oct: (cx, cy) in the octant's authalic-
 * barycentric frame (BRAW then the authalic warp, h9_warp_inv == Python's
 * AuthalicWarp.undo — the b_raw → b_oct step), plus oid 0..7. This is exactly
 * the (cx,cy) the descent in uuid_from_cxcy_full operates on (minus the
 * descent-internal preamble nudge), so an embedder that owns its own descent
 * can bin it to ANY depth — the L29 cap is in the discretiser, not here.
 * Respects hex9_set_use_warp: with the warp off this returns the unwarped
 * BRAW coordinate (== Python b_oct with no warp), verified bit-identical.
 *
 * This is the integration seam for embedders with a native continuous octant
 * coordinate (e.g. hhg9's b_oct = 2×float64 + oid). See
 * tools/boct_backend_design.md and tools/boct_parity_probe.py.
 *
 *   - (lon, lat) are degrees, WGS84, lon/lat order (NOT [lat,lon]).
 *   - cx/cy/oid are written per item; oid bit-packing is the encoder's own
 *     (bit0=u<0, bit1=v<0, bit2=w<0) — embedders must map it to their face id.
 *   - The batch form releases the GIL in the Python bindings and parallelises
 *     with OpenMP when built with it (independent per-point work).
 * Return 0 on success.
 */
int  hex9_project(double lon, double lat, double *cx, double *cy, int *oid);
int  hex9_project_many(const double *lon, const double *lat, size_t n,
                       double *cx, double *cy, int *oid);

/* Inverse of hex9_project: WARPED octant coordinate b_oct (cx, cy, oid) →
 * (lon, lat) degrees, WGS84, lon/lat order. Un-warps (AuthalicWarp.do) then
 * runs the geometric inverse; honours hex9_set_use_warp identically to
 * hex9_project, so the two are exact inverses in both warp modes (to the warp's
 * Newton tolerance, ~1e-9). Returns 0 on success, non-zero if the result is
 * not finite. The _many batch form parallelises with OpenMP. */
int  hex9_unproject(double cx, double cy, int oid, double *lon, double *lat);
int  hex9_unproject_many(const double *cx, const double *cy, const int *oid,
                         size_t n, double *lon, double *lat);

/* ── Labels ────────────────────────────────────────────────────────────────
 * Write a NUL-terminated label for the UUID at `layer` into buf. Returns the
 * string length (excluding NUL), or -1 on error / insufficient buffer.
 * _key appends ".<key_tail>". A 40-byte buffer is always sufficient.
 *
 * Labels name the CANONICAL bin (== hex9_bin), spelled for humans. They are
 * COSMETIC, not load-bearing: the body shows the full-hex digit and drops the
 * half-hex/mode, so labels are NOT cross-layer prefix-nested at split-hex
 * (6/7/8). Common-prefix locality lives in the full UUID (the load-bearer) and
 * is a downstream concern. The full UUID is resolved to its canonical bin
 * first, so label(full) == label(hex9_bin(full)).
 */
int  hex9_label(const uint8_t uuid[16], int layer, char *buf, size_t buflen);
int  hex9_label_key(const uint8_t uuid[16], int layer, char *buf, size_t buflen);

/* ── Label parsing / mesh prefix ops ───────────────────────────────────────
 *
 * Labels are the human/SQL-facing cell names produced by hex9_label (body
 * digits, one char per layer 0..L) and hex9_label_key (".<k>" appended). They
 * name the canonical bin, so parse(label_key(u, L)) == hex9_bin(u, L) (the
 * keyed tail is the canonical bin tail; parse rebuilds the bin exactly). NOTE:
 * canonical labels are NOT cross-layer prefix-nested at split-hex — the layer-l
 * label is NOT in general a prefix of deeper labels (the mode-0 home can switch
 * ancestry). Do not coarsen by cutting a label; re-bin from the full UUID (see
 * test/bin_prefix_guard).
 *
 * hex9_parse_label accepts either form. BARE labels search the six possible
 * tails and verify by canonical re-encode, returning the canonical bin.
 * Returns the layer (0..29), or -1 on error. Note layer-0 labels (single char)
 * name L0 cells, which the cell/grid APIs do not otherwise reach (layer >= 1).
 *
 * FOSSIL CAVEAT (docs/addressing-doctrine.md, F1): bare labels are NOT unique
 * at split-hex (6/7/8) bodies — the same body names two cells and the parse
 * silently returns one of them. Labels are names, not addresses; only the full
 * UUID is guaranteed recoverable.
 */
int  hex9_parse_label(const char *label, uint8_t out_uuid[16]);

/* Geographic centroid of the labelled cell (same convention as the grid
 * enumerator's per-cell centroid, including the half-hex 4-vertex mean).
 * Returns 0 on success. */
int  hex9_label_centroid(const char *label, double *lon, double *lat);

/* Deepest common ancestor (in the address hierarchy) of n cells given as
 * UUIDs (full, or bin at >= `layer`), all treated at `layer`. Writes the
 * ancestor's label into buf and, when out_uuid is non-NULL, its bin UUID.
 * Returns the ancestor layer (0..layer), or -1 when there is none (cells
 * span L0 hexes) / on error. A 32-byte buf is always sufficient.
 *
 * NOTE: this is ancestry in the ADDRESS tree (descent containment), the
 * usual DGGS caveat applies — hexagon children straddle their parent's
 * geometric boundary, so it is not exact geometric containment. The
 * common prefix lets a mesh be stored as one ancestor label + per-cell
 * suffixes (label chars ancestor_layer+1 .. layer).
 *
 * FOSSIL CAVEAT (docs/addressing-doctrine.md, F1): the returned UUID is
 * recovered by bare-prefix parse, so it is only trustworthy when the
 * ancestor body avoids split-hex ambiguity. The label and layer are
 * always right. */
int  hex9_common_ancestor(const uint8_t *uuids, size_t n, int layer,
                          char *buf, size_t buflen, uint8_t *out_uuid);

/* ── Single-cell geometry ──────────────────────────────────────────────────
 *
 * A cell hexagon is returned as a closed ring of interleaved (lon, lat) pairs.
 * `densify` (>=0) subdivides each of the 6 edges into 3^densify segments:
 *   point count = hex9_ring_npoints(densify) = 6 * 3^densify + 1.
 * The caller sizes out_lonlat to hold 2 * npoints doubles.
 *
 * hex9_cell_ring accepts full UUIDs (from hex9_encode) and bin UUIDs (from
 * hex9_bin); bin UUIDs are lossy (decode to a representative cell). Returns
 * the number of (lon,lat) points written, or -1 on error.
 */
int  hex9_ring_npoints(int densify);
int  hex9_cell_ring(const uint8_t uuid[16], int layer, int densify,
                    double *out_lonlat, int max_points);

/* ── Grid enumeration (set-returning) ───────────────────────────────────────
 *
 * Enumerate every cell at `layer` (1..29) whose geographic centroid lies
 * within the lon/lat bounding box. The glue then filters to the real bounds
 * geometry and builds polygons.
 *
 * hex9_grid_create allocates an opaque handle the glue iterates and must free
 * with hex9_grid_destroy. `max_cells` caps the enumeration (the
 * hex9.grid_max_cells GUC); if the estimate exceeds it, the call returns NULL
 * and writes a message to errbuf. `densify` is stored as the default ring
 * resolution but each ring is built on demand (see hex9_grid_cell_ring), so
 * the per-cell geometry cost is paid during the glue's per-row iteration.
 */
typedef struct hex9_grid hex9_grid;

hex9_grid *hex9_grid_create(double lon_min, double lat_min,
                            double lon_max, double lat_max,
                            int layer, int densify,
                            int64_t max_cells,
                            char *errbuf, size_t errlen);

int        hex9_grid_count(const hex9_grid *g);
void       hex9_grid_cell_uuid(const hex9_grid *g, int i, uint8_t out_uuid[16]);
/* Full L29 identity UUID of grid cell i — the reversible address; h9_bin(id,L')
 * is exact at every layer L' <= the grid layer (unlike re-binning a bin). */
void       hex9_grid_cell_id(const hex9_grid *g, int i, uint8_t out_uuid[16]);
void       hex9_grid_cell_centroid(const hex9_grid *g, int i, double *lon, double *lat);
/* Build cell i's ring at `densify` into out_lonlat; returns point count or -1. */
int        hex9_grid_cell_ring(const hex9_grid *g, int i, int densify,
                               double *out_lonlat, int max_points);
void       hex9_grid_destroy(hex9_grid *g);

/* ── Adaptive multi-layer grid (population digest) ──────────────────────────
 *
 * Aggregate weighted points into a mixed-layer cell set where cell values
 * respect a population ceiling, by bottom-up digestion: points are binned
 * at max_layer; a cell whose accumulated value reaches `floor_` EMITS
 * itself, consuming whole points (first-fit, in deterministic order) until
 * the next would push past `ceiling` — always at least one point, so a
 * single point heavier than the ceiling still lands; unconsumed excess
 * re-bins to the parent layer and the pass repeats, down to min_layer,
 * which emits everything that remains. Dense areas resolve into fine
 * cells, sparse areas aggregate upward, and the sample is captured
 * exactly: emitted values sum to the input weight total.
 *
 * Notes:
 *   - `weight` may be NULL (all points weigh 1). Values accumulate by
 *     weight ("the number at each address"), not by point count.
 *   - Emitted cells of different layers overlap geometrically — a parent
 *     holds only what its descendants did not digest. The cell set
 *     partitions the SAMPLE, not the surface.
 *   - Cells whose value exceeds the ceiling can only arise from a single
 *     overweight point, or at min_layer (which must absorb remainders).
 *   - min_layer 0..29 (0 = the twelve L0 cells), max_layer >= min_layer.
 *
 * Cells are reported in digestion order: max_layer first, UUID-sorted
 * within each layer.
 */
typedef struct hex9_adaptive hex9_adaptive;

/* INPUT IS n FULL UUIDs (n×16 bytes, from hex9_encode) — bin input is
 * REJECTED, like the k-family: the digest re-bins across layers, which is
 * only guaranteed from the full uuid (doctrine F3); and addresses, not
 * coordinates, are the digest's natural input (Ben's ruling 2026-06-12 —
 * pairs directly with full-uuid sample tables, and spares the per-point
 * encode the old geometry-coupled form paid internally). */
hex9_adaptive *hex9_adaptive_create(const uint8_t *uuids,
                                    const double *weight, size_t n,
                                    int min_layer, int max_layer,
                                    double ceiling, double floor_,
                                    char *errbuf, size_t errlen);
int   hex9_adaptive_count(const hex9_adaptive *a);
void  hex9_adaptive_cell(const hex9_adaptive *a, int i, uint8_t out_uuid[16],
                         int *layer, double *value, int64_t *npoints);
/* Representative FULL uuid of a point cell i digested. hex9_adaptive_cell's
 * out_uuid is the layer-scoped bin KEY; render geometry from THIS instead, so
 * decode takes identity_from_uuid's guaranteed full path, not the bin fossil. */
void  hex9_adaptive_cell_full(const hex9_adaptive *a, int i, uint8_t out_full[16]);
/* Per input point: the index of the emitted cell that digested it. */
void  hex9_adaptive_assign(const hex9_adaptive *a, int64_t *out);
void  hex9_adaptive_destroy(hex9_adaptive *a);

/* ── Neighbours / k-ring / k-disk ──────────────────────────────────────────
 *
 * Symbolic adjacency on the H9 mesh — exact integer lattice arithmetic, no
 * floating point (core/h9_kring.h; algebra validated against the geometric
 * mesh at L1..L4 by tools/kring_probe.cpp, and against the encoder oracle
 * by test/gc_kring.c).
 *
 * INPUT IS A FULL UUID (from hex9_encode) ONLY — bin input is REJECTED
 * (returns -1). Ben's ruling (2026-06-11): a bin is a layer-scoped key,
 * not an address; resolving one is approximate (the key tail cannot carry
 * the leaf parent mode), and the approximation risks mis-use and failed
 * conclusions from an unknown underlying error. Output cells are written
 * as canonical bin UUIDs at `layer`, sorted — keys for joining (e.g. to
 * h9_grid), not addresses for further traversal; to walk onward, encode a
 * point and pass the full UUID.
 *
 * Every cell has 6 neighbours except the 12 half-hex cells per layer that
 * make up the 6 octahedron-vertex hexagons; those have 5 (4 sides + the
 * partner half-hex, which shares the same hexagonal outline).
 */

/* The (up to 6) edge-adjacent cells. out_uuids holds 6*16 bytes.
 * Returns the neighbour count (5 or 6), or -1 on error. */
int     hex9_neighbors(const uint8_t uuid[16], int layer, uint8_t *out_uuids);

/* Nominal cell count of a k-disk, 1 + 3k(k+1) — an upper bound (the true
 * count is smaller when the disk covers an octahedron-vertex hexagon).
 * Use it to size out_uuids. Returns -1 for k < 0. */
int64_t hex9_disk_ncells(int k);

/* Cells at graph distance exactly k (ring) / at most k including the
 * centre (disk). Writes up to max_cells UUIDs into out_uuids; returns the
 * count written, or -1 on error (bad args or buffer too small). */
int64_t hex9_k_ring(const uint8_t uuid[16], int layer, int k,
                    uint8_t *out_uuids, int64_t max_cells);
int64_t hex9_k_disk(const uint8_t uuid[16], int layer, int k,
                    uint8_t *out_uuids, int64_t max_cells);

/* ── Diagnostics (optional; mirrors h9_diag) ───────────────────────────────
 * Writes the BRAW (pre-warp) and BARY (post-warp) descent coordinates for a
 * point into buf. Internal aid, not part of the stable surface.
 */
int  hex9_diag(double lon, double lat, char *buf, size_t buflen);

#ifdef __cplusplus
}  /* extern "C" */
#endif

#endif /* HEX9_C_H */
