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
 * Part of the Hex9 (H9) Project. Licensed under the Apache License 2.0
 * (see LICENSE and COPYRIGHT — the canonical statements for the work).
 * Apache-2.0 is one-way GPL-compatible, so the GPL'd PostGIS extension may
 * link this core freely; an earlier note here claiming GPL "to match
 * PostGIS" predated that understanding and was an error.
 **********************************************************************/

#ifndef HEX9_C_H
#define HEX9_C_H

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ── Version ───────────────────────────────────────────────────────────────
 *
 * HEX9_VERSION is a COMPILE-TIME constant: it records the version of the
 * header a consumer built against. hex9_version() is a RUNTIME call: it
 * reports the version of the shared library actually loaded. Comparing the
 * two is the only way a consumer can detect that it is running against a
 * different libhex9 from the one it was compiled for.
 *
 * That check matters more than it looks. The shared library carries an
 * SOVERSION, so a major bump installs under a NEW filename and does not
 * replace its predecessor. A consumer built against 1.x therefore keeps
 * loading 1.x after an upgrade, quite happily — and since 2.0.0 changed what
 * addresses the chain produces, it would go on emitting the old regime's
 * addresses while every version string around it said 2.0.0. Silent, and
 * invisible in the data. Consumers that persist addresses SHOULD refuse to
 * start on a mismatch: see postgis_hex9's _PG_init for the pattern.
 */
#define HEX9_VERSION       "2.2.1"
#define HEX9_VERSION_MAJOR 2
#define HEX9_VERSION_MINOR 2
#define HEX9_VERSION_PATCH 1

/* ── Lifecycle / configuration ─────────────────────────────────────────────
 *
 * hex9_init() builds the addressing chain once: the authalic latitude
 * series and the Sphere-L6 warp field (~1 s — the embedded v4 fund blob ships
 * its own gradients). It is idempotent and should be called from the
 * extension's _PG_init(). Returns 0 on success; on failure writes a message
 * into errbuf (if errlen > 0) and returns non-zero. On failure the library
 * degrades to an unwarped identity field: it will still produce addresses,
 * but they are NOT Hex9 addresses. Treat a non-zero return as fatal.
 * (hex9_warp_init is the pre-2.1.0 name, kept forever as an alias; the
 * function was never only about the warp.) One init serves BOTH datums —
 * the sphere entry points below need no series, but they do need the warp.
 *
 * There is ONE addressing regime. Until 2.0.0 a WGS84-trained field could be
 * selected at runtime, which meant a single point had two possible addresses
 * with nothing in the 16 bytes to tell them apart. Both the selector and that
 * field are gone. Data addressed by 1.x must be re-derived from its source
 * geometry — decode-and-re-encode silently invalidates anything derived from
 * it. See docs/warp-regimes.md.
 *
 * TWO DATUMS, ONE REGIME (2.1.0). Every address has always been minted on the
 * unit sphere; WGS84 enters only as the authalic latitude reduction at the
 * lon/lat boundary. The *_sphere twins below run the identical chain minus
 * that reduction: their (lon, lat) are already-spherical degrees. They exist
 * for callers that own their datum — another body's authalic frame, celestial
 * RA/dec — and make libhex9 a pure spherical addressing engine. The datum is
 * part of the function's IDENTITY, never process state: both families coexist
 * in one process, and each function stays immutable. But the 16 bytes do not
 * record the datum, and a WGS84-minted and sphere-minted address for the same
 * numeric (lon, lat) differ below ~layer 5 — the datum is DATASET metadata,
 * owned by the caller, exactly like the body the data belongs to. Never mix
 * datums within one dataset; never surface the choice as ambient state (no
 * GUC, no setting) — only as distinct immutable functions.
 */
const char *hex9_version(void);                       /* libhex9 build/version string */
int         hex9_lmax(void);                          /* deepest addressable layer (29 legacy / 30) */
int         hex9_init(char *errbuf, size_t errlen);
int         hex9_warp_init(char *errbuf, size_t errlen);   /* historic alias of hex9_init */

/* ── Dev / A-B research toggles ────────────────────────────────────────────
 *
 * Neither is wired to SQL, and neither should be: the addressing functions are
 * IMMUTABLE, so their result must not depend on mutable process state. An
 * index built under one setting would silently fail to match its own table
 * under the other. These exist for the hhg9 parity harness.
 */

/* 0 = identity (NOT Hex9 addresses — the raw, non-equal-area lattice),
 * 1 = apply the warp (default). The hhg9 b_oct backend's warp/no-warp A-B. */
void        hex9_set_use_warp(int on);

/* Encoder selection: 0 = nearest-neighbour (embedder/A-B research only),
 * 1 = containment (grid-canonical; the default). */
void        hex9_set_encoder(int mode);

/* ── Point addressing ──────────────────────────────────────────────────────*/

/* Encode (lon, lat) to its layer-29 cell UUID. Returns 0 on success.
 * (lon, lat) are degrees, WGS84. The _sphere twin takes already-spherical
 * degrees instead (see the datum doctrine above); everything else is
 * identical, twin by twin, for every *_sphere function in this header. */
int  hex9_encode(double lon, double lat, uint8_t out_uuid[16]);
int  hex9_encode_sphere(double lon, double lat, uint8_t out_uuid[16]);

/* Decode a UUID to the representative (lon, lat). Returns 0 on success.
 * GUARANTEED for full UUIDs only. Bins are layer-scoped keys, not
 * addresses: decoding one resolves it through the identity machinery,
 * which mis-locates meta-bearing (split-hex / seam-flavour) cells —
 * docs/addressing-doctrine.md, F2. */
int  hex9_decode(const uint8_t uuid[16], double *lon, double *lat);
int  hex9_decode_sphere(const uint8_t uuid[16], double *lon, double *lat);

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
int  hex9_encode_many_sphere(const double *lon, const double *lat, size_t n,
                             uint8_t *out_uuid);
int  hex9_decode_many(const uint8_t *uuid, size_t n,
                      double *lon, double *lat);
int  hex9_decode_many_sphere(const uint8_t *uuid, size_t n,
                             double *lon, double *lat);
int  hex9_bin_many(const uint8_t *uuid, int layer, size_t n,
                   uint8_t *out_uuid);

/* ── Canonical cell parent / ancestor (mode-0 convention) ──────────────────
 *
 * The CELL-level roll-up: which single layer-(L-1) cell is the canonical
 * parent of a layer-L cell (the parent x_cell containing the cell's mode-0
 * d_cell). Distinct from hex9_bin, which answers the POINT question and is
 * only guaranteed from a FULL uuid. Input is a layer-L bin UUID (L >= 1 for
 * parent); ancestor = iterated one-level parent down to `layer`, pass-through
 * when already there. Every parent has exactly 9 canonical children.
 * Return 0 on success; 1 on an L0 parent request, an input coarser than the
 * target layer, or a malformed uuid (batch forms fail whole-batch). */
int  hex9_cell_parent(const uint8_t uuid[16], uint8_t out_uuid[16]);
int  hex9_cell_parent_many(const uint8_t *uuid, size_t n, uint8_t *out_uuid);
int  hex9_cell_ancestor(const uint8_t uuid[16], int layer, uint8_t out_uuid[16]);
int  hex9_cell_ancestor_many(const uint8_t *uuid, int layer, size_t n,
                             uint8_t *out_uuid);

/* ── Hamiltonian curve addressing (space-filling curve) ────────────────────
 *
 * The H9 space-filling curve visits every layer-L cell exactly once, in an
 * order that is edge-adjacent between consecutive indices (local
 * Hamiltonicity) and REFINES: a cell's 9 lineage children occupy curve
 * indices index*9 .. index*9+8. Emitted by a closed 36-state transducer
 * walked down the cell's lineage chain (core/h9_curve.h; tables verbatim
 * from hhg9, machine-verified closure). Its two jobs are SORTING (curve
 * order is locality order — ORDER BY / CLUSTER / point-cloud ordering) and
 * GENERATION (enumerate a region's cells in locality order).
 *
 * The packed CURVE-UUID is the sortable form: nibble 0 = 0xC (type marker —
 * an h9-uuid's nibble 0 is a root digit 0..11, so the two kinds coexist in
 * one column; curve-uuids sort after all h9-uuids), nibble 1 = axiom slot,
 * then one base-9 rank nibble per layer, 0xF-padded. Bytewise uuid order at
 * a fixed layer IS curve order. Unlike an h9-uuid body, a curve-uuid
 * truncates EXACTLY: dropping rank nibbles gives the lineage ancestor's
 * curve address (hex9_curve_bin — no fold, no geometry). NOTE in
 * mixed-layer collections the 0xF padding makes an ancestor sort AFTER its
 * descendants (post-order); at one layer, or grouped by layer, order is the
 * curve itself.
 *
 * hex9_curve accepts a canonical bin uuid at any layer (from hex9_bin /
 * hex9_grid / hex9_adaptive) or a full uuid (re-binned canonically at its
 * own depth first); a curve-uuid passes through unchanged. The curve is
 * defined on canonical hexagon bins. Encode is pure address arithmetic —
 * exact at any depth. Returns 0 on success, non-zero on malformed input
 * (a T1 transducer miss means a NON-CANONICAL address, never a table gap).
 */
int  hex9_curve(const uint8_t uuid[16], uint8_t out_curve[16]);
int  hex9_curve_many(const uint8_t *uuid, size_t n, uint8_t *out_curve);

/* Constructive inverse: curve-uuid -> the cell's canonical bin uuid at the
 * curve address's layer. Forward-fits from the root, selecting each child
 * from the parent's 9 canonical children (one bounded child-oracle call
 * per level). An h9-uuid input passes through unchanged. Returns 0 on
 * success. */
int  hex9_curve_decode(const uint8_t curve[16], uint8_t out_uuid[16]);
int  hex9_curve_decode_many(const uint8_t *curve, size_t n, uint8_t *out_uuid);

/* 1 if uuid[16] is a packed curve-uuid (nibble 0 == 0xC — positional test
 * only), else 0. */
int  hex9_is_curve(const uint8_t uuid[16]);

/* Layer of a packed curve-uuid (0..lmax), or -1 if not a curve-uuid. */
int  hex9_curve_layer(const uint8_t curve[16]);

/* Layer-`layer` curve ancestor — PURE PREFIX TRUNCATION of the rank
 * nibbles, exact on the curve/lineage tree (contrast hex9_bin's fossil
 * caveats: the curve-uuid has no such defect). Input must be a curve-uuid
 * at >= layer. Returns 0 on success. */
int  hex9_curve_bin(const uint8_t curve[16], int layer, uint8_t out_curve[16]);

/* Curve index as a DECIMAL STRING (indices reach 12*9^L - 1, overflowing
 * int64 above L18; the string feeds SQL numeric / Python int exactly).
 * Accepts an h9 uuid (encoded first) or a curve-uuid (pure arithmetic).
 * Writes the NUL-terminated numeral into buf (40 bytes always suffice);
 * returns the string length, or -1 on error. */
int  hex9_curve_index(const uint8_t uuid[16], char *buf, size_t buflen);

/* Curve-uuid from a decimal index string at a given layer (the inverse of
 * hex9_curve_index on the representation; a bare index does not
 * self-describe its layer). Returns 0 on success, non-zero for a malformed
 * numeral or an index out of range for the layer. */
int  hex9_curve_pack(const char *index_dec, int layer, uint8_t out_curve[16]);

/* Human-readable curve label 'c<slot hex char><base-9 rank chars>' (e.g.
 * 'c112504'; purely positional, length carries the layer). Accepts an h9
 * uuid (encoded first) or a curve-uuid. Returns the string length, or -1
 * on error / insufficient buffer (40 bytes always suffice). */
int  hex9_curve_label(const uint8_t uuid[16], char *buf, size_t buflen);

/* Parse a curve label back to a curve-uuid. Returns 0 on success. */
int  hex9_curve_parse_label(const char *label, uint8_t out_curve[16]);

/* Number of layer-`to_layer` descendants of a layer-`from_layer` cell:
 * 9^(to_layer - from_layer). Returns -1 on bad layers or int64 overflow
 * (depth difference > 19). */
int64_t hex9_curve_ncells(int from_layer, int to_layer);

/* Enumerate every layer-`layer` LINEAGE descendant of a cell IN CURVE
 * ORDER (the generation primitive: emit order == ascending curve index).
 * The lineage relation — iterated one-generation canonical parents, the
 * curve's own tree — is transitive and exactly 9^d per cell, but its
 * descendant sets displace 1/6 of the ancestor's area (hexagon-band
 * spelling); for the geometrically-bounded partition use
 * hex9_owned_cells. `uuid` is an h9 bin/full uuid or a curve-uuid
 * (decoded first); `layer` >= the cell's own layer. Writes bin uuids to
 * out_bins (n*16 bytes, required) and, when out_curves is non-NULL, the
 * matching curve-uuids. Both arrays hold
 * hex9_curve_ncells(cell_layer, layer) rows; the call fails if that
 * exceeds max_cells. Returns the count, or -1 on error. */
int64_t hex9_curve_cells(const uint8_t uuid[16], int layer,
                         uint8_t *out_bins, uint8_t *out_curves,
                         int64_t max_cells);

/* ── Lineage vs ownership: the two exact hierarchy relations ───────────────
 *
 * Every cell has exactly 9 canonical one-generation children (the two
 * relations coincide there), but beyond one generation they diverge on
 * hexagon-band cells (~1/9 of cells per extra level):
 *
 *   LINEAGE   — iterated one-generation parents (hex9_cell_parent, the
 *               curve tree). Transitive; curve-prefix truncation IS the
 *               lineage ancestor; 1/6 of descendant area displaced.
 *   OWNERSHIP — the direct deep fold (hex9_cell_ancestor). NOT
 *               transitive; geometrically bounded (only rim splits
 *               protrude, by their far half); owned sets PARTITION every
 *               layer at exactly 9^d per zone.
 *
 * This is the "owned sub-zone" relation proposed for OGC API — DGGS
 * (opengeospatial/ogcapi-discrete-global-grid-systems#108): sub-zones
 * give coverage, identifiers give lineage, aggregation needs ownership.
 */

/* The 9 canonical children of a cell (one generation — lineage ==
 * ownership), IN CURVE-RANK ORDER (the deterministic order the curve
 * induces). Input is an h9 bin/full uuid or a curve-uuid. out_uuids holds
 * 9*16 bytes. Returns 0 on success; 1 for a cell at lmax (no children)
 * or malformed input. */
int  hex9_cell_children(const uint8_t uuid[16], uint8_t *out_uuids);

/* Enumerate every layer-`layer` OWNED sub-zone of a cell — the cells
 * whose hex9_cell_ancestor at the cell's layer IS the cell — exactly
 * 9^(layer - cell_layer), CURVE-SORTED. The owned sets over all zones at
 * one layer partition the globe, so aggregation by owned sub-zone is
 * loss-free and double-count-free. Same signature conventions as
 * hex9_curve_cells (h9 bin default in out_bins; curve-uuids optional).
 * Returns the count, or -1 on error. */
int64_t hex9_owned_cells(const uint8_t uuid[16], int layer,
                         uint8_t *out_bins, uint8_t *out_curves,
                         int64_t max_cells);

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
int  hex9_project_sphere(double lon, double lat, double *cx, double *cy, int *oid);
int  hex9_project_many(const double *lon, const double *lat, size_t n,
                       double *cx, double *cy, int *oid);

/* Inverse of hex9_project: WARPED octant coordinate b_oct (cx, cy, oid) →
 * (lon, lat) degrees, WGS84, lon/lat order. Un-warps (AuthalicWarp.do) then
 * runs the geometric inverse; honours hex9_set_use_warp identically to
 * hex9_project, so the two are exact inverses in both warp modes (to the warp's
 * Newton tolerance, ~1e-9). Returns 0 on success, non-zero if the result is
 * not finite. The _many batch form parallelises with OpenMP. */
int  hex9_unproject(double cx, double cy, int oid, double *lon, double *lat);
int  hex9_unproject_sphere(double cx, double cy, int oid, double *lon, double *lat);
int  hex9_unproject_many(const double *cx, const double *cy, const int *oid,
                         size_t n, double *lon, double *lat);
int  hex9_unproject_many_sphere(const double *cx, const double *cy, const int *oid,
                                size_t n, double *lon, double *lat);
int  hex9_project_many_sphere(const double *lon, const double *lat, size_t n,
                              double *cx, double *cy, int *oid);

/* NOTE (datum doctrine): the w_oct family below has NO _sphere twins on
 * purpose — b_oct <-> w_oct is a pure rotation, so sphere callers compose
 * hex9_project_sphere/boct_to_woct (and back). Likewise label centroids
 * compose hex9_parse_label + hex9_decode_sphere, and grid centroids compose
 * hex9_grid_cell_id + hex9_decode_sphere. Only functions where the reduction
 * is internal and non-composable get a twin. */

/* ── Warped octahedral cartesian CRS (w_oct): the 3D storage baseline ───────
 *
 * A seamless 3D coordinate for the whole sphere — the octahedral analogue of
 * ECEF. xyz lies on the unit octahedron (|x|+|y|+|z| == 1), and the octant id
 * is implicit in the sign octant: oid == ((z<0)<<2)|((y<0)<<1)|(x<0). So a 3D
 * octree rooted at the origin splits into the 8 octants at level 1, and the
 * point's position IS its address (bin it with the encoder to any layer). There
 * is NO seam and NO cone-point handling: a 2-sphere embeds in 3D losslessly, so
 * the flattening obstruction that forces a cut in any 2D map simply does not
 * arise here. (The 2D octahedral net is a rendering device, not this.)
 *
 * This is the POST-WARP octahedral cartesian: the pure per-oid rotation lift of
 * the WARPED b_oct that hex9_project emits. b_oct <-> xyz is a rotation only;
 * the authalic warp lives entirely upstream, in lon/lat <-> b_oct. So the
 * addressing pipeline for a stored point is xyz -> b_oct (hex9_woct_to_boct) ->
 * descent (hex9_bin / hex9_encode), with NO warp re-applied. It is a distinct
 * domain from the pre-warp geometric octahedral cartesian (c_oct, which lifts
 * b_raw off the WGS84 ellipsoid) — same lift geometry, warped source. Named
 * w_oct to match hhg9, where c_oct is load-bearing.
 *
 * ALTITUDE is orthogonal: store a straight WGS84 altitude (reference units)
 * alongside these coordinates; it is never folded into the octahedral radius.
 *
 *   - hex9_to_woct / hex9_from_woct         : WGS84 lon/lat <-> xyz (cross the warp).
 *   - hex9_boct_to_woct / hex9_woct_to_boct : b_oct (cx,cy,oid) <-> xyz (rotation only).
 *   xyz is 3 doubles; the _many forms are row-major (item i at xyz + i*3) and
 *   parallelise with OpenMP. Return 0 on success.
 */
int  hex9_to_woct(double lon, double lat, double xyz[3]);
int  hex9_from_woct(const double xyz[3], double *lon, double *lat);
int  hex9_boct_to_woct(double cx, double cy, int oid, double xyz[3]);
int  hex9_woct_to_boct(const double xyz[3], double *cx, double *cy, int *oid);
int  hex9_to_woct_many(const double *lon, const double *lat, size_t n, double *xyz);
int  hex9_from_woct_many(const double *xyz, size_t n, double *lon, double *lat);
int  hex9_boct_to_woct_many(const double *cx, const double *cy, const int *oid,
                            size_t n, double *xyz);
int  hex9_woct_to_boct_many(const double *xyz, size_t n,
                            double *cx, double *cy, int *oid);

/* ── Face-coordinate addressing — bring your own ⟨polyhedron⟩ ───────────────
 * The ⟨face⟩⟨cell⟩ layers exposed directly at the (cx, cy, oid) chart, for
 * callers who supply their OWN sphere→octahedron projection (research
 * projections, AKW-comparison studies, coordinates from no sphere at all).
 * Descent is projection-blind; every downstream operation (bins, lineage,
 * k-rings, curve) works unchanged.
 *
 * DOCTRINE — non-canonical addresses. A uuid minted here means "this cell
 * under the CALLER's projection", not AKW's, yet is bit-indistinguishable
 * from a canonical address (no spare bits mark it — deliberately). The
 * projection identity is dataset metadata, exactly as the sphere/WGS84 datum
 * is: never mix projections within one dataset. See docs/projection.md.
 *
 * The (cx, cy) contract is hex9's face chart — the values hex9_project /
 * hex9_woct_to_boct emit. Guarantees split at this seam: ⟨face⟩⟨cell⟩ is
 * deterministic (universality); everything upstream is the caller's.
 *
 * encode: always the grid-canonical containment descent (hex9_set_encoder
 * does not apply). decode: the cell's lattice centroid at the uuid's bin
 * layer, in the chart of the returned oid — the same centroid hex9_decode
 * unprojects, emitted before unprojection. cell_ring: (cx, cy, oid) per
 * vertex, per-vertex frames (each vertex in the chart of ITS oid); closed
 * ring of hex9_ring_npoints(densify) points. Returns: 0 on success (ring:
 * point count, -1 on error); _many forms OR the per-item results. */
int  hex9_encode_boct(double cx, double cy, int oid, uint8_t out_uuid[16]);
int  hex9_encode_boct_many(const double *cx, const double *cy, const int *oid,
                           size_t n, uint8_t *out_uuid);
int  hex9_decode_boct(const uint8_t uuid[16],
                      double *cx, double *cy, int *oid);
int  hex9_decode_boct_many(const uint8_t *uuid, size_t n,
                           double *cx, double *cy, int *oid);
int  hex9_cell_ring_boct(const uint8_t uuid[16], int layer, int densify,
                         double *out_cx, double *out_cy, int *out_oid,
                         int max_points);

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
/* ── Cell lattice identity (integer UV) ─────────────────────────────────────
 *
 * The cell's EXACT integer lattice identity at `layer` — the address-side
 * view of its geometry, with no floating point and no datum: these values
 * sit upstream of b_oct, upstream of the warp, upstream of lon/lat, so a
 * WGS84-minted and a sphere-minted address for the same cell yield identical
 * keys, and no *_sphere twins exist or ever will.
 *
 * Centre key: (c_ia, c_ib) in octant c_oid's frame — the LATTICE HEX CENTRE
 * (for ext half-hexes, the octahedron-vertex hexagon's centre). This is the
 * mesh anchor, deliberately distinct from h9_decode's representative point
 * (the geometric centroid, which for ext cells is the 4-own-vertex mean):
 * anchor questions come here, point questions go to decode.
 *
 * Vertex keys: the 6 ring vertices in ring order, as CANONICAL pool keys —
 * out-of-face vertices via the exact-on-seam resolution (doctrine F5),
 * on-boundary vertices via the lexicographically smallest of their
 * equivalent representations across the shared edge(s). Per-point
 * deterministic, so the same physical vertex yields the SAME key from every
 * cell that touches it, including across octant seams. Pool by (ia, ib,
 * oid) equality to build shared-vertex meshes (cf. hhg9 HexMesh); across
 * layers, scale keys by 3^(L_fine - L) — every layer-L lattice point is a
 * finer-layer point.
 *
 * `ext` flags the seam-straddling chain cells (centred on an octant edge;
 * 12·3^layer per layer — the L0-descended chains): their identity lives in
 * the governing mode-0 frame and ring vertices 4-5 are the reflected
 * across-seam half. For these, this centre key and h9_decode's
 * representative point genuinely differ — decode returns the mode-0 half's
 * centroid (a point-in-cell guarantee), while c_ia/c_ib is the whole hex's
 * lattice centre, ON the seam (the mesh anchor).
 *
 * Conversion to b_oct is pure arithmetic:
 *     cx = ia * u1 / 3^layer,  cy = ib * v3 / 3^layer     (hex9_uv_units)
 *
 * `layer` may be a per-cell array in the _many form (digests mix layers).
 * Returns 0 on success; non-zero on a malformed uuid (batch: any row).
 */
int  hex9_cell_uv(const uint8_t uuid[16], int layer,
                  int64_t *c_ia, int64_t *c_ib, int *c_oid,
                  int64_t v_ia[6], int64_t v_ib[6], int v_oid[6],
                  int *ext);
int  hex9_cell_uv_many(const uint8_t *uuid, const int32_t *layer, size_t n,
                       int64_t *c_ia, int64_t *c_ib, int32_t *c_oid,
                       int64_t *v_ia, int64_t *v_ib, int32_t *v_oid,
                       int32_t *ext);
void hex9_uv_units(double *u1, double *v3);

int  hex9_ring_npoints(int densify);
int  hex9_cell_ring(const uint8_t uuid[16], int layer, int densify,
                    double *out_lonlat, int max_points);
int  hex9_cell_ring_sphere(const uint8_t uuid[16], int layer, int densify,
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

/* Sphere-datum twin: the bbox is spherical lon/lat, and the HANDLE REMEMBERS —
 * every lon/lat this grid subsequently emits (cell_centroid, cell_ring) is in
 * the datum it was created with. The datum is per-handle state, fixed at
 * create; the accessors below serve both datums unchanged. */
hex9_grid *hex9_grid_create_sphere(double lon_min, double lat_min,
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

/* ── Grid verbs: aim / walk_to / vision_cone ───────────────────────────────
 *
 * Cell-first field primitives over the validated neighbour algebra — the C
 * twins of the python wheel's hex9.verbs (ported from the hhg9
 * fox-and-rabbits PoC). Every verb takes and returns layer-scoped BIN KEYS
 * (from hex9_bin / hex9_grid); full-uuid re-derivation happens internally
 * (the guaranteed path — bins are not addresses). Bearings and distances
 * are advisory FP over cell centroids: they steer, they never mint — every
 * returned key comes from the canonical encode/bin/neighbour chain, so the
 * universality contract is untouched. Bearings are degrees clockwise from
 * north (great-circle initial bearing); the seam changes nothing.
 *
 *   aim         : the neighbour whose centroid bearing best matches.
 *   walk_to     : shortest hex path src → dest (A* over hex9_neighbors,
 *                 obstacle keys impassable; dest passable even if listed).
 *                 Writes the key list (first = src, last = dest); returns
 *                 the count, 0 when no path within max_expand node
 *                 expansions (max_expand <= 0 selects the default 5000),
 *                 -1 on error / out buffer too small.
 *   vision_cone : keys within k rings, within half_angle of bearing, and
 *                 not occluded by an obstacle key (sight lines sampled at
 *                 a third of the local hex pitch; an obstacle never
 *                 occludes itself — you can see the hedge, just not
 *                 through it). The source key is always included. Returns
 *                 the count (keys UUID-sorted), or -1 on error.
 */
int  hex9_aim(const uint8_t key[16], int layer, double bearing,
              uint8_t out_key[16]);
int64_t hex9_walk_to(const uint8_t src_key[16], const uint8_t dest_key[16],
                     int layer, const uint8_t *obstacles, size_t n_obstacles,
                     int64_t max_expand, uint8_t *out_keys, int64_t max_cells);
int64_t hex9_vision_cone(const uint8_t src_key[16], int layer,
                         double bearing, double half_angle, int k,
                         const uint8_t *obstacles, size_t n_obstacles,
                         uint8_t *out_keys, int64_t max_cells);

/* ── E4H: the aperture-4 structural tail (address extension) ───────────────
 *
 * An E4H address extends a hex9 HOST bin at attach layer L with an 0xE break
 * marker, one HALF nibble (0/1: which state-cut trapezoid of the host
 * hexagon), and up to (lmax-2-L) TAIL nibbles from {0, 1..5} descending the
 * aperture-4 half-hex rep-4 carrier — exact nesting, straight edges, and
 * suffix-local truncation (truncation = binning, unlike the a9 digits):
 *
 *     nibbles: [h9 body 0..L] [0xE] [half] [d1..dB] [0xF pad] [key_tail]
 *     label  : <h9-label>E<half><d1..dB>            e.g. 0031586E1213
 *
 * E4H uuids are ADDRESSES (ruling 2026-08-04), under the same universality
 * contract as the h9 chain: the classifier is EXACT — a frozen det-math FP
 * seed, one snap at 2^-46, then pure integer arithmetic in Z[1/2, sqrt3]
 * (core/h9_e4h.h; tables frozen from the hhg9 reference by
 * tools/gen_e4h_tables.py) — so the same point mints the bit-identical E4H
 * uuid on every platform at every depth. Digit semantics: 0 = centre child;
 * 1..5 = the five-symbol enumeration (the digit of an edge child names the
 * octahedral axis of the neighbour octant its class points at); matched
 * pairs are global — the two halves of every fine hexagon share their final
 * digit (hex9_e4h_partner walks the pair). Normative reference:
 * hhg9/h9/e4h.py and docs/dggs-transport-tilings.md §4b-§4d; conformance
 * corpus test_data/e4h_pin.tsv.
 *
 * 0xE cannot occur in any h9 or curve uuid, so hex9_is_e4h is decisive; all
 * OTHER uuid-consuming entry points in this header REJECT E4H input (the
 * marker-guard doctrine — E4H tails have their own operations; the h9
 * machinery must not silently mis-read them).
 *
 * encode: layer 0..lmax-2, depth 0..lmax-2-layer (the nibble budget).
 * decode/partner: representative lon/lat (leaf trapezoid centroid / the
 * mirror half across the leaf's long side). bin: TAIL truncation to
 * `depth` digits (suffix-local, exact). split: host bin + half + digits.
 * Errors: 0 ok; 1 bad args; 2 grammar (not E4H / invalid nibbles); 3 point
 * two seams from host (unreachable for canonically binned points — census
 * hhg9 experimental/e4h/e4h_symbolic.py). Batch forms OR the item results.
 */
int  hex9_e4h_encode(double lon, double lat, int layer, int depth,
                     uint8_t out_uuid[16]);
int  hex9_e4h_encode_sphere(double lon, double lat, int layer, int depth,
                            uint8_t out_uuid[16]);
int  hex9_e4h_encode_many(const double *lon, const double *lat, size_t n,
                          int layer, int depth, uint8_t *out_uuid);
int  hex9_e4h_encode_many_sphere(const double *lon, const double *lat, size_t n,
                                 int layer, int depth, uint8_t *out_uuid);
int  hex9_e4h_decode(const uint8_t uuid[16], double *lon, double *lat);
int  hex9_e4h_decode_sphere(const uint8_t uuid[16], double *lon, double *lat);
int  hex9_e4h_decode_many(const uint8_t *uuid, size_t n,
                          double *lon, double *lat);
int  hex9_e4h_decode_many_sphere(const uint8_t *uuid, size_t n,
                                 double *lon, double *lat);
int  hex9_e4h_partner(const uint8_t uuid[16], double *lon, double *lat);
int  hex9_e4h_partner_sphere(const uint8_t uuid[16], double *lon, double *lat);
int  hex9_e4h_split(const uint8_t uuid[16], uint8_t out_host[16],
                    int *half, uint8_t digits[28], int *ndigits);
int  hex9_e4h_bin(const uint8_t uuid[16], int depth, uint8_t out_uuid[16]);
int  hex9_e4h_depth(const uint8_t uuid[16]);   /* tail digit count, -1 if not E4H */
int  hex9_is_e4h(const uint8_t uuid[16]);      /* 1 iff any nibble == 0xE */
int  hex9_e4h_label(const uint8_t uuid[16], char *buf, size_t buflen);
int  hex9_e4h_parse_label(const char *label, uint8_t out_uuid[16]);

/* ── Diagnostics (optional; mirrors h9_diag) ───────────────────────────────
 * Writes the BRAW (pre-warp) and BARY (post-warp) descent coordinates for a
 * point into buf. Internal aid, not part of the stable surface.
 */
int  hex9_diag(double lon, double lat, char *buf, size_t buflen);

#ifdef __cplusplus
}  /* extern "C" */
#endif

#endif /* HEX9_C_H */
