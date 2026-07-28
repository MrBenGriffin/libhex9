# The AKW projection

AKW — the **Anders Kaseorg Warp** projection — is the map at the heart of
hex9: an equal-area projection of the sphere onto the regular octahedron.
Every hex9 address is minted through it, but the projection is a first-class
artifact in its own right, independent of cells, and this document is its
public definition: the layer architecture, the canonical interfaces, what is
deliberately *not* canonical, and how the projection can be admitted into
rendering toolchains that have never heard of it.

## Name and provenance

The projection began as a question about smoothness. The first construction
was numerical: a smooth triangulation of WGS84 over the triangular hex9
transport (at L4), solved by spherical Kamada–Kawai relaxation. It worked —
the curvature was demonstrably smooth — but the result was non-analytical: a
mesh, not a map. A render of that curvature was posted to Mathematics Stack
Exchange ([question 5016695](https://math.stackexchange.com/questions/5016695/),
posed 2024-12-28) asking for an analytical equivalent, and **Anders Kaseorg**
answered (2025) with the closed-form model that remains the analytic core of
the forward projection to this day: the tangent-substitution map with
coupling parameter α ≈ 3.2278 — recognisable in `core/h9_math.h` as the
`tan`/fourth-root forward with its cross-coupling term. (Stack Exchange
content is CC BY-SA; this attribution, and the citation in the paper, are
part of the licence's terms, not just its spirit.) The `AK` codenames in
`core/h9_math.h` (the AK inverse; the AJ analytic-Jacobian forward) record
this lineage. AKW = that analytical model (**AK**) plus the trained
authalizing residual field (**W**, the warp — see `docs/warp-regimes.md`).

## The layer architecture

The full chain from geodetic coordinates to a cell address factorises into
orthogonal stages:

    (⟨auth⟩) ⟨sphere⟩ ⟨polyhedron⟩ ⟨face⟩ ⟨cell⟩

| Stage | What it does | Continuous? | ABI |
|---|---|---|---|
| ⟨auth⟩ | WGS84 → authalic sphere (the only ellipsoid-aware stage; optional — callers with their own datum skip it via the `*_sphere` twins) | yes | inside `hex9_project`; bypassed by `hex9_project_sphere` |
| ⟨sphere⟩ | unit-sphere point | yes | — |
| ⟨polyhedron⟩ | AK analytic map + authalizing warp → point **on the unit octahedron surface**, `[ox, oy, oz]`, \|ox\|+\|oy\|+\|oz\| = 1 | yes — seamless | `hex9_to_woct`, `hex9_from_woct` |
| ⟨face⟩ | face chart: octant id + in-face coordinates `(oid, cx, cy)`; `oid = sign(ox,oy,oz)` | per-face | `hex9_project`/`hex9_unproject` (+`_sphere`, +`_many`), `hex9_woct_to_boct`, `hex9_boct_to_woct` |
| ⟨cell⟩ | aperture-9 descent → uuid / bins / curve | discrete | `hex9_encode` and everything above it |

Two boundaries in this table carry doctrine:

- **PROJ's world ends after ⟨polyhedron⟩.** Everything up to and including
  `[ox, oy, oz]` is a continuous coordinate operation — the kind of thing a
  PROJ pipeline can carry (the sphere→surface map is a homeomorphism; the
  fold at octahedron edges bends the *embedding*, not the map). Everything
  from ⟨face⟩ down involves a discrete component and belongs to the DGGS
  library, not to a coordinate engine. A face-indexed tuple is not a
  coordinate system, and no apology is needed for that — it is an address
  in the making.
- **The DGGS is married to the octahedron; the recipe may not be.** The AKW
  *construction* — analytic base map onto a deltahedron plus a trained
  fundamental-domain warp — plausibly generalises to other triangle-faced
  polyhedra (face-id derivation is an argmax over face normals; the
  octahedron's `sign(xyz)` is simply the case where the argmax is free).
  That generalisation is a separate research project. hex9 addresses are
  octahedral, and nothing in this document promises otherwise.

## Canonical interfaces

**`w_oct` `[ox, oy, oz]` is the canonical continuous form** of an AKW-projected
point: seamless, single-valued, no discrete component, face recoverable as
`sign`. It is already the storage CRS for point clouds (PDAL); it is equally
the correct interchange form for any consumer who wants the projection
without the cells. Altitude is carried orthogonally (WGS84 height), never
baked into the surface point.

**`b_oct` `(oid, cx, cy)` is the canonical face-indexed form**, one chart per
octant. It is the input to descent and the last stop before addressing.

Both forms are **frozen**. The projection's outputs are part of libhex9's
deterministic floating-point program (see CHANGELOG 2.1.0, "universality"):
the same input yields the bit-identical result on every conforming platform,
enforced by `regime_pin` and the via-sphere oracle across the CI matrix. A
consumer of the bare projection inherits that guarantee with no extra work —
and, symmetrically, any change to the projection's arithmetic is a regime
change: deliberate, versioned, goldens re-derived. The projection cannot be
"improved" quietly, because addresses are downstream of every bit of it.

## Nets are presentation, not projection

There is no canonical flat AKW map. Any 2D unfolding of the octahedron
chooses a net — butterfly, rhombus, and kin — and every net pays with seams
placed by editorial decision, not geometry (a seamless flattening is
forbidden outright: Gauss–Bonnet). hex9 therefore ships **no canonical net**.
Nets exist as renderers — presentation-layer assemblies of the eight face
charts — and every rendered map should be understood as *a* net of AKW, not
*the* AKW map. (The hhg9 Python package provides net rendering; the C
library deliberately does not.)

## Bring your own ⟨polyhedron⟩

The layering is not just descriptive — it is an API. Since 2.2.0 the
⟨face⟩⟨cell⟩ layers accept and emit face coordinates directly
(`hex9_encode_boct`, `hex9_decode_boct`, `hex9_cell_ring_boct`, and their
Python twins), so a researcher with her **own** sphere→octahedron map — a
different warp, Snyder-style equal-area, a relaxation mesh, coordinates from
no sphere at all — gets the entire cell machinery: descent, bins, lineage,
k-rings, the curve. It fires both ways: alternative octahedral projections
can be studied *within* the hex9 frame, and AKW can be measured against them
on identical addressing infrastructure.

The contract, precisely:

- **The chart is hex9's warped face frame** — the (cx, cy, oid) values
  `hex9_project` and `hex9_woct_to_boct` emit. Feeding that chart back in
  reproduces the canonical chain bit-for-bit (`project → encode_boct` ==
  `encode`); feeding your own projection's coordinates mints your cells.
- **Non-canonical addresses.** A uuid minted from a non-hex9 projection
  means "this cell under *that* projection" yet is bit-indistinguishable
  from a canonical address — no spare bits mark it, deliberately (bins,
  labels and the curve consume the full layout). The projection identity is
  dataset metadata, exactly as the sphere/WGS84 datum is: never mix
  projections within one dataset. The mixing failure is silent semantic
  collision, which is why this paragraph exists.
- **Guarantees split at the seam.** ⟨face⟩⟨cell⟩ retains the universality
  regime — deterministic descent, bit-identical everywhere. Equal-area
  properties, and the determinism of everything upstream of the chart, are
  the caller's projection's to provide.
- **The warp must be live** (`hex9_init`): the chart is the warped frame,
  and an uninitialised library collapses it onto the raw frame.

## Admitting AKW to rendering toolchains (backdrops, QGIS, cartopy)

Three separable things are commonly conflated here:

- **WKT2 (ISO 19162) describes.** A fully self-contained WKT2 definition of
  an AKW net CRS can be written today, with no authority and no
  registration. Software executes it only if its engine implements the
  method it names.
- **Authorities name.** EPSG registration of a new method is a petition to
  IOGP that historically *follows* adoption (Equal Earth: paper 2018 → PROJ
  within weeks → EPSG code after). The OGC Naming Authority can mint
  definition URIs on a lighter process. Identity, not execution.
- **PROJ executes.** QGIS, GDAL, rasterio, and cartopy all delegate
  coordinate operations to PROJ. For rendering, PROJ is the only door that
  matters.

The practical admission path, in three stages, each superseding the last:

1. **Now, stock PROJ, no upstream contribution:** rendering forces a net
   anyway, and for a chosen net the whole lon/lat → net-plane map can be
   shipped as a **`+proj=tinshift`** file — a dense triangulation emitted
   from the known forward map, evaluated piecewise-linearly by unmodified
   PROJ. Accuracy is mesh-density-bounded (irrelevant for backdrops); seams
   sit on triangulation edges, which tinshift handles natively. Wrapped in
   WKT2 via PROJ's *PROJ-based operation method* (a `DerivedProjectedCRS`
   embedding the pipeline string), it becomes a single shareable `.wkt`
   that QGIS loads as a custom CRS, and the same pipeline serves
   pyproj/cartopy.
2. **With the paper: a native `+proj=` implementation** — exact, fast, the
   warp field as a data-backed operation (PROJ has precedent via
   proj-data grids).
3. **After adoption, if wanted: EPSG/OGC identifiers**, Equal-Earth-style.

## See also

- `docs/warp-regimes.md` — the warp field's own history and doctrine.
- `docs/addressing-doctrine.md` — what happens after ⟨face⟩.
- `CHANGELOG.md` 2.1.0 — the determinism (universality) regime the
  projection's outputs live under.
