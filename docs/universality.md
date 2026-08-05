# Universality — the bit-exact address contract

Referenced from `core/h9_det_math.h`, `CMakeLists.txt`, `test/regime_pin.c`,
`test/wheel_pin.py`, and the CI workflows. This document is the normative
statement of hex9's universality contract: what is promised, what a
conforming platform is, how conformance is checked, and what a regime change
means. The doctrine already lives, enforced, in those files; this is the one
place it is stated whole.

## The promise

The same (lon, lat) yields the **bit-identical** 128-bit cell uuid — and the
bit-identical curve uuid — on every conforming platform. Not "agrees to
within a tolerance": identical bytes.

The promise matters because addresses are *stored*. They are minted once,
written into databases, joined against other datasets, exchanged between
systems, and re-derived years later on hardware that did not exist when they
were minted. An address that depends on which libm the minting machine linked
is not an identifier; it is a measurement. Hex9 addresses are identifiers.

Accuracy cannot buy this promise. A point may sit arbitrarily close to a cell
boundary — distance-to-boundary has no floor — so *any* platform-dependent
last-ulp difference eventually flips a nibble for some point. The first
multi-platform CI run (2026-07-27) measured exactly that: Linux and Windows
flipped deep nibbles against macOS-generated goldens for boundary-adjacent
points. The primary source was compiler FMA contraction on arm64; libm
variance (worst: `pow`) was the residual. The conclusion is structural, not
incidental: **universality is bought with bit-identical arithmetic, or not at
all.**

## The regime: the floating-point program is the definition

The encode chain's floating-point program — every operation, in order, at
IEEE-754 binary64 — **is** the definition of the address. There is no
abstract "true" cell against which implementations are measured for accuracy;
the program's output is the truth, and a conforming implementation reproduces
it bit-for-bit. `hex9_c.h` states the same doctrine from the C side: there is
ONE addressing regime (and since 2.1.0, two datums within it — sphere and
WGS84 — minting on the same regime).

What that program consists of, and how each nondeterminism class was closed:

- **Transcendentals are owned, not linked.** `core/h9_det_math.h` vendors
  `sin`/`cos`/`tan`/`atan2`/`hypot` verbatim from musl 1.2.5 (fdlibm
  lineage — coefficients, evaluation order, branch structure unchanged),
  header-only and `h9_`-prefixed. Apple, glibc and mingw libm disagree in
  the last ulps; the vendored kernels are identical everywhere.
- **`pow` is eliminated.** Every `pow(x, 0.25)` in the chain became the
  quarter-root `sqrt(sqrt(x))` — two IEEE-required correctly-rounded
  operations: deterministic on every conforming platform, and more accurate
  than libm `pow`.
- **`sqrt` and `fmod` remain libm, deliberately.** IEEE 754 requires them
  correctly rounded, so every conforming platform returns the same bits.
- **FMA contraction is pinned off, with both belts.** `-ffp-contract=off` in
  the build (GNU/Clang) and `#pragma STDC FP_CONTRACT OFF` in the header.
  arm64 clang fuses `a*b+c` by default; x86-64 cannot; a fused
  multiply-add changes addresses.
- **Excess precision is excluded at compile time.** `h9_det_math.h` raises
  `#error` unless `FLT_EVAL_METHOD == 0` — doubles evaluated as doubles.
  This is what rules 32-bit x87 builds out of the platform contract.
- **The warp field is an embedded artifact**, not a computation: the trained
  authalizing residual (the W of AKW — `docs/warp-regimes.md`) ships as
  frozen bytes, identical in every build.

Nothing in this list is a style choice. Reordering one addition in a kernel
changes addresses; so does "optimising" the quarter-root back to `pow`, or
relaxing contraction per-target "for speed". Any such edit is a **regime
change** (below), never a refactor.

## The conforming platform

A platform conforms if all of the following hold. Each is checked
mechanically — none is taken on trust:

1. **IEEE-754 binary64 arithmetic**, doubles evaluated in double precision
   (`FLT_EVAL_METHOD == 0`). Enforced by `#error` at compile time. In
   practice: x86-64 (SSE2), arm64, and wasm conform; 32-bit x87 does not
   and cannot.
2. **No FMA contraction** in the chain's translation units. Pinned by the
   build; a build that re-enables it (e.g. `-ffast-math`) is nonconforming
   and will fail the pins.
3. **Correctly rounded `sqrt` and `fmod`** — IEEE-required, so any
   conforming C library qualifies.

Proven in CI today: macOS (arm64), Linux (x86-64), Windows MinGW — the
Windows path reproduces the goldens bit-for-bit (MSVC can never build
libhex9 — GNU-as `.incbin` warp embed — so MinGW *is* the native Windows
path; `docs/building-on-windows.md`). Every published wheel platform is
additionally gated on the pins at build time (below).

## The instruments

Two distinct roles, kept separate on purpose — one validates, one pins:

- **`test/via_sphere` — the validator.** 4,010 rows against the frozen hhg9
  Python reference: the independent oracle that says the chain computes the
  *right* thing. Goldens are only ever generated from a build where this
  passes.
- **`test/regime_pin.c` — the pin.** Asserts the chain still produces the
  frozen addresses in `test_data/regime_pin.tsv`: 237 rows of
  `name  lon  lat  uuid  curve` covering named landmarks, poles, the
  antimeridian, octant corners and seams (where a regime change shows up
  first), plus a deterministic equal-area LCG sample — including `s005` and
  `s030`, the two points that historically flipped under FMA variance.
  Comparison is bit-exact; a failure reports the first differing nibble, so
  the depth of divergence tells you whether you are looking at a projection
  change (shallow) or arithmetic drift (deep tail). If regime_pin and
  via_sphere both fail, believe via_sphere.
- **`test/boct_io.c` — the seam.** Pins `project → encode_boct` byte-equal
  to `encode`, so the bring-your-own-projection surface cannot drift off
  the regime.
- **`test/wheel_pin.py` — the shipping gate.** Runs against the *installed*
  wheel, on every platform a wheel is built for, before anything is
  published (`.github/workflows/wheels.yml`). A failure is a regime
  violation or a broken build — never a tolerance issue.
- **`hex9.selftest()` — the runtime check.** The installed wheel carries a
  subset of the pins and re-mints them on demand, so an embedder can verify
  *their* environment conforms — catching a nonconforming rebuild at import
  time rather than as silent address drift in production.

- **`test/e4h_parity.c` + `test_data/e4h_pin.tsv` — the E4H corpus.** The
  aperture-4 structural tail (2.3.0) joins the frozen chain with its own
  validator-and-pin pair: the corpus rows were minted by the normative
  reference (`hhg9/h9/e4h.py`, generator `tools/gen_e4h_pin.py`) across all
  octants, seams, cone points and depths to the full nibble budget, and the
  C classifier must reproduce every uuid byte-for-byte. The E4H address
  definition is: the frozen det-math seed (project → host lattice-centre
  frame → seam unfold), ONE snap at `2^-46`, then pure integer descent in
  ℤ[½, √3] (`core/h9_e4h.h`; constants frozen from the reference by
  `tools/gen_e4h_tables.py`). Because every classification decision after
  the snap is exact integer arithmetic, E4H tails cannot drift with depth,
  platform, or optimisation level — the all-depths question is closed by
  construction. Knife-edge points (classification margin under `1e-5` in
  the residual frame) are excluded from the corpus: on a cut line either
  side is a valid address, exactly as in the reference's own census.

The uuid marker space is a registry: nibble sequences containing `0xE` are
E4H addresses, nibble 0 = `0xC` is a curve uuid, `0xF` is padding; neither
marker can occur in the other kinds, so the tests are decisive and every
h9/curve entry point rejects marked input it does not speak for.

Third parties claim conformance the same way the wheels do: reproduce
`test_data/regime_pin.tsv` bit-for-bit (and `test_data/e4h_pin.tsv` for the
E4H surface). The corpus is plain TSV,
regenerable from source (`tools/gen_regime_pin.c`) on any conforming
platform, and is the interchange artifact — an independent implementation
(another language, another runtime) that reproduces all 237 rows speaks the
regime; one that "agrees closely" does not.

A note on the wider hex9 ecosystem: only implementations that transliterate
the program can make the bit-exact promise. `hhg9` (the Python research
implementation) computes through numpy and platform libm; it agrees
*spatially* to deep levels but is not, and does not claim to be, bit-exact
at full depth. A pure-Python conforming implementation is possible — CPython
neither fuses nor carries excess precision — but only by porting the `h9_*`
kernels, not by calling `math.sin`.

## Regime changes

A regime change is any edit that moves any address: a kernel edit, a
reordering, a change to the projection or the warp field, a change to
contraction or evaluation settings. The procedure is deliberate and loud:

1. Decide it. A regime change moves every stored address in every
   downstream database. It is a product decision, not an implementation
   detail.
2. Validate the new chain against `test/via_sphere` (regenerating that
   reference first if the change is upstream of it).
3. Regenerate the pins — `tools/gen_regime_pin` — only from a build where
   via_sphere passes. The generator is a recorder, not a validator.
4. Bump the **major version**. The major version is the regime's name:
   all 2.x builds mint identical addresses.
5. Record it in the CHANGELOG as a regime change, as 2.0.0
   (WGS84-trained → via-sphere) and 2.1.0 (deterministic chain) are.

An unexplained pin diff is the other case: a bug. There is no third case.

E4H tails are addresses under the same rule: introducing them (2.3.0) moved
no existing address and is a minor; any later change to the E4H program —
the seed, the snap scale, the frozen tables (`core/h9_e4h_tables.h` carries
its generating hhg9 commit), the descent — moves stored E4H addresses and
is a regime change: major version, corpus regenerated deliberately, never
to "make it pass".

## History, briefly

- **2.0.0** exposed the gap: a regime change moved addresses from layer 7
  down and the whole suite passed, because the one pinned anchor happened to
  agree across regimes. Westminster's layer-8 bin moved with no test to say
  so. `regime_pin` exists so that can never happen silently again.
- **2.1.0 (2026-07-27)** closed the platform gap: first multi-platform CI
  proved libm/FMA variance flips deep nibbles; `h9_det_math.h` and the
  pinned build made the arithmetic identical everywhere, and the pins went
  strict — bit-exact on every CI platform, every wheel.
