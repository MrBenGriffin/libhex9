cmake_minimum_required(VERSION 3.16)

# ── macOS deployment target ─────────────────────────────────────────────────
# macOS only, and only as a default the caller can override.
#
# Without this, libhex9 builds against whatever SDK the host has while
# PostgreSQL's PGXS compiles the extension with its own (older) target — on a
# Postgres.app build that is -mmacosx-version-min=11.0. The link then warns
#     "building for macOS-11.0, but linking with dylib ... built for newer
#      version 26.5"
# and the resulting module CLAIMS to support macOS 11 while depending on a
# library that does not. Harmless when both sit on the same machine; a lie the
# moment anything is distributed.
#
# 11.0 (Big Sur) matches what the common PostgreSQL distributions target.
# Must be set BEFORE project() — CMake consumes it while configuring the
# compiler, so setting it later has no effect on the flags actually used.
if(APPLE AND NOT DEFINED CMAKE_OSX_DEPLOYMENT_TARGET)
  set(CMAKE_OSX_DEPLOYMENT_TARGET "11.0" CACHE STRING
      "Minimum macOS version libhex9 targets (match your PostgreSQL's)")
endif()

project(libhex9 VERSION 2.2.1 LANGUAGES CXX C)

# ── Version consistency ─────────────────────────────────────────────────────
# HEX9_VERSION in the public header is hand-maintained (the header is shipped
# and hand-written, not generated). Assert it matches project() so the two
# cannot drift: hex9_version() reported "0.1.0" for months precisely because
# nothing checked. A consumer's compile-time/runtime version guard is only as
# trustworthy as this equality.
file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/hex9_c.h" _h9_ver_line
     REGEX "^#define[ \t]+HEX9_VERSION[ \t]+\"")
string(REGEX REPLACE ".*\"([^\"]+)\".*" "\\1" _h9_hdr_ver "${_h9_ver_line}")
if(NOT _h9_hdr_ver STREQUAL PROJECT_VERSION)
  message(FATAL_ERROR
    "libhex9: version mismatch between build and public header.\n"
    "    CMakeLists.txt project(): ${PROJECT_VERSION}\n"
    "    hex9_c.h HEX9_VERSION:    ${_h9_hdr_ver}\n"
    "Update both together — consumers compare HEX9_VERSION against "
    "hex9_version() to detect a stale shared library.")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# ── Determinism: pinned FP contraction ──────────────────────────────────────
# The encode chain's floating-point program IS the address definition — the
# same lon/lat must yield the same uuid on every platform (universality,
# 2026-07-27). Transcendentals are already owned (core/h9_det_math.h); the
# remaining nondeterminism is compiler FMA contraction: arm64 clang fuses
# a*b+c by default, x86-64 cannot, and regime_pin caught the difference as
# platform-dependent deep-nibble flips. Pinning contraction OFF makes the
# arithmetic identical everywhere. Do NOT relax per-target "for speed" — a
# fused multiply-add changes addresses, which is a regime change.
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  add_compile_options(-ffp-contract=off)
endif()
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
include(GNUInstallDirs)

# ── Addressing layout: L30 (default, reclaims the h_term nibble) vs L29 (legacy)
# L30 is canonical: the h_term nibble is reclaimed as body[30], a full UUID IS the
# max-depth bin, and the deepest tail is the 3-bit canonical key (r_mo|p_c2). Build
# with -DHEX9_USE_L29=ON for the legacy 29-layer on-disk format (e.g. to reproduce
# test/golden_l29.py). The L30 invariant suite (test/l30_invariants.py) gates this.
option(HEX9_USE_L29 "Use the legacy 29-layer UUID layout (OFF = 30-layer, default)" OFF)
if(HEX9_USE_L29)
  add_compile_definitions(HEX9_USE_L29)
endif()

# ── Embedded warp blob ──────────────────────────────────────────────────────
# core/h9_warp_embedded.cpp .incbin's this file at compile time. It is the
# Sphere-L6 fundamental-domain (v4 wedge-fold) sidecar generated on the hhg9
# python side — checked into core/, so a normal build needs nothing external.
# If it is missing, the assembler error from .incbin is opaque; fail early with
# a pointer to the regeneration step instead.
#
# ONE blob since 2.0.0. The WGS84-trained F6 v3 field that used to sit beside
# it was a second addressing regime — same point, different address, nothing
# in the 16 bytes to distinguish them. Removed with the regime.
set(H9_WARP_BLOB_PATH
  "${CMAKE_CURRENT_SOURCE_DIR}/core/Sphere_l6_fund.f64g.h9warp")
if(NOT EXISTS "${H9_WARP_BLOB_PATH}")
  message(FATAL_ERROR
    "libhex9: embedded warp blob not found:\n"
    "    ${H9_WARP_BLOB_PATH}\n"
    "The library .incbin's this v4 fund (delta+gradient) sidecar at compile time.\n"
    "It is normally checked into core/. Regenerate it from a python env with\n"
    "hhg9 importable (input: Sphere_l6_fund_warp_data.npz):\n"
    "    python tools/support/export_warp_fund_v4.py")
endif()

# ── Warp embed: .incbin (default) vs generated TU (MSVC / opt-in) ───────────
# core/h9_warp_embedded.cpp embeds the blob with GNU-as `.incbin`, which MSVC
# cannot compile. Swap in a generated TU carrying the same h9:: symbols
# (array style for MSVC — its string literals cap at 64 KB). Same mechanism
# as the DuckDB extension tree; the generator script is shared. Runtime is
# identical either way: bytes in .rodata, no file load.
option(HEX9_EMBED_PORTABLE "Use the generated (no-.incbin) warp embed TU" OFF)
set(H9_WARP_EMBED_SRC core/h9_warp_embedded.cpp)
if(MSVC OR HEX9_EMBED_PORTABLE)
  set(H9_WARP_EMBED_SRC "${CMAKE_CURRENT_BINARY_DIR}/h9_warp_embedded_gen.cpp")
  if(NOT EXISTS "${H9_WARP_EMBED_SRC}")
    find_package(Python3 REQUIRED COMPONENTS Interpreter)
    if(MSVC)
      set(_h9_embed_style "--array")
    else()
      set(_h9_embed_style "")
    endif()
    execute_process(
      COMMAND "${Python3_EXECUTABLE}"
              "${CMAKE_CURRENT_SOURCE_DIR}/duckdb/scripts/gen_warp_embed.py"
              "${H9_WARP_BLOB_PATH}" "${H9_WARP_EMBED_SRC}" ${_h9_embed_style}
      RESULT_VARIABLE _h9_embed_rc)
    if(NOT _h9_embed_rc EQUAL 0)
      message(FATAL_ERROR "libhex9: gen_warp_embed.py failed (${_h9_embed_rc})")
    endif()
  endif()
endif()

# ── libhex9 shared library ──────────────────────────────────────────────────
add_library(hex9 SHARED
  hex9_c.cpp
  core/h9_warp_runtime.cpp
  ${H9_WARP_EMBED_SRC})

target_include_directories(hex9
  PUBLIC  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
          $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/core)      # the C++ core
target_compile_definitions(hex9 PRIVATE
  H9_WARP_ENABLE=1
  H9_WARP_BLOB="${H9_WARP_BLOB_PATH}")
set_target_properties(hex9 PROPERTIES
  VERSION ${PROJECT_VERSION} SOVERSION ${PROJECT_VERSION_MAJOR}
  OUTPUT_NAME hex9
  PUBLIC_HEADER hex9_c.h)

# ── OpenMP (optional) ───────────────────────────────────────────────────────
# Serial-correct without it; lights up the batch loops when present. Apple
# clang needs hints to brew's keg-only libomp — pass -DOpenMP_ROOT=$(brew --prefix libomp).
find_package(OpenMP QUIET)
if(OpenMP_CXX_FOUND)
  target_link_libraries(hex9 PRIVATE OpenMP::OpenMP_CXX)
  message(STATUS "libhex9: OpenMP enabled")
else()
  message(STATUS "libhex9: OpenMP NOT found — serial build")
endif()

# ── libhex9 static library ──────────────────────────────────────────────────
# For static linkers — notably the Rust `-sys` crate at geoplegma/ (consumed by
# GeoPlegma as a DGGRS backend). Same sources/headers/defs as the shared lib;
# OUTPUT_NAME hex9 → libhex9.a. OpenMP is intentionally NOT linked here so the
# archive has no external runtime dependency (serial-correct; the batch pragmas
# are no-ops without -fopenmp), keeping the -sys link line trivial.
add_library(hex9Static STATIC
  hex9_c.cpp
  core/h9_warp_runtime.cpp
  ${H9_WARP_EMBED_SRC})
target_include_directories(hex9Static
  PUBLIC  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
          $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/core)
target_compile_definitions(hex9Static PRIVATE
  H9_WARP_ENABLE=1
  H9_WARP_BLOB="${H9_WARP_BLOB_PATH}")
set_target_properties(hex9Static PROPERTIES OUTPUT_NAME hex9 POSITION_INDEPENDENT_CODE ON)

if(MSVC)
  # M_PI et al. hide behind this on MSVC (no-op elsewhere).
  target_compile_definitions(hex9 PRIVATE _USE_MATH_DEFINES)
  target_compile_definitions(hex9Static PRIVATE _USE_MATH_DEFINES)
endif()

# ── Wheel mode ──────────────────────────────────────────────────────────────
# Driven by pyproject.toml (scikit-build-core): builds ONLY the nanobind
# module as hex9/_core, statically linked (self-contained — no shared lib to
# bundle, no OpenMP: serial-correct, the archive's doctrine), and installs
# nothing else. The `hex9` python package sources live in python/hex9/.
option(HEX9_WHEEL "Build for a Python wheel (hex9._core; no lib/header installs)" OFF)

# ── install (cross-platform) ────────────────────────────────────────────────
if(NOT HEX9_WHEEL)
  install(TARGETS hex9 hex9Static
    LIBRARY        DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE        DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME        DESTINATION ${CMAKE_INSTALL_BINDIR}
    PUBLIC_HEADER  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
endif()

# ── tests ────────────────────────────────────────────────────────────────────
enable_testing()

# ── PDAL plugin (opt-in) ─────────────────────────────────────────────────────
# filters.hex9 — projects WGS84 lon/lat to Hex9 w_oct. Off by default so the
# PostGIS/CLI build is unaffected; enable with -DHEX9_PDAL=ON (needs PDAL).
# After enable_testing() so the plugin's guarded ctest registers.
option(HEX9_PDAL "Build the PDAL filters.hex9 plugin" OFF)
if (HEX9_PDAL)
  find_package(PDAL REQUIRED)
  add_subdirectory(pdal)
endif()

# ── DuckDB extension (opt-in) ────────────────────────────────────────────────
# duckdb/ is a self-contained extension-template tree (its own Makefile drives
# DuckDB's build, so it is NOT an add_subdirectory — see duckdb/README.md for
# the one-time clones it needs). This target is a thin driver over that flow.
option(HEX9_DUCKDB "Build the DuckDB h9 extension (drives duckdb/Makefile)" OFF)
if (HEX9_DUCKDB)
  add_custom_target(duckdb_h9 ALL
    COMMAND ${CMAKE_COMMAND} -E env GEN=ninja make -C ${CMAKE_SOURCE_DIR}/duckdb release
    COMMENT "Building DuckDB h9 extension via duckdb/Makefile (extension-template flow)"
    VERBATIM)
endif()

add_executable(smoke test/smoke.c)
target_link_libraries(smoke PRIVATE hex9)
add_test(NAME smoke COMMAND smoke)

add_executable(batch test/batch.c)
target_link_libraries(batch PRIVATE hex9)
add_test(NAME batch COMMAND batch)

add_executable(grid test/grid.c)
target_link_libraries(grid PRIVATE hex9)
add_test(NAME grid COMMAND grid)

add_executable(kring test/kring.c)
target_link_libraries(kring PRIVATE hex9)
add_test(NAME kring COMMAND kring)

add_executable(labels test/labels.c)
target_link_libraries(labels PRIVATE hex9 m)
add_test(NAME labels COMMAND labels)

add_executable(adaptive test/adaptive.c)
target_link_libraries(adaptive PRIVATE hex9)
add_test(NAME adaptive COMMAND adaptive)

add_executable(cellgeom test/cellgeom.c)
target_link_libraries(cellgeom PRIVATE hex9)
add_test(NAME cellgeom COMMAND cellgeom)

# Standing guard against the "prefix-cut bins to coarsen" lapse (split-hex
# 6/7/8 canonical leaf == mode-0 parent; re-bin from the FULL uuid).
add_executable(bin_prefix_guard test/bin_prefix_guard.c)
target_link_libraries(bin_prefix_guard PRIVATE hex9)
add_test(NAME bin_prefix_guard COMMAND bin_prefix_guard)

add_executable(curve test/curve.c)
target_link_libraries(curve PRIVATE hex9)
add_test(NAME curve COMMAND curve)

# Authalic latitude series (via-sphere front-end) vs the hhg9 oracle.
# Header-only — needs no hex9 link beyond libm.
add_executable(authalic test/authalic.c)
target_link_libraries(authalic PRIVATE m)
add_test(NAME authalic COMMAND authalic)

# Sphere-datum twins: parity with the WGS84 chain (minus reduction),
# divergence between datums, and the grid handle's datum memory.
add_executable(sphere_mode test/sphere_mode.c)
target_link_libraries(sphere_mode PRIVATE hex9 m)
add_test(NAME sphere_mode COMMAND sphere_mode)

# Integer lattice identity surface: arithmetic parity with the geometric
# chain, canonical shared vertex keys, ext census.
add_executable(cell_uv test/cell_uv.c)
target_link_libraries(cell_uv PRIVATE hex9 m)
add_test(NAME cell_uv COMMAND cell_uv)

# v4 wedge-fold warp vs the hhg9 fold build (frozen Python reference).
# Header-only internals; reads the fund blob + parity reference from the
# source tree.
add_executable(fund_warp test/fund_warp.cpp)
target_compile_definitions(fund_warp PRIVATE
  H9_REPO_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
add_test(NAME fund_warp COMMAND fund_warp)

# End-to-end via-sphere chain vs the hhg9 Python reference (public ABI).
# THE ORACLE: every frozen expectation in the suite is downstream of this one.
add_executable(via_sphere test/via_sphere.c)
target_link_libraries(via_sphere PRIVATE hex9 m)
target_compile_definitions(via_sphere PRIVATE
  H9_REPO_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
add_test(NAME via_sphere COMMAND via_sphere)

# The same chain at its worst-conditioned points: polar approach, the poles
# from every seam meridian, straddled seams into the polar band, and the
# octant corners. Covers the 88..90 deg pole x seam band that the hhg9 pytests
# reach from neither side. Reference: tools/support/gen_pole_seam_ref.py.
add_executable(pole_seam test/pole_seam.c)
target_link_libraries(pole_seam PRIVATE hex9 m)
target_compile_definitions(pole_seam PRIVATE
  H9_REPO_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
add_test(NAME pole_seam COMMAND pole_seam)

# Regime pin: the frozen lon/lat -> uuid -> curve goldens and the ownership
# ladder. Catches ANY projection/warp change — the guard that was missing when
# the 2.0.0 regime change slipped past the whole suite. Regenerate goldens
# deliberately with the gen_regime_pin tool below, never to "make it pass".
add_executable(regime_pin test/regime_pin.c)
target_link_libraries(regime_pin PRIVATE hex9)
target_compile_definitions(regime_pin PRIVATE
  H9_REPO_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
add_test(NAME regime_pin COMMAND regime_pin)

# Producer for the regime-pin goldens. Not a test — run by hand, and ONLY from
# a build where via_sphere passes (it records, it does not validate).
add_executable(gen_regime_pin EXCLUDE_FROM_ALL tools/gen_regime_pin.c)
target_link_libraries(gen_regime_pin PRIVATE hex9)
target_compile_definitions(gen_regime_pin PRIVATE
  H9_REPO_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

add_executable(full_id test/full_id.c)
target_link_libraries(full_id PRIVATE hex9 m)
add_test(NAME full_id COMMAND full_id)

# Face-coordinate (bring-your-own-projection) surface: the boct seam is
# byte-identical to the canonical chain, both 2D and 3D routes; centroids
# and rings cohere with their lon/lat twins through unproject.
add_executable(boct_io test/boct_io.c)
target_link_libraries(boct_io PRIVATE hex9 m)
add_test(NAME boct_io COMMAND boct_io)

# w_oct 3D storage CRS: seamless round-trip, unit-octahedron invariant,
# oid==sign(xyz), and b_oct<->xyz rotation round-trip.
add_executable(woct test/woct.c)
target_link_libraries(woct PRIVATE hex9 m)
add_test(NAME woct COMMAND woct)

# Great-circle adjacency validation against the encoder oracle (Ben's
# "demonstrably safe ground" harness). Default 400 pts/circle keeps ctest
# fast; soak-run by hand with a larger count, e.g. `gc_kring 5000`.
add_executable(gc_kring test/gc_kring.c)
target_link_libraries(gc_kring PRIVATE hex9 m)
add_test(NAME gc_kring COMMAND gc_kring)

# Dev tool: derives/validates the k-ring seam algebra against the geometric
# mesh (not part of the default build or test suite — run by hand).
add_executable(kring_probe EXCLUDE_FROM_ALL tools/diagnosis/kring_probe.cpp)
target_include_directories(kring_probe PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/core)

# Dev probe: claim-1 measurement (origin-descent full UUID decode vs centroid)
# and the h9_bin(full,L)==bin acceptance check. Run by hand.
add_executable(full_uuid_probe EXCLUDE_FROM_ALL tools/diagnosis/full_uuid_probe.cpp)
target_include_directories(full_uuid_probe PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/core)
target_link_libraries(full_uuid_probe PRIVATE hex9)

# Dev probe: re-measure NN vs containment encoder disagreement. Run by hand.
add_executable(encoder_ab_probe EXCLUDE_FROM_ALL tools/diagnosis/encoder_ab_probe.cpp)
target_include_directories(encoder_ab_probe PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(encoder_ab_probe PRIVATE hex9)

# Dev probe: PROTOTYPE 2D-Newton encoder inversion vs the beam. Run by hand.
add_executable(newton_probe EXCLUDE_FROM_ALL tools/diagnosis/newton_probe.cpp)
target_include_directories(newton_probe PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/core)
target_link_libraries(newton_probe PRIVATE hex9)

# Dev probe: stress the 2D-Newton inversion at the octant hot zones. Run by hand.
add_executable(newton_hotzones EXCLUDE_FROM_ALL tools/diagnosis/newton_hotzones.cpp)
target_include_directories(newton_hotzones PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/core)
target_link_libraries(newton_hotzones PRIVATE hex9)

# Dev probe: verify + benchmark the analytic-Jacobian inversion. Run by hand.
add_executable(newton_aj_probe EXCLUDE_FROM_ALL tools/diagnosis/newton_aj_probe.cpp)
target_include_directories(newton_aj_probe PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/core)
target_link_libraries(newton_aj_probe PRIVATE hex9)

# ── Python module (optional: built only when nanobind is available) ──────────
option(HEX9_PYTHON "Build the nanobind Python module" ON)
if(HEX9_PYTHON)
  find_package(Python 3.9 COMPONENTS Interpreter Development.Module QUIET)
  if(Python_FOUND)
    execute_process(
      COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
      OUTPUT_VARIABLE NB_CMAKE_DIR OUTPUT_STRIP_TRAILING_WHITESPACE
      RESULT_VARIABLE NB_RC ERROR_QUIET)
    if(NB_RC EQUAL 0)
      list(APPEND CMAKE_PREFIX_PATH "${NB_CMAKE_DIR}")
      find_package(nanobind CONFIG QUIET)
    endif()
  endif()
  if(nanobind_FOUND)
    if(HEX9_WHEEL)
      # Wheel target: import name hex9._core, static core, self-contained.
      nanobind_add_module(_core STABLE_ABI python/hex9_ext.cpp)
      target_compile_definitions(_core PRIVATE HEX9_PY_MODULE=_core)
      if(MSVC)
        target_compile_definitions(_core PRIVATE _USE_MATH_DEFINES)
      endif()
      target_link_libraries(_core PRIVATE hex9Static)
      install(TARGETS _core LIBRARY DESTINATION hex9)
      message(STATUS "libhex9: wheel module hex9._core enabled")
    else()
      nanobind_add_module(hex9_ext STABLE_ABI python/hex9_ext.cpp)
      target_link_libraries(hex9_ext PRIVATE hex9)
      message(STATUS "libhex9: nanobind Python module enabled")
      add_test(NAME smoke_py
               COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test/smoke.py)
      set_tests_properties(smoke_py PROPERTIES
               ENVIRONMENT "PYTHONPATH=$<TARGET_FILE_DIR:hex9_ext>")
    endif()
  else()
    message(STATUS "libhex9: nanobind not found — skipping Python module "
                   "(pip install nanobind to enable)")
  endif()
endif()
