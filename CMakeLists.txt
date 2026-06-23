cmake_minimum_required(VERSION 3.16)
project(libhex9 LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
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
# F6 v3 (delta+gradient) sidecar generated from the sinkhorn warp on the hhg9
# python side — checked into core/, so a normal build needs nothing external.
# If it is missing, the assembler error from .incbin is opaque; fail early with
# a pointer to the regeneration step instead.
set(H9_WARP_BLOB_PATH
  "${CMAKE_CURRENT_SOURCE_DIR}/core/WGS84_l5_warp_f6.full.f64g.h9warp")
if(NOT EXISTS "${H9_WARP_BLOB_PATH}")
  message(FATAL_ERROR
    "libhex9: embedded warp blob not found:\n"
    "    ${H9_WARP_BLOB_PATH}\n"
    "The library .incbin's this F6 v3 (delta+gradient) sidecar at compile time.\n"
    "It is normally checked into core/. To regenerate it, run the producers in\n"
    "tools/support/ from a python env with hhg9 importable (see that README):\n"
    "    python tools/support/export_warp_grads.py   # -> WGS84_l5_warp_grads.npz\n"
    "    python tools/support/export_warp_v3.py       # -> the .h9warp blob above\n"
    "Inputs come from the hhg9/hex9 python repo (WGS84_l5_warp_data.npz).")
endif()

# ── libhex9 shared library ──────────────────────────────────────────────────
add_library(hex9 SHARED
  hex9_c.cpp
  core/h9_warp_runtime.cpp
  core/h9_warp_embedded.cpp)

target_include_directories(hex9
  PUBLIC  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
          $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/core)      # the C++ core
target_compile_definitions(hex9 PRIVATE
  H9_WARP_ENABLE=1
  H9_WARP_BLOB="${H9_WARP_BLOB_PATH}")
set_target_properties(hex9 PROPERTIES
  VERSION 0.1.0 SOVERSION 0
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
  core/h9_warp_embedded.cpp)
target_include_directories(hex9Static
  PUBLIC  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
          $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/core)
target_compile_definitions(hex9Static PRIVATE
  H9_WARP_ENABLE=1
  H9_WARP_BLOB="${H9_WARP_BLOB_PATH}")
set_target_properties(hex9Static PROPERTIES OUTPUT_NAME hex9 POSITION_INDEPENDENT_CODE ON)

# ── install (cross-platform) ────────────────────────────────────────────────
install(TARGETS hex9 hex9Static
  LIBRARY        DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE        DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME        DESTINATION ${CMAKE_INSTALL_BINDIR}
  PUBLIC_HEADER  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# ── tests ────────────────────────────────────────────────────────────────────
enable_testing()
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
target_link_libraries(labels PRIVATE hex9)
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

add_executable(full_id test/full_id.c)
target_link_libraries(full_id PRIVATE hex9)
add_test(NAME full_id COMMAND full_id)

# Great-circle adjacency validation against the encoder oracle (Ben's
# "demonstrably safe ground" harness). Default 400 pts/circle keeps ctest
# fast; soak-run by hand with a larger count, e.g. `gc_kring 5000`.
add_executable(gc_kring test/gc_kring.c)
target_link_libraries(gc_kring PRIVATE hex9)
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
    nanobind_add_module(hex9_ext STABLE_ABI python/hex9_ext.cpp)
    target_link_libraries(hex9_ext PRIVATE hex9)
    message(STATUS "libhex9: nanobind Python module enabled")
    add_test(NAME smoke_py
             COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test/smoke.py)
    set_tests_properties(smoke_py PROPERTIES
             ENVIRONMENT "PYTHONPATH=$<TARGET_FILE_DIR:hex9_ext>")
  else()
    message(STATUS "libhex9: nanobind not found — skipping Python module "
                   "(pip install nanobind to enable)")
  endif()
endif()
