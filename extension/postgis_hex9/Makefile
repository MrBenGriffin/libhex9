# Makefile — PGXS build for postgis_hex9, the PostGIS Hex9 (H9) extension.
#
# Lives inside the libhex9 repository (extension/postgis_hex9/) and compiles
# the H9 core from ../../core directly — one repo installs both the library
# and the extension; the core cannot drift from the SQL surface. (PostGIS
# upstream will not take this into core, so it ships as an independent
# installable, h3-pg style.)
#
# Usage (from extension/postgis_hex9/):
#   make
#   sudo make install
#   make installcheck   # requires PostgreSQL with postgis already loaded

PG_CONFIG   ?= pg_config
PGXS        := $(shell $(PG_CONFIG) --pgxs)

# ── Extension identity ────────────────────────────────────────────────────────
EXTENSION   = postgis_hex9
MODULE_big  = postgis_hex9-3
# Objects live under ./build/ so Xcode (which scans the source dir for .o
# inputs) doesn't pick them up and double-link the same symbols. Both the
# pattern rules below and PGXS's link step see these paths consistently.
OBJS        = build/lwgeom_hex9.o build/hex9_c.o build/h9_warp_runtime.o build/h9_warp_embedded.o
DATA        = postgis_hex9--1.4.0.sql \
              postgis_hex9--1.0.0--1.1.0.sql \
              postgis_hex9--1.1.0--1.2.0.sql \
              postgis_hex9--1.2.0--1.3.0.sql \
              postgis_hex9--1.3.0--1.4.0.sql

# ── Regression tests ──────────────────────────────────────────────────────────
REGRESS      = hex9
REGRESS_OPTS = --load-extension=postgis --load-extension=postgis_hex9

# ── H9 core (shared with libhex9 / the Python module) ─────────────────────────
# The glue (lwgeom_hex9.cpp) sees ONLY the C ABI in hex9_c.h; the core C++
# headers are compiled solely into hex9_c.o.
H9_ROOT     := ../..
H9_CORE     := $(H9_ROOT)/core
H9_WARP_BIN := $(abspath $(H9_CORE)/WGS84_l5_warp_f6.full.f64g.h9warp)

# ── PostGIS source headers ────────────────────────────────────────────────────
# liblwgeom headers are not installed by packagers; point POSTGIS_SRC at an
# unpacked PostGIS source tree matching the loaded postgis version.
POSTGIS_SRC ?= /Users/ben/Documents/Projects/postgis-3.6.2
PG_CPPFLAGS  = \
    -I$(POSTGIS_SRC)/liblwgeom \
    -I$(POSTGIS_SRC)/libpgcommon \
    -I$(POSTGIS_SRC) \
    -I$(H9_CORE) \
    -I$(H9_ROOT) \
    -I.

# ── macOS: postgis-3.dylib is a Mach-O bundle; cannot link against it. ────────
# liblwgeom symbols are resolved at runtime in _PG_init via h9_lwgeom_resolve()
# (see h9_lwgeom_shim.h) — load_file($libdir/postgis-3) + dlsym — so no LOAD is
# needed in hex9.sql.in and the module loads in any backend order. The module
# carries no undefined liblwgeom symbols, but -undefined dynamic_lookup is kept
# for the remaining PG-core symbols resolved via -bundle_loader postgres.
ifeq ($(shell uname -s),Darwin)
SHLIB_LINK = -undefined dynamic_lookup
endif

# ── OpenMP (optional; parallelizes the batch encode/decode loops) ─────────────
# hex9_c.cpp drives hex9_encode_many (and the decode/bin batch loops) with
# #pragma omp parallel for. PGXS does not enable OpenMP, so without this the
# extension's batch paths run serial — correct results, just no thread speedup.
# Mirrors libhex9's CMake find_package(OpenMP); serial fallback when libomp is
# absent (no behaviour change). Force serial with H9_OPENMP=0.
H9_OPENMP ?= 1
ifeq ($(H9_OPENMP),1)
  ifeq ($(shell uname -s),Darwin)
    # Apple clang needs Homebrew's keg-only libomp (the hint CMake documents).
    LIBOMP := $(shell brew --prefix libomp 2>/dev/null)
    ifneq ($(wildcard $(LIBOMP)/include/omp.h),)
      OPENMP_CFLAGS := -Xpreprocessor -fopenmp -I$(LIBOMP)/include
      OPENMP_LIBS   := -L$(LIBOMP)/lib -lomp
      # Homebrew libomp is single-arch, but Postgres.app's pg_config emits a
      # universal (x86_64 + arm64) build. macOS is dropping Intel, so collapse
      # to the host arch rather than ship a slice with unresolved OpenMP symbols.
      # (arch flags live in CFLAGS/LDFLAGS from PGXS — normalized post-include.)
      H9_HOST_ARCH := $(shell uname -m)
    endif
  else
    # gcc and Linux clang accept -fopenmp directly for compile and link.
    OPENMP_CFLAGS := -fopenmp
    OPENMP_LIBS   := -fopenmp
  endif
endif
SHLIB_LINK += $(OPENMP_LIBS)

include $(PGXS)

# ── macOS: fix stale SDK path baked into pg_config --cppflags after Xcode upgrade
ifeq ($(shell uname -s),Darwin)
  _REAL_SDK := $(shell xcrun --sdk macosx --show-sdk-path 2>/dev/null)
  ifneq ($(_REAL_SDK),)
    override CFLAGS   := $(shell echo "$(CFLAGS)"   | sed "s|-isysroot [^ ]*|-isysroot $(_REAL_SDK)|g")
    override CPPFLAGS := $(shell echo "$(CPPFLAGS)" | sed "s|-isysroot [^ ]*|-isysroot $(_REAL_SDK)|g")
    override LDFLAGS  := $(shell echo "$(LDFLAGS)"  | sed "s|-isysroot [^ ]*|-isysroot $(_REAL_SDK)|g")
  endif
  # OpenMP-on-macOS (see above): collapse the universal arch flags to host-only
  # so the single-arch libomp links. Strips every "-arch X" pair, re-adds host.
  ifneq ($(H9_HOST_ARCH),)
    $(info postgis_hex9: OpenMP on macOS — building $(H9_HOST_ARCH)-only (libomp is single-arch))
    override CFLAGS  := $(shell echo "$(CFLAGS)"  | sed -E 's/-arch [^ ]+//g') -arch $(H9_HOST_ARCH)
    override LDFLAGS := $(shell echo "$(LDFLAGS)" | sed -E 's/-arch [^ ]+//g') -arch $(H9_HOST_ARCH)
  endif
endif

# ── C++ compilation ───────────────────────────────────────────────────────────
# PGXS only knows about .c; add an explicit pattern rule for .cpp.
# Output goes into build/ — see OBJS comment.
# H9_GIT_REV stamps every build so h9_version() can be matched against the
# loaded dylib (`-dirty` suffix when working tree has uncommitted edits).
H9_GIT_REV := $(shell git -C "$(CURDIR)" describe --always --dirty 2>/dev/null || echo unknown)
H9_CXXFLAGS = $(CFLAGS) -std=c++17 -fno-exceptions $(CPPFLAGS) $(PG_CPPFLAGS) \
              $(OPENMP_CFLAGS) -DH9_GIT_REV='"$(H9_GIT_REV)"' -MMD -MP
build/%.o: %.cpp | build
	$(CXX) $(H9_CXXFLAGS) -c -o $@ $<

# Core sources compiled straight from the libhex9 tree. h9_warp_embedded.cpp
# .incbin's the 2 MB warp blob at compile time (self-contained dylib, no
# runtime file load). The define carries its own quotes for the .incbin
# string literal, matching how libhex9's CMake passes it.
build/%.o: $(H9_CORE)/%.cpp | build
	$(CXX) $(H9_CXXFLAGS) -DH9_WARP_BLOB='"$(H9_WARP_BIN)"' -c -o $@ $<

# The C ABI implementation (repo root) — the only object that compiles the
# core C++ headers; the glue links against it through hex9_c.h.
build/hex9_c.o: $(H9_ROOT)/hex9_c.cpp | build
	$(CXX) $(H9_CXXFLAGS) -c -o $@ $<

# Header dependency tracking: -MMD -MP emits build/*.d listing each object's
# header prerequisites, so editing a header (e.g. h9_ct.h, h9_grid.h) triggers
# a rebuild. Without this, header-only changes silently produced a stale .dylib.
-include $(OBJS:.o=.d)

build:
	mkdir -p build

clean: clean-build
.PHONY: clean-build
clean-build:
	rm -rf build

# ── SQL file generation ───────────────────────────────────────────────────────
# MODULE_PATHNAME is left literal: PostgreSQL substitutes it at script-run
# time from the control file's module_pathname ($libdir/postgis_hex9-3) —
# portable across installations, and the same mechanism the shipped upgrade
# scripts rely on. (An absolute build-time path here once made the regression
# suite silently load a stale dylib from another PostgreSQL.)
postgis_hex9--1.4.0.sql: hex9.sql.in
	cp $< $@

all: postgis_hex9--1.4.0.sql

clean: clean-hex9-generated
.PHONY: clean-hex9-generated
clean-hex9-generated:
	rm -f postgis_hex9--1.4.0.sql
