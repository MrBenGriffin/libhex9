# Makefile — PGXS build for postgis_hex9, the PostGIS Hex9 (H9) extension.
#
# Lives inside the libhex9 repository (extension/postgis_hex9/) and links the
# CMake-installed libhex9 for the H9 core: only the PostGIS glue
# (lwgeom_hex9.cpp) is compiled here — it sees ONLY the C ABI in hex9_c.h — and
# the core (C ABI impl + warp runtime + embedded warp blob, incl. its OpenMP
# batch loops) comes from the shared libhex9 produced by the top-level build:
#
#     cmake -S ../.. -B ../../build && cmake --build ../../build
#     cmake --install ../../build            # -> $(HEX9_PREFIX)/{lib,include}
#
# This keeps the H9 core a single source of truth — the extension can no longer
# drift from the library — and inherits libhex9's OpenMP-accelerated batch
# loops. libhex9.{so,dylib} must be present on the DB host at runtime; the rpath
# below points dyld/ld.so at the install prefix.
#
# Usage (from extension/postgis_hex9/):
#   make [HEX9_PREFIX=/your/cmake/install/prefix]
#   sudo make install
#   make installcheck   # requires PostgreSQL with postgis already loaded

PG_CONFIG   ?= pg_config
PGXS        := $(shell $(PG_CONFIG) --pgxs)

# ── Extension identity ────────────────────────────────────────────────────────
EXTENSION   = postgis_hex9
MODULE_big  = postgis_hex9-3
# Only the PostGIS glue is compiled here now; the H9 core is linked from the
# installed libhex9 (see below). Object lives under ./build/ so Xcode (which
# scans the source dir for .o inputs) doesn't pick it up and double-link.
OBJS        = build/lwgeom_hex9.o
DATA        = postgis_hex9--1.4.0.sql \
              postgis_hex9--1.0.0--1.1.0.sql \
              postgis_hex9--1.1.0--1.2.0.sql \
              postgis_hex9--1.2.0--1.3.0.sql \
              postgis_hex9--1.3.0--1.4.0.sql

# ── Regression tests ──────────────────────────────────────────────────────────
REGRESS      = hex9
REGRESS_OPTS = --load-extension=postgis --load-extension=postgis_hex9

# ── libhex9 (CMake-installed) ─────────────────────────────────────────────────
# The glue includes only hex9_c.h (the C ABI) and links the shared library.
# HEX9_PREFIX is the CMAKE_INSTALL_PREFIX used for `cmake --install` (default
# /usr/local). Override it if you installed libhex9 elsewhere.
HEX9_PREFIX ?= /usr/local
H9_INCDIR   ?= $(HEX9_PREFIX)/include
H9_LIBDIR   ?= $(HEX9_PREFIX)/lib

# ── PostGIS source headers ────────────────────────────────────────────────────
# liblwgeom headers are not installed by packagers; point POSTGIS_SRC at an
# unpacked PostGIS source tree matching the loaded postgis version.
POSTGIS_SRC ?= /Users/ben/Documents/Projects/postgis-3.6.2
PG_CPPFLAGS  = \
    -I$(POSTGIS_SRC)/liblwgeom \
    -I$(POSTGIS_SRC)/libpgcommon \
    -I$(POSTGIS_SRC) \
    -I$(H9_INCDIR) \
    -I.

# ── Link the installed libhex9 ────────────────────────────────────────────────
# -lhex9 from the install libdir; rpath so the module finds it at load time.
SHLIB_LINK  = -L$(H9_LIBDIR) -lhex9 -Wl,-rpath,$(H9_LIBDIR)

# ── macOS: postgis-3.dylib is a Mach-O bundle; cannot link against it. ────────
# liblwgeom symbols are resolved at runtime in _PG_init via h9_lwgeom_resolve()
# (see h9_lwgeom_shim.h) — load_file($libdir/postgis-3) + dlsym — so no LOAD is
# needed in hex9.sql.in and the module loads in any backend order. The module
# carries no undefined liblwgeom symbols, but -undefined dynamic_lookup is kept
# for the remaining PG-core symbols resolved via -bundle_loader postgres.
ifeq ($(shell uname -s),Darwin)
SHLIB_LINK += -undefined dynamic_lookup
# A plain CMake build emits a host-arch libhex9.dylib, but Postgres.app's
# pg_config emits universal (x86_64 + arm64) arch flags. Collapse the module to
# the host arch so it links against — and matches — the single-arch libhex9.
# (arch flags live in CFLAGS/LDFLAGS from PGXS — normalized post-include below.)
H9_HOST_ARCH := $(shell uname -m)
endif

include $(PGXS)

# ── macOS: fix stale SDK path baked into pg_config --cppflags after Xcode upgrade
ifeq ($(shell uname -s),Darwin)
  _REAL_SDK := $(shell xcrun --sdk macosx --show-sdk-path 2>/dev/null)
  ifneq ($(_REAL_SDK),)
    override CFLAGS   := $(shell echo "$(CFLAGS)"   | sed "s|-isysroot [^ ]*|-isysroot $(_REAL_SDK)|g")
    override CPPFLAGS := $(shell echo "$(CPPFLAGS)" | sed "s|-isysroot [^ ]*|-isysroot $(_REAL_SDK)|g")
    override LDFLAGS  := $(shell echo "$(LDFLAGS)"  | sed "s|-isysroot [^ ]*|-isysroot $(_REAL_SDK)|g")
  endif
  # Collapse the universal arch flags to host-only so the host-arch libhex9
  # links. Strips every "-arch X" pair, re-adds host.
  ifneq ($(H9_HOST_ARCH),)
    $(info postgis_hex9: macOS — building $(H9_HOST_ARCH)-only to match the installed libhex9)
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
              -DH9_GIT_REV='"$(H9_GIT_REV)"' -MMD -MP
build/%.o: %.cpp | build
	$(CXX) $(H9_CXXFLAGS) -c -o $@ $<

# Header dependency tracking: -MMD -MP emits build/*.d listing each object's
# header prerequisites, so editing a header (e.g. h9_lwgeom_shim.h) triggers a
# rebuild. Without this, header-only changes silently produced a stale .dylib.
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
