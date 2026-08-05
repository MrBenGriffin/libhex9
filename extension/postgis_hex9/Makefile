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

# ── Machine-local settings ───────────────────────────────────────────────────
# Optional, gitignored. The two variables below (PG_CONFIG, POSTGIS_SRC) are
# machine-specific and easy to get silently wrong; pin them here once:
#
#     # Makefile.local
#     PG_CONFIG   = /Applications/Postgres.app/Contents/Versions/18/bin/pg_config
#     POSTGIS_SRC = /path/to/postgis-3.6.3
#
# Included first so it wins over the defaults below.
-include Makefile.local

PSQL        ?= psql
PG_CONFIG   ?= pg_config
PGXS        := $(shell $(PG_CONFIG) --pgxs)

# ── Extension identity ────────────────────────────────────────────────────────
EXTENSION   = postgis_hex9
MODULE_big  = postgis_hex9-3
# Only the PostGIS glue is compiled here now; the H9 core is linked from the
# installed libhex9 (see below). Object lives under ./build/ so Xcode (which
# scans the source dir for .o inputs) doesn't pick it up and double-link.
OBJS        = build/lwgeom_hex9.o
DATA        = postgis_hex9--2.3.0.sql \
              postgis_hex9--2.1.0--2.3.0.sql \
              postgis_hex9--2.1.0.sql \
              postgis_hex9--2.0.0.sql \
              postgis_hex9--1.0.0--1.1.0.sql \
              postgis_hex9--1.1.0--1.2.0.sql \
              postgis_hex9--1.2.0--1.3.0.sql \
              postgis_hex9--1.3.0--1.4.0.sql \
              postgis_hex9--1.4.0--1.5.0.sql \
              postgis_hex9--1.5.0--2.0.0.sql \
              postgis_hex9--2.0.0--2.1.0.sql

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
# liblwgeom headers are installed by NO packager (not Homebrew, not Postgres.app,
# not the PGDG builds), so this glue needs an unpacked PostGIS SOURCE tree — and
# one whose MINOR version matches the PostGIS the server actually loaded, since
# LWGEOM/GSERIALIZED layouts are only stable within a minor release.
#
# There is deliberately NO default path. A wrong default is worse than none: it
# produces include errors that point at someone else's directory instead of
# naming the variable you need to set.
#
# If unset, try to find one — ask the running server which PostGIS it loaded,
# then look for a matching source tree in the usual places. The version-matched
# candidate is preferred; whatever is chosen is PRINTED, never assumed silently.
ifeq ($(origin POSTGIS_SRC), undefined)
  H9_GIS_VER := $(shell for db in postgres contrib_regression template1; do \
      v=`$(PSQL) -X -tA -c "SELECT postgis_lib_version()" $$db 2>/dev/null`; \
      [ -n "$$v" ] && { echo $$v; break; }; done)
  H9_GIS_ROOTS := $(CURDIR)/../.. $(CURDIR)/../../.. $(HOME)/src \
                  $(HOME)/Documents/Projects /usr/local/src /usr/src /opt
  H9_GIS_HITS := $(wildcard \
      $(foreach d,$(H9_GIS_ROOTS),$(d)/postgis-$(H9_GIS_VER)/liblwgeom/liblwgeom.h) \
      $(foreach d,$(H9_GIS_ROOTS),$(d)/postgis-*/liblwgeom/liblwgeom.h) \
      $(foreach d,$(H9_GIS_ROOTS),$(d)/postgis/liblwgeom/liblwgeom.h))
  POSTGIS_SRC := $(patsubst %/liblwgeom/liblwgeom.h,%,$(firstword $(H9_GIS_HITS)))
  ifneq ($(POSTGIS_SRC),)
    $(info postgis_hex9: found PostGIS source at $(POSTGIS_SRC)$(if $(H9_GIS_VER), (server runs $(H9_GIS_VER))))
  endif
endif
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
# Single source of truth for the extension version: the control file. Keeps
# h9_version() from drifting (it reported "1.2.0" for three releases).
H9_EXT_VERSION := $(shell sed -n "s/^default_version[ 	]*=[ 	]*'\(.*\)'/\1/p" $(CURDIR)/postgis_hex9.control)
H9_CXXFLAGS = $(CFLAGS) -std=c++17 -fno-exceptions $(CPPFLAGS) $(PG_CPPFLAGS) \
              -DH9_GIT_REV='"$(H9_GIT_REV)"' \
              -DH9_EXT_VERSION='"$(H9_EXT_VERSION)"' -MMD -MP
build/%.o: %.cpp | build check-postgis-src
	$(CXX) $(H9_CXXFLAGS) -c -o $@ $<

# Named prerequisite of compilation only, so `make clean` still works with
# POSTGIS_SRC unset.
.PHONY: check-postgis-src
check-postgis-src:
	@if [ -z "$(POSTGIS_SRC)" ] || [ ! -f "$(POSTGIS_SRC)/liblwgeom/liblwgeom.h" ]; then \
	  echo "postgis_hex9: POSTGIS_SRC is not set to a PostGIS source tree."; \
	  echo ""; \
	  echo "  liblwgeom's headers are not shipped by any PostGIS package, so this"; \
	  echo "  extension must be built against an unpacked PostGIS SOURCE tree."; \
	  echo ""; \
	  echo "  Either set the environment variable:"; \
	  echo "      POSTGIS_SRC=/path/to/postgis-source make"; \
	  echo "  or pass it to make:"; \
	  echo "      make POSTGIS_SRC=/path/to/postgis-source"; \
	  echo "  or pin it once in Makefile.local (gitignored):"; \
	  echo "      echo 'POSTGIS_SRC = /path/to/postgis-source' >> Makefile.local"; \
	  echo ""; \
	  if [ -n "$(POSTGIS_SRC)" ]; then \
	    echo "  (tried \"$(POSTGIS_SRC)\" — no liblwgeom/liblwgeom.h there)"; \
	  fi; \
	  echo "  It must match the MINOR version of the PostGIS your server loaded"; \
	  echo "  (SELECT postgis_lib_version();) — struct layouts change between"; \
	  echo "  minor releases and this module resolves them by dlsym, so a"; \
	  echo "  mismatch corrupts memory rather than failing to load."; \
	  exit 1; \
	fi

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
postgis_hex9--2.3.0.sql: hex9.sql.in
	cp $< $@

all: postgis_hex9--2.3.0.sql

clean: clean-hex9-generated
.PHONY: clean-hex9-generated
clean-hex9-generated:
	rm -f postgis_hex9--2.3.0.sql

# ── Which PostgreSQL, and which PostGIS? ─────────────────────────────────────
# macOS routinely has two PostgreSQL installs: a Homebrew/`/usr/local` one that
# owns the `pg_config` on PATH — and therefore the headers and libs dev work
# reaches for — and a Postgres.app one that is usually the server actually
# running. Build against the first and install, and everything SUCCEEDS: the
# files land in a prefix the live server never reads, so you go on testing a
# stale module. Same silent-mismatch family as loading a superseded libhex9.
#
# PostGIS is a third: liblwgeom headers are not installed by packagers, so
# POSTGIS_SRC points at an unpacked source tree that nothing keeps in step with
# the PostGIS the server loaded. Within a patch release the structs are stable;
# across a minor release LWGEOM/GSERIALIZED layout changes would corrupt memory
# rather than fail cleanly — and because the module resolves liblwgeom by dlsym
# (h9_lwgeom_shim.h), there is no link-time check to catch it.
#
# So: ask the RUNNING server what it actually uses, and compare. Skipped
# silently when no server is reachable (a build need not have one).
PSQL ?= psql

.PHONY: check-pg-match
check-pg-match:
	@srv_lib=`$(PSQL) -X -tA -c "SELECT setting FROM pg_config() WHERE name='PKGLIBDIR'" postgres 2>/dev/null`; \
	if [ -z "$$srv_lib" ]; then \
	  echo "postgis_hex9: no reachable server — skipping install-target check."; \
	else \
	  our_lib=`$(PG_CONFIG) --pkglibdir`; \
	  if [ "$$srv_lib" != "$$our_lib" ]; then \
	    echo "postgis_hex9: ERROR — PG_CONFIG points at a different PostgreSQL"; \
	    echo "    this build installs to : $$our_lib"; \
	    echo "    running server loads from: $$srv_lib"; \
	    echo "  Installing here would succeed and change nothing the server sees."; \
	    echo "  Rebuild with the server's own pg_config, e.g.:"; \
	    echo "      make PG_CONFIG=`$(PSQL) -X -tA -c \"SELECT setting FROM pg_config() WHERE name='BINDIR'\" postgres 2>/dev/null`/pg_config"; \
	    exit 1; \
	  fi; \
	  srv_gis=""; \
	  for db in postgres contrib_regression template1; do \
	    srv_gis=`$(PSQL) -X -tA -c "SELECT postgis_lib_version()" $$db 2>/dev/null`; \
	    [ -n "$$srv_gis" ] && break; \
	  done; \
	  src_gis=`sed -n 's/^#define POSTGIS_LIB_VERSION "\(.*\)"/\1/p' $(POSTGIS_SRC)/postgis_config.h 2>/dev/null`; \
	  if [ -z "$$src_gis" ]; then \
	    echo "postgis_hex9: WARNING — cannot read POSTGIS_LIB_VERSION from"; \
	    echo "    $(POSTGIS_SRC)/postgis_config.h"; \
	    echo "  POSTGIS_SRC may be wrong or unconfigured; header/runtime match UNVERIFIED."; \
	  elif [ -z "$$srv_gis" ]; then \
	    echo "postgis_hex9: WARNING — no database found with PostGIS installed;"; \
	    echo "  cannot verify POSTGIS_SRC ($$src_gis) against the server. UNVERIFIED."; \
	  elif [ "$$srv_gis" != "$$src_gis" ]; then \
	    srv_mm=`echo $$srv_gis | cut -d. -f1,2`; src_mm=`echo $$src_gis | cut -d. -f1,2`; \
	    if [ "$$srv_mm" != "$$src_mm" ]; then \
	      echo "postgis_hex9: ERROR — POSTGIS_SRC is a different MINOR release"; \
	      echo "    headers : $$src_gis  ($(POSTGIS_SRC))"; \
	      echo "    server  : $$srv_gis"; \
	      echo "  liblwgeom structs are not stable across minor releases and this"; \
	      echo "  module resolves them by dlsym, so a mismatch corrupts memory"; \
	      echo "  rather than failing to load. Point POSTGIS_SRC at $$srv_gis source."; \
	      exit 1; \
	    fi; \
	    echo "postgis_hex9: NOTE — POSTGIS_SRC $$src_gis vs server $$srv_gis"; \
	    echo "  (patch-level difference; struct layouts are stable within $$srv_mm)"; \
	  fi; \
	fi

installcheck: check-pg-match
