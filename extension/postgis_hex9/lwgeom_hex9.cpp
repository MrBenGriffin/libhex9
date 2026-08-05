/*
 * lwgeom_hex9.cpp — PostgreSQL/PostGIS glue for the Hex9 (H9) DGGS.
 *
 * Part of the Hex9 (H9) Project
 * Copyright ©2025, Ben Griffin
 * Licensed under the Apache License, Version 2.0
 *
 * SQL functions (all use h9_ prefix to match existing Python reference):
 *
 *	 h9_encode(geometry)			→ uuid		 encode point to self-contained UUID
 *	 h9_decode(uuid)				→ geometry	 POINT(lon lat) SRID 4326
 *	 h9_bin(uuid, integer)		→ uuid		 IMMUTABLE — bin key at layer L (0..h9_lmax())
 *	 h9_cell(uuid, integer)		 → geometry	 cell polygon SRID 4326 (layer 1..h9_lmax())
 *	 h9_label(uuid, integer)		→ text		 human label e.g. '478232778'
 *	 h9_label_key(uuid, integer)	→ text		 label with key_tail e.g. '478232778.9'
 *	 h9_grid(geometry, integer)	 → TABLE(hex9 uuid, geom geometry, centroid geometry)
 *
 * Design boundary (see ../../hex9_c.h)
 * ─────────────────────────────────────
 *	 ALL Hex9 math lives behind the geometry-free C ABI in hex9_c.h: lon/lat
 *	 degrees and 16-byte UUIDs in, lon/lat rings and UUIDs out. This file is
 *	 pure glue — PostgreSQL types (Datum, pg_uuid_t, SRF state) and liblwgeom
 *	 geometry assembly (GSERIALIZED points/polygons, containment tests).
 *	 Nothing here touches core/h9_*.h; the same libhex9 surface serves the
 *	 standalone CLI and the Python module.
 *
 *	 The PG_FUNCTION_INFO_V1 macros and Datum prototypes must sit in
 *	 extern "C" blocks so the PostgreSQL loader can find them by name.
 */

extern "C" {
#include "postgres.h"
#include "fmgr.h"
#include "utils/uuid.h"		 /* pg_uuid_t, UUID_LEN */
#include "utils/builtins.h"	 /* CStringGetTextDatum etc. */
#include "utils/guc.h"		 /* DefineCustomIntVariable */
#include "catalog/pg_type.h"	/* UUIDOID, INT2OID */
#include "funcapi.h"			/* ReturnSetInfo, SRF_* macros */
#include "access/htup_details.h"
#include "utils/array.h"		 /* deconstruct_array (h9_common_ancestor) */
#include "utils/numeric.h"		 /* numeric_in/out (h9_curve_index/pack) */
#include "utils/lsyscache.h"	/* get_typlenbyvalalign (h9_encode_many) */

/* PostGIS */
#include "liblwgeom.h"
#include "lwgeom_pg.h"

PG_MODULE_MAGIC;

} /* extern "C" */

/* Runtime liblwgeom resolution — see header. Must follow liblwgeom.h (types)
 * and precede the redirect macros further down. */
#include "h9_lwgeom_shim.h"

#include "hex9_c.h"

/* ── GUC variables ────────────────────────────────────────────────────────
 *
 * Only resource limits are GUCs here. NOTHING that can change an address is,
 * and nothing that can should ever become one.
 *
 * The addressing functions are IMMUTABLE — they are used in functional
 * indexes and generated columns, so PostgreSQL is entitled to compute a value
 * once and reuse it forever. A session-settable input to those functions
 * means two sessions disagree about what a stored index entry means, and an
 * index built under one setting silently fails to match its own table under
 * the other. That is data corruption with no error raised.
 *
 * libhex9 keeps two such toggles for the hhg9 A-B parity harness
 * (hex9_set_use_warp, hex9_set_encoder). Both are deliberately unreachable
 * from SQL. Until 2.0.0, `hex9.use_warp` WAS a PGC_USERSET GUC wired to the
 * first of them — a latent instance of exactly this bug. It is gone; the core
 * default (warp on) is simply left in place.
 *
 * The sphere datum (2.1.0) is the counter-example done right: it changes what
 * an address means, so it is surfaced as DISTINCT immutable functions
 * (h9_encode_sphere, ...), never as a setting. The datum is part of the
 * function's identity — index-safe, session-independent — and dataset
 * metadata is the caller's job. See "Two datums, one regime" in
 * docs/warp-regimes.md. */
static int h9_grid_max_cells;

/* Resolve liblwgeom entry points (see h9_lwgeom_shim.h). Defined here, ahead of
 * the redirect macros, so the real symbol names reach dlsym. */
extern "C" void h9_lwgeom_resolve(void) {
	/* load_file runs PostGIS's _PG_init (wires liblwgeom's palloc/pfree
	 * handlers) and loads it RTLD_GLOBAL; no-op if already resident. */
	load_file("$libdir/postgis-3", false);

	#define H9_RESOLVE(fn)                                                     \
		do {                                                                   \
			h9lw_##fn = (h9lw_fp_##fn) dlsym(RTLD_DEFAULT, #fn);               \
			if (!h9lw_##fn)                                                    \
				ereport(ERROR,                                                 \
					(errcode(ERRCODE_UNDEFINED_FUNCTION),                      \
					 errmsg("postgis_hex9: cannot resolve liblwgeom symbol "   \
					        "\"%s\"", #fn),                                    \
					 errhint("A compatible PostGIS (3.x) must be installed."))); \
		} while (0)

	H9_RESOLVE(getPoint4d_p);
	H9_RESOLVE(gserialized_from_lwgeom);
	H9_RESOLVE(lwgeom_calculate_gbox);
	H9_RESOLVE(lwgeom_free);
	H9_RESOLVE(lwgeom_from_gserialized);
	H9_RESOLVE(lwpoint_free);
	H9_RESOLVE(lwpoint_make2d);
	H9_RESOLVE(lwpoly_construct);
	H9_RESOLVE(lwpoly_contains_point);
	H9_RESOLVE(lwpoly_free);
	H9_RESOLVE(ptarray_construct);
	H9_RESOLVE(ptarray_set_point4d);

	#undef H9_RESOLVE
}

/* Refuse to run against a libhex9 we were not built for.
 *
 * The shared library carries an SOVERSION, so a major upgrade installs under
 * a NEW filename (libhex9.2.dylib) and does NOT replace its predecessor
 * (libhex9.0.dylib). A module built against 1.x therefore keeps resolving and
 * loading 1.x indefinitely after an upgrade — no error, no missing symbol.
 *
 * That is benign for most libraries and dangerous for this one: 2.0.0 changed
 * the addressing chain, so such a module would keep emitting the PREVIOUS
 * regime's addresses while the control file, h9_version() and every release
 * note said 2.0.0. The divergence starts around layer 7 and is invisible in
 * the data — nothing in a 16-byte address records which regime made it.
 *
 * So: compare the header we compiled against (HEX9_VERSION) with the library
 * actually loaded (hex9_version()), and ERROR rather than proceed. A refusal
 * to load is recoverable in one command; silently mis-addressed data is not.
 */
static void h9_check_lib_version(void) {
	const char *runtime = hex9_version();   /* "libhex9 <ver> (<date> <time>)" */
	const char *p = runtime;

	/* Skip the leading "libhex9 " so we compare version to version. */
	const char *sp = strchr(runtime, ' ');
	if (sp) p = sp + 1;

	const size_t want = strlen(HEX9_VERSION);
	/* Match the version token exactly: equal prefix AND a terminator after it,
	 * so "2.0.0" does not satisfy a check for "2.0.0-rc1" or vice versa. */
	if (strncmp(p, HEX9_VERSION, want) != 0 ||
	    (p[want] != '\0' && p[want] != ' ')) {
		ereport(ERROR,
			(errmsg("postgis_hex9: libhex9 version mismatch"),
			 errdetail("Built against libhex9 %s, but loaded \"%s\".",
			           HEX9_VERSION, runtime),
			 errhint("A major libhex9 upgrade installs under a new SOVERSION "
			         "filename and leaves the old library in place, so this "
			         "module can keep loading the previous one. Rebuild and "
			         "reinstall the extension against the installed libhex9, "
			         "and remove superseded libhex9 shared libraries. This is "
			         "not cosmetic: libhex9 2.0.0 changed the addressing "
			         "chain, so a stale library produces different Hex9 "
			         "addresses for the same point.")));
	}
}

extern "C" PGDLLEXPORT void _PG_init(void);
extern "C" PGDLLEXPORT void _PG_init(void) {
	/* Resolve liblwgeom first — every geometry-producing function depends on
	 * it, and this also ensures PostGIS is loaded for the rest of init. */
	h9_lwgeom_resolve();

	/* Before anything can compute an address, prove we are talking to the
	 * libhex9 we were compiled against. */
	h9_check_lib_version();

	/* Warp init runs once per backend (idempotent in the core). On failure
	 * the core falls back to an identity field, which does NOT produce Hex9
	 * addresses — so this is an ERROR, not a warning. Emitting plausible
	 * wrong addresses is worse than refusing to start. */
	char errbuf[256];
	if (hex9_init(errbuf, sizeof(errbuf)) != 0) {
		ereport(ERROR,
			(errmsg("postgis_hex9: warp init failed: %s", errbuf),
			 errhint("Without the warp field the core falls back to an "
			         "identity projection, whose output is NOT a Hex9 "
			         "address. Refusing to start rather than emit them.")));
	}

	/* hex9.grid_max_cells — operator-tunable safety cap for h9_grid.
	 * Default 708,588 = 12 × 9⁵ = global L5 cell count. Bump it (e.g.
	 *   SET hex9.grid_max_cells = 5000000;
	 * ) when validating coarser global meshes or running area-uniformity
	 * tests at deeper layers. */
	DefineCustomIntVariable(
		"hex9.grid_max_cells",
		"Maximum estimated cells h9_grid() will produce before erroring.",
		"Estimated from bbox area / mean cell area at the requested layer. "
		"Raise via SET if your bbox legitimately exceeds the default.",
		&h9_grid_max_cells,
		708588,        /* boot */
		1000,          /* min */
		100000000,     /* max — 100M, plenty of headroom */
		PGC_USERSET,
		0,
		NULL, NULL, NULL);

	/* No addressing GUCs. See the note by h9_grid_max_cells: the core's warp
	 * and encoder toggles stay at their defaults here and are not reachable
	 * from SQL, because IMMUTABLE functions may not depend on session state. */
}

/* Redirect the natural liblwgeom names onto the runtime-resolved pointers
 * (h9_lwgeom_shim.h). Defined here — after every liblwgeom prototype and after
 * h9_lwgeom_resolve()'s body — so neither the library headers nor the resolver
 * are rewritten. Every call site below this point binds through the pointer,
 * leaving our module with no undefined liblwgeom symbols. */
#define getPoint4d_p            (*h9lw_getPoint4d_p)
#define gserialized_from_lwgeom (*h9lw_gserialized_from_lwgeom)
#define lwgeom_calculate_gbox   (*h9lw_lwgeom_calculate_gbox)
#define lwgeom_free             (*h9lw_lwgeom_free)
#define lwgeom_from_gserialized (*h9lw_lwgeom_from_gserialized)
#define lwpoint_free            (*h9lw_lwpoint_free)
#define lwpoint_make2d          (*h9lw_lwpoint_make2d)
#define lwpoly_construct        (*h9lw_lwpoly_construct)
#define lwpoly_contains_point   (*h9lw_lwpoly_contains_point)
#define lwpoly_free             (*h9lw_lwpoly_free)
#define ptarray_construct       (*h9lw_ptarray_construct)
#define ptarray_set_point4d     (*h9lw_ptarray_set_point4d)

/* ── h9_version() → text ─────────────────────────────────────────────────── */

extern "C" {
PG_FUNCTION_INFO_V1(h9_version);
/* H9_GIT_REV is injected by the Makefile via -DH9_GIT_REV=… so every
 * build is uniquely identifiable. The fallback "dev" applies when
 * compiled outside the Makefile (e.g. Xcode-only CLI builds). */
#ifndef H9_GIT_REV
#define H9_GIT_REV "dev"
#endif

/* Set by the Makefile from postgis_hex9.control's default_version, so this
 * string cannot drift from the extension's real version — it said "1.2.0"
 * through three releases because it was a hand-edited literal. */
#ifndef H9_EXT_VERSION
#define H9_EXT_VERSION "unknown"
#endif

Datum h9_version(PG_FUNCTION_ARGS) {
	char buf[256];
	snprintf(buf, sizeof(buf),
	         "postgis_hex9 " H9_EXT_VERSION "+" H9_GIT_REV
	         " built " __DATE__ " " __TIME__
	         " (%s)", hex9_version());
	PG_RETURN_TEXT_P(cstring_to_text(buf));
}
} /* extern "C" */

/* ── h9_lmax() → integer ─────────────────────────────────────────────────── */

extern "C" {
PG_FUNCTION_INFO_V1(h9_lmax);
/* Deepest addressable layer of the loaded libhex9: 30 for the default
 * (reclaimed) layout, 29 on a legacy HEX9_USE_L29 build. Authoritative — the
 * value is compiled into the linked library (hex9_lmax() from the C ABI), so it
 * reports the layout the extension is actually running, not what the SQL was
 * generated against. */
Datum h9_lmax(PG_FUNCTION_ARGS) {
	PG_RETURN_INT32(hex9_lmax());
}
} /* extern "C" */

/* ── Internal helpers ────────────────────────────────────────────────────── */

/*
 * Extract lon/lat degrees from a PostGIS geometry.
 * Returns false (and sets *err) if the geometry is not a POINT.
 */
static bool geom_to_lonlat(GSERIALIZED *g, double *lon, double *lat,
							const char **err) {
	LWGEOM	*lwg = lwgeom_from_gserialized(g);
	if (lwg->type != POINTTYPE) {
		*err = "h9 functions require a POINT geometry";
		lwgeom_free(lwg);
		return false;
	}
	LWPOINT *pt	= lwgeom_as_lwpoint(lwg);
	POINT4D	p;
	getPoint4d_p(pt->point, 0, &p);
	*lon = p.x;
	*lat = p.y;
	lwgeom_free(lwg);
	return true;
}

/*
 * Build a 2D POINT geometry from lon/lat degrees. SRID 4326 for the WGS84
 * functions; the sphere twins pass SRID_UNKNOWN (0) — spherical degrees have
 * no EPSG identity here, the datum/body is dataset metadata (never claim
 * 4326 for coordinates that are not WGS84: downstream ST_Transform would
 * silently corrupt).
 */
static GSERIALIZED *lonlat_to_geom(double lon, double lat, int32_t srid = 4326) {
	LWPOINT	*pt	= lwpoint_make2d(srid, lon, lat);
	GSERIALIZED *g	= gserialized_from_lwgeom((LWGEOM *)pt, NULL);
	lwpoint_free(pt);
	return g;
}

/*
 * Copy a raw 16-byte UUID value into a palloc'd pg_uuid_t.
 */
static pg_uuid_t *make_pg_uuid(const uint8_t bytes[UUID_LEN]) {
	pg_uuid_t *u = (pg_uuid_t *) palloc(sizeof(pg_uuid_t));
	memcpy(u->data, bytes, UUID_LEN);
	return u;
}

/* Normalise a closed ring of interleaved (lon, lat) pairs — as returned by
 * hex9_cell_ring / hex9_grid_cell_ring — across the antimeridian / octant
 * seam and assemble it into a GSERIALIZED POLYGON (SRID 4326). Shifts each
 * vertex's lon by ±360° relative to its predecessor so no consecutive pair
 * differs by more than 180°, then re-centres the ring. Same stopgap as the
 * legacy corner path: the proper split (ST_Split etc.) is projection-
 * dependent and left to the caller. Mutates lonlat[] in place. */
static GSERIALIZED *h9_polygon_from_ring(double *lonlat, int n_ring,
                                         int32_t srid = 4326) {
	/* Pairwise-normalise only the unique vertices (exclude the duplicate close
	 * at n_ring-1): a pole-winding ring legitimately accumulates ±360° around
	 * the circuit, so shifting the closing vertex relative to its predecessor
	 * would leave it != vertex 0 and break ring closure. */
	for (int i = 1; i < n_ring - 1; ++i) {
		while (lonlat[2*i] - lonlat[2*(i-1)] >  180.0) lonlat[2*i] -= 360.0;
		while (lonlat[2*i] - lonlat[2*(i-1)] < -180.0) lonlat[2*i] += 360.0;
	}
	double mean_lon = 0.0;
	for (int i = 0; i < n_ring - 1; ++i) mean_lon += lonlat[2*i];  /* exclude duplicate close */
	mean_lon /= (double)(n_ring - 1);
	if      (mean_lon >  180.0) for (int i=0; i<n_ring; ++i) lonlat[2*i] -= 360.0;
	else if (mean_lon < -180.0) for (int i=0; i<n_ring; ++i) lonlat[2*i] += 360.0;
	/* Close the ring exactly on the normalised vertex 0, so the LinearRing
	 * is always closed. */
	lonlat[2*(n_ring-1)]     = lonlat[0];
	lonlat[2*(n_ring-1) + 1] = lonlat[1];

	POINTARRAY **rings = (POINTARRAY **) palloc(sizeof(POINTARRAY *));
	POINTARRAY  *pa    = ptarray_construct(false, false, (uint32_t)n_ring);
	rings[0] = pa;
	for (int i = 0; i < n_ring; ++i) {
		POINT4D pt = { lonlat[2*i], lonlat[2*i + 1], 0.0, 0.0 };
		ptarray_set_point4d(pa, (uint32_t)i, &pt);
	}

	LWPOLY      *poly = lwpoly_construct(srid, NULL, 1, rings);
	GSERIALIZED *g    = gserialized_from_lwgeom((LWGEOM *)poly, NULL);
	lwpoly_free(poly);
	return g;
}

/* ── h9_encode(geometry) → uuid ─────────────────────────────────────────── */

extern "C" {
static Datum h9_encode_common(PG_FUNCTION_ARGS, bool sphere) {
	if (PG_ARGISNULL(0)) PG_RETURN_NULL();

	GSERIALIZED *g = PG_GETARG_GSERIALIZED_P(0);
	double lon, lat;
	const char *err = NULL;
	if (!geom_to_lonlat(g, &lon, &lat, &err))
		ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
						errmsg("%s", err)));

	uint8_t uuid_bytes[UUID_LEN];
	const int rc = sphere ? hex9_encode_sphere(lon, lat, uuid_bytes)
	                      : hex9_encode(lon, lat, uuid_bytes);
	if (rc != 0)
		ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
						errmsg("h9_encode: could not encode point (%g, %g)",
						       lon, lat)));
	PG_RETURN_UUID_P(make_pg_uuid(uuid_bytes));
}

PG_FUNCTION_INFO_V1(h9_encode);
Datum h9_encode(PG_FUNCTION_ARGS) { return h9_encode_common(fcinfo, false); }

PG_FUNCTION_INFO_V1(h9_encode_sphere);
Datum h9_encode_sphere(PG_FUNCTION_ARGS) { return h9_encode_common(fcinfo, true); }
} /* extern "C" */

/* ── h9_encode_many(geometry[]) → uuid[] ────────────────────────────────────
 *
 * Batch encode: one OpenMP-parallel pass (hex9_encode_many) over an array of
 * POINTs, returning the max-depth UUIDs in input order. Same result as mapping
 * h9_encode over the array, but it crosses the SQL/C boundary once and runs the
 * independent point work in parallel — the fast path for encoding a whole
 * column via array_agg(geom). NULL elements yield NULL UUIDs (position
 * preserved); a non-POINT element raises an error, as h9_encode does. The
 * uuid[] result feeds directly into h9_adaptive(uuid[], …). */

extern "C" {
static Datum h9_encode_many_common(PG_FUNCTION_ARGS, bool sphere) {
	if (PG_ARGISNULL(0)) PG_RETURN_NULL();

	ArrayType *arr      = PG_GETARG_ARRAYTYPE_P(0);
	Oid        elemtype = ARR_ELEMTYPE(arr);
	int16      typlen;
	bool       typbyval;
	char       typalign;
	get_typlenbyvalalign(elemtype, &typlen, &typbyval, &typalign);

	Datum *elems;
	bool  *nulls;
	int    n;
	deconstruct_array(arr, elemtype, typlen, typbyval, typalign,
	                  &elems, &nulls, &n);

	if (n == 0)
		PG_RETURN_ARRAYTYPE_P(construct_empty_array(UUIDOID));

	/* Gather coordinates; NULL elements get a placeholder (0,0) that is
	 * encoded but discarded — the output null mask drops it. */
	double *lon = (double *) palloc((size_t)n * sizeof(double));
	double *lat = (double *) palloc((size_t)n * sizeof(double));
	for (int i = 0; i < n; i++) {
		if (nulls[i]) { lon[i] = 0.0; lat[i] = 0.0; continue; }
		GSERIALIZED *g   = (GSERIALIZED *) PG_DETOAST_DATUM(elems[i]);
		const char  *err = NULL;
		if (!geom_to_lonlat(g, &lon[i], &lat[i], &err))
			ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
							errmsg("h9_encode_many: %s (element %d)", err, i + 1)));
	}

	uint8_t *out = (uint8_t *) palloc((size_t)n * UUID_LEN);
	if (sphere) hex9_encode_many_sphere(lon, lat, (size_t)n, out);
	else        hex9_encode_many(lon, lat, (size_t)n, out);
	pfree(lon);
	pfree(lat);

	Datum *res = (Datum *) palloc((size_t)n * sizeof(Datum));
	for (int i = 0; i < n; i++)
		res[i] = nulls[i] ? (Datum) 0
		                  : UUIDPGetDatum(make_pg_uuid(out + (size_t)i * UUID_LEN));

	int dims[1] = { n };
	int lbs[1]  = { 1 };
	ArrayType *result = construct_md_array(res, nulls, 1, dims, lbs,
	                                       UUIDOID, UUID_LEN, false, TYPALIGN_CHAR);
	PG_RETURN_ARRAYTYPE_P(result);
}

PG_FUNCTION_INFO_V1(h9_encode_many);
Datum h9_encode_many(PG_FUNCTION_ARGS) {
	return h9_encode_many_common(fcinfo, false);
}

PG_FUNCTION_INFO_V1(h9_encode_many_sphere);
Datum h9_encode_many_sphere(PG_FUNCTION_ARGS) {
	return h9_encode_many_common(fcinfo, true);
}
} /* extern "C" */

/* ── h9_decode(uuid) → geometry ─────────────────────────────────────────── */
/*																			 */
/*	 Decode a UUID back to a POINT geometry (SRID 4326).						*/
/*	 The UUID is self-contained: no companion byte is needed.				 */
/*	 Bin UUIDs decode to the cell's exact geographic centroid (identity		 */
/*	 path — coherent with h9_grid/h9_cell for every encoding flavour);		 */
/*	 full UUIDs keep the beam backward walk (representative point).			*/

extern "C" {
static Datum h9_decode_common(PG_FUNCTION_ARGS, bool sphere) {
	if (PG_ARGISNULL(0)) PG_RETURN_NULL();

	pg_uuid_t *u = PG_GETARG_UUID_P(0);

	double lon, lat;
	const int rc = sphere ? hex9_decode_sphere(u->data, &lon, &lat)
	                      : hex9_decode(u->data, &lon, &lat);
	if (rc != 0)
		ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
						errmsg("h9_decode: could not decode UUID")));

	GSERIALIZED *g = lonlat_to_geom(lon, lat, sphere ? 0 : 4326);
	PG_RETURN_POINTER(g);
}

PG_FUNCTION_INFO_V1(h9_decode);
Datum h9_decode(PG_FUNCTION_ARGS) { return h9_decode_common(fcinfo, false); }

PG_FUNCTION_INFO_V1(h9_decode_sphere);
Datum h9_decode_sphere(PG_FUNCTION_ARGS) { return h9_decode_common(fcinfo, true); }
} /* extern "C" */

/* ── h9_bin(uuid, integer) → uuid ───────────────────────────────────────── */
/*																			 */
/*	 Returns the CANONICAL bin-key UUID at the given layer (0..h9_lmax()) — the		*/
/*	 identity coarsening h9_grid enumerates with, so the cell geometrically	*/
/*	 contains the point and h9_bin == h9_kdisk(.,0) == h9_grid.h9_bin (JOIN	*/
/*	 your bins straight to the grid). At a split-hex (6/7/8) leaf it resolves */
/*	 to the canonical mode-0 parent, not the address's mode-1 sibling.		*/
/*	 All points in the same cell at layer L produce the same bin UUID.		*/
/*	 Declared IMMUTABLE and STRICT in SQL — safe for GENERATED STORED cols	*/
/*	 and functional indexes. (Bin VALUES changed for ~20% of split-hex cells */
/*	 vs the pre-canonical h9_bin — rebuild any persisted bin column/index.)	*/

extern "C" {
PG_FUNCTION_INFO_V1(h9_bin);
Datum h9_bin(PG_FUNCTION_ARGS) {
	if (PG_ARGISNULL(0) || PG_ARGISNULL(1)) PG_RETURN_NULL();

	pg_uuid_t *u	 = PG_GETARG_UUID_P(0);
	int32		layer = PG_GETARG_INT32(1);

	if (layer < 0 || layer > hex9_lmax())
		ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
						errmsg("h9_bin: layer must be 0..%d, got %d", hex9_lmax(), layer)));

	uint8_t out_bytes[UUID_LEN];
	if (hex9_bin(u->data, layer, out_bytes) != 0)
		ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
						errmsg("h9_bin: could not bin UUID at layer %d", layer)));

	PG_RETURN_UUID_P(make_pg_uuid(out_bytes));
}
} /* extern "C" */

/* ── h9_cell(uuid, layer, densify) → geometry ─────────────────────────────── */
/*                                                                             */
/*   Returns the hexagonal cell polygon (SRID 4326) for the H9 UUID at the    */
/*   given layer (1..h9_lmax()). The third arg `densify` is an optional non-negative  */
/*   offset (default 0 via SQL): each of the 6 hex edges is subdivided into    */
/*   3^densify segments. Output ring size: 6·3^densify + 1 points.             */
/*                                                                             */
/*       densify = 0 → 7 points (corners only — same as old h9_cell).          */
/*       densify = 1 → 19 points; densify = 2 → 55; densify = 3 → 163.         */
/*                                                                             */
/*   Use case for densify > 0: rendering large hexes (layer ≤ 3) where         */
/*   straight-line edges in (lon, lat) visibly diverge from the underlying     */
/*   great-circle / lattice geometry. For layer ≥ 5 it's rarely needed.        */
/*                                                                             */
/*   The ring is built by the core's identity path (h9_cell_geom.h via         */
/*   hex9_cell_ring) — the exact construction the grid enumerator uses, so     */
/*   h9_cell == h9_grid for every UUID flavour.                                */

extern "C" {
static Datum h9_cell_common(PG_FUNCTION_ARGS, bool sphere) {
	if (PG_ARGISNULL(0) || PG_ARGISNULL(1) || PG_ARGISNULL(2)) PG_RETURN_NULL();
	pg_uuid_t *u       = PG_GETARG_UUID_P(0);
	int32      layer   = PG_GETARG_INT32(1);
	int32      densify = PG_GETARG_INT32(2);
	if (layer < 1 || layer > hex9_lmax())
		ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
						errmsg("h9_cell: layer must be 1..%d, got %d", hex9_lmax(), layer)));
	if (densify < 0)
		ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
						errmsg("h9_cell: densify must be >= 0, got %d", densify)));
	if (layer + densify > hex9_lmax())
		ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
						errmsg("h9_cell: layer + densify must be <= %d "
						       "(layer=%d, densify=%d -> %d)",
						       hex9_lmax(), layer, densify, layer + densify)));
	/* Hard cap to keep ring sizes sane (6·3^densify+1 = 118099 at densify=9). */
	if (densify > 9)
		ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
						errmsg("h9_cell: densify must be <= 9, got %d", densify)));

	const int n_ring = hex9_ring_npoints(densify);
	double *lonlat = (double *) palloc((size_t)n_ring * 2 * sizeof(double));
	const int got = sphere
	    ? hex9_cell_ring_sphere(u->data, layer, densify, lonlat, n_ring)
	    : hex9_cell_ring(u->data, layer, densify, lonlat, n_ring);
	if (got != n_ring)
		ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
						errmsg("h9_cell: not a valid H9 UUID at layer %d", layer)));
	GSERIALIZED *g = h9_polygon_from_ring(lonlat, n_ring, sphere ? 0 : 4326);
	pfree(lonlat);
	PG_RETURN_POINTER(g);
}

PG_FUNCTION_INFO_V1(h9_cell);
Datum h9_cell(PG_FUNCTION_ARGS) { return h9_cell_common(fcinfo, false); }

PG_FUNCTION_INFO_V1(h9_cell_sphere);
Datum h9_cell_sphere(PG_FUNCTION_ARGS) { return h9_cell_common(fcinfo, true); }
} /* extern "C" */

/* ── h9_label(uuid, integer) → text ─────────────────────────────────────── */
/*																			 */
/*	 Returns the body nibbles at the given layer as a plain text label.		 */
/*	 Example: h9_label(h9_bin(h9_encode(pt), 8), 8) → '478232778'			 */

extern "C" {
PG_FUNCTION_INFO_V1(h9_label);
Datum h9_label(PG_FUNCTION_ARGS) {
	if (PG_ARGISNULL(0) || PG_ARGISNULL(1)) PG_RETURN_NULL();

	pg_uuid_t *u	 = PG_GETARG_UUID_P(0);
	int32		layer = PG_GETARG_INT32(1);

	if (layer < 0 || layer > hex9_lmax())
		ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
						errmsg("h9_label: layer must be 0..%d, got %d", hex9_lmax(), layer)));

	char buf[40];
	int	len = hex9_label(u->data, layer, buf, sizeof(buf));
	if (len < 0)
		ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
						errmsg("h9_label: could not label UUID at layer %d", layer)));
	PG_RETURN_TEXT_P(cstring_to_text_with_len(buf, len));
}
} /* extern "C" */

/* ── h9_label_key(uuid, integer) → text ─────────────────────────────────── */
/*																			 */
/*	 Returns the label with the key_tail nibble appended as '.X'.			 */
/*	 Example: h9_label_key(h9_encode(pt), 8) → '478232778.9'				*/

extern "C" {
PG_FUNCTION_INFO_V1(h9_label_key);
Datum h9_label_key(PG_FUNCTION_ARGS) {
	if (PG_ARGISNULL(0) || PG_ARGISNULL(1)) PG_RETURN_NULL();

	pg_uuid_t *u	 = PG_GETARG_UUID_P(0);
	int32		layer = PG_GETARG_INT32(1);

	if (layer < 0 || layer > hex9_lmax())
		ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
						errmsg("h9_label_key: layer must be 0..%d, got %d", hex9_lmax(), layer)));

	char buf[40];
	int	len = hex9_label_key(u->data, layer, buf, sizeof(buf));
	if (len < 0)
		ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
						errmsg("h9_label_key: could not label UUID at layer %d", layer)));
	PG_RETURN_TEXT_P(cstring_to_text_with_len(buf, len));
}
} /* extern "C" */

/* ── h9_grid(bounds, layer, densify) → TABLE(hex9 uuid, geom, centroid) ─── */
/*																			 */
/*	 Enumerates every H9 cell at layer L whose centre lies within bounds.	 */
/*																			 */
/*	 The core (hex9_grid_create) enumerates by bounding box — integer-UV	   */
/*	 supercell BFS with strict octant-region pruning + v4/v5 seam			   */
/*	 reflection, deduplicated. The glue then filters each cell's centroid	   */
/*	 against the actual bounds geometry (lwpoly_contains_point) and builds   */
/*	 the polygon ring on demand per row.									   */
/*																			 */
/*	 Safety cap: hex9.grid_max_cells GUC.  Default 708 588 = 12 × 9⁵        */
/*	 (exactly one global L5 mesh — generous for ordinary local queries,	   */
/*	 lets a coarse global pass without override).							 */

/* Point-in-polygon test using liblwgeom (no GEOS symbols needed).
 * Handles POLYGON and MULTIPOLYGON; passes through for other types. */
static bool h9_lwgeom_contains_pt(const LWGEOM *geom, const POINT2D *pt) {
	switch (geom->type) {
	case POLYGONTYPE:
		/* LW_INSIDE=1, LW_BOUNDARY=0, LW_OUTSIDE=-1.
		 * '>= LW_FALSE(=0)' passes INSIDE and BOUNDARY, rejects OUTSIDE. */
		return lwpoly_contains_point((const LWPOLY *)geom, pt) >= LW_FALSE;
	case MULTIPOLYGONTYPE:
	case COLLECTIONTYPE: {
		const LWCOLLECTION *col = (const LWCOLLECTION *)geom;
		for (uint32_t i = 0; i < col->ngeoms; i++)
			if (h9_lwgeom_contains_pt(col->geoms[i], pt)) return true;
		return false;
	}
	default:
		return true;	 /* non-polygon bounds: let through, bbox already filtered */
	}
}

struct H9GridState {
	int          layer;
	int          densify;        /* per-row polygon densify offset (0 = legacy) */
	int32_t      srid;           /* 4326, or 0 for the sphere datum */
	int          count;
	int          idx;
	hex9_grid   *grid;           /* core handle — freed by the reset callback */
	double      *ringbuf;        /* palloc'd, 2 * ring points doubles */
	int          n_ring;
	LWGEOM      *bounds_lwg;
};

/* The hex9_grid handle is heap-allocated by the core (not palloc), so tie its
 * lifetime to the SRF's multi_call_memory_ctx: the callback fires on context
 * reset/delete — normal completion AND early abort (e.g. LIMIT). */
static void h9_grid_state_release(void *arg) {
	hex9_grid_destroy((hex9_grid *)arg);
}

extern "C" {
static Datum h9_grid_common(PG_FUNCTION_ARGS, bool sphere) {
	FuncCallContext *funcctx;
	H9GridState     *state;

	if (SRF_IS_FIRSTCALL()) {
		MemoryContext oldctx;
		funcctx = SRF_FIRSTCALL_INIT();
		oldctx  = MemoryContextSwitchTo(funcctx->multi_call_memory_ctx);

		/* Argument order: (bounds geometry, layer integer, densify integer).
		 * densify always present at the C level — SQL CREATE FUNCTION uses
		 * DEFAULT 0 so callers can write either h9_grid(g, L) or
		 * h9_grid(g, L, k). */
		GSERIALIZED *bounds_raw = PG_GETARG_GSERIALIZED_P(0);
		int32        layer      = PG_GETARG_INT32(1);
		int32        densify    = PG_GETARG_INT32(2);

		if (layer < 1 || layer > hex9_lmax())
			ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
							errmsg("h9_grid: layer must be 1..%d, got %d", hex9_lmax(), layer)));
		if (densify < 0)
			ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
							errmsg("h9_grid: densify must be >= 0, got %d", densify)));
		if (layer + densify > hex9_lmax())
			ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
							errmsg("h9_grid: layer + densify must be <= %d "
							       "(layer=%d, densify=%d -> %d)",
							       hex9_lmax(), layer, densify, layer + densify)));
		if (densify > 9)
			ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
							errmsg("h9_grid: densify must be <= 9, got %d", densify)));

		GSERIALIZED *bounds = (GSERIALIZED *) palloc(VARSIZE(bounds_raw));
		memcpy(bounds, bounds_raw, VARSIZE(bounds_raw));

		/* Get the bbox from the LWGEOM's actual float8 coords, NOT the
		 * cached float32 header bbox.  The cached bbox has ULP ~ |coord|·2^-23,
		 * which at lat 51.5° is ~6 µdeg — way larger than micro-bbox queries.
		 * E.g., 2e-9° span at L24 rounds OUTWARD to ~6 µdeg, inflating the
		 * area estimate by ~3000× and triggering false overflow rejections
		 * (saw "estimated 169048 cells at layer 24" for an actual-111-cell bbox). */
		LWGEOM *bounds_lwg = lwgeom_from_gserialized(bounds);
		GBOX   gbox;
		if (lwgeom_calculate_gbox(bounds_lwg, &gbox) != LW_SUCCESS)
			ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
							errmsg("h9_grid: could not determine bounds extent")));

		/* Safety cap: spherical-zone area for the bbox over mean H9 cell
		 * area at this layer, plus a perimeter term for partial-coverage
		 * cells at the bbox edge.  Checked BEFORE enumeration so a runaway
		 * request fails fast instead of materialising millions of cells
		 * (the core re-checks the actual count against the same cap).
		 *
		 * The zone formula  A = R²·(sinφ₂ − sinφ₁)·(λ₂ − λ₁)  is exact for
		 * spherical zones and reduces to the flat formula for tiny bboxes.
		 *
		 * Verified against actual CLI counts (no change at small bboxes):
		 *     L23 micro-bbox (2e-9° × 3e-9°): predicts ≈22, actual 10
		 *     L24 same bbox:                  predicts ≈124, actual 111
		 *     L13 Greenwich (0.018° × 0.011°): predicts ≈115k, actual 109k
		 *     L5  global (-180..180,-89..89):  predicts ≈712k, actual 708,588
		 *                                      (was 1,119,576 under flat). */
		const double EARTH_R_M     = 6378137.0;
		const double EARTH_AREA_M2 = 5.10072e14;
		const double D2R           = M_PI / 180.0;
		const double sin_lat_hi    = sin(gbox.ymax * D2R);
		const double sin_lat_lo    = sin(gbox.ymin * D2R);
		const double dlon_rad      = (gbox.xmax - gbox.xmin) * D2R;
		const double bbox_area_m2  = EARTH_R_M * EARTH_R_M
		                             * fabs(sin_lat_hi - sin_lat_lo)
		                             * fabs(dlon_rad);
		const double cell_area_m2  = EARTH_AREA_M2 / (12.0 * pow(9.0, (double)layer));
		const long   est_interior  = (long)(bbox_area_m2 / cell_area_m2);
		const long   est_perimeter = (long)ceil(4.0 * sqrt(bbox_area_m2) / sqrt(cell_area_m2));
		const long   est_cells     = est_interior + est_perimeter + 1;
		if (est_cells > h9_grid_max_cells)
			ereport(ERROR,
				(errcode(ERRCODE_PROGRAM_LIMIT_EXCEEDED),
				 errmsg("h9_grid: estimated %ld cells at layer %d exceeds "
						"limit %d; use a coarser layer or smaller bounds",
						est_cells, layer, h9_grid_max_cells)));

		/* Hand off to the core enumerator.  Same algorithm Python
		 * HexMesh.create_clipped uses; UUIDs derived via the containment-
		 * based xy_regions encoder. */
		char errbuf[256];
		hex9_grid *grid = sphere
		    ? hex9_grid_create_sphere(gbox.xmin, gbox.ymin,
		                              gbox.xmax, gbox.ymax,
		                              layer, densify,
		                              (int64_t)h9_grid_max_cells,
		                              errbuf, sizeof(errbuf))
		    : hex9_grid_create(gbox.xmin, gbox.ymin,
		                       gbox.xmax, gbox.ymax,
		                       layer, densify,
		                       (int64_t)h9_grid_max_cells,
		                       errbuf, sizeof(errbuf));
		if (!grid)
			ereport(ERROR, (errcode(ERRCODE_PROGRAM_LIMIT_EXCEEDED),
							errmsg("h9_grid: %s", errbuf)));

		/* Free the core handle whenever the SRF context goes away. */
		MemoryContextCallback *cb =
			(MemoryContextCallback *) palloc(sizeof(MemoryContextCallback));
		cb->func = h9_grid_state_release;
		cb->arg  = grid;
		MemoryContextRegisterResetCallback(funcctx->multi_call_memory_ctx, cb);

		const int n_ring = hex9_ring_npoints(densify);

		state              = (H9GridState *) palloc(sizeof(H9GridState));
		state->layer       = layer;
		state->densify     = densify;
		state->srid        = sphere ? 0 : 4326;
		state->count       = hex9_grid_count(grid);
		state->idx         = 0;
		state->grid        = grid;
		state->ringbuf     = (double *) palloc((size_t)n_ring * 2 * sizeof(double));
		state->n_ring      = n_ring;
		state->bounds_lwg  = bounds_lwg;
		funcctx->user_fctx = state;

		get_call_result_type(fcinfo, NULL, &funcctx->tuple_desc);
		BlessTupleDesc(funcctx->tuple_desc);

		MemoryContextSwitchTo(oldctx);
	}

	funcctx = SRF_PERCALL_SETUP();
	state   = (H9GridState *) funcctx->user_fctx;

	while (state->idx < state->count) {
		const int i = state->idx++;

		/* Polygon-containment filter on the centroid (matches Python's
		 * matplotlib Path test in create_clipped).  Skipped when bounds is
		 * not a polygon — the bbox prune inside the core already applied. */
		double clon, clat;
		hex9_grid_cell_centroid(state->grid, i, &clon, &clat);
		POINT2D pt = { clon, clat };
		if (!h9_lwgeom_contains_pt(state->bounds_lwg, &pt))
			continue;

		uint8_t id_uuid[UUID_LEN], bin_uuid[UUID_LEN];
		hex9_grid_cell_id(state->grid, i, id_uuid);
		hex9_grid_cell_uuid(state->grid, i, bin_uuid);
		if (hex9_grid_cell_ring(state->grid, i, state->densify,
		                        state->ringbuf, state->n_ring) != state->n_ring)
			ereport(ERROR, (errcode(ERRCODE_INTERNAL_ERROR),
							errmsg("h9_grid: could not build ring for cell %d", i)));

		Datum values[4];
		bool  isnull[4] = {false, false, false, false};
		/* h9_id: full reversible identity; h9_bin: layer-scoped grouping key. */
		values[0] = UUIDPGetDatum(make_pg_uuid(id_uuid));
		values[1] = UUIDPGetDatum(make_pg_uuid(bin_uuid));
		values[2] = PointerGetDatum(h9_polygon_from_ring(state->ringbuf, state->n_ring,
		                                                 state->srid));
		/* Cell centroid (POINT, in the handle's datum SRID) — computed by the
		 * core in the cell's own frame, not the mean of the rendered polygon's
		 * vertices. */
		values[3] = PointerGetDatum(lonlat_to_geom(clon, clat, state->srid));
		HeapTuple tup = heap_form_tuple(funcctx->tuple_desc, values, isnull);
		SRF_RETURN_NEXT(funcctx, HeapTupleGetDatum(tup));
	}

	if (state->bounds_lwg) lwgeom_free(state->bounds_lwg);
	SRF_RETURN_DONE(funcctx);
}

PG_FUNCTION_INFO_V1(h9_grid);
Datum h9_grid(PG_FUNCTION_ARGS) { return h9_grid_common(fcinfo, false); }

PG_FUNCTION_INFO_V1(h9_grid_sphere);
Datum h9_grid_sphere(PG_FUNCTION_ARGS) { return h9_grid_common(fcinfo, true); }
} /* extern "C" */

/* ── Neighbours / k-ring / k-disk (SETOF uuid) ──────────────────────────── */
/*																			 */
/*	 Symbolic adjacency on the H9 mesh — exact integer lattice arithmetic	   */
/*	 in the core (hex9_c.h), no floating point. Input may be a full UUID	   */
/*	 (from h9_encode) or a bin UUID at a layer >= `layer`; output cells are  */
/*	 bin UUIDs at `layer`, sorted. Every cell has 6 neighbours except the	   */
/*	 12 half-hex cells per layer at the octahedron vertices, which have 5.   */

struct H9UuidSetState {
	int64_t  count;
	int64_t  idx;
	uint8_t *uuids;              /* palloc'd, count * 16 bytes */
};

/* Shared per-call iteration for the three set-returning adjacency functions.
 * fill() runs once in the SRF's multi-call context and returns the cell
 * count (its uuids buffer must be palloc'd there too). */
template <typename Fill>
static Datum h9_uuid_set_srf(FunctionCallInfo fcinfo, Fill fill) {
	FuncCallContext *funcctx;

	if (SRF_IS_FIRSTCALL()) {
		funcctx = SRF_FIRSTCALL_INIT();
		MemoryContext oldctx = MemoryContextSwitchTo(funcctx->multi_call_memory_ctx);

		H9UuidSetState *state = (H9UuidSetState *) palloc(sizeof(H9UuidSetState));
		state->idx   = 0;
		state->count = fill(state);
		funcctx->user_fctx = state;

		MemoryContextSwitchTo(oldctx);
	}

	funcctx = SRF_PERCALL_SETUP();
	H9UuidSetState *state = (H9UuidSetState *) funcctx->user_fctx;

	if (state->idx < state->count) {
		const uint8_t *u = state->uuids + 16 * state->idx++;
		SRF_RETURN_NEXT(funcctx, UUIDPGetDatum(make_pg_uuid(u)));
	}
	SRF_RETURN_DONE(funcctx);
}

/* Adjacency input must be a full UUID from h9_encode — bins are layer-scoped
 * keys, not addresses (their key tail cannot carry the meta the resolution
 * needs); give the doctrine error rather than a generic one. */
static void h9_reject_bin_input(const pg_uuid_t *u, const char *fname) {
	if ((u->data[15] >> 4) == 0x0Fu)
		ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
						errmsg("%s: bin UUIDs are layer keys, not addresses — "
						       "pass the full UUID from h9_encode", fname),
						errhint("Re-derive adjacency from the original point; "
						        "output bins are for joining (e.g. to h9_grid), "
						        "not for further traversal.")));
}

extern "C" {
PG_FUNCTION_INFO_V1(h9_neighbors);
Datum h9_neighbors(PG_FUNCTION_ARGS) {
	return h9_uuid_set_srf(fcinfo, [&](H9UuidSetState *state) -> int64_t {
		pg_uuid_t *u     = PG_GETARG_UUID_P(0);
		int32      layer = PG_GETARG_INT32(1);
		if (layer < 1 || layer > hex9_lmax())
			ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
							errmsg("h9_neighbors: layer must be 1..%d, got %d", hex9_lmax(), layer)));
		h9_reject_bin_input(u, "h9_neighbors");
		state->uuids = (uint8_t *) palloc(6 * 16);
		int n = hex9_neighbors(u->data, layer, state->uuids);
		if (n < 0)
			ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
							errmsg("h9_neighbors: not a valid H9 UUID at layer %d", layer)));
		return n;
	});
}
} /* extern "C" */

/* Common body for h9_kring / h9_kdisk: size by the nominal disk count,
 * fill via the core, error on bad args. */
typedef int64_t (*h9_kfill_fn)(const uint8_t uuid[16], int layer, int k,
                               uint8_t *out_uuids, int64_t max_cells);

static int64_t h9_kring_fill(PG_FUNCTION_ARGS, H9UuidSetState *state,
                             const char *fname, h9_kfill_fn fn) {
	pg_uuid_t *u     = PG_GETARG_UUID_P(0);
	int32      layer = PG_GETARG_INT32(1);
	int32      k     = PG_GETARG_INT32(2);
	if (layer < 1 || layer > hex9_lmax())
		ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
						errmsg("%s: layer must be 1..%d, got %d", fname, hex9_lmax(), layer)));
	if (k < 0)
		ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
						errmsg("%s: k must be >= 0, got %d", fname, k)));
	h9_reject_bin_input(u, fname);
	/* Nominal disk size 1 + 3k(k+1) is an upper bound for both ring and
	 * disk. Cap the buffer below MaxAllocSize (1 GB / 16 ≈ 67M cells). */
	const int64_t ncells = hex9_disk_ncells(k);
	if (ncells > 60000000)
		ereport(ERROR, (errcode(ERRCODE_PROGRAM_LIMIT_EXCEEDED),
						errmsg("%s: k=%d implies up to %lld cells; use a smaller k",
						       fname, k, (long long)ncells)));
	state->uuids = (uint8_t *) palloc((Size)ncells * 16);
	int64_t n = fn(u->data, layer, k, state->uuids, ncells);
	if (n < 0)
		ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
						errmsg("%s: not a valid H9 UUID at layer %d", fname, layer)));
	return n;
}

extern "C" {
PG_FUNCTION_INFO_V1(h9_kring);
Datum h9_kring(PG_FUNCTION_ARGS) {
	return h9_uuid_set_srf(fcinfo, [&](H9UuidSetState *state) {
		return h9_kring_fill(fcinfo, state, "h9_kring", hex9_k_ring);
	});
}

PG_FUNCTION_INFO_V1(h9_kdisk);
Datum h9_kdisk(PG_FUNCTION_ARGS) {
	return h9_uuid_set_srf(fcinfo, [&](H9UuidSetState *state) {
		return h9_kring_fill(fcinfo, state, "h9_kdisk", hex9_k_disk);
	});
}
} /* extern "C" */

/* ── h9_parse_label(text) → uuid ────────────────────────────────────────── */
/*																			 */
/*	 Recover the canonical bin UUID from a label. Accepts both forms:		   */
/*	 bare labels from h9_label (canonical-only — unique per layer, verified  */
/*	 by re-encode) and keyed labels from h9_label_key (flavour-blind). The   */
/*	 layer is implicit in the label length (body chars = layer + 1).		   */

extern "C" {
PG_FUNCTION_INFO_V1(h9_parse_label);
Datum h9_parse_label(PG_FUNCTION_ARGS) {
	if (PG_ARGISNULL(0)) PG_RETURN_NULL();

	char *label = text_to_cstring(PG_GETARG_TEXT_PP(0));
	uint8_t uuid[UUID_LEN];
	if (hex9_parse_label(label, uuid) < 0)
		ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
						errmsg("h9_parse_label: not a valid H9 label: \"%s\"",
						       label)));
	PG_RETURN_UUID_P(make_pg_uuid(uuid));
}
} /* extern "C" */

/* ── h9_label_centroid(text) → geometry ─────────────────────────────────── */
/*																			 */
/*	 Geographic centroid (POINT, SRID 4326) of the labelled cell — same	   */
/*	 convention as h9_grid's centroid column and h9_decode on a bin UUID.	   */

extern "C" {
PG_FUNCTION_INFO_V1(h9_label_centroid);
Datum h9_label_centroid(PG_FUNCTION_ARGS) {
	if (PG_ARGISNULL(0)) PG_RETURN_NULL();

	char *label = text_to_cstring(PG_GETARG_TEXT_PP(0));
	double lon, lat;
	if (hex9_label_centroid(label, &lon, &lat) != 0)
		ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
						errmsg("h9_label_centroid: not a valid H9 label: \"%s\"",
						       label)));
	PG_RETURN_POINTER(lonlat_to_geom(lon, lat));
}
} /* extern "C" */

/* ── h9_common_ancestor(uuid[], integer) → (label, hex9, layer) ─────────── */
/*																			 */
/*	 Deepest common ancestor (in the ADDRESS hierarchy) of the given cells,  */
/*	 all treated at `layer`. Returns its label, bin UUID, and layer; NULL	   */
/*	 when the cells span L0 hexes (no common ancestor). The usual DGGS		   */
/*	 caveat: descent containment, not exact geometric containment.			   */

extern "C" {
PG_FUNCTION_INFO_V1(h9_common_ancestor);
Datum h9_common_ancestor(PG_FUNCTION_ARGS) {
	if (PG_ARGISNULL(0) || PG_ARGISNULL(1)) PG_RETURN_NULL();

	ArrayType *arr   = PG_GETARG_ARRAYTYPE_P(0);
	int32      layer = PG_GETARG_INT32(1);

	if (layer < 0 || layer > hex9_lmax())
		ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
						errmsg("h9_common_ancestor: layer must be 0..%d, got %d",
						       hex9_lmax(), layer)));

	Datum *elems;
	bool  *nulls;
	int    n;
	/* uuid: typlen 16, by-reference, char alignment. */
	deconstruct_array(arr, UUIDOID, UUID_LEN, false, TYPALIGN_CHAR,
	                  &elems, &nulls, &n);
	if (n == 0)
		ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
						errmsg("h9_common_ancestor: uuid array must not be empty")));

	uint8_t *uuids = (uint8_t *) palloc((size_t)n * UUID_LEN);
	for (int i = 0; i < n; i++) {
		if (nulls[i])
			ereport(ERROR, (errcode(ERRCODE_NULL_VALUE_NOT_ALLOWED),
							errmsg("h9_common_ancestor: uuid array must not "
							       "contain NULLs")));
		memcpy(uuids + (size_t)i * UUID_LEN, DatumGetUUIDP(elems[i])->data,
		       UUID_LEN);
	}

	char    label[40];
	uint8_t anc_uuid[UUID_LEN];
	int anc_layer = hex9_common_ancestor(uuids, (size_t)n, layer,
	                                     label, sizeof(label), anc_uuid);
	pfree(uuids);
	if (anc_layer < 0)
		PG_RETURN_NULL();   /* cells span L0 hexes — no common ancestor */

	TupleDesc tupdesc;
	if (get_call_result_type(fcinfo, NULL, &tupdesc) != TYPEFUNC_COMPOSITE)
		ereport(ERROR, (errcode(ERRCODE_FEATURE_NOT_SUPPORTED),
						errmsg("h9_common_ancestor: composite return type "
						       "expected")));
	BlessTupleDesc(tupdesc);

	Datum values[3];
	bool  isnull[3] = {false, false, false};
	values[0] = CStringGetTextDatum(label);
	values[1] = UUIDPGetDatum(make_pg_uuid(anc_uuid));
	values[2] = Int32GetDatum(anc_layer);
	PG_RETURN_DATUM(HeapTupleGetDatum(heap_form_tuple(tupdesc, values, isnull)));
}
} /* extern "C" */

/* ── h9_adaptive(uuids, weights, …) → TABLE(hex9, layer, value, npoints) ── */
/*																			 */
/*	 Population-digest multi-layer grid: aggregate weighted addresses into   */
/*	 a mixed-layer cell set whose cell values respect a population ceiling,  */
/*	 by bottom-up digestion (see hex9_c.h for the full model). Dense areas   */
/*	 resolve into fine cells, sparse areas aggregate upward; emitted values  */
/*	 sum exactly to the input weight total. `weights` may be NULL (all 1).   */
/*	 Cells of different layers overlap geometrically — the set partitions	   */
/*	 the SAMPLE, not the surface.											   */
/*																			 */
/*	 INPUT IS FULL UUIDs (h9_encode) — bin input is rejected (the digest	   */
/*	 re-bins across layers, guaranteed only from the full uuid; Ben's		   */
/*	 ruling 2026-06-12: addresses, not coordinates, are the digest's		   */
/*	 input — start from geometry via array_agg(h9_encode(geom))).			   */

struct H9AdaptiveCell {
	uint8_t uuid[UUID_LEN];   /* layer-scoped bin key -> hex9 column */
	uint8_t full[UUID_LEN];   /* representative full uuid -> geometry render */
	int     layer;
	double  value;
	int64_t npoints;
};

struct H9AdaptiveState {
	int             count;
	int             idx;
	bool            sphere;   /* sphere twin: geom SRID 0, density /steradian */
	H9AdaptiveCell *cells;
};

extern "C" {
static Datum h9_adaptive_common(PG_FUNCTION_ARGS, bool sphere) {
	FuncCallContext *funcctx;

	if (SRF_IS_FIRSTCALL()) {
		funcctx = SRF_FIRSTCALL_INIT();
		MemoryContext oldctx = MemoryContextSwitchTo(funcctx->multi_call_memory_ctx);

		if (PG_ARGISNULL(0) || PG_ARGISNULL(2) || PG_ARGISNULL(3) ||
		    PG_ARGISNULL(4) || PG_ARGISNULL(5))
			ereport(ERROR, (errcode(ERRCODE_NULL_VALUE_NOT_ALLOWED),
							errmsg("h9_adaptive: only weights may be NULL")));

		ArrayType *uu_arr    = PG_GETARG_ARRAYTYPE_P(0);
		int32      min_layer = PG_GETARG_INT32(2);
		int32      max_layer = PG_GETARG_INT32(3);
		double     ceiling_  = PG_GETARG_FLOAT8(4);
		double     floor_    = PG_GETARG_FLOAT8(5);

		Datum *uu_elems;
		bool  *uu_nulls;
		int    n;
		deconstruct_array(uu_arr, UUIDOID, UUID_LEN, false, TYPALIGN_CHAR,
		                  &uu_elems, &uu_nulls, &n);
		if (n == 0)
			ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
							errmsg("h9_adaptive: uuids array must not be empty")));

		uint8_t *uuids = (uint8_t *) palloc((size_t)n * UUID_LEN);
		for (int i = 0; i < n; i++) {
			if (uu_nulls[i])
				ereport(ERROR, (errcode(ERRCODE_NULL_VALUE_NOT_ALLOWED),
								errmsg("h9_adaptive: uuids array must not "
								       "contain NULLs")));
			memcpy(uuids + (size_t)i * UUID_LEN,
			       DatumGetUUIDP(uu_elems[i])->data, UUID_LEN);
		}

		/* Optional weights: float8[] of the same cardinality, no NULLs. */
		double *weight = NULL;
		if (!PG_ARGISNULL(1)) {
			ArrayType *w_arr = PG_GETARG_ARRAYTYPE_P(1);
			Datum *w_elems;
			bool  *w_nulls;
			int    wn;
			deconstruct_array(w_arr, FLOAT8OID, sizeof(float8), FLOAT8PASSBYVAL,
			                  TYPALIGN_DOUBLE, &w_elems, &w_nulls, &wn);
			if (wn != n)
				ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
								errmsg("h9_adaptive: weights has %d elements, "
								       "points has %d", wn, n)));
			weight = (double *) palloc((size_t)n * sizeof(double));
			for (int i = 0; i < n; i++) {
				if (w_nulls[i])
					ereport(ERROR, (errcode(ERRCODE_NULL_VALUE_NOT_ALLOWED),
									errmsg("h9_adaptive: weights array must "
									       "not contain NULLs")));
				weight[i] = DatumGetFloat8(w_elems[i]);
			}
		}

		char errbuf[256];
		hex9_adaptive *a = hex9_adaptive_create(uuids, weight, (size_t)n,
		                                        min_layer, max_layer,
		                                        ceiling_, floor_,
		                                        errbuf, sizeof(errbuf));
		if (!a)
			ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
							errmsg("h9_adaptive: %s", errbuf)));

		/* Copy the digest out and free the core handle immediately — the
		 * SRF then iterates plain palloc'd memory. */
		const int count = hex9_adaptive_count(a);
		H9AdaptiveCell *cells = (H9AdaptiveCell *)
			palloc((size_t)Max(count, 1) * sizeof(H9AdaptiveCell));
		for (int i = 0; i < count; i++) {
			hex9_adaptive_cell(a, i, cells[i].uuid, &cells[i].layer,
			                   &cells[i].value, &cells[i].npoints);
			hex9_adaptive_cell_full(a, i, cells[i].full);
		}
		hex9_adaptive_destroy(a);
		pfree(uuids);
		if (weight) pfree(weight);

		H9AdaptiveState *state = (H9AdaptiveState *) palloc(sizeof(H9AdaptiveState));
		state->count  = count;
		state->idx    = 0;
		state->sphere = sphere;
		state->cells  = cells;
		funcctx->user_fctx = state;

		get_call_result_type(fcinfo, NULL, &funcctx->tuple_desc);
		BlessTupleDesc(funcctx->tuple_desc);

		MemoryContextSwitchTo(oldctx);
	}

	funcctx = SRF_PERCALL_SETUP();
	H9AdaptiveState *state = (H9AdaptiveState *) funcctx->user_fctx;

	if (state->idx < state->count) {
		const H9AdaptiveCell &c = state->cells[state->idx++];
		Datum values[7];
		bool  isnull[7] = {false, false, false, false, false, false, false};
		values[0] = UUIDPGetDatum(make_pg_uuid(c.uuid));
		values[1] = Int32GetDatum(c.layer);
		values[2] = Float8GetDatum(c.value);
		values[3] = Int64GetDatum(c.npoints);
		/* density: exact for the digest — cells are equal-area per layer,
		 * 12·9^layer cells tile the surface. WGS84: persons/km² over Earth's
		 * 510065622 km². Sphere twin: value per STERADIAN — the unit sphere's
		 * cell area (4π/(12·9^L) sr) is intrinsic, so no body radius is
		 * needed (the sphere datum carries none); per-km² on a body is
		 * density × 4π / body_area, caller-side. grade (log₉ graduation;
		 * +1 ⇒ 9× denser) is datum-free. Canonical so callers don't
		 * re-derive (and re-invert) the formula. See the API README for the
		 * interpretation caveats (npoints=1 = point-mass reading, source
		 * quantisation, floor/ceiling binding). */
		values[4] = Float8GetDatum(c.value * 12.0 * pow(9.0, (double)c.layer)
		                           / (state->sphere ? 4.0 * M_PI : 510065622.0));
		if (c.value > 0.0)
			values[5] = Float8GetDatum((double)c.layer + log(c.value) / log(9.0));
		else
			isnull[5] = true;   /* grade undefined for non-positive value */
		/* hexagon polygon, corners only — same identity render as h9_cell,
		 * so no companion h9_grid is needed to display the digest (and the
		 * displayable layer is no longer bounded by grid enumeration) */
		{
			double lonlat[2 * 7];
			const int got = state->sphere
			    ? hex9_cell_ring_sphere(c.full, c.layer, 0, lonlat, 7)
			    : hex9_cell_ring(c.full, c.layer, 0, lonlat, 7);
			if (got != 7)
				ereport(ERROR, (errcode(ERRCODE_INTERNAL_ERROR),
								errmsg("h9_adaptive: cell ring failed at layer %d",
								       c.layer)));
			values[6] = PointerGetDatum(h9_polygon_from_ring(lonlat, 7,
			                                                 state->sphere ? 0 : 4326));
		}
		HeapTuple tup = heap_form_tuple(funcctx->tuple_desc, values, isnull);
		SRF_RETURN_NEXT(funcctx, HeapTupleGetDatum(tup));
	}
	SRF_RETURN_DONE(funcctx);
}

PG_FUNCTION_INFO_V1(h9_adaptive);
Datum h9_adaptive(PG_FUNCTION_ARGS) { return h9_adaptive_common(fcinfo, false); }

PG_FUNCTION_INFO_V1(h9_adaptive_sphere);
Datum h9_adaptive_sphere(PG_FUNCTION_ARGS) { return h9_adaptive_common(fcinfo, true); }
} /* extern "C" */

/* ── Hamiltonian curve addressing (h9_curve family) ───────────────────────
 *
 * The H9 space-filling curve: every layer-L cell visited once, consecutive
 * indices edge-adjacent, and the order REFINES (a cell's 9 lineage children
 * occupy indices index*9 .. index*9+8). Two jobs: SORTING (curve order is
 * locality order — ORDER BY h9_curve(bin), CLUSTER, tile streaming) and
 * GENERATION (h9_curve_cells enumerates a cell's descendants in locality
 * order). The packed curve-uuid (nibble 0 = 0xC) is the sortable form:
 * native uuid btree order at a fixed layer IS curve order, and curve-uuids
 * coexist with h9-uuids in one column (they sort after all of them). In
 * mixed-layer collections the 0xF padding makes an ancestor sort AFTER its
 * descendants (post-order); group by h9_curve_layer when that matters.
 * Unlike an h9-uuid body a curve-uuid truncates EXACTLY (h9_curve_bin).
 * See hex9_c.h and core/h9_curve.h for the transducer doctrine. */

extern "C" {
/* Canonical cell parent / ancestor (mode-0 convention) — the CELL-level
 * roll-up of a bin, distinct from h9_bin's POINT question. h9_cell_parent
 * is the curve's generation step: the curve/lineage tree is ITERATED
 * one-generation parents. h9_cell_ancestor is the direct deep fold — the
 * two differ on hexagon-band cells beyond one generation. */
PG_FUNCTION_INFO_V1(h9_cell_parent);
Datum h9_cell_parent(PG_FUNCTION_ARGS) {
	pg_uuid_t *u = PG_GETARG_UUID_P(0);
	uint8_t out[UUID_LEN];
	if (hex9_cell_parent(u->data, out) != 0)
		ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
						errmsg("h9_cell_parent: input must be a bin UUID at "
						       "layer >= 1 (L0 cells have no parent)")));
	PG_RETURN_UUID_P(make_pg_uuid(out));
}

PG_FUNCTION_INFO_V1(h9_cell_ancestor);
Datum h9_cell_ancestor(PG_FUNCTION_ARGS) {
	pg_uuid_t *u	 = PG_GETARG_UUID_P(0);
	int32		layer = PG_GETARG_INT32(1);
	uint8_t out[UUID_LEN];
	if (layer < 0 || layer > hex9_lmax())
		ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
						errmsg("h9_cell_ancestor: layer must be 0..%d, got %d",
						       hex9_lmax(), layer)));
	if (hex9_cell_ancestor(u->data, layer, out) != 0)
		ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
						errmsg("h9_cell_ancestor: input must be a bin UUID at "
						       "a layer >= %d", layer)));
	PG_RETURN_UUID_P(make_pg_uuid(out));
}

PG_FUNCTION_INFO_V1(h9_curve);
Datum h9_curve(PG_FUNCTION_ARGS) {
	pg_uuid_t *u = PG_GETARG_UUID_P(0);
	uint8_t out[UUID_LEN];
	if (hex9_curve(u->data, out) != 0)
		ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
						errmsg("h9_curve: not a valid H9 UUID (a T1 transducer "
						       "miss means a non-canonical address)")));
	PG_RETURN_UUID_P(make_pg_uuid(out));
}

PG_FUNCTION_INFO_V1(h9_curve_decode);
Datum h9_curve_decode(PG_FUNCTION_ARGS) {
	pg_uuid_t *u = PG_GETARG_UUID_P(0);
	uint8_t out[UUID_LEN];
	if (hex9_curve_decode(u->data, out) != 0)
		ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
						errmsg("h9_curve_decode: not a valid curve UUID")));
	PG_RETURN_UUID_P(make_pg_uuid(out));
}

PG_FUNCTION_INFO_V1(h9_is_curve);
Datum h9_is_curve(PG_FUNCTION_ARGS) {
	pg_uuid_t *u = PG_GETARG_UUID_P(0);
	PG_RETURN_BOOL(hex9_is_curve(u->data) == 1);
}

PG_FUNCTION_INFO_V1(h9_curve_layer);
Datum h9_curve_layer(PG_FUNCTION_ARGS) {
	pg_uuid_t *u = PG_GETARG_UUID_P(0);
	const int layer = hex9_curve_layer(u->data);
	if (layer < 0)
		ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
						errmsg("h9_curve_layer: not a curve UUID (nibble 0 != c)")));
	PG_RETURN_INT32(layer);
}

PG_FUNCTION_INFO_V1(h9_curve_bin);
Datum h9_curve_bin(PG_FUNCTION_ARGS) {
	pg_uuid_t *u	 = PG_GETARG_UUID_P(0);
	int32		layer = PG_GETARG_INT32(1);
	uint8_t out[UUID_LEN];
	if (layer < 0 || layer > hex9_lmax())
		ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
						errmsg("h9_curve_bin: layer must be 0..%d, got %d",
						       hex9_lmax(), layer)));
	if (hex9_curve_bin(u->data, layer, out) != 0)
		ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
						errmsg("h9_curve_bin: input must be a curve UUID at a "
						       "layer >= %d", layer)));
	PG_RETURN_UUID_P(make_pg_uuid(out));
}

/* numeric because indices reach 12·9^L − 1: beyond bigint above layer 18. */
PG_FUNCTION_INFO_V1(h9_curve_index);
Datum h9_curve_index(PG_FUNCTION_ARGS) {
	pg_uuid_t *u = PG_GETARG_UUID_P(0);
	char buf[64];
	if (hex9_curve_index(u->data, buf, sizeof buf) < 0)
		ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
						errmsg("h9_curve_index: not a valid H9 or curve UUID")));
	return DirectFunctionCall3(numeric_in, CStringGetDatum(buf),
	                           ObjectIdGetDatum(InvalidOid), Int32GetDatum(-1));
}

PG_FUNCTION_INFO_V1(h9_curve_pack);
Datum h9_curve_pack(PG_FUNCTION_ARGS) {
	Datum num	 = PG_GETARG_DATUM(0);
	int32 layer  = PG_GETARG_INT32(1);
	if (layer < 0 || layer > hex9_lmax())
		ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
						errmsg("h9_curve_pack: layer must be 0..%d, got %d",
						       hex9_lmax(), layer)));
	char *s = DatumGetCString(DirectFunctionCall1(numeric_out, num));
	uint8_t out[UUID_LEN];
	if (hex9_curve_pack(s, layer, out) != 0)
		ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
						errmsg("h9_curve_pack: '%s' is not a valid curve index "
						       "at layer %d (integer 0..12*9^layer-1)", s, layer)));
	pfree(s);
	PG_RETURN_UUID_P(make_pg_uuid(out));
}

PG_FUNCTION_INFO_V1(h9_curve_label);
Datum h9_curve_label(PG_FUNCTION_ARGS) {
	pg_uuid_t *u = PG_GETARG_UUID_P(0);
	char buf[64];
	const int len = hex9_curve_label(u->data, buf, sizeof buf);
	if (len < 0)
		ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
						errmsg("h9_curve_label: not a valid H9 or curve UUID")));
	PG_RETURN_TEXT_P(cstring_to_text_with_len(buf, len));
}

PG_FUNCTION_INFO_V1(h9_curve_from_label);
Datum h9_curve_from_label(PG_FUNCTION_ARGS) {
	text *t = PG_GETARG_TEXT_PP(0);
	char *label = text_to_cstring(t);
	uint8_t out[UUID_LEN];
	if (hex9_curve_parse_label(label, out) != 0)
		ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
						errmsg("h9_curve_from_label: '%s' is not a curve label "
						       "('c' + slot hex char + base-9 rank chars)", label)));
	pfree(label);
	PG_RETURN_UUID_P(make_pg_uuid(out));
}
} /* extern "C" */

/* ── Descendant enumeration SRFs ──────────────────────────────────────────
 *
 * h9_curve_cells — every layer-L LINEAGE descendant (the curve tree), in
 * curve order. h9_owned_cells — every layer-L OWNED sub-zone (cells whose
 * h9_cell_ancestor IS the zone; exactly 9^d, the sets partition the
 * layer — the aggregation relation of OGC API-DGGS issue #108),
 * curve-sorted. Both return (h9_bin, h9_curve) rows — the h9 bin is the
 * primary key form (selection over format: curve->h9 is the expensive
 * direction, h9->curve the cheap one); input is an h9 bin/full uuid or a
 * curve-uuid. Capped by the hex9.grid_max_cells GUC like h9_grid. */

struct H9CurveCellsState {
	int64_t  count;
	int64_t  idx;
	uint8_t *bins;               /* count * 16 */
	uint8_t *curves;             /* count * 16 */
};

typedef int64_t (*h9_enum_fn)(const uint8_t uuid[16], int layer,
                              uint8_t *out_bins, uint8_t *out_curves,
                              int64_t max_cells);

static Datum h9_cells_srf(FunctionCallInfo fcinfo, const char *fname,
                          h9_enum_fn fn) {
	FuncCallContext *funcctx;

	if (SRF_IS_FIRSTCALL()) {
		funcctx = SRF_FIRSTCALL_INIT();
		MemoryContext oldctx = MemoryContextSwitchTo(funcctx->multi_call_memory_ctx);

		pg_uuid_t *u	 = PG_GETARG_UUID_P(0);
		int32		layer = PG_GETARG_INT32(1);
		if (layer < 0 || layer > hex9_lmax())
			ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
							errmsg("%s: layer must be 0..%d, got %d",
							       fname, hex9_lmax(), layer)));

		H9CurveCellsState *state = (H9CurveCellsState *) palloc(sizeof *state);
		state->idx = 0;

		/* Size from the cell's own layer; the ABI re-derives and caps. */
		uint8_t curve[UUID_LEN];
		if (hex9_curve(u->data, curve) != 0)
			ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
							errmsg("%s: not a valid H9 or curve UUID", fname)));
		const int from = hex9_curve_layer(curve);
		const int64_t want = hex9_curve_ncells(from, layer);
		if (want < 0)
			ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
							errmsg("%s: layer %d is above the cell's "
							       "own layer %d", fname, layer, from)));
		if (want > (int64_t)h9_grid_max_cells)
			ereport(ERROR, (errcode(ERRCODE_PROGRAM_LIMIT_EXCEEDED),
							errmsg("%s: %lld cells at layer %d exceeds "
							       "limit %d; use a shallower layer or raise "
							       "hex9.grid_max_cells",
							       fname, (long long)want, layer, h9_grid_max_cells)));

		state->bins   = (uint8_t *) palloc((Size)want * 16);
		state->curves = (uint8_t *) palloc((Size)want * 16);
		state->count  = fn(u->data, layer, state->bins, state->curves, want);
		if (state->count < 0)
			ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
							errmsg("%s: enumeration failed at layer %d",
							       fname, layer)));
		funcctx->user_fctx = state;

		get_call_result_type(fcinfo, NULL, &funcctx->tuple_desc);
		BlessTupleDesc(funcctx->tuple_desc);
		MemoryContextSwitchTo(oldctx);
	}

	funcctx = SRF_PERCALL_SETUP();
	H9CurveCellsState *state = (H9CurveCellsState *) funcctx->user_fctx;

	if (state->idx < state->count) {
		const int64_t i = state->idx++;
		Datum values[2];
		bool  isnull[2] = {false, false};
		values[0] = UUIDPGetDatum(make_pg_uuid(state->bins + i * 16));
		values[1] = UUIDPGetDatum(make_pg_uuid(state->curves + i * 16));
		HeapTuple tup = heap_form_tuple(funcctx->tuple_desc, values, isnull);
		SRF_RETURN_NEXT(funcctx, HeapTupleGetDatum(tup));
	}
	SRF_RETURN_DONE(funcctx);
}

extern "C" {
PG_FUNCTION_INFO_V1(h9_curve_cells);
Datum h9_curve_cells(PG_FUNCTION_ARGS) {
	return h9_cells_srf(fcinfo, "h9_curve_cells", hex9_curve_cells);
}

PG_FUNCTION_INFO_V1(h9_owned_cells);
Datum h9_owned_cells(PG_FUNCTION_ARGS) {
	return h9_cells_srf(fcinfo, "h9_owned_cells", hex9_owned_cells);
}

/* h9_cell_children(uuid) → SETOF uuid — the 9 canonical children in
 * curve-rank order (one generation: lineage == ownership). */
PG_FUNCTION_INFO_V1(h9_cell_children);
Datum h9_cell_children(PG_FUNCTION_ARGS) {
	return h9_uuid_set_srf(fcinfo, [&](H9UuidSetState *state) -> int64_t {
		pg_uuid_t *u = PG_GETARG_UUID_P(0);
		state->uuids = (uint8_t *) palloc(9 * 16);
		if (hex9_cell_children(u->data, state->uuids) != 0)
			ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
							errmsg("h9_cell_children: not a valid H9 UUID, or "
							       "the cell is at the deepest layer %d",
							       hex9_lmax())));
		return 9;
	});
}
} /* extern "C" */

/* ═══ E4H: the aperture-4 structural tail (Hex9 2.3.0) ══════════════════════
 * SQL twins over the hex9_e4h_* ABI — see hex9_c.h's E4H doctrine block and
 * libhex9 docs/universality.md: E4H uuids are ADDRESSES minted by the exact
 * classifier, and the h9_* machinery above rejects them (the 0xE marker
 * guard); these are their own operations. */

/* deconstruct a uuid[] argument into a packed 16-byte-row buffer */
static size_t h9_uuid_array_bytes(ArrayType *arr, uint8_t **out) {
	Datum *elems;
	bool  *nulls;
	int    n;
	deconstruct_array(arr, UUIDOID, UUID_LEN, false, TYPALIGN_CHAR,
	                  &elems, &nulls, &n);
	uint8_t *buf = (uint8_t *) palloc((size_t)(n > 0 ? n : 1) * UUID_LEN);
	size_t m = 0;
	for (int i = 0; i < n; i++) {
		if (nulls[i]) continue;
		memcpy(buf + m * UUID_LEN, DatumGetUUIDP(elems[i])->data, UUID_LEN);
		m++;
	}
	*out = buf;
	return m;
}

extern "C" {
static Datum h9e_encode_common(PG_FUNCTION_ARGS, bool sphere) {
	if (PG_ARGISNULL(0) || PG_ARGISNULL(1) || PG_ARGISNULL(2)) PG_RETURN_NULL();
	GSERIALIZED *g     = PG_GETARG_GSERIALIZED_P(0);
	int32        layer = PG_GETARG_INT32(1);
	int32        depth = PG_GETARG_INT32(2);
	double lon, lat;
	const char *err = NULL;
	if (!geom_to_lonlat(g, &lon, &lat, &err))
		ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
						errmsg("%s", err)));
	if (layer < 0 || depth < 0 || layer + 2 + depth > hex9_lmax())
		ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
						errmsg("h9e_encode: need layer >= 0, depth >= 0 and "
						       "layer + 2 + depth <= %d (the nibble budget); "
						       "got layer %d depth %d",
						       hex9_lmax(), layer, depth)));
	uint8_t u[UUID_LEN];
	const int rc = sphere ? hex9_e4h_encode_sphere(lon, lat, layer, depth, u)
	                      : hex9_e4h_encode(lon, lat, layer, depth, u);
	if (rc != 0)
		ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
						errmsg("h9e_encode: could not encode point (%g, %g)"
						       "%s", lon, lat,
						       rc == 3 ? " — point two seams from its host"
						               : "")));
	PG_RETURN_UUID_P(make_pg_uuid(u));
}

PG_FUNCTION_INFO_V1(h9e_encode);
Datum h9e_encode(PG_FUNCTION_ARGS) { return h9e_encode_common(fcinfo, false); }

PG_FUNCTION_INFO_V1(h9e_encode_sphere);
Datum h9e_encode_sphere(PG_FUNCTION_ARGS) { return h9e_encode_common(fcinfo, true); }

static Datum h9e_point_common(PG_FUNCTION_ARGS, bool sphere, bool partner) {
	if (PG_ARGISNULL(0)) PG_RETURN_NULL();
	pg_uuid_t *u = PG_GETARG_UUID_P(0);
	double lon, lat;
	int rc;
	if (partner) rc = sphere ? hex9_e4h_partner_sphere(u->data, &lon, &lat)
	                         : hex9_e4h_partner(u->data, &lon, &lat);
	else         rc = sphere ? hex9_e4h_decode_sphere(u->data, &lon, &lat)
	                         : hex9_e4h_decode(u->data, &lon, &lat);
	if (rc != 0)
		ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
						errmsg("%s: not an E4H uuid (or bad tail grammar)",
						       partner ? "h9e_partner" : "h9e_decode")));
	PG_RETURN_POINTER(lonlat_to_geom(lon, lat, sphere ? 0 : 4326));
}

PG_FUNCTION_INFO_V1(h9e_decode);
Datum h9e_decode(PG_FUNCTION_ARGS) { return h9e_point_common(fcinfo, false, false); }

PG_FUNCTION_INFO_V1(h9e_decode_sphere);
Datum h9e_decode_sphere(PG_FUNCTION_ARGS) { return h9e_point_common(fcinfo, true, false); }

PG_FUNCTION_INFO_V1(h9e_partner);
Datum h9e_partner(PG_FUNCTION_ARGS) { return h9e_point_common(fcinfo, false, true); }

PG_FUNCTION_INFO_V1(h9e_partner_sphere);
Datum h9e_partner_sphere(PG_FUNCTION_ARGS) { return h9e_point_common(fcinfo, true, true); }

PG_FUNCTION_INFO_V1(h9e_bin);
Datum h9e_bin(PG_FUNCTION_ARGS) {
	if (PG_ARGISNULL(0) || PG_ARGISNULL(1)) PG_RETURN_NULL();
	pg_uuid_t *u     = PG_GETARG_UUID_P(0);
	int32      depth = PG_GETARG_INT32(1);
	uint8_t out[UUID_LEN];
	if (hex9_e4h_bin(u->data, depth, out) != 0)
		ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
						errmsg("h9e_bin: not an E4H uuid, or depth %d exceeds "
						       "the tail (truncation only — the tail cannot "
						       "be deepened)", depth)));
	PG_RETURN_UUID_P(make_pg_uuid(out));
}

PG_FUNCTION_INFO_V1(h9e_depth);
Datum h9e_depth(PG_FUNCTION_ARGS) {
	if (PG_ARGISNULL(0)) PG_RETURN_NULL();
	pg_uuid_t *u = PG_GETARG_UUID_P(0);
	PG_RETURN_INT32(hex9_e4h_depth(u->data));
}

PG_FUNCTION_INFO_V1(h9_is_e4h);
Datum h9_is_e4h(PG_FUNCTION_ARGS) {
	if (PG_ARGISNULL(0)) PG_RETURN_NULL();
	pg_uuid_t *u = PG_GETARG_UUID_P(0);
	PG_RETURN_BOOL(hex9_is_e4h(u->data) != 0);
}

PG_FUNCTION_INFO_V1(h9e_host);
Datum h9e_host(PG_FUNCTION_ARGS) {
	if (PG_ARGISNULL(0)) PG_RETURN_NULL();
	pg_uuid_t *u = PG_GETARG_UUID_P(0);
	uint8_t host[UUID_LEN], digits[28];
	int half, nd;
	if (hex9_e4h_split(u->data, host, &half, digits, &nd) != 0)
		ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
						errmsg("h9e_host: not an E4H uuid")));
	PG_RETURN_UUID_P(make_pg_uuid(host));
}

PG_FUNCTION_INFO_V1(h9e_label);
Datum h9e_label(PG_FUNCTION_ARGS) {
	if (PG_ARGISNULL(0)) PG_RETURN_NULL();
	pg_uuid_t *u = PG_GETARG_UUID_P(0);
	char buf[64];
	const int len = hex9_e4h_label(u->data, buf, sizeof buf);
	if (len < 0)
		ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
						errmsg("h9e_label: not an E4H uuid")));
	PG_RETURN_TEXT_P(cstring_to_text_with_len(buf, len));
}

PG_FUNCTION_INFO_V1(h9e_parse_label);
Datum h9e_parse_label(PG_FUNCTION_ARGS) {
	if (PG_ARGISNULL(0)) PG_RETURN_NULL();
	char *label = text_to_cstring(PG_GETARG_TEXT_PP(0));
	uint8_t u[UUID_LEN];
	if (hex9_e4h_parse_label(label, u) != 0)
		ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
						errmsg("h9e_parse_label: not a valid E4H label: \"%s\"",
						       label)));
	PG_RETURN_UUID_P(make_pg_uuid(u));
}
} /* extern "C" */

/* ═══ Grid verbs: h9_aim / h9_walk_to / h9_vision_cone (Hex9 2.3.0) ═════════
 * Cell-first field primitives (the SQL twins of the python wheel's
 * hex9.verbs): bin KEYS in and out — the full-uuid re-derivation is internal
 * (deliberately unlike the k-family, whose outputs are join keys; a verb's
 * output is the next verb's input). Bearings are degrees clockwise from
 * north, advisory FP over centroids — they steer, they never mint. */

extern "C" {
PG_FUNCTION_INFO_V1(h9_aim);
Datum h9_aim(PG_FUNCTION_ARGS) {
	if (PG_ARGISNULL(0) || PG_ARGISNULL(1) || PG_ARGISNULL(2)) PG_RETURN_NULL();
	pg_uuid_t *u       = PG_GETARG_UUID_P(0);
	int32      layer   = PG_GETARG_INT32(1);
	float8     bearing = PG_GETARG_FLOAT8(2);
	uint8_t out[UUID_LEN];
	if (hex9_aim(u->data, layer, bearing, out) != 0)
		ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
						errmsg("h9_aim: not a valid H9 key at layer %d "
						       "(or non-finite bearing)", layer)));
	PG_RETURN_UUID_P(make_pg_uuid(out));
}

PG_FUNCTION_INFO_V1(h9_walk_to);
Datum h9_walk_to(PG_FUNCTION_ARGS) {
	if (PG_ARGISNULL(0) || PG_ARGISNULL(1) || PG_ARGISNULL(2)) PG_RETURN_NULL();
	pg_uuid_t *src   = PG_GETARG_UUID_P(0);
	pg_uuid_t *dest  = PG_GETARG_UUID_P(1);
	int32      layer = PG_GETARG_INT32(2);
	uint8_t   *obs   = NULL;
	size_t     nobs  = 0;
	if (!PG_ARGISNULL(3))
		nobs = h9_uuid_array_bytes(PG_GETARG_ARRAYTYPE_P(3), &obs);
	int64_t max_expand = PG_ARGISNULL(4) ? 0 : PG_GETARG_INT32(4);
	if (max_expand <= 0) max_expand = 5000;

	const int64_t cap = max_expand + 2;
	uint8_t *path = (uint8_t *) palloc((Size)cap * UUID_LEN);
	const int64_t n = hex9_walk_to(src->data, dest->data, layer, obs, nobs,
	                               max_expand, path, cap);
	if (n < 0)
		ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
						errmsg("h9_walk_to: invalid keys or layer %d", layer)));
	if (n == 0) PG_RETURN_NULL();          /* no path within max_expand */

	Datum *res = (Datum *) palloc((size_t)n * sizeof(Datum));
	for (int64_t i = 0; i < n; i++)
		res[i] = UUIDPGetDatum(make_pg_uuid(path + (size_t)i * UUID_LEN));
	int dims[1] = { (int)n };
	int lbs[1]  = { 1 };
	PG_RETURN_ARRAYTYPE_P(construct_md_array(res, NULL, 1, dims, lbs,
	                                         UUIDOID, UUID_LEN, false,
	                                         TYPALIGN_CHAR));
}

PG_FUNCTION_INFO_V1(h9_vision_cone);
Datum h9_vision_cone(PG_FUNCTION_ARGS) {
	return h9_uuid_set_srf(fcinfo, [&](H9UuidSetState *state) -> int64_t {
		pg_uuid_t *u       = PG_GETARG_UUID_P(0);
		int32      layer   = PG_GETARG_INT32(1);
		float8     bearing = PG_GETARG_FLOAT8(2);
		float8     half    = PG_GETARG_FLOAT8(3);
		int32      k       = PG_GETARG_INT32(4);
		uint8_t   *obs     = NULL;
		size_t     nobs    = 0;
		if (!PG_ARGISNULL(5))
			nobs = h9_uuid_array_bytes(PG_GETARG_ARRAYTYPE_P(5), &obs);
		if (k < 0 || hex9_disk_ncells(k) > 60000000)
			ereport(ERROR, (errcode(ERRCODE_PROGRAM_LIMIT_EXCEEDED),
							errmsg("h9_vision_cone: k=%d out of range", k)));
		const int64_t cap = hex9_disk_ncells(k);
		state->uuids = (uint8_t *) palloc((Size)cap * UUID_LEN);
		const int64_t n = hex9_vision_cone(u->data, layer, bearing, half, k,
		                                   obs, nobs, state->uuids, cap);
		if (n < 0)
			ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
							errmsg("h9_vision_cone: invalid key, layer %d, "
							       "or non-finite angles", layer)));
		return n;
	});
}
} /* extern "C" */

/* ── True-bearing correction (Hex9 2.3.0) ───────────────────────────────────
 * A physically measured true-north heading is a WGS84 geodesic azimuth; the
 * verbs' bearings are spherical trig over geodetic coordinates. These convert
 * between the two at the latitude where the heading holds (the source cell):
 * tan(b_verb) = (M/N)·tan(b_true). Identity at poles and cardinal bearings;
 * largest ~0.19° at equatorial diagonals. Magnetic headings must already be
 * reduced to true north (declination is epoch/model-dependent — dataset
 * metadata). Meaningless on the sphere datum. */
extern "C" {
PG_FUNCTION_INFO_V1(h9_bearing_from_true);
Datum h9_bearing_from_true(PG_FUNCTION_ARGS) {
	PG_RETURN_FLOAT8(hex9_bearing_from_true(PG_GETARG_FLOAT8(0),
	                                        PG_GETARG_FLOAT8(1)));
}

PG_FUNCTION_INFO_V1(h9_bearing_to_true);
Datum h9_bearing_to_true(PG_FUNCTION_ARGS) {
	PG_RETURN_FLOAT8(hex9_bearing_to_true(PG_GETARG_FLOAT8(0),
	                                      PG_GETARG_FLOAT8(1)));
}
} /* extern "C" */
