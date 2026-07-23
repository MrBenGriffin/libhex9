/* h9_extension.cpp — DuckDB extension entry point for libhex9.
 *
 * Load policy (Ben's ruling, 2026-07-21): warp init failure is a HARD LOAD
 * ERROR. On failure the core degrades to an identity field whose output is
 * NOT Hex9 addresses; an analytics engine silently producing wrong-regime
 * addresses would poison downstream data (same policy as postgis_hex9's
 * _PG_init, deliberately not PDAL's warn-and-continue).
 *
 * The version guard mirrors postgis_hex9: HEX9_VERSION is what this shim was
 * compiled against, hex9_version() is the core actually linked. Here the core
 * is compiled into the extension so they cannot drift, but the guard is kept
 * so a future dynamically-linked build fails loudly rather than emitting a
 * different regime's addresses (see hex9_c.h's version notes).
 */
#include "hex9_extension.hpp"

#include "duckdb.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

#include "hex9_c.h"

#include <cstring>

namespace duckdb {

void RegisterH9Scalar(ExtensionLoader &loader);
void RegisterH9Curve(ExtensionLoader &loader);
void RegisterH9Table(ExtensionLoader &loader);

static void LoadInternal(ExtensionLoader &loader) {
	char errbuf[256] = {0};
	if (hex9_init(errbuf, sizeof(errbuf)) != 0) {
		throw IOException("hex9: warp init failed — refusing to load (the degraded field would not "
		                  "produce Hex9 addresses): %s",
		                  errbuf);
	}
	/* hex9_version() is "libhex9 <ver> (<date> <time>)"; match the version
	 * token exactly (equal prefix AND a terminator after it), as postgis_hex9
	 * does in h9_check_lib_version. */
	const char *runtime = hex9_version();
	const char *p = runtime;
	if (const char *sp = std::strchr(runtime, ' ')) {
		p = sp + 1;
	}
	const size_t want = std::strlen(HEX9_VERSION);
	if (std::strncmp(p, HEX9_VERSION, want) != 0 || (p[want] != '\0' && p[want] != ' ')) {
		throw IOException("hex9: extension built against libhex9 %s but running %s — refusing to load "
		                  "(the addressing regime may differ)",
		                  HEX9_VERSION, runtime);
	}
	RegisterH9Scalar(loader);
	RegisterH9Curve(loader);
	RegisterH9Table(loader);
}

void Hex9Extension::Load(ExtensionLoader &loader) {
	LoadInternal(loader);
}

std::string Hex9Extension::Name() {
	return "hex9";
}

std::string Hex9Extension::Version() const {
#ifdef EXT_VERSION_HEX9
	return EXT_VERSION_HEX9;
#else
	return "";
#endif
}

} // namespace duckdb

extern "C" {
DUCKDB_CPP_EXTENSION_ENTRY(hex9, loader) {
	duckdb::LoadInternal(loader);
}
}
