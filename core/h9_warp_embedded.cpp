/* h9_warp_embedded.cpp — embed the warp blob via .incbin (compile-time).
 *
 * Replaces the former ~13 MB ASCII-array header. Runtime is identical: the
 * assembler bakes the binary's bytes into the library's read-only data section,
 * so there is NO runtime file load and the library stays self-contained.
 * H9_WARP_BLOB is the absolute path to the .h9warp binary, set by CMake.
 */
#include "h9_warp_embedded.h"
#include <cstddef>

#ifndef H9_WARP_BLOB
#error "H9_WARP_BLOB (absolute path to the .h9warp binary) must be defined by the build"
#endif

/* Stringize the build-provided path into a quoted literal for .incbin
 * (robust to however the build system passes the -D). */
#define H9_STR2(x) #x
#define H9_STR(x)  H9_STR2(x)

__asm__(
#if defined(__APPLE__)
    ".const_data\n"
    ".globl _h9_warp_blob\n"
    ".p2align 4\n"
    "_h9_warp_blob:\n"
    ".incbin " H9_STR(H9_WARP_BLOB) "\n"
    ".globl _h9_warp_blob_end\n"
    "_h9_warp_blob_end:\n"
#else
    ".section .rodata\n"
    ".globl h9_warp_blob\n"
    ".balign 16\n"
    "h9_warp_blob:\n"
    ".incbin " H9_STR(H9_WARP_BLOB) "\n"
    ".globl h9_warp_blob_end\n"
    "h9_warp_blob_end:\n"
#endif
);

extern "C" const unsigned char h9_warp_blob[];
extern "C" const unsigned char h9_warp_blob_end[];

namespace h9 {
const unsigned char *const EMBEDDED_WARP_DATA = h9_warp_blob;
const std::size_t          EMBEDDED_WARP_SIZE =
    static_cast<std::size_t>(h9_warp_blob_end - h9_warp_blob);
}
