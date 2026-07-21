cmake_minimum_required(VERSION 3.5...3.29)

set(TARGET_NAME hex9)
set(EXTENSION_NAME ${TARGET_NAME}_extension)
set(LOADABLE_EXTENSION_NAME ${TARGET_NAME}_loadable_extension)
project(${TARGET_NAME})

# Root of the libhex9 source tree. In-tree this is the parent directory; in the
# community-extensions wrapper repo it points at the libhex9 submodule instead —
# this is the only line that changes on lift-out.
set(HEX9_ROOT ${CMAKE_CURRENT_LIST_DIR}/..)

set(H9_WARP_BLOB_PATH "${HEX9_ROOT}/core/Sphere_l6_fund.f64g.h9warp")
if(NOT EXISTS "${H9_WARP_BLOB_PATH}")
  message(FATAL_ERROR "hex9 warp blob not found:\n    ${H9_WARP_BLOB_PATH}")
endif()

include_directories(src/include "${HEX9_ROOT}" "${HEX9_ROOT}/core")

# The three core sources are compiled directly (rather than linking the parent
# build's hex9Static) so this extension builds from its own tree under the
# extension-template CI, including the wasm targets.
set(EXTENSION_SOURCES
    src/h9_extension.cpp
    src/h9_scalar.cpp
    src/h9_curve.cpp
    src/h9_table.cpp
    "${HEX9_ROOT}/hex9_c.cpp"
    "${HEX9_ROOT}/core/h9_warp_runtime.cpp"
    "${HEX9_ROOT}/core/h9_warp_embedded.cpp")

# The core's `.incbin` inline asm doesn't survive every toolchain: Emscripten's
# wasm backend rejects it, and MSVC x64 has no __asm__ at all. Swap in a
# generated TU carrying the same h9:: symbols (array style for MSVC, whose
# string literals cap at 64 KB).
if(EMSCRIPTEN OR MSVC)
  set(H9_WARP_EMBED_TU "${CMAKE_CURRENT_BINARY_DIR}/h9_warp_embedded_gen.cpp")
  if(MSVC)
    set(H9_EMBED_STYLE "--array")
  else()
    set(H9_EMBED_STYLE "")
  endif()
  if(NOT EXISTS "${H9_WARP_EMBED_TU}")
    find_package(Python3 REQUIRED COMPONENTS Interpreter)
    execute_process(
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_LIST_DIR}/scripts/gen_warp_embed.py"
              "${H9_WARP_BLOB_PATH}" "${H9_WARP_EMBED_TU}"
      RESULT_VARIABLE H9_EMBED_RC)
    if(NOT H9_EMBED_RC EQUAL 0)
      message(FATAL_ERROR "hex9: gen_warp_embed.py failed (${H9_EMBED_RC})")
    endif()
  endif()
  list(REMOVE_ITEM EXTENSION_SOURCES "${HEX9_ROOT}/core/h9_warp_embedded.cpp")
  list(APPEND EXTENSION_SOURCES "${H9_WARP_EMBED_TU}")
endif()

build_static_extension(${TARGET_NAME} ${EXTENSION_SOURCES})
build_loadable_extension(${TARGET_NAME} " " ${EXTENSION_SOURCES})

foreach(tgt ${EXTENSION_NAME} ${LOADABLE_EXTENSION_NAME})
  # _USE_MATH_DEFINES: MSVC/mingw hide M_PI without it (no-op elsewhere).
  target_compile_definitions(${tgt} PRIVATE
    H9_WARP_BLOB="${H9_WARP_BLOB_PATH}"
    _USE_MATH_DEFINES)
  set_target_properties(${tgt} PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON)
endforeach()

install(
  TARGETS ${EXTENSION_NAME}
  EXPORT "${DUCKDB_EXPORT_SET}"
  LIBRARY DESTINATION "${INSTALL_LIB_DIR}"
  ARCHIVE DESTINATION "${INSTALL_LIB_DIR}")
