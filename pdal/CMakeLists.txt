# filters.hex9 — PDAL plugin projecting WGS84 lon/lat to Hex9 w_oct.
# Loaded by PDAL via PDAL_DRIVER_PATH; links the in-tree libhex9 C ABI.
add_library(pdal_plugin_filter_hex9 SHARED Hex9Filter.cpp)

target_include_directories(pdal_plugin_filter_hex9 PRIVATE
    ${CMAKE_SOURCE_DIR}          # hex9_c.h
    ${PDAL_INCLUDE_DIRS})

target_link_libraries(pdal_plugin_filter_hex9 PRIVATE
    hex9                         # the in-tree libhex9 target
    ${PDAL_LIBRARIES})

# PDAL discovers plugins by filename: libpdal_plugin_filter_<name>.{so,dylib} on
# Linux/macOS, pdal_plugin_filter_<name>.dll (NO "lib") on Windows. Let CMake use
# the platform-default prefix so the name is right on every OS.
set_target_properties(pdal_plugin_filter_hex9 PROPERTIES
    CXX_STANDARD 17
    OUTPUT_NAME "pdal_plugin_filter_hex9")

# filters.hex9bin — aggregate a point cloud into Hex9 cells (adaptive digest).
add_library(pdal_plugin_filter_hex9bin SHARED Hex9BinFilter.cpp)
target_include_directories(pdal_plugin_filter_hex9bin PRIVATE
    ${CMAKE_SOURCE_DIR} ${PDAL_INCLUDE_DIRS})
target_link_libraries(pdal_plugin_filter_hex9bin PRIVATE hex9 ${PDAL_LIBRARIES})
set_target_properties(pdal_plugin_filter_hex9bin PROPERTIES
    CXX_STANDARD 17
    OUTPUT_NAME "pdal_plugin_filter_hex9bin")

# ── install ──────────────────────────────────────────────────────────────────
# PDAL discovers plugins on PDAL_DRIVER_PATH, or in its own lib dir. Point this
# at your PDAL plugin dir (e.g. `$(pdal-config --prefix)/lib`) for zero-config
# discovery, or leave it at the default lib dir and set PDAL_DRIVER_PATH to it.
# INSTALL_RPATH lets the installed plugin find the installed libhex9.
include(GNUInstallDirs)
set(HEX9_PDAL_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}" CACHE STRING
    "Install directory for the Hex9 PDAL plugins")
set_target_properties(pdal_plugin_filter_hex9 pdal_plugin_filter_hex9bin PROPERTIES
    INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR};${PDAL_LIBRARY_DIRS}")
install(TARGETS pdal_plugin_filter_hex9 pdal_plugin_filter_hex9bin
        LIBRARY DESTINATION "${HEX9_PDAL_INSTALL_DIR}")

# ── guarded ctest: run a real pipeline through the plugin, verify vs the C ABI ─
# Needs the `pdal` CLI on PATH (present when PDAL was found). Drives
# cmake/run_pdal_test.cmake so the plugin is exercised exactly as a user would.
find_program(PDAL_EXECUTABLE pdal)
if(PDAL_EXECUTABLE)
    add_executable(pdal_verify ${CMAKE_SOURCE_DIR}/test/pdal_verify.c)
    target_include_directories(pdal_verify PRIVATE ${CMAKE_SOURCE_DIR})
    target_link_libraries(pdal_verify PRIVATE hex9)
    add_test(NAME pdal_filter COMMAND ${CMAKE_COMMAND}
        -DPDAL_EXE=${PDAL_EXECUTABLE}
        -DDRIVER_DIR=$<TARGET_FILE_DIR:pdal_plugin_filter_hex9>
        -DPTS=${CMAKE_SOURCE_DIR}/test/pdal_pts.csv
        -DVERIFY_EXE=$<TARGET_FILE:pdal_verify>
        -DWORKDIR=${CMAKE_CURRENT_BINARY_DIR}
        -DLIBDIR=$<TARGET_FILE_DIR:hex9>
        -P ${CMAKE_SOURCE_DIR}/cmake/run_pdal_test.cmake)

    add_executable(pdal_bin_verify ${CMAKE_SOURCE_DIR}/test/pdal_bin_verify.c)
    add_test(NAME pdal_hex9bin COMMAND ${CMAKE_COMMAND}
        -DPDAL_EXE=${PDAL_EXECUTABLE}
        -DDRIVER_DIR=$<TARGET_FILE_DIR:pdal_plugin_filter_hex9bin>
        -DPTS=${CMAKE_SOURCE_DIR}/test/pdal_bin_pts.csv
        -DVERIFY_EXE=$<TARGET_FILE:pdal_bin_verify>
        -DWORKDIR=${CMAKE_CURRENT_BINARY_DIR}
        -DLIBDIR=$<TARGET_FILE_DIR:hex9>
        -DNPOINTS=24
        -P ${CMAKE_SOURCE_DIR}/cmake/run_pdal_bin_test.cmake)
endif()
