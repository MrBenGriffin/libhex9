# run_pdal_bin_test.cmake — ctest driver for the filters.hex9bin plugin.
# Args (via -D): PDAL_EXE DRIVER_DIR PTS VERIFY_EXE WORKDIR LIBDIR NPOINTS
set(pipe "${WORKDIR}/pdal_bin_pipeline.json")
set(out  "${WORKDIR}/pdal_bin_out.csv")

file(WRITE "${pipe}"
"[
  {\"type\":\"readers.text\",\"filename\":\"${PTS}\",\"spatialreference\":\"EPSG:4326\"},
  {\"type\":\"filters.hex9bin\",\"min_layer\":4,\"max_layer\":9,\"ceiling\":25,\"floor\":4,\"aggregates\":\"Z:mean,Z:max,Z:min,Classification:mode,Intensity:mean\"},
  {\"type\":\"writers.text\",\"filename\":\"${out}\",\"order\":\"H9Layer,H9Count,H9Value,Z,ZMax,ZMin,Classification,Intensity\",\"keep_unspecified\":\"false\",\"precision\":6}
]")

set(ENV{PDAL_DRIVER_PATH} "${DRIVER_DIR}")
set(ENV{DYLD_LIBRARY_PATH} "${LIBDIR}")
set(ENV{LD_LIBRARY_PATH} "${LIBDIR}")
if(WIN32)
    set(ENV{PATH} "${LIBDIR};$ENV{PATH}")
endif()

execute_process(COMMAND "${PDAL_EXE}" pipeline "${pipe}" RESULT_VARIABLE rc)
if(NOT rc EQUAL 0)
    message(FATAL_ERROR "pdal pipeline (hex9bin) failed (exit ${rc})")
endif()

execute_process(COMMAND "${VERIFY_EXE}" "${out}" "${NPOINTS}" RESULT_VARIABLE vrc)
if(NOT vrc EQUAL 0)
    message(FATAL_ERROR "pdal_bin_verify failed (exit ${vrc})")
endif()
