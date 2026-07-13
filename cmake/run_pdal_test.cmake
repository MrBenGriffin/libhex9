# run_pdal_test.cmake — ctest driver for the filters.hex9 PDAL plugin.
# Writes a pipeline, runs `pdal pipeline` with the built plugin on PDAL_DRIVER_PATH,
# then verifies the output against the libhex9 C ABI.
# Args (via -D): PDAL_EXE DRIVER_DIR PTS VERIFY_EXE WORKDIR LIBDIR
set(pipe "${WORKDIR}/pdal_pipeline.json")
set(out  "${WORKDIR}/pdal_out.csv")

file(WRITE "${pipe}"
"[
  {\"type\":\"readers.text\",\"filename\":\"${PTS}\",\"spatialreference\":\"EPSG:4326\"},
  {\"type\":\"filters.hex9\"},
  {\"type\":\"writers.text\",\"filename\":\"${out}\",\"order\":\"X,Y,Z,WoctX,WoctY,WoctZ\",\"keep_unspecified\":\"false\",\"precision\":15}
]")

set(ENV{PDAL_DRIVER_PATH} "${DRIVER_DIR}")
set(ENV{DYLD_LIBRARY_PATH} "${LIBDIR}")   # macOS: locate libhex9.dylib for pdal & verifier
set(ENV{LD_LIBRARY_PATH} "${LIBDIR}")     # Linux
if(WIN32)
    set(ENV{PATH} "${LIBDIR};$ENV{PATH}")  # Windows: no rpath, hex9.dll must be on PATH
endif()

execute_process(COMMAND "${PDAL_EXE}" pipeline "${pipe}" RESULT_VARIABLE rc)
if(NOT rc EQUAL 0)
    message(FATAL_ERROR "pdal pipeline failed (exit ${rc})")
endif()

execute_process(COMMAND "${VERIFY_EXE}" "${out}" RESULT_VARIABLE vrc)
if(NOT vrc EQUAL 0)
    message(FATAL_ERROR "pdal_verify failed (exit ${vrc})")
endif()
