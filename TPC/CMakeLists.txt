# -*- mode: cmake -*-
# Create a library called "lib<name>" which includes the source files given in
# the array .
# The extension is already found.  Any number of sources could be listed here.

set(INCLUDE_DIRECTORIES
${CMAKE_SOURCE_DIR}/TPC
${CMAKE_SOURCE_DIR}/STEER
${CMAKE_BINARY_DIR}/STEER
${CMAKE_SOURCE_DIR}/STAT
${CMAKE_SOURCE_DIR}/RAW
${CMAKE_SOURCE_DIR}/ANALYSIS
${ROOT_INCLUDE_DIR}
)

include_directories( ${INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
${ROOT_LIBRARY_DIR}
) 

link_directories( ${LINK_DIRECTORIES})

SetModule()

include(CMake_libTPCbase.txt)

include(CMake_libTPCcalib.txt)

include(CMake_libTPCfast.txt)

include(CMake_libTPCrec.txt)

include(CMake_libTPCmon.txt)

include(CMake_libTPCsim.txt)


