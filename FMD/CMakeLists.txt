# -*- mode: cmake -*-
# Create a library called "lib<name>" which includes the source files given in
# the array .
# The extension is already found.  Any number of sources could be listed here.

set(INCLUDE_DIRECTORIES
${CMAKE_SOURCE_DIR}/FMD
${CMAKE_SOURCE_DIR}/FMD/flow
${CMAKE_SOURCE_DIR}/FMD/analysis
${CMAKE_SOURCE_DIR}/STEER
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

include(CMake_libFMDbase.txt)

include(CMake_libFMDrec.txt)

include(CMake_libFMDutil.txt)

include(CMake_libFMDflow.txt)

include(CMake_libFMDsim.txt)

include(CMake_libFMDanalysis.txt)

