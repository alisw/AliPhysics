# -*- mode: cmake -*-
# Create a library called "lib<name>" which includes the source files given in
# the array .
# The extension is already found.  Any number of sources could be listed here.

set(INCLUDE_DIRECTORIES
${CMAKE_SOURCE_DIR}/TRD
${CMAKE_SOURCE_DIR}/TRD/Cal
${CMAKE_SOURCE_DIR}/TRD/qaGui
${CMAKE_SOURCE_DIR}/TRD/qaAnalysis
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

include(CMake_libTRDgui.txt)

include(CMake_libTRDqaAnalysis.txt)

include(CMake_libTRDrec.txt)

include(CMake_libTRDbase.txt)

include(CMake_libTRDsim.txt)

include(CMake_libTRDfast.txt)



