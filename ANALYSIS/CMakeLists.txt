# -*- mode: cmake -*-
# Create a library called "lib<name>" which includes the source files given in
# the array .
# The extension is already found.  Any number of sources could be listed here.

set(INCLUDE_DIRECTORIES
${CMAKE_SOURCE_DIR}/ANALYSIS
${CMAKE_SOURCE_DIR}/ANALYSIS/Aliengui
${CMAKE_SOURCE_DIR}/STEER
${CMAKE_SOURCE_DIR}/TPC
${CMAKE_SOURCE_DIR}/ITS
${ROOT_INCLUDE_DIR}
)

include_directories( ${INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
${ROOT_LIBRARY_DIR}
) 

link_directories( ${LINK_DIRECTORIES})

SetModule()

include(CMake_libANALYSISalice.txt)

include(CMake_libANALYSIS.txt)

include(CMake_libANALYSISRL.txt)

include(CMake_libAliengui.txt)

include(CMake_binaliengui.txt)



