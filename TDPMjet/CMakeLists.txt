# -*- mode: cmake -*-

# Create a library called "lib<name>" which includes the source files given in
# the array .
# The extension is already found.  Any number of sources could be listed here.

set(INCLUDE_DIRECTORIES
${CMAKE_SOURCE_DIR}/TDPMjet
${CMAKE_SOURCE_DIR}/STEER
${CMAKE_SOURCE_DIR}/EVGEN
${CMAKE_SOURCE_DIR}/DPMJET
${ROOT_INCLUDE_DIR}
)

include_directories( ${INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
${ROOT_LIBRARY_DIR}
) 

link_directories( ${LINK_DIRECTORIES})

SetModule()

include(CMake_libTDPMjet.txt)


