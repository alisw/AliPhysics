# -*- mode: cmake -*-
# Create a library called "lib<name>" which includes the source files given in
# the array .
# The extension is already found.  Any number of sources could be listed here.

set(INCLUDE_DIRECTORIES
${CMAKE_SOURCE_DIR}/MUON
${CMAKE_SOURCE_DIR}/MUON/mapping
${CMAKE_SOURCE_DIR}/STEER 
${CMAKE_SOURCE_DIR}/RAW
${ROOT_INCLUDE_DIR}
)

include_directories( ${INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
${ROOT_LIBRARY_DIR}
) 

link_directories( ${LINK_DIRECTORIES})

SetModule()

include(CMake_libMUONraw.txt)

include(CMake_libMUONcalib.txt)

include(CMake_libMUONrec.txt)

include(CMake_libMUONgraphics.txt)

include(CMake_libMUONsim.txt)

include(CMake_libMUONtrigger.txt)

include(CMake_libMUONcore.txt)

include(CMake_libMUONmapping.txt)

include(CMake_libMUONbase.txt)

include(CMake_libMUONgeometry.txt)

include(CMake_libMUONevaluation.txt)

include(CMake_libMUONshuttle.txt)

include(CMake_binmchview.txt)

