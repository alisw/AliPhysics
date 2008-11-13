# -*- mode: cmake -*-
# Create a library called "lib<name>" which includes the source files given in
# the array .
# The extension is already found.  Any number of sources could be listed here.

set(INCLUDE_DIRECTORIES
${CMAKE_SOURCE_DIR}/ANALYSIS
${CMAKE_SOURCE_DIR}/EMCAL
${CMAKE_SOURCE_DIR}/EVE
${CMAKE_SOURCE_DIR}/EVE/EveBase
${CMAKE_SOURCE_DIR}/EVE/EveDet
${CMAKE_SOURCE_DIR}/EVE/EveHLT
${CMAKE_SOURCE_DIR}/STEER
${CMAKE_SOURCE_DIR}/RAW
${CMAKE_SOURCE_DIR}/MUON
${CMAKE_SOURCE_DIR}/MUON/mapping
${CMAKE_SOURCE_DIR}/TRD
${CMAKE_SOURCE_DIR}/TOF
${CMAKE_SOURCE_DIR}/ITS
${CMAKE_SOURCE_DIR}/TPC
${CMAKE_SOURCE_DIR}/PMD
${CMAKE_SOURCE_DIR}/T0
${CMAKE_SOURCE_DIR}/VZERO
${CMAKE_SOURCE_DIR}/HLT/BASE
${CMAKE_SOURCE_DIR}/HLT/BASE/HOMER
${ROOT_INCLUDE_DIR}
)

include_directories( ${INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
${ROOT_LIBRARY_DIR}
) 

link_directories( ${LINK_DIRECTORIES})

SetModule()

include(CMake_libEveBase.txt)

include(CMake_libEveDet.txt)

include(CMake_libEveHLT.txt)

include(CMake_binalieve.txt)

