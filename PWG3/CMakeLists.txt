# -*- mode: cmake -*-
# Create a library called "lib<name>" which includes the source files given in
# the array .
# The extension is already found.  Any number of sources could be listed here.

set(INCLUDE_DIRECTORIES
${CMAKE_SOURCE_DIR}/PWG3
${CMAKE_SOURCE_DIR}/PWG3/base
${CMAKE_SOURCE_DIR}/PWG3/muon
${CMAKE_SOURCE_DIR}/PWG3/vertexingHF
${CMAKE_SOURCE_DIR}/PWG3/vertexingOld
${CMAKE_SOURCE_DIR}/STEER 
${CMAKE_SOURCE_DIR}/ANALYSIS
${CMAKE_SOURCE_DIR}/CORRFW
${ROOT_INCLUDE_DIR}
)

include_directories( ${INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
${ROOT_LIBRARY_DIR}
) 

link_directories( ${LINK_DIRECTORIES})

SetModule()

include(CMake_libPWG3base.txt)

include(CMake_libPWG3.txt)

include(CMake_libPWG3vertexingHF.txt)

include(CMake_libPWG3vertexingOld.txt)

include(CMake_libPWG3muon.txt)
