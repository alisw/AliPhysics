# -*- mode: cmake -*-

# Create a library called "lib<name>" which includes the source files given in
# the array .
# The extension is already found.  Any number of sources could be listed here.

set(INCLUDE_DIRECTORIES
${CMAKE_SOURCE_DIR}/PWG4
${CMAKE_SOURCE_DIR}/PWG4/PartCorrBase
${CMAKE_SOURCE_DIR}/PWG4/PartCorrDep
${CMAKE_SOURCE_DIR}/PWG4/JetTasks
${CMAKE_SOURCE_DIR}/PWG4/GammaConv
${CMAKE_SOURCE_DIR}/PHOS
${CMAKE_SOURCE_DIR}/STEER
${CMAKE_SOURCE_DIR}/ANALYSIS 
${CMAKE_SOURCE_DIR}/JETAN 
${ROOT_INCLUDE_DIR}
)

include_directories( ${INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
${ROOT_LIBRARY_DIR}
) 

link_directories( ${LINK_DIRECTORIES})

SetModule()

include(CMake_libPWG4PartCorrBase.txt)
include(CMake_libPWG4PartCorrDep.txt)
include(CMake_libPWG4JetTasks.txt)
include(CMake_libPWG4GammaConv.txt)
