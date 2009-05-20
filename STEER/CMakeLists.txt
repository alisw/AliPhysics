# -*- mode: cmake -*-
# Create all libraries of STEER subdirectory 

set(INCLUDE_DIRECTORIES
${CMAKE_SOURCE_DIR}/STEER
${CMAKE_BINARY_DIR}/STEER
${CMAKE_SOURCE_DIR}/RAW
${CMAKE_SOURCE_DIR}/HLT/sim
${ROOT_INCLUDE_DIR}
)

include_directories( ${INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
${ROOT_LIBRARY_DIR}
) 

link_directories( ${LINK_DIRECTORIES})

SetModule()

# include files with the descriptions how to build the
# different libraries of this subdirectory

include(CMake_libESD.txt)

include(CMake_libCDB.txt)

include(CMake_libSTEERBase.txt)

include(CMake_libSTEER.txt)

include(CMake_libAOD.txt)
