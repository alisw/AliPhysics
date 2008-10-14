# -*- mode: cmake -*-
# Create a library called "lib<name>" which includes the source files given in
# the array .
# The extension is already found.  Any number of sources could be listed here.

set(INCLUDE_DIRECTORIES
${CMAKE_SOURCE_DIR}/EMCAL
${CMAKE_SOURCE_DIR}/EMCAL/Calib
${CMAKE_SOURCE_DIR}/EMCAL/SMcalib
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

include(CMake_libEMCALrec.txt)

include(CMake_libEMCALbase.txt)

include(CMake_libEMCALsim.txt)

include(CMake_libEMCALcalib.txt)


