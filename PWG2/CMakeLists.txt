# -*- mode: cmake -*-
# Create a library called "lib<name>" which includes the source files given in
# the array .
# The extension is already found.  Any number of sources could be listed here.

set(INCLUDE_DIRECTORIES
${CMAKE_SOURCE_DIR}/ANALYSIS
${CMAKE_SOURCE_DIR}/PWG2
${CMAKE_SOURCE_DIR}/PWG2/AOD
${CMAKE_SOURCE_DIR}/PWG2/EBYE
${CMAKE_SOURCE_DIR}/PWG2/EBYE/LRC
${CMAKE_SOURCE_DIR}/PWG2/FEMTOSCOPY
${CMAKE_SOURCE_DIR}/PWG2/FEMTOSCOPY/AliFemto
${CMAKE_SOURCE_DIR}/PWG2/FEMTOSCOPY/AliFemtoUser
${CMAKE_SOURCE_DIR}/PWG2/FLOW
${CMAKE_SOURCE_DIR}/PWG2/SPECTRA
${CMAKE_SOURCE_DIR}/STEER
${CMAKE_SOURCE_DIR}/CORRFW
${ALIROOT_INSTALL_DIR}/include 
${ROOT_INCLUDE_DIR}
)

include_directories( ${INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
${ROOT_LIBRARY_DIR}
) 

link_directories( ${LINK_DIRECTORIES})

SetModule()

include(CMake_libPWG2spectra.txt)

include(CMake_libPWG2flow.txt)

include(CMake_libPWG2femtoscopyUser.txt)

include(CMake_libPWG2femtoscopy.txt)

include(CMake_libPWG2ebye.txt)

include(CMake_libPWG2AOD.txt)



