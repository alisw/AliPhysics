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
${CMAKE_SOURCE_DIR}/PWG2/EVCHAR
${CMAKE_SOURCE_DIR}/PWG2/FEMTOSCOPY
${CMAKE_SOURCE_DIR}/PWG2/FEMTOSCOPY/AliFemto
${CMAKE_SOURCE_DIR}/PWG2/FEMTOSCOPY/AliFemtoUser
${CMAKE_SOURCE_DIR}/PWG2/FLOW/AliFlowCommon
${CMAKE_SOURCE_DIR}/PWG2/FLOW/AliFlowTasks
${CMAKE_SOURCE_DIR}/PWG2/KINK
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

include(CMake_libPWG2kink.txt)

include(CMake_libPWG2flowCommon.txt)

include(CMake_libPWG2flowTasks.txt)

include(CMake_libPWG2femtoscopyUser.txt)

include(CMake_libPWG2femtoscopy.txt)

include(CMake_libPWG2ebye.txt)

include(CMake_libPWG2AOD.txt)

include(CMake_libPWG2evchar.txt)


