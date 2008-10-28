# -*- mode: cmake -*-
# Create a library called "lib<name>" which includes the source files given in
# the array .
# The extension is already found.  Any number of sources could be listed here.

set(INCLUDE_DIRECTORIES
${CMAKE_SOURCE_DIR}/PHOS
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

include(CMake_libPHOSUtils.txt)

include(CMake_libPHOSrec.txt)

include(CMake_libPHOSbase.txt) 

include(CMake_libPHOSsim.txt) 

include(CMake_libPHOSshuttle.txt) 


