# -*- mode: cmake -*-
# Create a library called "lib<name>" which includes the source files given in
# the array .
# The extension is already found.  Any number of sources could be listed here.

set(INCLUDE_DIRECTORIES
RAW
${CMAKE_SOURCE_DIR}/RAW
${CMAKE_SOURCE_DIR}/STEER
${ROOT_INCLUDE_DIR}
)

include_directories( ${INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
${ROOT_LIBRARY_DIR}
) 

link_directories( ${LINK_DIRECTORIES})

SetModule()

include(CMake_libRAWDatasim.txt)

include(CMake_libMDC.txt)

include(CMake_libMDCFilter.txt)

include(CMake_libRAWDatarec.txt)

include(CMake_libRAWDatabase.txt)

include(CMake_libRAWDatarecOnline.txt)

include(CMake_binalimdc.txt)

include(CMake_binalirawdump.txt)

include(CMake_bindateStream.txt)

