# AliRoot Build System Module to find and configure GEANT3
#
# Author: Johny Jose (johny.jose@cern.ch)
#         Port of previous Makefile build to cmake

cmake_minimum_required(VERSION 2.8.4 FATAL_ERROR)

set(GEANT3DIR $ENV{GEANT3DIR})
ALICE_RootConfig(GEANT3_PLATFORM --arch)

find_path(GEANT3_SYSTEM_DIR NAMES TGeant3 PATHS ${ALICE}/geant3 ${GEANT3DIR})
message(STATUS "GEANT3 system directory: ${GEANT3_SYSTEM_DIR}")

find_path(GEANT3_INCLUDE_DIR NAMES TGeant3.h PATHS ${ALICE}/geant3/TGeant3 ${GEANT3DIR}/TGeant3)
message(STATUS "GEANT3 include directory: ${GEANT3_INCLUDE_DIR}")
set(GEANT3_ROOT ${GEANT3_INCLUDE_DIR})

find_path(GEANT3_LIBRARY_DIR NAMES libgeant321.so PATHS ${ALICE}/geant3/lib/tgt_${GEANT3_PLATFORM} ${GEANT3DIR}/lib/tgt_${GEANT3_PLATFORM})
message(STATUS "GEANT3 library directory: ${GEANT3_LIBRARY_DIR}")

find_library(GEANT3_LIB NAMES geant321 PATHS ${ALICE}/geant3/lib/tgt_${GEANT3_PLATFORM} ${GEANT3DIR}/lib/tgt_${GEANT3_PLATFORM})

if (GEANT3_INCLUDE_DIR AND GEANT3_LIB AND GEANT3_SYSTEM_DIR)
   set(GEANT3_FOUND TRUE)
endif (GEANT3_INCLUDE_DIR AND GEANT3_LIB AND GEANT3_SYSTEM_DIR)

if (GEANT3_FOUND)
  message(STATUS "GEANT3 found")
else()
  message(STATUS "GEANT3 not found")
endif (GEANT3_FOUND)
