# **************************************************************************
# * Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
# *                                                                        *
# * Author: The ALICE Off-line Project.                                    *
# * Contributors are mentioned in the code where appropriate.              *
# *                                                                        *
# * Permission to use, copy, modify and distribute this software and its   *
# * documentation strictly for non-commercial purposes is hereby granted   *
# * without fee, provided that the above copyright notice appears in all   *
# * copies and that both the copyright notice and this permission notice   *
# * appear in the supporting documentation. The authors make no claims     *
# * about the suitability of this software for any purpose. It is          *
# * provided "as is" without express or implied warranty.                  *
# **************************************************************************

# Treat warnings as errors
# with Clang
if(${CMAKE_CXX_COMPILER_ID} MATCHES Clang)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror")
endif(${CMAKE_CXX_COMPILER_ID} MATCHES Clang)
# with gcc
if(CMAKE_COMPILER_IS_GNUCXX)
  if (NOT ${CMAKE_SYSTEM_VERSION} MATCHES el5)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror")
  endif (NOT ${CMAKE_SYSTEM_VERSION} MATCHES el5)
endif(CMAKE_COMPILER_IS_GNUCXX)

add_subdirectory(ANALYSIS)
add_subdirectory(ANALYSISalice)
add_subdirectory(ANALYSISaliceBase)
add_subdirectory(ANALYSIScalib)
add_subdirectory(ESDfilter)

# Installing the macros
install(DIRECTORY macros DESTINATION ANALYSIS)
install(DIRECTORY BadChunkFilter DESTINATION ANALYSIS)

message(STATUS "ANALYSIS enabled")
