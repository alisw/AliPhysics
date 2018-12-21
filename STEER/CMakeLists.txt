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
  if(NOT ${CMAKE_SYSTEM_VERSION} MATCHES el5 AND NOT ${ROOT_VERSION_MAJOR} LESS 6)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror")
  endif()
endif(CMAKE_COMPILER_IS_GNUCXX)

# Include the libraries
add_subdirectory (STEERBase)
add_subdirectory (ESD)
add_subdirectory (AOD)
add_subdirectory (CDB)
add_subdirectory (STEER)
add_subdirectory (Utilities)

# Macros
install(DIRECTORY macros DESTINATION STEER)

message(STATUS "STEER enabled")
