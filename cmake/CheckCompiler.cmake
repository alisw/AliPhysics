
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

# Compiler major and minor version
#       - CLANG_MAJOR.CLANG_MINOR or
#       - GCC_MAJOR.GCC_MINOR.GCC_PATCH

message(STATUS "Found ${CMAKE_CXX_COMPILER_ID} compiler, version ${CMAKE_CXX_COMPILER_VERSION}")

# Clang compiler
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    string(REGEX REPLACE "^.*[ ]([0-9]+)\\.[0-9].*$" "\\1" CLANG_MAJOR "${CMAKE_CXX_COMPILER_VERSION}")
    string(REGEX REPLACE "^.*[ ][0-9]+\\.([0-9]).*$" "\\1" CLANG_MINOR "${CMAKE_CXX_COMPILER_VERSION}")
    message(STATUS "Compiler MAJOR ${CLANG_MAJOR}, MINOR ${CLANG_MINOR}")
endif()

# GNU compiler
if(CMAKE_COMPILER_IS_GNUCXX)
    string(REGEX REPLACE "^([0-9]+).*$"                   "\\1" GCC_MAJOR "${CMAKE_CXX_COMPILER_VERSION}")
    string(REGEX REPLACE "^[0-9]+\\.([0-9]+).*$"          "\\1" GCC_MINOR "${CMAKE_CXX_COMPILER_VERSION}")
    string(REGEX REPLACE "^[0-9]+\\.[0-9]+\\.([0-9]+).*$" "\\1" GCC_PATCH "${CMAKE_CXX_COMPILER_VERSION}")
    message(STATUS "Compiler major ${GCC_MAJOR}, minor ${GCC_MINOR}, patch ${GCC_PATCH}")
endif()
