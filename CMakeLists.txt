#--------------------------------------------------------------------------#
# Set Basic CMake Configuration                                            #
#--------------------------------------------------------------------------#

cmake_minimum_required(VERSION 2.8.11 FATAL_ERROR)

project(AliPhysics CXX C)

message(STATUS "CMake platform: ${CMAKE_SYSTEM}")
message(STATUS "Build folder: ${AliPhysics_BINARY_DIR}")
message(STATUS "Source folder: ${AliPhysics_SOURCE_DIR}")
message(STATUS "Installation folder: ${CMAKE_INSTALL_PREFIX}")

# You can change the build type using 
# cmake -DCMAKE_BUILD_TYPE=DEBUG | RELEASE | RELWITHDEBINFO | MINSIZEREL ...
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RELWITHDEBINFO)
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")


if(CMAKE_INSTALL_PREFIX STREQUAL "${AliPhysics_SOURCE_DIR}")
    message(FATAL_ERROR "Please choose a different installation point than the source tree!")
endif()

# Path to additonal modules
set(CMAKE_MODULE_PATH "${AliPhysics_SOURCE_DIR}/cmake")

# AliRoot version extracted from Git repository
# Sets 4 git variables
#  - GIT_REFSPEC - complete name of the current reference
#  - ALIROOT_BRANCH - name of the branch or tag extracted from the current reference
#  - GIT_SHA1 - current hash in the long format
#  - GIT_SHORT_SHA1 - current hash in the short format
#  - ALIROOT_VERSION = ALIROOT_BRANCH
#  - ALIROOT_REVISION = GIT_SHORT_SHA1
include(CheckGitVersion)

#       - CLANG_MAJOR.CLANG_MINOR or
#       - GCC_MAJOR.GCC_MINOR.GCC_PATCH
include(CheckCompiler)

# Shared library suffix
if (NOT CMAKE_SYSTEM_NAME STREQUAL Windows)
  set(CMAKE_SHARED_LIBRARY_SUFFIX .so)
endif (NOT CMAKE_SYSTEM_NAME STREQUAL Windows)

# ROOT configuration mandatory
if(ROOTSYS)
    find_package(ROOT REQUIRED)

    # ROOT must be build with XML2 support
    if(NOT ROOT_HASXML)
        message(FATAL_ERROR "ROOT was not build with xml2 support. Please reinstall or rebuild ROOT with xml2 support")
    endif(NOT ROOT_HASXML)
else()
    message(FATAL_ERROR "ROOT installation not found!\nPlease point to the ROOT installation using -DROOTSYS=ROOT_INSTALL_DIR")
endif(ROOTSYS)

# ROOT dictionaries and maps
include(CMakeALICE)

# Checking first for DIM, DATE, AMORE and daqDA 
# in case ROOT Extra static library needs to be enabled

# DATE
# date-config needs DIMDIR and ODIR set
if(DATE_CONFIG)
    if(DIMDIR AND ODIR)
        find_package(DATE)
    else()
        message(FATAL_ERROR "DATE enabled but no DIMDIR and ODIR set. Please set DIMDIR and ODIR")
    endif()
endif()

# If no Fortran, i.e on Windows
# We need to specify ROOT fortran
# (f95 comes before gfortran in default module)
include(CheckLanguage)

# Optional: Doxygen for reference doc generation
find_package(Doxygen)
if(DOXYGEN_FOUND AND DOXYGEN_DOT_FOUND)
    message(STATUS "Doxygen ${DOXYGEN_VERSION} and Graphviz found: ${DOXYGEN_EXECUTABLE}, ${DOXYGEN_DOT_EXECUTABLE}")
    set(CAN_GENERATE_DOC TRUE)
else()
    message(STATUS "Doxygen and Graphviz not found. Disabling support for documentation generation")
endif(DOXYGEN_FOUND AND DOXYGEN_DOT_FOUND)

# FastJet
find_package(FASTJET)

# ZEROMQ
find_package(ZeroMQ)

# General flags -> Should be moved into a configuration file
set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)
set(CMAKE_NO_SYSTEM_FROM_IMPORTED TRUE)

# AliRoot modules
add_subdirectory(CORRFW)
add_subdirectory(EVEANALYSIS)
add_subdirectory(EVENTMIX)
add_subdirectory(HLTANALYSIS)
add_subdirectory(JETAN)
add_subdirectory(OADB)
add_subdirectory(TENDER)

# PWG libraries
add_subdirectory(PWG)
add_subdirectory(PWGCF)
add_subdirectory(PWGGA)
add_subdirectory(PWGDQ)
add_subdirectory(PWGHF)
add_subdirectory(PWGJE)
add_subdirectory(PWGLF)
add_subdirectory(PWGPP)
add_subdirectory(PWGUD)
