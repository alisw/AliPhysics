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

#--------------------------------------------------------------------------#
# Set Basic CMake Configuration                                            #
#--------------------------------------------------------------------------#

cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)

# AliRoot will refuse to compile against a wrong ROOT version
set(ROOT_VERSION_MIN_5 "5.34.30")
set(ROOT_VERSION_MIN_6 "6.04.00")

### CMP0025   Compiler id for Apple Clang is now AppleClang.
### CMP0042   MACOSX_RPATH is enabled by default.

foreach(p
  CMP0025 # CMake 3.0
  CMP0042 # CMake 3.0
  )
  if(POLICY ${p})
  cmake_policy(SET ${p} NEW)
  endif()
endforeach()

project(AliRoot CXX C)

message(STATUS "CMake platform: ${CMAKE_SYSTEM}")
message(STATUS "Build folder: ${AliRoot_BINARY_DIR}")
message(STATUS "Source folder: ${AliRoot_SOURCE_DIR}")
message(STATUS "Installation folder: ${CMAKE_INSTALL_PREFIX}")

# CMake supports different build types by default. We want the DEBUG build type
# to have "-g -O0" flags: by default it only has "-g"
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")

# You can change the build type using
# cmake -DCMAKE_BUILD_TYPE=DEBUG | RELEASE | RELWITHDEBINFO | MINSIZEREL ...
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RELWITHDEBINFO)
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build type: ${CMAKE_BUILD_TYPE} (${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}})")

if(CMAKE_INSTALL_PREFIX STREQUAL "${AliRoot_SOURCE_DIR}")
  message(FATAL_ERROR "Please choose a different installation point than the source tree!")
endif()

# Path to additonal modules
set(CMAKE_MODULE_PATH "${AliRoot_SOURCE_DIR}/cmake")

# AliRoot version extracted from Git repository
# Sets 4 git variables
#  - GIT_REFSPEC - complete name of the current reference
#  - ALIROOT_BRANCH - name of the branch or tag extracted from the current reference
#  - GIT_SHA1 - current hash in the long format
#  - GIT_SHORT_SHA1 - current hash in the short format
#  - ALIROOT_VERSION = ALIROOT_BRANCH
#  - ALIROOT_REVISION = GIT_SHORT_SHA1
include(CheckGitVersion)

# - CLANG_MAJOR.CLANG_MINOR or
# - GCC_MAJOR.GCC_MINOR.GCC_PATCH
include(CheckCompiler)

# Utility to generate PARfiles
include(cmake/GenParFile.cmake)

# Shared library suffix
if (NOT CMAKE_SYSTEM_NAME STREQUAL Windows)
  set(CMAKE_SHARED_LIBRARY_SUFFIX .so)
endif (NOT CMAKE_SYSTEM_NAME STREQUAL Windows)

# Optional: Doxygen for reference doc generation
find_package(Doxygen)
if(DOXYGEN_FOUND AND DOXYGEN_DOT_FOUND)
  message(STATUS "Doxygen ${DOXYGEN_VERSION} and Graphviz found: ${DOXYGEN_EXECUTABLE}, ${DOXYGEN_DOT_EXECUTABLE}")
  set(CAN_GENERATE_DOC TRUE)
elseif(DOXYGEN_ONLY)
  message(FATAL_ERROR "Doxygen and Graphviz not found. Cannot generate documentation")
else()
  message(STATUS "Doxygen and Graphviz not found. Disabling support for documentation generation")
endif()

# Doxygen documentation is always included
add_subdirectory(doxygen)

if(DEFINED DOXYGEN_ONLY AND DOXYGEN_ONLY)
  # In this case we do not prepare targets for the rest, but for Doxygen only.
  # This is useful because no checks are performed, and it means that Doxygen
  # documentation can be generated on a relatively "dumb" server that has as
  # requirements only cmake, doxygen and graphviz.
  message(WARNING "Preparing build for Doxygen documentation only as requested. AliRoot will not be built!")
else()

  # Standard CMake inclusions and checks if we are actually building AliRoot,
  # i.e. not only the documentation

  # Checking first for DIM, DATE, AMORE and daqDA 
  # in case ROOT Extra static library needs to be enabled

  # DATE
  # date-config needs DIMDIR and ODIR set
  # Setting variables in cache to be accesible by ccmake
  set(DATE_CONFIG CACHE STRING "date-config script location")
  set(DIMDIR CACHE STRING "dim installation folder")
  set(ODIR CACHE STRING "dim plaform, usually \"linux\"")

  if(DIMDIR AND ODIR)
    find_package(DIM)
  endif()

  if(DATE_CONFIG)
    if(DIM_FOUND)
      find_package(DATE)
    else()
      message(FATAL_ERROR "DATE enabled but no DIMDIR and ODIR set. Please set DIMDIR and ODIR")
    endif()
  endif()

  # daqDA
  # Setting variables in cache to be accesible by ccmake
  set(daqDA CACHE STRING "daqDA library folder installation")
  if(daqDA)
    find_package(daqDA)
  endif(daqDA)

  # AMORE -> requires ROOT to be installed
  # Setting variables in cache to be accesible by ccmake
  set(AMORE_CONFIG CACHE STRING "amore-config script location")
  if(AMORE_CONFIG)
    if(ROOTSYS)
      find_package(AMORE)
    else()
      message(FATAL_ERROR "AMORE enabled but no ROOTSYS defined")
    endif()
  endif(AMORE_CONFIG)

  # id DARPM we enable DA if it is not set
  # Setting variables in cache to be accesible by ccmake
  set(DARPM CACHE STRING "Enable DARPM: ON")
  set(DA CACHE STRING "Enable DA: ON")
  if(DARPM)
    if(NOT DA)
      set(DA ON)
    endif(NOT DA)
  endif(DARPM)

  # DA is enabled
  # Setting variables in cache to be accesible by ccmake
  set(ALIROOT_STATIC CACHE STRING "ENABLE static building of AliRoot: ON")
  if(DA)
    if(NOT DIM_FOUND)
      set(DA FALSE)
      message(FATAL_ERROR "DAs enabled but no DIMDIR and ODIR set. Please set DIMDIR to DIM installation and ODIR to platform (default linux)")
    endif()

    if(NOT DATE_FOUND)
      set(DA FALSE)
      message(FATAL_ERROR "DAs enabled but no DATE support found. Please point to your date installation using \"DATE_CONFIG\" variable")
    endif()

    if(NOT daqDA_FOUND)
      set(DA FALSE)
      message(FATAL_ERROR "DAs enabled but no daqDA support found. Please point to your daqDA installation using \"daqDA\" variable")
    endif()

    if(NOT AMORE_FOUND)
      set(DA FALSE)
      message(FATAL_ERROR "DAs enabled but no AMORE support found. Please point to your AMORE installation using \"AMORE_CONFIG\" variable")
    endif()
      
    # Enable static libraries
    set(ALIROOT_STATIC TRUE)
    message(STATUS "DAs enabled")
  endif(DA)

  # MDC rpm creation enables the static build
  # Setting variables in cache to be accesible by ccmake
  set(MDCRPM CACHE STRING "Enable mdc rpm build")
  if(MDCRPM)
    set(ALIROOT_STATIC TRUE)
    message(STATUS "AliMDC RPM enabled. AliRoot static build enabled")
  endif(MDCRPM)

  # ROOT configuration mandatory
  # FindRoot has to be after the enabling of DA due to libRootExtra
  # Do not put it above DA
  # Setting variables in cache to be accesible by ccmake
  set(ROOTSYS CACHE STRING "ROOT installation folder")
  if(ROOTSYS)
    find_package(ROOT REQUIRED)

    if(ROOT_VERSION_MAJOR EQUAL 5)
      set(ROOT_VERSION_MIN "${ROOT_VERSION_MIN_5}")
    elseif(ROOT_VERSION_MAJOR EQUAL 6)
      set(ROOT_VERSION_MIN "${ROOT_VERSION_MIN_6}")
    else()
      message(FATAL_ERROR "AliRoot works with ROOT 5 or 6 only")
    endif()

    if(ROOT_VERSION_NORM VERSION_LESS "${ROOT_VERSION_MIN}")
      message(FATAL_ERROR "AliRoot requires at least ROOT ${ROOT_VERSION_MIN} (you have ${ROOT_VERSION_NORM})")
    endif()

    # ROOT must be built with XML2 support
    if(NOT ROOT_HASXML)
      message(FATAL_ERROR "ROOT was not built with xml2 support. Please reinstall or rebuild ROOT with xml2 support")
    endif(NOT ROOT_HASXML)


  else(ROOTSYS)
    message(FATAL_ERROR "ROOT installation not found!\nPlease point to the ROOT installation using -DROOTSYS=ROOT_INSTALL_DIR")
  endif(ROOTSYS)

  # FTGL: prefer library from ROOT, look for both libFTGL and libftgl
  unset(LIB_FTGL)
  if(${ROOT_FEATURES} MATCHES "builtin_ftgl")
    message(STATUS "Using FTGL from ROOT")
    find_library(LIB_FTGL NAMES FTGL PATHS "${ROOT_LIBDIR}" NO_DEFAULT_PATH)
  else()
    message(STATUS "Using FTGL from the system")
    find_library(LIB_FTGL NAMES ftgl FTGL)
  endif()
  if(NOT LIB_FTGL)
    message(FATAL_ERROR "Required library FTGL not found. Please install it, or reconfigure ROOT with --enable-builtin-ftgl")
  endif()

  # ROOT dictionaries and maps
  include(CMakeALICE)

  # If no Fortran, i.e on Windows
  # We need to specify ROOT fortran
  # (f95 comes before gfortran in default module)
  include(CheckLanguage)

  if(ROOT_FORTRAN)
    message(STATUS "Using the Fortran compiler defined by ROOT configuration: ${ROOT_FORTRAN}")
    set(CMAKE_Fortran_COMPILER ${ROOT_FORTRAN})
  else()
    message(STATUS "Using default system Fortran compiler")
  endif(ROOT_FORTRAN)

  check_language(Fortran)
  if(CMAKE_Fortran_COMPILER)
    enable_language(Fortran OPTIONAL)
  else()
    message(STATUS "No Fortran support. Disabling LHAPDF, PHYTIA6, MICROCERN, etc.")
  endif()

  # FastJet
  find_package(FASTJET)

  # ZEROMQ
  find_package(ZeroMQ)

  # Generating the AliRoot-config.cmake file
  configure_file(${PROJECT_SOURCE_DIR}/cmake/AliRoot-config.cmake.in ${CMAKE_BINARY_DIR}/version/AliRoot-config.cmake @ONLY)
  install(FILES ${PROJECT_BINARY_DIR}/version/AliRoot-config.cmake DESTINATION etc)

  # General flags -> Should be moved into a configuration file
  set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)
  set(CMAKE_NO_SYSTEM_FROM_IMPORTED TRUE)

  # List of modules for which PARfiles are enabled -- use the cache as global scope and clean it
  set(ALIPARFILES "" CACHE INTERNAL "ALIPARFILES" FORCE)

  # AliRoot base modules
  add_subdirectory(STEER)
  add_subdirectory(RAW)
  add_subdirectory(ANALYSIS)

  # Include Vc own cmake
  include(Vc/Vc.cmake)
  if (NOT ROOT_HASVC)
    add_subdirectory(Vc)
  endif (NOT ROOT_HASVC)

  # AliRoot common macros
  add_subdirectory(macros)

  # AliRoot modules
  add_subdirectory(ACORDE)
  add_subdirectory(AD)
  add_subdirectory(BCM)
  add_subdirectory(EMCAL)
  add_subdirectory(EPOS)
  add_subdirectory(EVE)
  add_subdirectory(EVGEN)
  add_subdirectory(FASTSIM)
  add_subdirectory(FIT)
  add_subdirectory(FMD)
  add_subdirectory(HLT)
  add_subdirectory(HMPID)
  add_subdirectory(ITS)
  add_subdirectory(MFT)
  add_subdirectory(MONITOR)
  add_subdirectory(MUON)
  add_subdirectory(PHOS)
  add_subdirectory(PMD)
  add_subdirectory(PYTHIA8)
  add_subdirectory(STARLIGHT)
  add_subdirectory(STAT)
  add_subdirectory(STRUCT)
  add_subdirectory(T0)
  add_subdirectory(TDPMjet)
  add_subdirectory(TEvtGen)
  add_subdirectory(THerwig)
  add_subdirectory(TOF)
  add_subdirectory(TPC)
  add_subdirectory(TRD)
  add_subdirectory(TRIGGER)
  add_subdirectory(TTherminator)
  add_subdirectory(VZERO)
  add_subdirectory(ZDC)
  add_subdirectory(ALIGN)

  # Fortran modules
  if(CMAKE_Fortran_COMPILER)
    add_subdirectory(DIME)
    add_subdirectory(DPMJET)
    add_subdirectory(HERWIG)
    add_subdirectory(HIJING)
    add_subdirectory(LHAPDF)
    add_subdirectory(MICROCERN)
    add_subdirectory(PYTHIA6)
    add_subdirectory(TEPEMGEN)
    add_subdirectory(THbtp)
    add_subdirectory(THijing)
    add_subdirectory(THydjet)
    add_subdirectory(TPHIC)
    add_subdirectory(TUHKMgen)
    add_subdirectory(TAmpt)
  endif(CMAKE_Fortran_COMPILER)

  add_subdirectory(data)
  add_subdirectory(GRP)
  add_subdirectory(OCDB)
  add_subdirectory(QAref)

  # Enable SHUTTLE compilation
  # Check if DIMDIR and ODIR are set
  # # Setting variables in cache to be accesible by ccmake
  set(SHUTTLE CACHE STRING "Enable SHUTTLE build : ON")
  if(SHUTTLE)
    if(NOT ROOT_HASMONALISA)
        message(FATAL_ERROR "ROOT was not build with MonAlisa support")
    endif()

    if(NOT ROOT_HASLDAP)
        message(FATAL_ERROR "ROOT was not build with LDAP support")
    endif()

    find_package(ApMon)

    if(ROOT_HASALIEN STREQUAL "no")
      message(FATAL_ERROR "Shuttle needs ROOT build with AliEn support. Please build ROOT with AliEn support. Do not forget to set ALIEN to your AliEn installation")
    endif()

    if(DIM_FOUND AND ALIEN)
        add_subdirectory(SHUTTLE)
    else()
      message(FATAL_ERROR "SHUTTLE not enabled! Please specify DIMDIR, ODIR and ALIEN")
    endif()
  endif(SHUTTLE)

  # Enable SHUTTLE/TestShuttle compilation
  # independently from SHUTTLE (and from monalisa, ldap, dim)
  # Setting variables in cache to be accesible by ccmake
  set(TESTSHUTTLE CACHE STRING "Enable TestShuttle build : ON")
  if(TESTSHUTTLE)
    message(STATUS "TestShuttle enabled")
    add_subdirectory(SHUTTLE/TestShuttle)
  endif(TESTSHUTTLE)

  add_subdirectory(ALIROOT)

  # List modules with PARfiles
  string(REPLACE ";" " " ALIPARFILES_FLAT "${ALIPARFILES}")
  message(STATUS "PARfile target enabled for the following modules: ${ALIPARFILES_FLAT}")

  # Install Find macros for common packages to be further used by AliPhysics
  install(FILES
    cmake/CMakeALICE.cmake
    cmake/FindROOT.cmake
    cmake/FindFASTJET.cmake
    cmake/FindZeroMQ.cmake
    cmake/GetGitRevisionDescription.cmake
    cmake/GetGitRevisionDescription.cmake.in
    DESTINATION cmake)

endif(DEFINED DOXYGEN_ONLY AND DOXYGEN_ONLY)
