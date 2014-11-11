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

cmake_minimum_required(VERSION 2.8.8 FATAL_ERROR)

project(AliRoot CXX C)

# Path to additonal modules
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

#       - CLANG_MAJOR.CLANG_MINOR or
#       - GCC_MAJOR.GCC_MINOR.GCC_PATCH
include(CheckCompiler)

# ROOT dictionaries and maps
include(CMakeALICE)

# ROOT configuration mandatory
if(ROOTSYS)
    find_package(ROOT REQUIRED)

    # ROOT must be build with XML2 support
    if(ROOT_HASXML STREQUAL "no")
        message(FATAL_ERROR "ROOT was not build with xml2 support. Please reinstall or rebuild ROOT with xml2 support")
    endif(ROOT_HASXML STREQUAL "no")
else()
    message(FATAL_ERROR "ROOT installation not found!\nPlease point to the ROOT installation using -DROOTSYS=ROOT_INSTALL_DIR")
endif(ROOTSYS)

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
    message("No Fortran support. Disabling LHAPDF, PHYTIA6, MICROCERN")
endif()


# DATE
find_package(DATE)

# FastJet
find_package(FASTJET)

# AliRoot version extracted from Git repository
# Sets 4 git variables
#  - GIT_REFSPEC - complete name of the current reference
#  - ALIROOT_BRANCH - name of the branch or tag extracted from the current reference
#  - GIT_SHA1 - current hash in the long format
#  - GIT_SHORT_SHA1 - current hash in the short format
include(CheckGitVersion)

message(STATUS "CMake platform: ${CMAKE_SYSTEM}")
message(STATUS "Build folder: ${AliRoot_BINARY_DIR}")
message(STATUS "Source folder: ${AliRoot_SOURCE_DIR}")
message(STATUS "Installation folder: ${CMAKE_INSTALL_PREFIX}")


# AliRoot modules
add_subdirectory(STEER)
add_subdirectory(STAT)
add_subdirectory(OADB)
add_subdirectory(ANALYSIS)
add_subdirectory(RAW)
add_subdirectory(JETAN)
add_subdirectory(ITS)
add_subdirectory(TPC)
add_subdirectory(VZERO)
add_subdirectory(T0)
add_subdirectory(TOF)
add_subdirectory(TRD)
add_subdirectory(ZDC)
add_subdirectory(STRUCT)
add_subdirectory(HMPID)
add_subdirectory(PMD)
add_subdirectory(MUON)
add_subdirectory(FMD)
add_subdirectory(ACORDE)
add_subdirectory(EMCAL)
add_subdirectory(PHOS)
add_subdirectory(TRIGGER)
add_subdirectory(BCM)
add_subdirectory(CORRFW)
add_subdirectory(FASTSIM)
if(CMAKE_Fortran_COMPILER)
  add_subdirectory(LHAPDF)
  add_subdirectory(PYTHIA6)
  add_subdirectory(MICROCERN)
  add_subdirectory(TUHKMgen)
  add_subdirectory(THydjet)
  add_subdirectory(TPHIC)
  add_subdirectory(TEPEMGEN)
  add_subdirectory(THbtp)
endif(CMAKE_Fortran_COMPILER)
add_subdirectory(EVGEN)
add_subdirectory(TTherminator)
# Include Vc own cmake
include(Vc/Vc.cmake)
add_subdirectory(Vc)
add_subdirectory(HLT)
add_subdirectory(STARLIGHT)
add_subdirectory(EVE)
add_subdirectory(AD)
