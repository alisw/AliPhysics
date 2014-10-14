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

# If no Fortran, i.e on Windows
include(CheckLanguage)
check_language(Fortran)
if(CMAKE_Fortran_COMPILER)
    enable_language(Fortran OPTIONAL)
else()
    message("No Fortran support. Disabling LHAPDF, PHYTIA6, MICROCERN")
endif()

# Path to additonal modules
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
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

# DATE
find_package(DATE)

# Configure ARVerion.h using Git informatiion
# Sets 3 git variables
#  - GIT_REFSPEC - complete name of the current reference
#  - ALIROOT_BRANCH - name of the branch or tag extracted from the current reference
#  - GIT_SHA1 - current hash in the long format
#  - GIT_SHORT_SHA1 - current hash in the short format
if(EXISTS ${PROJECT_SOURCE_DIR}/.git/)
    include(GetGitRevisionDescription)
    
    find_package(Git)
    
    if(GIT_FOUND)
        get_git_head_revision(GIT_REFSPEC GIT_SHA1)

        # GIT_REFSPEC is empty for detached mode = tags in detached mode or checkout to specific hash

        # returns the closest reference to the current hash
        # name of the current tag or heads/branch in the case of branches
        git_describe(ALIROOT_GIT_TAG "--all" "--abbrev=0")

        STRING(REGEX REPLACE "^(.+/)(.+)/(.*)$" "\\2" BRANCH_TYPE "${GIT_REFSPEC}" )
        
        # the revision is not set in the case of a branch, it means we are doing development
        # and the revision will trigger a reconfiguration
        if(BRANCH_TYPE STREQUAL "heads")
            set(ALIROOT_REVISION "ThisIsaBranchNoRevisionProvided")
            STRING(REGEX REPLACE "^(.+/)(.+/)(.*)$" "\\3" SHORT_BRANCH "${GIT_REFSPEC}" )
        else()
            set(SHORT_BRANCH ${ALIROOT_GIT_TAG})
            set(ALIROOT_REVISION ${GIT_SHA1})
        endif()

        set(ALIROOT_BRANCH ${SHORT_BRANCH})
  
        # generate the short version of the revision hash
        execute_process(COMMAND git rev-parse --short ${GIT_SHA1} 
                          WORKING_DIRECTORY ${CMAKE_SOURCE_DIR} 
                          OUTPUT_STRIP_TRAILING_WHITESPACE 
                          RESULT_VARIABLE res
                          OUTPUT_VARIABLE GIT_SHORT_SHA1)

        # if the rev-parse fails we set the short sha to the long initial one
        if(NOT res EQUAL 0)
            set(GIT_SHORT_SHA1 ${GIT_SHA1})
        endif()
  
        message(STATUS "Aliroot branch/tag: \"${ALIROOT_BRANCH}\" - Revision:  \"${GIT_SHORT_SHA1}\" ")

    else()
        message(STATUS "Git not installed. I can't tell you which revision you are using!")
    endif(GIT_FOUND)
else()
    message("AliRoot sources not downloaded from a Version Control System. I can't tell which revision you are using!")
endif(EXISTS ${PROJECT_SOURCE_DIR}/.git/)

configure_file(${PROJECT_SOURCE_DIR}/cmake/ARVersion.h.tmp ${CMAKE_BINARY_DIR}/version/ARVersion.h @ONLY)
install(FILES ${PROJECT_BINARY_DIR}/version/ARVersion.h DESTINATION include)

string (REPLACE "-" "." ALIROOT_BRANCH_SPEC "${ALIROOT_BRANCH}")

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
add_subdirectory(ITS)
add_subdirectory(TPC)
add_subdirectory(VZERO)
add_subdirectory (T0)
add_subdirectory (TOF)
add_subdirectory (TRD)
add_subdirectory (ZDC)
#add_subdirectory (STRUCT)
#add_subdirectory (HMPID)
#add_subdirectory (PMD)
#add_subdirectory (MUON)
#add_subdirectory (FMD)
#add_subdirectory (ACORDE)
#add_subdirectory (EMCAL)
#add_subdirectory (PHOS)
#add_subdirectory (TRIGGER)
#add_subdirectory (BCM)
#add_subdirectory (CORRFW)
#add_subdirectory (FASTSIM)
#if(CMAKE_Fortran_COMPILER_WORKS)
#  add_subdirectory (LHAPDF)
#  add_subdirectory (PYTHIA6)
#  add_subdirectory (MICROCERN)
#endif(CMAKE_Fortran_COMPILER_WORKS)
#add_subdirectory (EVGEN)
#add_subdirectory (TUHKMgen)
#add_subdirectory (TTherminator)
#add_subdirectory (TPHIC)
#add_subdirectory (THydjet)
#add_subdirectory (PWG)
#add_subdirectory (Vc)
#add_subdirectory (HLT)
