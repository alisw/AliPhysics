# -*- mode: cmake -*-

# Check if cmake has the required version
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

# Set name of our project to AliRoot. To be done after check of cmake version
project(AliRoot NONE)

# Set AliRoot variables
set(CMAKE_VERBOSE_MAKEFILES OFF)
#fix for MACOSX
SET(CMAKE_INCLUDE_SYSTEM_FLAG_C "-isystem")
SET(CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem") 

set(ALICE_INSTALL $ENV{ALICE_INSTALL})
set(ALICE_ROOT $ENV{ALICE_ROOT})
set(ALICE_STATIC_BUILD OFF CACHE BOOL "Build also statics libs")
execute_process(COMMAND root-config --arch OUTPUT_VARIABLE ALICE_TARGET)
string (REGEX REPLACE "\n" "" ALICE_TARGET ${ALICE_TARGET})

if(NOT ALICE_INSTALL)
message(FATAL_ERROR "Please set environment variable ALICE_INSTALL to the AliRoot installation directory")
endif(NOT ALICE_INSTALL)
if(NOT ALICE_ROOT)
message(FATAL_ERROR "Please set environment variable ALICE_ROOT to the AliRoot source directory")
endif(NOT ALICE_ROOT)


message("ALICE_ROOT: ${ALICE_ROOT}")


# where to look first for cmake modules, before ${CMAKE_ROOT}/Modules/ 
set(PROJECT_BINARY_DIR ${ALICE_INSTALL})
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin/tgt_${ALICE_TARGET})
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib/tgt_${ALICE_TARGET})
set(CMAKE_MODULE_PATH "${ALICE_ROOT}/cmake/modules")

# To fix a hack that has been introduced for the Mac
if(APPLE)
  set(CMAKE_INSTALL_NAME_TOOL CMAKE_INSTALL_NAME_TOOL-NOTFOUND)
endif(APPLE)

# Setup ROOT and C, C++ , Fortran Compilers
find_package(ROOT 5.0.0 REQUIRED)
set(CMAKE_C_COMPILER $ENV{CC})
enable_language(C)
set(CMAKE_CXX_COMPILER $ENV{CXX})
enable_language(CXX)
set(CMAKE_Fortran_COMPILER $ENV{F77})
enable_language(Fortran)

# Load some basic macros which are needed later on
#Include(WriteConfigFile)

include (CTest)
include(ALICEMacros)

# Check if the user wants to build the project in the source directory
Check_out_of_Source_Build()

# searches for needed packages

Set(ALIROOT_INSTALL_DIR ${ALICE_INSTALL}/${ALICE_TARGET})

# Setup GEANT3
find_package(GEANT3)

# Setup RuleChecker
find_package(RuleChecker)

# Check if DATE is installed
find_package(DATE)

# Setup system dependent flags 
SetupSystem()

Include(CreateSubCMakeFiles)
CreateSubCMakeFiles()


# Set the library version in the main CMakeLists.txt

set(ALIROOT_MAJOR_VERSION 0)
set(ALIROOT_MINOR_VERSION 0)
set(ALIROOT_PATCH_VERSION 0)
set(ALIROOT_VERSION "${ALIROOT_MAJOR_VERSION}.${ALIROOT_MINOR_VERSION}.${ALIROOT_PATCH_VERSION}")

set(ALIROOT_LIBRARY_PROPERTIES ${ALIROOT_LIBRARY_PROPERTIES}
    VERSION "${ALIROOT_VERSION}"
    SOVERSION "${ALIROOT_MAJOR_VERSION}"
    SUFFIX ".so"
)

# Set version variables
find_program( READLINK_PROGRAM readlink)

if(READLINK_PROGRAM) 
  execute_process(COMMAND readlink -f ${ALICE_ROOT} 
    COMMAND xargs svn info
    OUTPUT_VARIABLE _svn_out)
else(READLINK_PROGRAM)
  execute_process(COMMAND svn info ${ALICE_ROOT} 
    OUTPUT_VARIABLE _svn_out)
endif(READLINK_PROGRAM)

string(REGEX REPLACE "^.*Revision: ([^\n]*).*$" "\\1" 
  ALIROOT_SVN_REVISION ${_svn_out})
string(REGEX REPLACE "^.*URL: ([^\n]*).*$" "\\1"
  ALIROOT_SVN_BRANCH ${_svn_out})
string(REGEX REPLACE "^.*AliRoot/([^\n]*)" "\\1" ALIROOT_SVN_BRANCH ${ALIROOT_SVN_BRANCH})


# Additional targets
add_custom_target(check-all COMMENT "Checking coding conventions")
message(STATUS "REVISION: ${ALIROOT_SVN_REVISION}")
message(STATUS "BRANCH: ${ALIROOT_SVN_BRANCH}")

# Delete and recreate header ARVersion Header file 
execute_process(COMMAND cmake -E remove "${CMAKE_CURRENT_BINARY_DIR}/STEER/ARVersion.h")
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/STEER/ARVersion.h.in ${CMAKE_CURRENT_BINARY_DIR}/STEER/ARVersion.h @ONLY)

# Get list of enabled root features
execute_process(COMMAND root-config --features OUTPUT_VARIABLE FEATURES)
string(REGEX MATCH "pythia6" PYTHIA6 ${FEATURES})
string(REGEX MATCH "pythia8" PYTHIA8 ${FEATURES})
string(REGEX MATCH "opengl" OPENGL ${FEATURES})
if(PYTHIA6)
  message(STATUS "Pythia6 enabled")
endif(PYTHIA6)
if(PYTHIA8)
  message(STATUS "Pythia8 enabled")
endif(PYTHIA8)
if(OPENGL)
  message(STATUS "OpenGL enabled")	
endif(OPENGL)	
option (TFLUKA "Build TFluka" OFF)
option (THYDJET "Build THydjet" OFF)
option (SHUTTLE "Build SHUTTLE" OFF)


add_custom_target(STEER-all)
add_custom_target(PHOS-all)	
add_custom_target(TRD-all) 
add_custom_target(TPC-all)
add_custom_target(ZDC-all) 
add_custom_target(MUON-all) 
add_custom_target(PMD-all) 
add_custom_target(FMD-all)
add_custom_target(TOF-all) 
add_custom_target(ITS-all) 
add_custom_target(ACORDE-all) 
add_custom_target(HMPID-all) 
add_custom_target(T0-all) 
add_custom_target(BCM-all) 
add_custom_target(STRUCT-all) 
add_custom_target(EVGEN-all) 
add_custom_target(RALICE-all) 
add_custom_target(VZERO-all) 
add_custom_target(THijing-all)
add_custom_target(THbtp-all) 
add_custom_target(EMCAL-all) 
add_custom_target(THerwig-all) 
add_custom_target(TEPEMGEN-all) 
add_custom_target(FASTSIM-all) 
add_custom_target(TPHIC-all) 
add_custom_target(RAW-all) 
add_custom_target(MONITOR-all) 
add_custom_target(ANALYSIS-all) 
add_custom_target(JETAN-all)
add_custom_target(HLT-all) 
add_custom_target(LHC-all) 
add_custom_target(ESDCheck-all) 
add_custom_target(STAT-all) 
add_custom_target(TTherminator-all) 
add_custom_target(CORRFW-all) 
add_custom_target(DPMJET-all) 
add_custom_target(TDPMjet-all) 
add_custom_target(PWG0-all) 
add_custom_target(PWG1-all) 
add_custom_target(PWG2-all) 
add_custom_target(PWG3-all) 
add_custom_target(PWG4-all)
add_custom_target(TRIGGER-all) 
add_custom_target(TUHKMgen-all) 
add_custom_target(EPOS-all) 
add_custom_target(PYTHIA8-all) 
add_custom_target(EVE-all) 
add_custom_target(TFluka-all)
add_custom_target(THydjet-all)
add_custom_target(SHUTTLE-all)
add_custom_target(LHAPDF-all) 
add_custom_target(HIJING-all) 
add_custom_target(MICROCERN-all) 
add_custom_target(HERWIG-all) 
add_custom_target(PYTHIA6-all)
add_custom_target(ALIROOT-all)



# Recurse into the given subdirectories.  This does not actually
# cause another cmake executable to run.  The same process will walk through
# the project's entire directory structure.
Add_subdirectory(STEER)
Add_subdirectory(PHOS) 
Add_subdirectory(TRD) 
Add_subdirectory(TPC)
Add_subdirectory(ZDC) 
Add_subdirectory(MUON) 
Add_subdirectory(PMD) 
Add_subdirectory(FMD)
Add_subdirectory(TOF) 
Add_subdirectory(ITS) 
Add_subdirectory(ACORDE) 
Add_subdirectory(HMPID) 
Add_subdirectory(T0) 
Add_subdirectory(BCM) 
Add_subdirectory(STRUCT) 
Add_subdirectory(EVGEN) 
Add_subdirectory(RALICE) 
Add_subdirectory(VZERO) 
Add_subdirectory(THijing)
Add_subdirectory(THbtp) 
Add_subdirectory(EMCAL) 
Add_subdirectory(THerwig) 
Add_subdirectory(TEPEMGEN) 
Add_subdirectory(FASTSIM) 
Add_subdirectory(TPHIC) 
Add_subdirectory(RAW) 
Add_subdirectory(MONITOR) 
Add_subdirectory(ANALYSIS) 
Add_subdirectory(JETAN)
Add_subdirectory(HLT) 
Add_subdirectory(LHC) 
Add_subdirectory(ESDCheck) 
Add_subdirectory(STAT) 
Add_subdirectory(TTherminator) 
Add_subdirectory(CORRFW) 
Add_subdirectory(DPMJET) 
Add_subdirectory(TDPMjet) 
Add_subdirectory(PWG0) 
Add_subdirectory(PWG1) 
Add_subdirectory(PWG2) 
Add_subdirectory(PWG3) 
Add_subdirectory(PWG4)
Add_subdirectory(TRIGGER) 
Add_subdirectory(TUHKMgen) 
Add_subdirectory(EPOS) 
if(PYTHIA8)
Add_subdirectory(PYTHIA8) 
endif(PYTHIA8)
if(OPENGL)
Add_subdirectory(EVE) 
endif(OPENGL)
if(TFLUKA)
Add_subdirectory(TFluka)
endif(TFLUKA)
if(THYDJET)	
Add_subdirectory(THydjet)
endif(THYDJET)
if(SHUTTLE)
Add_subdirectory(SHUTTLE)
endif(SHUTTLE)	
Add_subdirectory(LHAPDF) 
Add_subdirectory(HIJING) 
Add_subdirectory(MICROCERN) 
Add_subdirectory(HERWIG) 
if(PYTHIA6)
Add_subdirectory(PYTHIA6)
endif(PYTHIA6)
Add_subdirectory(ALIROOT)
add_subdirectory(test)
add_test ("Aliroot_executable" aliroot -q) 
#copy testing files
configure_file("CTestCustom.ctest" "${PROJECT_BINARY_DIRECTORY}" COPYONLY)
