# -*- mode: cmake -*-

# Check if cmake has the required version
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

# Set name of our project to AliRoot. To be done after check of cmake version
project(AliRoot NONE)

# Set AliRoot variables
set(CMAKE_VERBOSE_MAKEFILES OFF)
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

	
# Recurse into the given subdirectories.  This does not actually
# cause another cmake executable to run.  The same process will walk through
# the project's entire directory structure.
Add_subdirectory(ACORDE) 
Add_subdirectory(ANALYSIS) 
Add_subdirectory(BCM) 
Add_subdirectory(CORRFW) 
Add_subdirectory(DPMJET) 
Add_subdirectory(EMCAL) 
Add_subdirectory(EPOS) 
Add_subdirectory(ESDCheck) 
if(OPENGL)
Add_subdirectory(EVE) 
endif(OPENGL)
Add_subdirectory(EVGEN) 
Add_subdirectory(FASTSIM) 
Add_subdirectory(FMD)
Add_subdirectory(HERWIG) 
Add_subdirectory(HIJING) 
Add_subdirectory(HLT) 
Add_subdirectory(HMPID) 
Add_subdirectory(ITS) 
Add_subdirectory(JETAN)
Add_subdirectory(LHAPDF) 
Add_subdirectory(LHC) 
Add_subdirectory(MICROCERN) 
Add_subdirectory(MONITOR) 
Add_subdirectory(MUON) 
Add_subdirectory(PHOS) 
Add_subdirectory(PMD) 
Add_subdirectory(PWG0) 
Add_subdirectory(PWG1) 
Add_subdirectory(PWG2) 
Add_subdirectory(PWG3) 
Add_subdirectory(PWG4)
if(PYTHIA6)
Add_subdirectory(PYTHIA6)
Add_subdirectory(RALICE) 
Add_subdirectory(TPHIC) 
endif(PYTHIA6)
if(PYTHIA8)
Add_subdirectory(PYTHIA8) 
endif(PYTHIA8)
Add_subdirectory(RAW) 
if(SHUTTLE)
Add_subdirectory(SHUTTLE)
endif(SHUTTLE)	
Add_subdirectory(STAT) 
Add_subdirectory(STEER)
Add_subdirectory(STRUCT) 
Add_subdirectory(T0) 
Add_subdirectory(TDPMjet) 
Add_subdirectory(TEPEMGEN) 
Add_subdirectory(THbtp) 
Add_subdirectory(THerwig) 
Add_subdirectory(THijing)
if(TFLUKA)
Add_subdirectory(TFluka)
endif(TFLUKA)
if(THYDJET)	
Add_subdirectory(THydjet)
endif(THYDJET)
Add_subdirectory(TOF) 
Add_subdirectory(TPC)
Add_subdirectory(TRD) 
Add_subdirectory(TRIGGER) 
Add_subdirectory(TTherminator) 
Add_subdirectory(TUHKMgen) 
Add_subdirectory(VZERO) 
Add_subdirectory(ZDC) 
Add_subdirectory(ALIROOT)
add_subdirectory(test)
add_test ("Aliroot_executable" aliroot -q) 
#copy testing files
if(EXISTS "${PROJECT_SOURCE_DIRECTORY}/CTestCustom.cmake")
configure_file("CTestCustom.cmake" "${PROJECT_BINARY_DIRECTORY}" COPYONLY)
endif(EXISTS "${PROJECT_SOURCE_DIRECTORY}/CTestCustom.cmake")
execute_process(COMMAND "cmake" "-E" "copy_directory" "${ALICE_ROOT}/test" "${ALICE_INSTALL}/test")

message("Note : In order for the warnings test to execute successfully you need to log your make build output to make.log in the ALICE_INSTALL Directory using this command \n make -k > make.log 2>&1")
