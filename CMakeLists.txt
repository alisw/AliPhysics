# -*- mode: cmake -*-

# Check if cmake has the required version
Cmake_Minimum_Required(VERSION 2.6 FATAL_ERROR)

# Set name of our project to AliRoot. To be done after check of cmake version
Project(AliRoot NONE)

# Set AliRoot variables
Set(ALICE_ROOT $ENV{ALICE_ROOT})
Set(ALICE      $ENV{ALICE})
Set(ALICE_STATIC_BUILD OFF CACHE BOOL "Build also statics libs")

MESSAGE("ALICE_ROOT: ${ALICE_ROOT}")
Set(CMAKE_VERBOSE_MAKEFILE ON)

# where to look first for cmake modules, before ${CMAKE_ROOT}/Modules/ 
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
Set(CMAKE_MODULE_PATH "${ALICE_ROOT}/cmake/modules")

# To fix a hack that has been introduced for the Mac
If(APPLE)
  Set(CMAKE_INSTALL_NAME_TOOL CMAKE_INSTALL_NAME_TOOL-NOTFOUND)
Endif(APPLE)

# Setup ROOT
Find_package(ROOT 5.0.0 REQUIRED)

Set(CMAKE_C_COMPILER $ENV{CC})
Enable_Language(C)

Set(CMAKE_CXX_COMPILER $ENV{CXX})
Enable_Language(CXX)

Set(CMAKE_Fortran_COMPILER $ENV{F77})
Enable_Language(Fortran)

# Load some basic macros which are needed later on
#Include(WriteConfigFile)
Include(Dart)
Include(ALICEMacros)

# Check if the user wants to build the project in the source directory
Check_out_of_Source_Build()

# searches for needed packages

Set(ALIROOT_INSTALL_DIR ${ALICE_ROOT}/${ALICE_TARGET})

# Setup GEANT3
Find_package(GEANT3 REQUIRED)

# Setup RuleChecker
Find_Package(RuleChecker)

# Check if DATE is installed
Find_Package(DATE)

# Setup system dependent flags 
SetupSystem()


Option(GENERATE_FILES "Generate CMake input files from package files" OFF)
if(GENERATE_FILES)
  Include(CreateSubCMakeFiles)
  CreateSubCMakeFiles()
endif(GENERATE_FILES)

# Set the library version in the main CMakeLists.txt
SET(ALIROOT_MAJOR_VERSION 0)
SET(ALIROOT_MINOR_VERSION 0)
SET(ALIROOT_PATCH_VERSION 0)
SET(ALIROOT_VERSION "${ALIROOT_MAJOR_VERSION}.${ALIROOT_MINOR_VERSION}.${ALIROOT_PATCH_VERSION}")
SET(ALIROOT_LIBRARY_PROPERTIES ${ALIROOT_LIBRARY_PROPERTIES}
    VERSION "${ALIROOT_VERSION}"
    SOVERSION "${ALIROOT_MAJOR_VERSION}"
    SUFFIX ".so"
)

# Set version variables
Find_program( READLINK_PROGRAM readlink)

If(READLINK_PROGRAM) 
  Execute_process(COMMAND readlink -f ${ALICE_ROOT} 
    COMMAND xargs svn info
    OUTPUT_VARIABLE _svn_out)
Else(READLINK_PROGRAM)
  Execute_process(COMMAND svn info ${ALICE_ROOT} 
    OUTPUT_VARIABLE _svn_out)
Endif(READLINK_PROGRAM)

String(REGEX REPLACE "^.*Revision: ([^\n]*).*$" "\\1" 
  ALIROOT_SVN_REVISION ${_svn_out})
String(REGEX REPLACE "^.*URL: https://alisoft.cern.ch/AliRoot/([^\n]*).*$" "w \\1"
  ALIROOT_SVN_BRANCH ${_svn_out})
String(REPLACE "/" "\\\\/"
  ALIROOT_SVN_BRANCH ${ALIROOT_SVN_BRANCH})

# Additional targets
Add_Custom_Target(check-all COMMENT "Checking coding conventions")

# Configure special include file
# Only if there is  a change in the file
If(UNIX)
  Set(_in ${CMAKE_CURRENT_SOURCE_DIR}/STEER/ARVersion.h.in)
  Set(_out ${CMAKE_CURRENT_BINARY_DIR}/STEER/ARVersion.h)
  If(NOT EXISTS ${_out})
    Add_custom_target(ARversion ALL
      COMMAND sed < ${_in}
                -e "s/\\@ALIROOT_SVN_REVISION\\@/${ALIROOT_SVN_REVISION}/" 
                -e "s/\\@ALIROOT_SVN_BRANCH\\@/${ALIROOT_SVN_BRANCH}/" 
                -e "s/\\@ALIROOT_VERSION\\@/${ALIROOT_VERSION}/" >  ${_out})
    Install(FILES ${_out} DESTINATION ${ALIROOT_INSTALL_DIR}/include)
  Else(NOT EXISTS ${_out})
    FILE(READ "${_out}" contents)
    STRING(REGEX REPLACE ";" "\\\\;" contents "${contents}")
    STRING(REGEX REPLACE "\n" ";" contents "${contents}")
    Foreach(line ${contents})
      If(line MATCHES "^.*ALIROOT_SVN_REVISION.*$")
        STRING(REGEX REPLACE "^.*ALIROOT_SVN_REVISION[ ](.*)$" "\\1" _result_svn_rev "${line}")
      EndIf(line MATCHES "^.*ALIROOT_SVN_REVISION.*$")
      If(line MATCHES "^.*ALIROOT_SVN_BRANCH.*$")
        STRING(REGEX REPLACE "^.*ALIROOT_SVN_BRANCH[ ]\"(.*)\"$" "\\1" _result_svn_branch "${line}")
      EndIf(line MATCHES "^.*ALIROOT_SVN_BRANCH.*$")
      If(line MATCHES "^.*ALIROOT_VERSION.*$")
        STRING(REGEX REPLACE "^.*ALIROOT_VERSION[ ](.*)$" "\\1" _result_aliroot_version "${line}")
      EndIf(line MATCHES "^.*ALIROOT_VERSION.*$")
    EndForeach(line ${contents})
    If(_result_aliroot_version STREQUAL ${ALIROOT_VERSION} AND _result_svn_rev STREQUAL ${ALIROOT_SVN_REVISION} AND _result_svn_branch STREQUAL ${ALIROOT_SVN_BRANCH})
#      MESSAGE("VERSION: ${_result_aliroot_version}")
#      MESSAGE("VERSION: ${ALIROOT_VERSION}")
#      MESSAGE("REVISION: ${_result_svn_rev}")
#      MESSAGE("REVISION: ${ALIROOT_SVN_REVISION}")
#      MESSAGE("BRANCH: ${_result_svn_branch}")
#      MESSAGE("BRANCH: ${ALIROOT_SVN_BRANCH}")
    Else(_result_aliroot_version STREQUAL ${ALIROOT_VERSION} AND _result_svn_rev STREQUAL ${ALIROOT_SVN_REVISION} AND _result_svn_branch STREQUAL ${ALIROOT_SVN_BRANCH})
    Add_custom_target(ARversion ALL
      COMMAND sed < ${_in}
                -e "s/\\@ALIROOT_SVN_REVISION\\@/${ALIROOT_SVN_REVISION}/" 
                -e "s/\\@ALIROOT_SVN_BRANCH\\@/${ALIROOT_SVN_BRANCH}/" 
                -e "s/\\@ALIROOT_VERSION\\@/${ALIROOT_VERSION}/" >  ${_out})
      Install(FILES ${_out} DESTINATION ${ALIROOT_INSTALL_DIR}/include)
    EndIf(_result_aliroot_version STREQUAL ${ALIROOT_VERSION} AND _result_svn_rev STREQUAL ${ALIROOT_SVN_REVISION} AND _result_svn_branch STREQUAL ${ALIROOT_SVN_BRANCH})
  EndIf(NOT EXISTS ${_out})

Endif(UNIX)

# This would work, but just once... we want this to happen at every make
#  Configure_file(${CMAKE_CURRENT_SOURCE_DIR}/STEER/ARVersion.h.in 
#${CMAKE_CURRENT_BINARY_DIR}/STEER/ARVersion.h @ONLY)

# Recurse into the given subdirectories.  This does not actually
# cause another cmake executable to run.  The same process will walk through
# the project's entire directory structure.

#Add_subdirectory(PYTHIA8) 
#Add_subdirectory(TFluka) 
#Add_subdirectory(THydjet) 

Add_subdirectory(ACORDE) 
Add_subdirectory(ANALYSIS) 
Add_subdirectory(BCM) 
Add_subdirectory(CORRFW) 
Add_subdirectory(DPMJET) 
Add_subdirectory(EMCAL) 
Add_subdirectory(EPOS) 
Add_subdirectory(ESDCheck) 
Add_subdirectory(EVE) 
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
Add_subdirectory(PYTHIA6)
Add_subdirectory(RALICE) 
Add_subdirectory(RAW) 
Add_subdirectory(STAT) 
Add_subdirectory(STEER)
Add_subdirectory(STRUCT) 
Add_subdirectory(T0) 
Add_subdirectory(TDPMjet) 
Add_subdirectory(TEPEMGEN) 
Add_subdirectory(THbtp) 
Add_subdirectory(THerwig) 
Add_subdirectory(THijing)
#Add_subdirectory(TIsajet) 
Add_subdirectory(TOF) 
Add_subdirectory(TPC)
Add_subdirectory(TPHIC) 
Add_subdirectory(TRD) 
Add_subdirectory(TRIGGER) 
Add_subdirectory(TTherminator) 
Add_subdirectory(TUHKMgen) 
Add_subdirectory(VZERO) 
Add_subdirectory(ZDC) 
#
Add_subdirectory(ALIROOT)

