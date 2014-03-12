# -*- mode: cmake -*-

# Top level CMakeList.txt file for the AliRoot Build System
#
# Author: Johny Jose (johny.jose@cern.ch)
#         Port of previous Makefile build to cmake


cmake_minimum_required(VERSION 2.8.4 FATAL_ERROR)

#--------------------------------------------------------------------------------#
# Set Basic CMake Configuration                                                  #
#--------------------------------------------------------------------------------#

project(AliRoot CXX C Fortran)

# Environment variables
set(ALICE $ENV{ALICE})
if(NOT ALICE)
  message(WARNING "Environment variable ALICE is not set")
endif(NOT ALICE)

set(ALICE_ROOT $ENV{ALICE_ROOT})
if(NOT ALICE_ROOT)
  message(FATAL_ERROR "Please set ALICE_ROOT to the AliRoot source directory")
endif(NOT ALICE_ROOT)

set(ALICE_INSTALL $ENV{ALICE_INSTALL})
# Set AliRoot Installation directory 
if(NOT ALICE_INSTALL)
  message(STATUS "AliRoot files will be installed in ${ALICE_ROOT}. You can specify a different location by setting ALICE_INSTALL environment variable")
  set(ALICE_INSTALL ${ALICE_ROOT})
else()
  message(STATUS "AliRoot files will be installed in ${ALICE_INSTALL}")
endif(NOT ALICE_INSTALL)
if(${ALICE_INSTALL} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "Install dir should be different from build dir")
endif()

# If CMAKE_INSTALL_PREFIX was not set than we set it to ALICE_INSTALL. 
# If ALICE_INSTALL was not defined than it will default to source directory
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX ${ALICE_INSTALL})
endif()

# Environment for SHUTTLE
set(DIMDIR $ENV{DIMDIR})
set(ODIR $ENV{ODIR})

get_filename_component(PROJECT_SOURCE_DIR "${PROJECT_SOURCE_DIR}" REALPATH)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
set(CPACK_COMPONENTS_ALL MDC)

include(ALICE_CMake)

#--------------------------------------------------------------------------------#
# Find Packages                                                                  #
#--------------------------------------------------------------------------------#

# Subversion
# ------------------------------

find_package(Subversion REQUIRED)

# ROOT
# ------------------------------

find_package(ROOT REQUIRED)

# IRST Rule Checker 
# ------------------------------

find_package(RuleChecker)

# DATE
# ------------------------------ 

find_package(DATE)

# GEANT3
# ------------------------------

find_package(GEANT3)

# AMORE
# ------------------------------

find_package(AMORE)

# AliRoot CPack configuration                                                    
#-------------------------------

find_package(CPack)

# ZeroMQ
# ------------------------------
find_package(ZeroMQ)

# fastjet
# ------------------------------
find_package(FASTJET)

#--------------------------------------------------------------------------------#
# Environment/Compiler Configuration                                             #
#--------------------------------------------------------------------------------#

ALICE_CheckOutOfSourceBuild()
ALICE_RootConfig(CMAKE_CXX_COMPILER --cxx)
ALICE_RootConfig(CMAKE_C_COMPILER --cc)
ALICE_RootConfig(CMAKE_Fortran_COMPILER --f77)
ALICE_ConfigurePlatform()
ALICE_ConfigureCompiler()


#--------------------------------------------------------------------------------#
# AliRoot Path Settings                                                          #
#--------------------------------------------------------------------------------#

ALICE_DevFlagsOutput()
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin/tgt_${ALICE_TARGET}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib/tgt_${ALICE_TARGET}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib/tgt_${ALICE_TARGET}")
set(CMAKE_INCLUDE_EXPORT_DIRECTORY "${PROJECT_BINARY_DIR}/include")
file(MAKE_DIRECTORY ${CMAKE_INCLUDE_EXPORT_DIRECTORY})
set(EXPORTDIRS "${ALICE_INSTALL}/include")
set(BINLIBDIRS "-L${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")


#--------------------------------------------------------------------------------#
# AliRoot Modules                                                                #
#--------------------------------------------------------------------------------#

# Exclude Module from all target
# ------------------------------

#set(EXCLUDEMODULES THydjet SHUTTLE)


# Standard Modules
# ------------------------------
set(ALIROOTMODULES STEER PHOS TRD TPC ZDC MUON PMD FMD TOF ITS AD ACORDE HMPID T0 BCM STRUCT VZERO THijing THbtp EMCAL THerwig TEPEMGEN FASTSIM TPHIC RAW MONITOR ANALYSIS JETAN HLT STAT TTherminator CORRFW DPMJET TDPMjet STARLIGHT PWGPP PWG PWGCF PWGGA PWGHF PWGDQ PWGJE PWGLF PWGUD TRIGGER PYTHIA8 TAmpt OADB MFT TEvtGen FIT)

if(THydjet)
  list(APPEND ALIROOTMODULES THydjet)
endif(THydjet)

if(SHUTTLE)
  list(APPEND ALIROOTMODULES SHUTTLE)
endif(SHUTTLE)

# Additional Modules
# ------------------------------
list(APPEND ALIROOTMODULES TUHKMgen EPOS)

# EVE
# ------------------------------
ALICE_RootConfig(EVE --has-opengl)
if(NOT EVE STREQUAL "no")
  list(APPEND ALIROOTMODULES EVE)
endif(NOT EVE STREQUAL "no")

# CERN Modules
# ------------------------------
set(CERNMODULES LHAPDF HIJING MICROCERN HERWIG)

# Check for Pythia 6
if(EXISTS ${ROOTINCDIR}/TPythia6.h)
  set(PYTHIA6 TRUE)
  list(APPEND CERNMODULES PYTHIA6 EVGEN)
elseif(NOT NOPYTHIA)
  message(FATAL_ERROR "ROOT must have Pythia6 enabled !")
endif(EXISTS ${ROOTINCDIR}/TPythia6.h)

# Set Final list of modules 
set(MODULES ${ALIROOTMODULES} ${CERNMODULES} ALIROOT)

# Set Module directories
set(MODDIRS ${MODULES})

# Libraries for Linking
# ------------------------------
set(ALILIBS MUON TPC PMD TRD FMD TOF ITS PHOS AD ACORDE HMPID VZERO ZDC STRUCT T0 EVGEN STEER TRIGGER MFT)
set(LIBS ${ROOTCLIBS} ${ROOTPLIBS} ${SYSLIBS})

# Configure ARVerion.h
# ------------------------------
#If the sources were taken from Git
if(EXISTS ${PROJECT_SOURCE_DIR}/.git/)
 include(GetGitRevisionDescription)
 
 find_package(Git)
 
 if(GIT_FOUND)
  get_git_head_revision(GIT_REFSPEC GIT_SHA1)
  message(STATUS "Aliroot branch - ${GIT_REFSPEC} Revision - ${GIT_SHA1} ")

  STRING(REGEX REPLACE "^(.+/)(.+)/(.*)$" "\\2" BRANCH_TYPE "${GIT_REFSPEC}" )
  
  if(BRANCH_TYPE STREQUAL "heads")
   set(ALIROOT_REVISION "ThisIsaBranchNoRevisionProvided")
  else()
   set(ALIROOT_REVISION ${GIT_SHA1})
  endif()

  STRING(REGEX REPLACE "^(.+/)(.+/)(.*)$" "\\3" SHORT_BRANCH "${GIT_REFSPEC}" )

  set(ALIROOT_BRANCH ${SHORT_BRANCH})

 else()
  message("Git not installed. I can't tell you which revision you are using!")
  endif(GIT_FOUND)
else()
 message("AliRoot sources not downloaded from a Version Control System. I can't tell which revision you are using!")
endif(EXISTS ${PROJECT_SOURCE_DIR}/.git/)

configure_file(${PROJECT_SOURCE_DIR}/cmake/ARVersion.h.tmp ${CMAKE_INCLUDE_EXPORT_DIRECTORY}/ARVersion.h @ONLY)

string (REPLACE "-" "." ALIROOT_BRANCH_SPEC "${ALIROOT_BRANCH}")

#--------------------------------------------------------------------------------#
# Build Modules                                                                  #
#--------------------------------------------------------------------------------#

# Module Variables
# ------------------------------
set(ALLLIBS)
set(ALLEXECS)
set(INCLUDEFILES)
set(BINLIBS)
set(EXPORTFILES)

# Top Level Build Targets
# ------------------------------

add_custom_target(alilibs-static)
add_custom_target(include-headers ALL)
add_custom_target(DA-all)
add_custom_target(check-all) 
add_custom_target(check-hxml)
add_custom_target(smell-all) 
add_custom_target(htmldoc) 
add_custom_target(par-all)
add_custom_target(test-par-all)

#NY file(GLOB CFILES_LIST ${CMAKE_SOURCE_DIR}/macros/*.C)

add_custom_command(TARGET htmldoc
		   PRE_BUILD
		   COMMAND mkdir -p html
		   COMMAND rm -rf html/roothtml
		   COMMAND rm -f html/picts
		   COMMAND cd html && aliroot -q -b mkhtml.C\\\(0,1\\\)
		   WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

# Include Vc, which provides the static libVc and the relevant CMake macros
include(Vc/Vc.cmake)
add_subdirectory(Vc)

# Add Modules	
# ------------------------------
foreach(MODULE ${MODDIRS})
  set(${MODULE}ALIBS)
#RS  configure_file(${PROJECT_SOURCE_DIR}/cmake/Module.tmp ${PROJECT_SOURCE_DIR}/${MODULE}/CMakeLists.txt @ONLY)
  add_custom_target(${MODULE}-all ALL)
  add_custom_target(${MODULE}-par-all)
  add_custom_target(test-${MODULE}-par-all)
  add_custom_target(${MODULE}-all-static)
  add_dependencies(alilibs-static "${MODULE}-static")
  add_subdirectory(${MODULE})
  
#NY  set(dir "${CMAKE_SOURCE_DIR}/${MODULE}")
#NY  if (EXISTS "${dir}")
#NY      file(GLOB CFILES_MOD ${dir}/*.C)
#NY      list(LENGTH CFILES_MOD len)
#NY      if (${len} GREATER 0)
#NY        list(APPEND CFILES_LIST "${CFILES_MOD}")
#NY      endif (${len} GREATER 0)
#NY  endif(EXISTS "${dir}")
endforeach(MODULE)

#NY foreach(cfile ${CFILES_LIST})
#NY     add_custom_command(TARGET htmldoc
#NY 		       PRE_BUILD
#NY 		       COMMAND cd html && 
#NY 		    	       aliroot -b -q \"mkhtml.C(\"${cfile}\")\" > /dev/null
#NY  		       WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
#NY endforeach(cfile)

if(DATE_FOUND)
  include (CMakeDA)
endif(DATE_FOUND)

if(RULECHECKER_FOUND)
    add_custom_target(factfile DEPENDS ${FACTFILE})
    add_custom_command( OUTPUT ${FACTFILE}
                        COMMAND ${CMAKE_COMMAND} -E echo "Re-generating factfile"
                        COMMAND ${JAVA_RUNTIME}  -Xmx1024M -jar ${FACTEXTRACTOR_JAR} ${CMAKE_BINARY_DIR} ${CMAKE_BINARY_DIR}			
			DEPENDS check-hxml ${CMAKE_BINARY_DIR}/check-hxml-touchfile ${_factfile_deps} 
                        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
endif(RULECHECKER_FOUND)

#--------------------------------------------------------------------------------#
# Additional Targets                                                             #
#--------------------------------------------------------------------------------#

# AliMdc RPM
# ------------------------------

execute_process(
COMMAND root-config --libdir 
OUTPUT_VARIABLE ROOTALIBDIR)
string(REPLACE "\n" "" ROOTALIBDIR ${ROOTALIBDIR}) 

execute_process(
COMMAND root-config --arch 
OUTPUT_VARIABLE ARCHITECTURE)


set(AEXT a)
set(LIBPATH ${CMAKE_BINARY_DIR}/lib/tgt_${ALICE_TARGET})
set(BINPATH ${CMAKE_BINARY_DIR}/bin/tgt_${ALICE_TARGET})
set(EXPORTDIR ${ALICE_ROOT}/include)
set(BINLIBDIRS -L${ALICE_ROOT}/${LIBPATH})
set(RAWDIRO ${CMAKE_BINARY_DIR}/RAW)
set(ALIMDCSPECFILE ${RAWDIRO}/alimdc.spec)
set(RAWDatabaseALIB ${LIBPATH}/libRAWDatabase.${AEXT})
set(ESDALIB ${LIBPATH}/libESD.${AEXT})
set(STEERBaseALIB ${LIBPATH}/libSTEERBase.${AEXT})
set(MDCALIB ${LIBPATH}/libMDC.${AEXT})
set(MODDIRO RAW/tgt_${ALICE_TARGET})


if( ${ARCHITECTURE} MATCHES "linuxx8664gcc" )
set(ALIMDCARCHDIR x86_64)
elseif( ${ARCHITECTURE} MATCHES "linux" )
set(ALIMDCARCHDIR i386)
elseif( ${ARCHITECTURE} MATCHES "linuxia64gcc" )
set(ALIMDCARCHDIR ia64)
elseif(NOT DEFINED ${ARCHITECTURE})
message("alimdc-rpm: Unknown architecture" )
else( ${ARCHITECTURE} MATCHES "linuxx8664gcc" )
message("alimdc-rpm: Unknown architecture" )
endif( ${ARCHITECTURE} MATCHES "linuxx8664gcc" )


# Module version will be written in the generated header
configure_file ("RAW/alimdc.spec.in" "${CMAKE_CURRENT_BINARY_DIR}/RAW/alimdc.spec" @ONLY)

add_custom_target( alimdc-rpm
DEPENDS alimdc-static
)

add_custom_command(
TARGET alimdc-rpm
POST_BUILD 
COMMAND rm -rf ${CMAKE_CURRENT_BINARY_DIR}/alimdc-root
COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/alimdc-root/opt/alimdc/lib
COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/alimdc-root/opt/alimdc/include
COMMAND cp ${ALICE_ROOT}/RAW/mdc.h ${CMAKE_CURRENT_BINARY_DIR}/alimdc-root/opt/alimdc/include
COMMAND cp ${LIBPATH}/libAliMDC.a ${CMAKE_CURRENT_BINARY_DIR}/alimdc-root/opt/alimdc/lib
COMMAND cp ${ROOTALIBDIR}/libRoot.a ${CMAKE_CURRENT_BINARY_DIR}/alimdc-root/opt/alimdc/lib
COMMAND cp ${ROOTALIBDIR}/liblzma.a ${CMAKE_CURRENT_BINARY_DIR}/alimdc-root/opt/alimdc/lib
COMMAND cp ${ROOTALIBDIR}/libpcre.a ${CMAKE_CURRENT_BINARY_DIR}/alimdc-root/opt/alimdc/lib
COMMAND cp ${ROOTALIBDIR}/libfreetype.a ${CMAKE_CURRENT_BINARY_DIR}/alimdc-root/opt/alimdc/lib
COMMAND cp ${ROOTALIBDIR}/libpcre.a ${CMAKE_CURRENT_BINARY_DIR}/alimdc-root/opt/alimdc/lib
COMMAND rm -rf ${CMAKE_CURRENT_BINARY_DIR}/RPMS
COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/RPMS/${ALIMDCARCHDIR}/
COMMAND rpmbuild --verbose --define "_topdir ${CMAKE_CURRENT_BINARY_DIR}" --define "%buildroot ${CMAKE_CURRENT_BINARY_DIR}/alimdc-root" -bb ${ALIMDCSPECFILE} 
COMMAND cp -p ${CMAKE_CURRENT_BINARY_DIR}/RPMS/${ALIMDCARCHDIR}/alimdc-*.rpm .
COMMAND rm -rf ${CMAKE_CURRENT_BINARY_DIR}/alimdc-root
COMMAND rm -rf ${CMAKE_CURRENT_BINARY_DIR}/RPMS
COMMAND echo "***** alimdc RPM created and put in ${CMAKE_BINARY_DIR} folder *****"
)

add_custom_command(
TARGET ${RAWDIRO}
COMMAND mkdir-p ${RAWDIRO}
)


add_custom_target( alimdc-static
DEPENDS ${LIBPATH} ${BINPATH} RAW-all-static STEER-all-static )

add_custom_command( 
TARGET alimdc-static
POST_BUILD
COMMAND ar r ${LIBPATH}/alimdccxxo.a ${CMAKE_BINARY_DIR}/RAW/CMakeFiles/alimdc.dir/alimdc_main.cxx.o
COMMAND rm -rf ${LIBPATH}/libAliMDC.a
COMMAND ar x ${RAWDatabaseALIB}
COMMAND ar x ${MDCALIB}
COMMAND ar x ${ESDALIB}
COMMAND ar x ${STEERBaseALIB}
COMMAND ar r ${LIBPATH}/libAliMDC.a *.o
COMMAND ar t ${RAWDatabaseALIB} | xargs rm 
COMMAND ar t ${MDCALIB} | xargs rm 
COMMAND ar t ${ESDALIB} | xargs rm 
COMMAND ar t ${STEERBaseALIB} | xargs rm 
COMMAND ${CMAKE_CXX_COMPILER} ${LDFLAGS} -o ${BINPATH}/alimdca ${LIBPATH}/alimdccxxo.a ${LIBPATH}/libAliMDC.a ${ROOTALIBDIR}/libRoot.a ${ROOTALIBDIR}/libfreetype.a ${ROOTALIBDIR}/libpcre.a ${ROOTALIBDIR}/liblzma.a -pthread -lssl -ldl -lz
COMMAND @echo "Target alimdc-static made"
)

add_custom_command(
TARGET ${LIBPATH}
COMMAND mkdir -p ${LIBPATH}
)


#include (MakefileDA) 
include (InstallRequiredSystemLibraries)
include (CPack)

# Install data
# ------------------------------
install ( DIRECTORY data 
          DESTINATION . )
install ( DIRECTORY OCDB 
	  DESTINATION . 
	  PATTERN ".svn" EXCLUDE)
install ( DIRECTORY QAref 
	  DESTINATION . 
	  PATTERN ".svn" EXCLUDE)
install ( DIRECTORY GRP 
	  DESTINATION . 
	  PATTERN ".svn" EXCLUDE)

#
# EOF
#

