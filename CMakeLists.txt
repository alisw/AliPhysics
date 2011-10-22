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
set(ALICE_ROOT $ENV{ALICE_ROOT})
set(ALICE $ENV{ALICE})
set(ALICE_INSTALL $ENV{ALICE_INSTALL})
if(${ALICE_INSTALL} MATCHES ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "Install dir should be different from build dir")
endif(${ALICE_INSTALL} MATCHES ${CMAKE_BINARY_DIR})
# Environment for SHUTTLE
set(DIMDIR $ENV{DIMDIR})
set(ODIR $ENV{ODIR})

if(NOT ALICE)
  message(WARNING "Environment variable ALICE is not set")
endif(NOT ALICE)

if(NOT ALICE_ROOT)
  message(FATAL_ERROR "Please set ALICE_ROOT to the AliRoot source directory")
endif(NOT ALICE_ROOT)

# Set AliRoot Installation directory 
if(NOT ALICE_INSTALL)
  message(STATUS "AliRoot files will be installed in ${ALICE_ROOT}. You can specify a different location by setting ALICE_INSTALL environment variable")
  set(ALICE_INSTALL ${ALICE_ROOT})
else()
  message(STATUS "AliRoot files will be installed in ${ALICE_INSTALL}")
endif(NOT ALICE_INSTALL)
set(CMAKE_INSTALL_PREFIX ${ALICE_INSTALL})

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
set(ALIROOTMODULES STEER PHOS TRD TPC ZDC MUON PMD FMD TOF ITS ACORDE HMPID T0 BCM STRUCT VZERO THijing THbtp EMCAL THerwig TEPEMGEN FASTSIM TPHIC RAW MONITOR ANALYSIS JETAN HLT STAT TTherminator CORRFW DPMJET TDPMjet PWG0 PWG1 PWG2 PWG3 PWG4 TRIGGER PYTHIA8 TAmpt OADB MFT)

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
set(ALILIBS MUON TPC PMD TRD FMD TOF ITS PHOS ACORDE HMPID VZERO ZDC STRUCT T0 EVGEN STEER TRIGGER MFT)
set(LIBS ${ROOTCLIBS} ${ROOTPLIBS} ${SYSLIBS})

# Configure ARVerion.h
# ------------------------------
if(Subversion_FOUND)
  if(EXISTS ${PROJECT_SOURCE_DIR}/.svn/ )
  Subversion_WC_INFO(${PROJECT_SOURCE_DIR} PROJECT)
  set(ALIROOT_SVN_REVISION ${PROJECT_WC_REVISION})
  set(ALIROOT_SVN_BRANCH ${PROJECT_WC_URL})
  string(REGEX MATCH "[^/]+$" ALIROOT_SVN_BRANCH ${ALIROOT_SVN_BRANCH})

  elseif(EXISTS ${PROJECT_SOURCE_DIR}/.git/ )
    execute_process(WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} COMMAND git svn info OUTPUT_VARIABLE OUTPUTVAR)
    string(REPLACE "\n" ";" OUTPUTVAR ${OUTPUTVAR})
    foreach(f ${OUTPUTVAR})
      if(${f} MATCHES URL)
	string(REGEX MATCH "[^/]+$" ALIROOT_SVN_BRANCH ${f})
      elseif(${f} MATCHES Revision)
	string(REGEX MATCH "[^: ]+$" ALIROOT_SVN_REVISION ${f})
      endif(${f} MATCHES URL)   
    endforeach()
  endif(EXISTS ${PROJECT_SOURCE_DIR}/.svn/ )

  message(STATUS "Aliroot Revision - ${ALIROOT_SVN_REVISION} Branch - ${ALIROOT_SVN_BRANCH}")
  if(EXISTS ${CMAKE_INCLUDE_EXPORT_DIRECTORY}/ARVersion.h)
    file(REMOVE ${CMAKE_INCLUDE_EXPORT_DIRECTORY}/ARVersion.h)
  endif(EXISTS ${CMAKE_INCLUDE_EXPORT_DIRECTORY}/ARVersion.h)
  configure_file(${PROJECT_SOURCE_DIR}/cmake/ARVersion.h.tmp ${CMAKE_INCLUDE_EXPORT_DIRECTORY}/ARVersion.h @ONLY)
endif(Subversion_FOUND)
string (REPLACE "-" "." ALIROOT_SVN_BRANCH_SPEC "${ALIROOT_SVN_BRANCH}")

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
		   COMMAND cd html &&
		   	    aliroot -q -b \"mkhtml.C(0,1)\"
		   WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})


# Add Modules	
# ------------------------------
foreach(MODULE ${MODDIRS})
  set(${MODULE}ALIBS)
  configure_file(${PROJECT_SOURCE_DIR}/cmake/Module.tmp ${PROJECT_SOURCE_DIR}/${MODULE}/CMakeLists.txt @ONLY)
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

if(DATE_FOUND AND AMORE_FOUND)
  include (CMakeDA)
endif(DATE_FOUND AND AMORE_FOUND)

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


add_custom_target( alimdc-rpm
DEPENDS alimdc-static alimdc-specfile
)

add_custom_command(
TARGET alimdc-rpm
POST_BUILD 
COMMAND rm -rf ${ALICE_ROOT}/alimdc-root
COMMAND mkdir -p ${ALICE_ROOT}/alimdc-root/opt/alimdc/lib
COMMAND mkdir -p ${ALICE_ROOT}/alimdc-root/opt/alimdc/include
COMMAND cp ${ALICE_ROOT}/RAW/mdc.h ${ALICE_ROOT}/alimdc-root/opt/alimdc/include
COMMAND cp ${LIBPATH}/libAliMDC.a ${ALICE_ROOT}/alimdc-root/opt/alimdc/lib
COMMAND cp ${ROOTALIBDIR}/libRoot.a ${ALICE_ROOT}/alimdc-root/opt/alimdc/lib
COMMAND cp ${ROOTALIBDIR}/libfreetype.a ${ALICE_ROOT}/alimdc-root/opt/alimdc/lib
COMMAND cp ${ROOTALIBDIR}/libpcre.a ${ALICE_ROOT}/alimdc-root/opt/alimdc/lib
COMMAND rm -rf ${ALICE_ROOT}/RPMS
COMMAND mkdir -p ${ALICE_ROOT}/RPMS/${ALIMDCARCHDIR}/
COMMAND rpmbuild --verbose --define "_topdir ${ALICE_ROOT}" --define "_tmppath ${ALICE_ROOT}" -bb ${ALIMDCSPECFILE} 
COMMAND cp -p ${ALICE_ROOT}/RPMS/${ALIMDCARCHDIR}/alimdc-*.rpm .
COMMAND rm -rf ${ALICE_ROOT}/alimdc-root
COMMAND rm -rf ${ALICE_ROOT}/RPMS
COMMAND echo "***** alimdc RPM created and put in ${CMAKE_BINARY_DIR} folder *****"
)

add_custom_command(
TARGET ${RAWDIRO}
COMMAND mkdir-p ${RAWDIRO}
)


add_custom_target( alimdc-specfile
DEPENDS ${RAWDIRO}
)

add_custom_command( 
TARGET alimdc-specfile
POST_BUILD
COMMAND rm -rf ${ALIMDCSPECFILE}
COMMAND @echo "***** Making alimdc RPM spec-file ${ALIMDCSPECFILE} *****"
COMMAND @echo '\# RPM specfile for alimdc static libs' >> ${ALIMDCSPECFILE} 
COMMAND @echo '\# Package contains both ROOT and AliRoot' >> ${ALIMDCSPECFILE}
COMMAND @echo '\# static libs needed by mStreamRecorder' >> ${ALIMDCSPECFILE}
COMMAND @echo '\# in order to ROOT-ify the incoming raw' >> ${ALIMDCSPECFILE}
COMMAND @echo '\# data' >> ${ALIMDCSPECFILE}
COMMAND @echo '\# Example how-to build alimdc RPM:' >> ${ALIMDCSPECFILE}
COMMAND @echo '\# cd ${ALICE_ROOT}' >> ${ALIMDCSPECFILE}
COMMAND @echo '\# make alimdc-rpm' >> ${ALIMDCSPECFILE}
COMMAND @echo "" >> ${ALIMDCSPECFILE}
COMMAND @echo "Summary: AliMDC static libraries" >> ${ALIMDCSPECFILE}
COMMAND @echo "Name: alimdc" >> ${ALIMDCSPECFILE}
COMMAND @echo "Version:  ${ALIROOT_SVN_BRANCH_SPEC}" >> ${ALIMDCSPECFILE}
COMMAND @echo "Release: ${ALIROOT_SVN_REVISION}" >> ${ALIMDCSPECFILE}
COMMAND @echo '\# Copyright: CERN Alice Off-line' >> ${ALIMDCSPECFILE}
COMMAND @echo "License: CERN Alice Off-line" >> ${ALIMDCSPECFILE}
COMMAND @echo "Vendor: ALICE Core Off-line Group" >> ${ALIMDCSPECFILE}
COMMAND @echo "URL: http://aliceinfo.cern.ch" >> ${ALIMDCSPECFILE}
COMMAND @echo "Group: Applications/Alice" >> ${ALIMDCSPECFILE}
COMMAND @echo "Prefix: /opt/%{name}" >> ${ALIMDCSPECFILE}
COMMAND @echo "BuildRoot: %{_tmppath}/%{name}-root" >> ${ALIMDCSPECFILE}
COMMAND @echo "" >> ${ALIMDCSPECFILE}
COMMAND @echo '\# automatic dependencies' >> ${ALIMDCSPECFILE}
COMMAND @echo "AutoReqProv: yes" >> ${ALIMDCSPECFILE}
COMMAND @echo "" >> ${ALIMDCSPECFILE}
COMMAND @echo '\# list here required RPM packages for runtime' >> ${ALIMDCSPECFILE}
COMMAND @echo "Requires: glibc" >> ${ALIMDCSPECFILE}
COMMAND @echo "" >> ${ALIMDCSPECFILE}
COMMAND @echo "Provides: alimdc" >> ${ALIMDCSPECFILE}
COMMAND @echo "" >> ${ALIMDCSPECFILE}
COMMAND @echo '\# description of the package' >> ${ALIMDCSPECFILE}
COMMAND @echo "%description" >> ${ALIMDCSPECFILE}
COMMAND @echo "Package contains both ROOT and AliRoot" >> ${ALIMDCSPECFILE}
COMMAND @echo "static libs needed by mStreamRecorder" >> ${ALIMDCSPECFILE}
COMMAND @echo "in order to ROOT-ify the incoming raw" >> ${ALIMDCSPECFILE}
COMMAND @echo "data. The package version correspond to" >> ${ALIMDCSPECFILE}
COMMAND @echo "the AliRoot one." >> ${ALIMDCSPECFILE}
COMMAND @echo "" >> ${ALIMDCSPECFILE}
COMMAND @echo '\# list of files to be installed' >> ${ALIMDCSPECFILE}
COMMAND @echo "%files" >> ${ALIMDCSPECFILE}
COMMAND @echo '%defattr (-,root,root)' >> ${ALIMDCSPECFILE}
COMMAND @echo "%{prefix}/lib/libAliMDC.a" >> ${ALIMDCSPECFILE}
COMMAND @echo "%{prefix}/lib/libRoot.a" >> ${ALIMDCSPECFILE}
COMMAND @echo "%{prefix}/lib/libpcre.a" >> ${ALIMDCSPECFILE}
COMMAND @echo "%{prefix}/lib/libfreetype.a" >> ${ALIMDCSPECFILE}
COMMAND @echo "%{prefix}/include/mdc.h" >> ${ALIMDCSPECFILE}
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
COMMAND ${CMAKE_CXX_COMPILER} ${LDFLAGS} -o ${BINPATH}/alimdca ${LIBPATH}/alimdccxxo.a ${LIBPATH}/libAliMDC.a ${ROOTALIBDIR}/libRoot.a ${ROOTALIBDIR}/libfreetype.a ${ROOTALIBDIR}/libpcre.a -pthread -ldl -lz
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


#
# EOF
#

