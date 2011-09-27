#--------------------------------------------------------------------
# CMakeDA.cmake file for AliRoot Build System
#
# Author: Anshul Goel (anshul.goel@cern.ch)
#         Port of previous Makefile build to DAQ targets and RPM
#

cmake_minimum_required(VERSION 2.8.4 FATAL_ERROR)

# -----------Utilities--------------------

#list of detectors 
#set(ONLINEDETECTORS T00)
set(ONLINEDETECTORS SPD SDD SSD ACO GRP TST HLT EMC TRI T00 PMD CPV PHS FMD TPC TRD ZDC V00 MTR MCH HMP TOF)
function(expand output input)
    string(REGEX MATCH "\\\${[^}]*}" m "${input}")
    while(m)
        string(REGEX REPLACE "\\\${(.*)}" "\\1" v "${m}")
        string(REPLACE "\${${v}}" "${${v}}" input "${input}")
        string(REGEX MATCH "\\\${[^}]*}" m "${input}")
    endwhile()
    set("${output}" "${input}" PARENT_SCOPE)
endfunction()

#function to get module for detector
function (detector_module _module detector)
  #Map of online detectors to DA in pairs of ONLINEDETECTORNAME DAMODULE
  set (DETECTORMODULEMAP SPD ITS SDD ITS SSD ITS HMP HMPID PHS PHOS CPV PHOS MCH MUON MTR MUON T00 T0 V00 VZERO ACO ACORDE EMC EMCAL)
  list(FIND DETECTORMODULEMAP ${detector} _index)
  if(_index STREQUAL "-1")
    set(${_module} "${detector}" PARENT_SCOPE)
  else()
    math(EXPR _index "${_index}+1")
    list(GET DETECTORMODULEMAP ${_index} _index)
    set(${_module} ${_index} PARENT_SCOPE)
  endif(_index STREQUAL "-1")
endfunction()

#function to get subDA for detector
function (detector_subda _subda detector)
  #Map of online detectors to SUBDAMODULE in pairs of ONLINEDETECTORNAME SUBDAMODULE
  set (DETECTORSUBDAMAP SPD SPD SDD SDD SSD SSD CPV CPV MCH TRK MTR TRG)
  list(FIND DETECTORSUBDAMAP ${detector} _index)
  if(_index STREQUAL "-1")
    set(${_subda} "" PARENT_SCOPE)
  else()
    math(EXPR _index "${_index}+1")
    list(GET DETECTORSUBDAMAP ${_index} _index)
    set(${_subda} ${_index} PARENT_SCOPE)
  endif(_index STREQUAL "-1")
endfunction()

#function to extract info
function (getinfo _info pattern file)
  string(REGEX MATCH "${pattern}:[^\n]*" _match ${file})
  if(_match)
    string(REGEX REPLACE "${pattern}:[ ]*" "" _data ${_match})
    string(REPLACE "(" "\\(" _data "${_data}")
    string(REPLACE ")" "\\)" _data "${_data}")
    string(REPLACE "<" "\\<" _data "${_data}")
    string(REPLACE ">" "\\>" _data "${_data}")
  endif(_match)
  set(${_info} ${_data} PARENT_SCOPE)
endfunction()


foreach(oldetect ${ONLINEDETECTORS})
detector_module(h_module ${oldetect})
list(APPEND mod "-I${ALICE_ROOT}/${h_module}")
endforeach(oldetect ${ONLINEDETECTORS})
list(APPEND mod "-I${ALICE_ROOT}/include" "-I${ALICE_ROOT}/STEER" "-I${ALICE_ROOT}/ANALYSIS" "-I${ALICE_ROOT}/RAW" "-I${ALICE_ROOT}/STEER/STEER" "-I${ALICE_ROOT}/STEER/CDB" "-I${ROOTSYS}/include" "-I${ALICE_ROOT}/STEER/STEERBase" "-I${ALICE_ROOT}/STEER/AOD" "-I${ALICE_ROOT}/STEER/ESD" "-I${ALICE_ROOT}/MUON/mapping" "-I$ENV{AMORE}/include/amore")

# ----------Common stuff-------------------

file(GLOB_RECURSE _dafiles $ENV{ALICE_ROOT}/*da.cxx)
set (DAINSTALL "$ENV{ALICE_INSTALL}/DA")
if(DAINSTALL STREQUAL "/DA") 
  set(DAINSTALL "$ENV{ALICE_ROOT}/DA")
endif(DAINSTALL STREQUAL "/DA")


file(MAKE_DIRECTORY ${DAINSTALL})

string (REPLACE "-pedantic-errors" "" CXXFLAGS ${CXXFLAGS})

find_program(XML2 xml2-config)

if(AMORE_FOUND)
  #Set compiler flags
  set(CXXFLAGS "${CXXFLAGS} ${AMOREFLAGS}")
  set(CFLAGS "${CFLAGS} ${AMOREFLAGS}")
  set(CINTFLAGS "${CINTFLAGS} ${AMOREFLAGS}")
else()
  set(AMORELIBS)
endif(AMORE_FOUND)

execute_process(COMMAND ${XML2} --libs OUTPUT_VARIABLE _xml2libs OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND ${DATE_CONFIG} --monitorlibs=noshift OUTPUT_VARIABLE MONITORLIBS OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND $ENV{AMORE}/amore-config --ldflags-da-static OUTPUT_VARIABLE _amore1 OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND date-config --rcproxylibs OUTPUT_VARIABLE _amore2 OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND date-config --monitorlibs=noshift OUTPUT_VARIABLE _monitor1 OUTPUT_STRIP_TRAILING_WHITESPACE)
string(REPLACE "\n" " " _amore1 "${_amore1}") 
set(AMOREDALIBS "-static ${_amore1} ${_amore2}")
set(MONITORLIBS "${_monitor1}")

separate_arguments(MONITORLIBS)
separate_arguments(AMOREDALIBS)
set(SYSLIBS -ldl -lpthread ${_xml2libs})


set(EXTRAROOTLIB "libRootExtra.a")
set(ROOTLIB "RootExtra")
file(GLOB _extraroot "$ENV{ROOTSYS}/montercarlo/vmc/src/*.o" "$ENV{ROOTSYS}/tree/treeplayer/src/*.o" "$ENV{ROOTSYS}/io/xmlparser/src/*.o" "$ENV{ROOTSYS}/math/minuit2/src/*.o")

add_library(${ROOTLIB} STATIC ${_extraroot})      
set_target_properties(${ROOTLIB} PROPERTIES LINKER_LANGUAGE CXX)

set(DAQDALIB_PATH $ENV{DAQDALIB_PATH})	
if(DAQDALIB_PATH)
  set(DAQDADIR "${DAQDALIB_PATH}")
else()
  set(DAQDADIR "$ENV{ALICE}/daqDAlib")
endif(DAQDALIB_PATH)
set(DAQDALIB "${DAQDADIR}/libdaqDA.a")

include_directories(${DAQDADIR} ${ALICE_ROOT}/RAW ${ALICE_ROOT}/include ${ALICE_ROOT}/STEER)
include_directories(SYSTEM ${ROOTINCDIR})

# ----------Create All Valid targets---------
	  
foreach(detector ${ONLINEDETECTORS} )

  set(ONLINEDETECTORNAME ${detector})

  detector_module(DAMODULE ${ONLINEDETECTORNAME})

  detector_subda(SUBDAMODULE ${ONLINEDETECTORNAME})
  
  #ALIROOTALIBS

  set(ALIROOTALIBS)
	set(BASIC_TARGET "daqDA-${ONLINEDETECTORNAME}-all")
	add_custom_target(${BASIC_TARGET}
	WORKING_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
	) 
	set(BASIC_RPM "daqDA-${ONLINEDETECTORNAME}-all-rpm")
	add_custom_target(${BASIC_RPM}
	WORKING_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
        )
 list(APPEND ALIROOTALIBS RAWDatabase-static RAWDatarec-static RAWDatasim-static STEERBase-static STEER-static CDB-static ESD-static STAT-static AOD-static ANALYSIS-static ANALYSISalice-static )

  expand(ALIROOTALIBS2 "\${${DAMODULE}ALIBS}")
  expand(DAINCDIRS "\${${DAMODULE}INC}")
  list(APPEND ALIROOTALIBS ${ALIROOTALIBS2})
  
  include_directories(${DAMODULE} ${SUBDIR} ${DAINCDIRS})
#Get detector algorithms for this detector

  foreach(dafile ${_dafiles})

	string(REGEX MATCH "$ENV{ALICE_ROOT}/${DAMODULE}/${DAMODULE}${SUBDAMODULE}" match ${dafile})
#Found a valid target name
	if(match)
          string(REGEX REPLACE "$ENV{ALICE_ROOT}/${DAMODULE}/${DAMODULE}${SUBDAMODULE}(.*)da\\.cxx" "\\1" DANAME ${dafile})
#Check for default DA 
	  if(DANAME)
	    set(DATARGETNAME "daqDA-${ONLINEDETECTORNAME}-${DANAME}")
	else()
	    set(DATARGETNAME "daqDA-${ONLINEDETECTORNAME}")
	  endif(DANAME)
	add_dependencies(${BASIC_TARGET} ${DATARGETNAME})  
	add_dependencies(${BASIC_RPM} ${DATARGETNAME}-rpm) 
	 set(DATARGETDIR "${DAINSTALL}/${DAMODULE}/tgt_$ENV{ALICE_TARGET}")
	  file(MAKE_DIRECTORY ${DATARGETDIR})
	  set(DAOBJ "${DATARGETDIR}/${DAMODULE}${SUBDAMODULE}${DANAME}da.o")
	  set(DASRC "${DAMODULE}/${DAMODULE}${SUBDAMODULE}${DANAME}da.cxx")
	  set(DALIB "${DAMODULE}${SUBDAMODULE}${DANAME}DA")
	  set(DAEXE "${DAMODULE}${SUBDAMODULE}${DANAME}da.exe")
	  set(DADEP "${DATARGETDIR}/${DAMODULE}${SUBDAMODULE}${DANAME}da.d") 
	# DAVERSION
	  execute_process(COMMAND svn info $ENV{ALICE_ROOT}/${DASRC} OUTPUT_VARIABLE _daversion OUTPUT_STRIP_TRAILING_WHITESPACE)
	  string(REGEX REPLACE ".*Last Changed Rev: ([^\n]+)\n.*" "\\1" DAVERSION ${_daversion})

	# DAROOTRELEASE 
	  execute_process(COMMAND root-config --version OUTPUT_VARIABLE _darootrelease OUTPUT_STRIP_TRAILING_WHITESPACE)
	  string(REGEX REPLACE "/" "." DAROOTRELEASE ${_darootrelease})
	
	# DAALIROOTRELEASE
	  string(REGEX REPLACE ".*URL: .*/(.+)/${DASRC}.*$" "\\1" DAALIROOTRELEASE ${_daversion})
          string (REPLACE "-" "." DAALIROOTRELEASE "${DAALIROOTRELEASE}")
	
	  set(DAARCNAME "${DATARGETNAME}")
	  string(REPLACE "-" "" DAARCNAME "${DAARCNAME}")
	  set(DAARC "${DAARCNAME}-${DAVERSION}")
	  set(DATAR "${DATARGETDIR}/${DAARC}.src.tar.gz")
	  set(DASPECFILE "${DATARGETDIR}/${DAMODULE}${SUBDAMODULE}${DANAME}da.spec")
	  set(DAMAKEFILE "${DATARGETDIR}/${DAMODULE}${SUBDAMODULE}${DANAME}da.make")

  
	  file(READ "$ENV{ALICE_ROOT}/${DASRC}" _dasrc)
	  getinfo(DACONTACT "Contact" ${_dasrc})
	  getinfo(DALINKPAGE "Link" ${_dasrc})
	  getinfo(DAREFRUN "Reference Run" ${_dasrc})
	  getinfo(DARUNTYPE "Run Type" ${_dasrc})
	  getinfo(DATYPE "DA Type" ${_dasrc})
	  getinfo(DANUMBEROFEVENTS "Number of events needed" ${_dasrc})
  	  getinfo(DAINPUTFILES "Input Files" ${_dasrc})
	  getinfo(DAOUTPUTFILES "Output Files" ${_dasrc})
	  getinfo(DATRIGGERTYPE "Trigger types used" ${_dasrc})

	  if(NOT _dasrc)
	   message(FATAL_ERROR "SRC FILE MISSING")
	  endif(NOT _dasrc)

#----------- Targets ----------------------

	  set(CMAKE_CXX_FLAGS "${CXXFLAGS}")
	  set(CMAKE_C_FLAGS "${CFLAGS}")
	  set(CMAKE_Fortran_FLAGS ${FFLAGS})
	  set(CMAKE_SHARED_LINKER_FLAGS ${SOFLAGS}) 
	  set(CMAKE_MODULE_LINKER_FLAGS ${LDFLAGS})

add_dependencies(DA-all ${DATARGETNAME})
add_custom_target(${DATARGETNAME}
WORKING_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
)
add_dependencies(${DATARGETNAME} ${DAEXE})

set(ZIP)
foreach(_lib ${ALIROOTALIBS})
string(REGEX REPLACE "-static" "" _lib ${_lib})
list(APPEND ZIP && ar x "../lib${_lib}.a")
endforeach(_lib)
list (APPEND ZIP && ar r "../lib${DALIB}.a" "*.o")

add_custom_target(${DALIB}.a
COMMAND rm -rf junk${DAEXE} && mkdir -p junk${DAEXE}
COMMAND cd junk${DAEXE} ${ZIP}
COMMAND cd ../ && rm -rf junk${DAEXE}
DEPENDS ${ALIROOTALIBS} 
WORKING_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY} 
)

add_custom_command(
TARGET ${DATARGETNAME}-clean
COMMAND echo "***** Cleaning ${DAMODULE} ${SUBDAMODULE} ${DANAME} detector-algorithm *****"
COMMAND rm -f ${DASPECFILE}
COMMAND rm -f ${DATAR}
COMMAND rm -f ${DAEXE}
COMMAND rm -f ${DAOBJ}
COMMAND rm -f ${DADEP}
COMMAND rm -f ${DAMAKEFILE}
COMMAND rm -f ${DALIB}.a
COMMAND rm -f ${ALIROOTALIBS}
COMMAND rm -f ${EXTRAROOTLIB}
)

separate_arguments(SYSLIBS)
#execute_process(COMMAND ${AMORE}/amore-config --ldflags-da-static | tr "\n" " "  OUTPUT_VARIABLE _amore1 OUTPUT_STRIP_TRAILING_WHITESPACE)
#execute_process(COMMAND date-config --rcproxylibs OUTPUT_VARIABLE _amore2 OUTPUT_STRIP_TRAILING_WHITESPACE)
#execute_process(COMMAND date-config --monitorlibs=noshift OUTPUT_VARIABLE _monitor1 OUTPUT_STRIP_TRAILING_WHITESPACE)
  
#set(AMOREDALIBS "-static ${_amore1} ${_amore2}\\")
#set(MONITORLIBS "${_monitor1}\\")

add_custom_target(DAMAKEFILE_${DAEXE}_)
add_custom_command( 
TARGET DAMAKEFILE_${DAEXE}_
COMMAND rm -f ${DAMAKEFILE}
COMMAND @echo "***** Making makefile ${DAMAKEFILE} *****"
COMMAND @echo '\#***************************************************' >> ${DAMAKEFILE}
COMMAND @echo '\# Makefile for Detector Algorithm' >> ${DAMAKEFILE}
COMMAND @echo '\#' >> ${DAMAKEFILE}
COMMAND @echo '\# It is necessary to setup build environment before' >> ${DAMAKEFILE}
COMMAND @echo '\# using make:' >> ${DAMAKEFILE}
COMMAND @echo '\# - define path to daqDAlib (env. DAQDALIB_PATH)' >> ${DAMAKEFILE}
COMMAND @echo '\#' >> ${DAMAKEFILE}
COMMAND @echo '\#*****************************************************' >> ${DAMAKEFILE}
COMMAND @echo "" >> ${DAMAKEFILE}
COMMAND @echo "CXX=g++" >> ${DAMAKEFILE}
COMMAND @echo "LD=g++" >> ${DAMAKEFILE}
COMMAND @echo "" >> ${DAMAKEFILE}
COMMAND @echo "CXXFLAGS=${CXXFLAGS}" >> ${DAMAKEFILE}
COMMAND @echo "LDFLAGS=${LDFLAGS}" >> ${DAMAKEFILE}
COMMAND @echo "" >> ${DAMAKEFILE}
COMMAND @echo 'ifeq ($$(DAQDALIB_PATH),)' >> ${DAMAKEFILE}
COMMAND @echo "DAQDADIR=${ALICE}/daqDAlib" >> ${DAMAKEFILE}
COMMAND @echo "else" >> ${DAMAKEFILE}
COMMAND @echo 'DAQDADIR=$$(DAQDALIB_PATH)' >> ${DAMAKEFILE}
COMMAND @echo "endif" >> ${DAMAKEFILE}
COMMAND @echo "" >> ${DAMAKEFILE}
COMMAND @echo 'DAQDALIB=$$(DAQDADIR)/libdaqDA.a' >> ${DAMAKEFILE}
#COMMAND @echo 'AMOREDALIBS=-static $$(shell $$(AMORE)/amore-config --ldflags-da-static) $$(shell date-config --rcproxylibs)' >> ${DAMAKEFILE}
#COMMAND @echo 'MONITORLIBS=$$(shell date-config --monitorlibs=noshift)' >> ${DAMAKEFILE}
#COMMAND @echo 'AMOREDALIBS=-static ${_amore1} ${_amore2}' >> ${DAMAKEFILE}
#COMMAND @echo 'MONITORLIBS=${monitor1}' >> ${DAMAKEFILE}
COMMAND @echo "" >> ${DAMAKEFILE}
COMMAND @echo "${DAMODULE}${SUBDAMODULE}${DANAME}da.exe: ${DAMODULE}${SUBDAMODULE}${DANAME}da.o" >> ${DAMAKEFILE}
COMMAND @echo -e '\t$$(LD) $$(LDFLAGS) -o ''$$''@ ''$$''< ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/lib${DALIB}.a ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${EXTRAROOTLIB} $$(ROOTSYS)/lib/libRoot.a $$(ROOTSYS)/lib/libfreetype.a $$(ROOTSYS)/lib/libpcre.a ${SYSLIBS} $$(DAQDALIB) ${AMOREDALIBS} ${MONITORLIBS}' >> ${DAMAKEFILE}
COMMAND @echo "" >> ${DAMAKEFILE}
COMMAND @echo "${DAMODULE}${SUBDAMODULE}${DANAME}da.o: ${DAMODULE}${SUBDAMODULE}${DANAME}da.cxx" >> ${DAMAKEFILE}
COMMAND @echo -e '\t$$(CXX) -c $$(CXXFLAGS) -I$$(DAQDADIR) ${mod} ''$$''< -o ''$$''@' >> ${DAMAKEFILE}
COMMAND @echo "" >> ${DAMAKEFILE}
COMMAND @echo "clean:" >> ${DAMAKEFILE}
COMMAND @echo -e '\t@rm -f ${DAMODULE}${SUBDAMODULE}${DANAME}da.exe ${DAMODULE}${SUBDAMODULE}${DANAME}da.o' >> ${DAMAKEFILE}
COMMAND @echo "" >> ${DAMAKEFILE}
COMMAND @echo "install: ${DAMODULE}${SUBDAMODULE}${DANAME}da.exe" >> ${DAMAKEFILE}
COMMAND @echo -e '\tif [ $$(INSTALL_PATH) == '' ]; then \\' >> ${DAMAKEFILE}
COMMAND @echo -e '\techo Environment variable INSTALL_PATH undefined, cannot continue\; \\' >> ${DAMAKEFILE}
COMMAND @echo -e '\texit 1\; \\' >> ${DAMAKEFILE}
COMMAND @echo -e '\tfi' >> ${DAMAKEFILE}
COMMAND @echo -e '\t@cp -p ${DAMODULE}${SUBDAMODULE}${DANAME}da.exe $$(INSTALL_PATH)' >> ${DAMAKEFILE}
)

add_custom_target(DATAR_${DAEXE}_)
add_custom_command( 
TARGET DATAR_${DAEXE}_
COMMAND @echo "***** Making archive ${DATAR} *****"
COMMAND rm -rf ${DATAR}
COMMAND rm -rf junk
COMMAND mkdir junk && mkdir junk/${DAARC} 
COMMAND cp -a ${ALICE_ROOT}/${DASRC} junk/${DAARC} 
COMMAND cp -a ${DAMAKEFILE} junk/${DAARC}/Makefile 
COMMAND cp -a ${DASPECFILE} junk/${DAARC}/${DAMODULE}${SUBDAMODULE}${DANAME}da.spec
COMMAND cd junk && tar czf ${DATAR} * 
COMMAND cd .. && rm -rf junk
)
add_dependencies(DATAR_${DAEXE}_ DAMAKEFILE_${DAEXE}_ DASPECFILE_${DAEXE}_ ${DASRC}) 

add_custom_target(DASPECFILE_${DAEXE}_)

add_custom_command(
TARGET DASPECFILE_${DAEXE}_
COMMAND rm -f ${DASPECFILE}
COMMAND @echo "***** Making RPM spec-file ${DASPECFILE} *****"
COMMAND @echo '\# RPM specfile for $(DAMODULE)${SUBDAMODULE}$(DANAME) Detector Algorithm' >> ${DASPECFILE}
COMMAND @echo "Summary: ${ONLINEDETECTORNAME} Detector Algorithm" >> ${DASPECFILE}
COMMAND @echo "Name: ${DAARCNAME}" >> ${DASPECFILE}
COMMAND @echo "Version: ${DAVERSION}" >> ${DASPECFILE}
COMMAND @echo "Release: ${DAALIROOTRELEASE}" >> ${DASPECFILE}
COMMAND @echo "License: CERN Alice DAQ/Offine" >> ${DASPECFILE}
COMMAND @echo "Source: %{name}-%{version}.src.tar.gz" >> ${DASPECFILE}
COMMAND @echo "Group: Applications/Alice" >> ${DASPECFILE}
COMMAND @echo "Prefix: /opt/%{name}" >> ${DASPECFILE}
COMMAND @echo "BuildRoot: %{_tmppath}/%{name}-root" >> ${DASPECFILE}
COMMAND @echo "" >> ${DASPECFILE}
COMMAND @echo '\# automatic dependencies' >> ${DASPECFILE}
COMMAND @echo "AutoReqProv: yes" >> ${DASPECFILE}
COMMAND @echo "" >> ${DASPECFILE}
COMMAND @echo '\# list here required RPM packages for compilation' >> ${DASPECFILE}
COMMAND @echo "BuildRequires: date" >> ${DASPECFILE} 
COMMAND @echo "" >> ${DASPECFILE}
COMMAND @echo '\# list here required RPM packages for runtime' >> ${DASPECFILE}
COMMAND @echo "Requires: date" >> ${DASPECFILE}
COMMAND @echo "" >> ${DASPECFILE}    
COMMAND @echo '\# You can specify other dependencies in the description tag below.' >> ${DASPECFILE}
COMMAND @echo "" >> ${DASPECFILE}
COMMAND @echo '\# description of the package' >> ${DASPECFILE}
COMMAND @echo "%description" >> ${DASPECFILE}
COMMAND @echo "This is the ${ONLINEDETECTORNAME} ${DANAME} DA for online calibration." >> ${DASPECFILE}
COMMAND @echo "It uses data from ${DAMODULE} detectors at run time." >> ${DASPECFILE}
COMMAND @echo "Build requires: daqDAlib, date, AliRoot ${DAALIROOTRELEASE},ROOT ${DAROOTRELEASE}." >> ${DASPECFILE}
COMMAND @echo "Runtime requires: date." >> ${DASPECFILE}
COMMAND @echo "Contact: ${DACONTACT}" >> ${DASPECFILE}
COMMAND @echo "Link: ${DALINKPAGE}" >> ${DASPECFILE}
COMMAND @echo "Reference Run: ${DAREFRUN}" >> ${DASPECFILE}
COMMAND @echo "Run type: ${DARUNTYPE}" >> ${DASPECFILE}
COMMAND @echo "DA type: ${DATYPE}" >> ${DASPECFILE}
COMMAND @echo "Number of events needed: ${DANUMBEROFEVENTS}" >> ${DASPECFILE}
COMMAND @echo "Input files: ${DAINPUTFILES}" >> ${DASPECFILE}
COMMAND @echo "Output files: ${DAOUTPUTFILES}" >> ${DASPECFILE}
COMMAND @echo "Trigger types used: ${DATRIGGERTYPE}" >> ${DASPECFILE}
COMMAND @echo "" >> ${DASPECFILE}
COMMAND @echo '\#*****************************************************************' >> ${DASPECFILE}
COMMAND @echo '\# Do not modify following scripts' >> ${DASPECFILE}
COMMAND @echo '\#*****************************************************************' >> ${DASPECFILE}
COMMAND @echo "" >> ${DASPECFILE}
COMMAND @echo "%define debug_package %{nil}" >> ${DASPECFILE}
COMMAND @echo "%define __os_install_post %{nil}" >> ${DASPECFILE}
COMMAND @echo "" >> ${DASPECFILE}
COMMAND @echo '\# here is defined the installation root directory' >> ${DASPECFILE}
COMMAND @echo "%define pkgname %{name}-%{version}" >> ${DASPECFILE}
COMMAND @echo "%define destdir %{prefix}" >> ${DASPECFILE}
COMMAND @echo "" >> ${DASPECFILE}
COMMAND @echo '\# pre-compilation script: extract from tarball' >> ${DASPECFILE}
COMMAND @echo "%prep" >> ${DASPECFILE}
COMMAND @echo '\# extract archive' >> ${DASPECFILE}
COMMAND @echo "%setup -n %{pkgname}" >> ${DASPECFILE}
COMMAND @echo "" >> ${DASPECFILE}
COMMAND @echo '\# compile' >> ${DASPECFILE}
COMMAND @echo "%build" >> ${DASPECFILE}
COMMAND @echo "export DATE_SITE=" >> ${DASPECFILE}
COMMAND @echo ". /date/setup.sh" >> ${DASPECFILE}
COMMAND @echo "gmake" >> ${DASPECFILE}
COMMAND @echo "" >> ${DASPECFILE}
COMMAND @echo '\# install runtime files' >> ${DASPECFILE}
COMMAND @echo "%install" >> ${DASPECFILE}
COMMAND @echo '\#remove install dir if existing' >> ${DASPECFILE}
COMMAND echo  '[ -d ''$$''RPM_BUILD_ROOT ] && rm -rf ''$$''RPM_BUILD_ROOT' >> ${DASPECFILE}
COMMAND @echo '\#make install in install root directory' >> ${DASPECFILE}
COMMAND @echo 'mkdir -p ''$$''RPM_BUILD_ROOT%{prefix}' >> ${DASPECFILE}
COMMAND @echo 'make install INSTALL_PATH=''$$''RPM_BUILD_ROOT%{prefix}' >> ${DASPECFILE}
COMMAND @echo "" >> ${DASPECFILE}
COMMAND @echo '\# clean RPM build files' >> ${DASPECFILE}
COMMAND @echo "%clean" >> ${DASPECFILE}
COMMAND @echo '\# remove installed files' >> ${DASPECFILE}
COMMAND @echo 'rm -rf ''$$''RPM_BUILD_ROOT' >> ${DASPECFILE}
COMMAND @echo '\# remove source files' >> ${DASPECFILE}
COMMAND @echo 'rm -rf ''$$''RPM_BUILD_DIR/%{pkgname}' >> ${DASPECFILE}
COMMAND @echo "" >> ${DASPECFILE}
COMMAND @echo '\# pre-install script' >> ${DASPECFILE}
COMMAND @echo "%pre" >> ${DASPECFILE}
COMMAND @echo "" >> ${DASPECFILE}
COMMAND @echo '\# list of files to be installed' >> ${DASPECFILE}
COMMAND @echo "%files" >> ${DASPECFILE}
COMMAND @echo '%defattr (-,root,root)' >> ${DASPECFILE}
COMMAND @echo "%{destdir}" >> ${DASPECFILE}
COMMAND @echo "" >> ${DASPECFILE}
COMMAND @echo '\# post-install script' >> ${DASPECFILE}
COMMAND @echo '\# launched after rpm installed' >> ${DASPECFILE}
COMMAND @echo "%post" >> ${DASPECFILE}
COMMAND @echo "" >> ${DASPECFILE}
COMMAND @echo '\# post-uninstall script' >> ${DASPECFILE}
COMMAND @echo '\# launched after rpm removed' >> ${DASPECFILE}
COMMAND @echo "%postun" >> ${DASPECFILE}
)


add_custom_target( ${DATARGETNAME}-rpm
)
add_dependencies(${DATARGETNAME}-rpm DATAR_${DAEXE}_ DASPECFILE_${DAEXE}_ ${LIBPATH} ${DALIB}.a DADEP_${DAEXE}_ ${ROOTLIB}     
)
add_custom_command(TARGET ${DATARGETNAME}-rpm
COMMAND mkdir -p ${CMAKE_BINARY_DIR}/junk/SOURCES ${CMAKE_BINARY_DIR}/junk/SPECS ${CMAKE_BINARY_DIR}/junk/BUILD ${CMAKE_BINARY_DIR}/junk/RPMS ${CMAKE_BINARY_DIR}/junk/SRPMS
COMMAND cp ${DATAR} ${CMAKE_BINARY_DIR}/junk/SOURCES
COMMAND rpmbuild --verbose --define "_topdir ${CMAKE_BINARY_DIR}/junk" --nodeps -bb ${DASPECFILE}
COMMAND cp -p `find ${CMAKE_BINARY_DIR}/junk/ -name "${DAARC}-*.rpm"` . ;
COMMAND rm -rf junk
COMMAND echo "***** RPMS created and put ${CMAKE_BINARY_DIR} folder *****"
WORKING_DIRECTORY ${CMAKE_BINARY_DIR}  
)


add_custom_target(${DAEXE}
COMMAND echo "***** Making executable ${DAEXE} *****"
COMMAND g++ ${LDFLAGS} ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/lib${DALIB}.a ${EXTRAROOTLIB} ${ROOTSYS}/lib/libRoot.a ${ROOTSYS}/lib/libfreetype.a ${ROOTSYS}/lib/libpcre.a ${SYSLIBS} ${DAQDALIB} ${AMOREDALIBS} ${MONITORLIBS} -o ${DAEXE}
WORKING_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
)


#target_link_libraries(${DAEXE} "-L" "lib${DALIB}.a" ${DAOBJ} ${EXTRAROOTLIB} "${ROOTALIBDIR}/libRoot.a" "${ROOTALIBDIR}/libfreetype.a" "${ROOTALIBDIR}/libpcre.a" ${SYSLIBS} ${DAQDALIB} ${MONITORLIBS} ${AMOREDALIBS})

add_dependencies(${DAEXE} ${DASRC} DAOBJ_${DAEXE}_ ${BINPATH} ${LIBPATH} ${DALIB}.a ${DAQDALIB} ${ROOTLIB})

add_custom_target(DAOBJ_${DAEXE}_
)
add_custom_command(
TARGET DAOBJ_${DAEXE}_
COMMAND echo "***** Compiling ${DASRC} *****"
COMMAND g++ -c -DLinux -DDATE_SYS=Linux -Dlong32="int" -Dlong64="long long" -DdatePointer="long" -I/date/rorc -I/date/runControl -I/date/readList -I/date/eventBuilder -I/date/banksManager -I/date/bufferManager -I/date/db -I/date/commonDefs -I/date/monitoring -I/date/infoLogger -I/date/logbook -I${DAQDADIR} -I${ALICE_ROOT}/RAW -I${CMAKE_BINARY_DIR}/include -I$ENV{ROOTSYS}/include ${mod} ${date_head} ${ALICE_ROOT}/${DASRC} -o ${DAOBJ}
WORKING_DIRECTORY ${ALICE_ROOT} 
)
add_dependencies(DAOBJ_${DAEXE}_ DADEP_${DAEXE}_ )



add_custom_target(DADEP_${DAEXE}_)
add_custom_command(
TARGET DADEP_${DAEXE}_
COMMAND echo "***** Making detector-algorithm dependencies ${DADEP} *****"
COMMAND g++ -MM -DLinux -DDATE_SYS=Linux -Dlong32="int" -Dlong64="long long" -DdatePointer="long" -I/date/rorc -I/date/runControl -I/date/readList -I/date/eventBuilder -I/date/banksManager -I/date/bufferManager -I/date/db -I/date/commonDefs -I/date/monitoring -I/date/infoLogger -I/date/logbook -I${DAQDADIR} -I${ALICE_ROOT}/RAW -I${CMAKE_BINARY_DIR}/include -I$ENV{ROOTSYS}/include ${mod} ${date_head} ${ALICE_ROOT}/${DASRC} > ${DADEP}
WORKING_DIRECTORY ${ALICE_ROOT}
)


add_custom_command(TARGET clean
COMMAND rm -rf junk*.exe
WORKING_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
  

	endif(match)
  endforeach(dafile)
endforeach(detector)
