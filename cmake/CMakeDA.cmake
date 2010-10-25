cmake_minimum_required(VERSION 2.8.2 FATAL_ERROR)

# -----------Utilities--------------------

#list of detectors 
set(ONLINEDETECTORS SPD SDD SSD TPC TRD TOF HMP PHS CPV PMD MCH MTR FMD T00 V00 ZDC ACO TRI EMC HLT TST GRP)

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
  endif(_match)
  set(${_info} ${_data} PARENT_SCOPE)
endfunction()

# ----------Common stuff-------------------

file(GLOB_RECURSE _dafiles $ENV{ALICE_ROOT}/*da.cxx)
set(DAINSTALL "$ENV{ALICE_INSTALL}/DA")

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

separate_arguments(MONITORLIBS)

set(SYSLIBS -ldl -lpthread ${_xml2libs})


set(EXTRAROOTLIB "libRootExtra.a")

file(GLOB _extraroot "$ENV{ROOTSYS}/montercarlo/vmc/src/*.o" "$ENV{ROOTSYS}/tree/treeplayer/src/*.o" "$ENV{ROOTSYS}/io/xmlparser/src/*.o" "$ENV{ROOTSYS}/math/minuit2/src/*.o")

add_library(${EXTRAROOTLIB} STATIC ${_extraroot})	
set_target_properties(${EXTRAROOTLIB} PROPERTIES LINKER_LANGUAGE CXX)

set(DAQDALIB_PATH $ENV{DAQDALIB_PATH})	
if(DAQDALIB_PATH)
  set(DAQDADIR "${DAQDALIB_PATH}")
else()
  set(DAQDADIR "$ENV{ALICE}/daqDAlib")
endif(DAQDALIB_PATH)
set(DAQDALIB "${DAQDADIR}/libdaqDA.a")

include_directories(${DAQDADIR} RAW include STEER)
include_directories(SYSTEM ${ROOTINCDIR})

# ----------Create All Valid targets---------
	  
foreach(detector ${ONLINEDETECTORS} )

  set(ONLINEDETECTORNAME ${detector})

  detector_module(DAMODULE ${ONLINEDETECTORNAME})

  detector_subda(SUBDAMODULE ${ONLINEDETECTORNAME})
  
  #ALIROOTALIBS

  set(ALIROOTALIBS)

  list(APPEND ALIROOTALIBS RAWDatabase-static RAWDatarec-static RAWDatasim-static STEERBase-static STEER-static CDB-static ESD-static STAT-static AOD-static )

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

	  set(DATARGETDIR "${DAINSTALL}/${DAMODULE}/tgt_$ENV{ALICE_TARGET}")
	  file(MAKE_DIRECTORY ${DATARGETDIR})

	  set(DASRC "${DAMODULE}/${DAMODULE}${SUBDAMODULE}${DANAME}da.cxx")
	  set(DALIB "${DAMODULE}${SUBDAMODULE}${DANAME}DA")
	  set(DAEXE "${DAMODULE}${SUBDAMODULE}${DANAME}da.exe")
	  
	# DAVERSION
	  execute_process(COMMAND svn info $ENV{ALICE_ROOT}/${DASRC} OUTPUT_VARIABLE _daversion OUTPUT_STRIP_TRAILING_WHITESPACE)
	  string(REGEX REPLACE ".*Last Changed Rev: ([^\n]+)\n.*" "\\1" DAVERSION ${_daversion})

	# DAROOTRELEASE 
	  execute_process(COMMAND root-config --version OUTPUT_VARIABLE _darootrelease OUTPUT_STRIP_TRAILING_WHITESPACE)
	  string(REGEX REPLACE "/" "." DAROOTRELEASE ${_darootrelease})
	
	# DAALIROOTRELEASE
	  string(REGEX REPLACE ".*URL: .*/(.+)/${DASRC}.*$" "\\1" DAALIROOTRELEASE ${_daversion})
	
	  set(DAARCNAME "${DATARGETNAME}")
	  set(DAARC "${DAARCNAME}-${DAVERSION}")
	  set(DATAR "${DATARGETDIR}/${DAARC}.src.tar.gz")
	  set(DASPECFILE "${DATARGETDIR}/${DAMODULE}${SUBDAMODULE}${DANAME}da.spec")
	  set(DAMAKEFILE "${DATARGETDIR}/${DAMODULE}${SUBDAMODULE}${DANAME}da.make")


          
	  if(EXTRADAMODULE)
	    ##**set
	  endif(EXTRADAMODULE)

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

	  set(ZIP)
	  foreach(_lib ${ALIROOTALIBS})
	   string(REGEX REPLACE "-static" "" _lib ${_lib})
	   list(APPEND ZIP && ar x "../lib${_lib}.a")
	  endforeach(_lib)
 	  list (APPEND ZIP && ar r "../lib${DALIB}.a" "*.o")

	  add_custom_target( ${DALIB} COMMAND rm -rf junk${DAEXE} && mkdir -p junk${DAEXE} 
				COMMAND cd junk${DAEXE} ${ZIP}
				COMMAND cd ../ && rm -rf junk${DAEXE}
				DEPENDS ${ALIROOTALIBS}
				WORKING_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
	  add_custom_command( TARGET clean
	                     COMMAND rm -rf junk${DAEXE}
			     WORKING_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
	    

	  add_custom_target(${DATARGETNAME})
	  add_executable(${DAEXE} ${DASRC})
    	  set_property(TARGET ${DAEXE} PROPERTY EXCLUDE_FROM_ALL TRUE)
	  add_dependencies(${DAEXE} ${DALIB})
	  add_dependencies(DA-all ${DATARGETNAME})
	  target_link_libraries(${DAEXE} "-L${CMAKE_LIBRARY_OUTPUT_DIRECTORY}" "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/lib${DALIB}.a" ${EXTRAROOTLIB} "${ROOTLIBDIR}/libRoot.a" "${ROOTLIBDIR}/libfreetype.a" "${ROOTLIBDIR}/libpcre.a" ${SYSLIBS} ${DAQDALIB} ${MONITORLIBS} ${AMOREDALIBS})
	  add_dependencies(${DATARGETNAME} ${DAEXE})
	  

	endif(match)
  endforeach(dafile)
endforeach(detector)
