cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

# -----------Utilities--------------------

#list of detectors 
set(ONLINEDETECTORS SPD SDD SSD TPC TRD TOF HMP PHS CPV PMD MCH MTR FMD T00 V00 ZDC ACO TRI EMC HLT TST GRP)

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
set(DAINSTALL "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

file(MAKE_DIRECTORY ${DAINSTALL})

#**Set compiler flags


find_program(AMORE amore-config)
#temporary
find_program(DATE date-config)
find_program(ROOT_CONFIG root-config)
find_program(XML2 xml2-config)

if(AMORE)
  execute_process(COMMAND ${AMORE} --cflags --includes OUTPUT_VARIABLE AMOREFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
  #**Set compiler flags
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${AMOREFLAGS}")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${AMOREFLAGS}")
  
  execute_process(COMMAND ${AMORE} --ldflags-da-static OUTPUT_VARIABLE _lddaflags OUTPUT_STRIP_TRAILING_WHITESPACE)
  string(REGEX REPLACE "\n" " " _lddaflags ${_lddaflags})
  execute_process(COMMAND ${DATE} --rcproxylibs OUTPUT_VARIABLE _rcproxylib OUTPUT_STRIP_TRAILING_WHITESPACE)
  set(AMOREDALIBS "-static ${_lddaflags} ${_rcproxylib}") 	
else()
  set(AMORELIBS)
endif(AMORE)

execute_process(COMMAND ${ROOT_CONFIG} --libdir OUTPUT_VARIABLE ROOTLIBDIR OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND ${XML2} --libs OUTPUT_VARIABLE _xml2libs OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND ${DATE} --monitorlibs=noshift OUTPUT_VARIABLE MONITORLIBS OUTPUT_STRIP_TRAILING_WHITESPACE)


string(REGEX REPLACE " " ";" MONITORLIBS ${MONITORLIBS})
set(SYSLIBS -ldl -lpthread ${_xml2libs})
set(EXTRAROOTLIB "libRootExtra.a")
file(GLOB _extraroot "$ENV{ROOTSYS}/montercarlo/vmc/src/*.o" "$ENV{ROOTSYS}/tree/treeplayer/src/*.o" "$ENV{ROOTSYS}/io/xmlparser/src/*.o" "$ENV{ROOTSYS}/math/minuit2/src/*.o")
#message("${_extraroot}")


#add_custom_target(TARGET ${EXTRAROOTLIB} COMMAND ${CMAKE_AR} r $ENV{ALICE_INSTALL}/lib/tgt_$ENV{ALICE_TARGET}/${EXTRAROOTLIB} ${_extraroot} COMMAND pwd)

message("RAW SRCS ${RAWDatabase_SRC}")

# ----------Create All Valid targets---------
	  
foreach(detector ${ONLINEDETECTORS} )
  
  #ALIROOTALIBS
  set(ALIROOTALIBS)
  list(APPEND ALIROOTALIBS RAWDatabase_a RAWDatarec_a RAWDatasim_a STEERBase_a STEER_a CDB_a ESD_a STAT_a AOD_a )

  set(ONLINEDETECTORNAME ${detector})
  detector_module(DAMODULE ${ONLINEDETECTORNAME})
  detector_subda(SUBDAMODULE ${ONLINEDETECTORNAME})
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
          message("${DATARGETNAME}")

	  set(DATARGETDIR "${DAINSTALL}/${DAMODULE}/tgt_$ENV{ALICE_TARGET}")
	  file(MAKE_DIRECTORY ${DATARGETDIR})

	  set(DASRC "${DAMODULE}/${DAMODULE}${SUBDAMODULE}${DANAME}da.cxx")
	  set(DALIB "${DAMODULE}${SUBDAMODULE}${DANAME}DA")
	  set(DAEXE "${DAMODULE}${SUBDAMODULE}${DANAME}da")
	  
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
	  set(DAQDALIB_PATH $ENV{DAQDALIB_PATH})	
 	  if(DAQDALIB_PATH)
	    set(DAQDADIR "${DAQDALIB_PATH}")
	  else()
	    set(DAQDADIR "$ENV{ALICE}/daqDAlib")
	  endif(DAQDALIB_PATH)

	  set(DAQDALIB "${DAQDADIR}/libdaqDA.a")

#	  message(${DAVERSION})
#	  message(${DAALIROOTRELEASE})

	  ##**set(EXTRADAMODULE ALIROOTALIBS
#	  file(READ "$ENV{ALICE_ROOT}/${DAMODULE}/CMake_lib${DAMODULE}.txt" _modulesrc )
#	  message("${_modulesrc}")
#	  string(REGEX MATCHALL "[^ ]+\\.cxx" DAMODULE_SRC ${_modulesrc})
	  
#	  set(MODULEALIB ${DAMODULE}_SRC)
#	  
	  
#	  message("MODULE SRCS - ${DAMODULE_SRC}")
#	  string(REGEX MATCH "[^\n]*" 
#	  list(APPEND ALIROOTALIBS ${DAMODULE_SRC})

# Super Duper Hack :D
	  file(GLOB _damodule "$ENV{ALICE_ROOT}/${DAMODULE}/lib${DAMODULE}*.pkg" )
	  
	  message(${_damodule})
	  foreach(_submodule ${_damodule})
	    string(REGEX REPLACE ".*lib(${DAMODULE}.*)" "\\1" _submodule ${_submodule})
	    string(REGEX REPLACE "\\.pkg" "_a" _submodule ${_submodule})
	    list(APPEND ALIROOTALIBS "${_submodule}")	  
	    message("Adding ${_submodule}")
	  endforeach(_submodule)
	  
#	file(WRITE "$ENV{ALICE_INSTALL}/tmp" "list(APPEND ALIROOTALIBS ${DAMODULE}-all)")
#	  include("$ENV{ALICE_INSTALL}/tmp")

	  list(REMOVE_DUPLICATES ALIROOTALIBS)
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

#	  message("${_dasrc}")
#	  message("DACONTACT - ${DACONTACT}")
#	  message("DALINKPAGE - ${DALINKPAGE}")
#	  message("DAREFFUN - ${DAREFFUN}")
#	  message("DARUNTYPE - ${DARUNTYPE}")
#	  message("DATYPE - ${DATYPE}")
#	  message("DANUMBEROFEVENTS - ${DANUMBEROFEVENTS}")
#	  message("DAINPUTFILES - ${DAINPUTFILES}")
#	  message("DAOUTPUTFILES - ${DAOUTPUTFILES}")
#	  message("DATRIGGERTYPE - ${DATRIGGERTYPE}")

#----------- Targets ----------------------

  	  set(CMAKE_CXX_FLAGS "${CXXFLAGS}")
	  set(CMAKE_C_FLAGS "${CFLAGS}")
	  set(CMAKE_Fortran_FLAGS "${FFLAGS}")
	  set(CMAKE_SHARED_LINKER_FLAGS ${SOFLAGS}) 
	  set(CMAKE_MODULE_LINKER_FLAGS ${LDFLAGS})

	  include_directories(${MODULES} ${DAQDADIR} )
	  include_directories(SYSTEM ${ROOTINCDIR})
	  set(ZIP)
	  foreach(_lib ${ALIROOTALIBS})
	   string(REGEX REPLACE "-all" "_a" _lib ${_lib})
	   list(APPEND ZIP && ar x "../lib${_lib}.a")
	  endforeach(_lib)
 	  list (APPEND ZIP && ar r "../${DALIB}.a" "*.o")
 
	  add_custom_target("${DALIB}_x" rm -rf junk && mkdir -p junk 
				COMMAND cd junk ${ZIP}
				COMMAND cd ../ && rm -rf junk
				DEPENDS ${ALIROOTALIBS}
				WORKING_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

	  add_custom_target(${DATARGETNAME} DEPENDS ${DAEXE} )
	  add_executable(${DAEXE} ${DASRC} )
    	  set_property(TARGET ${DAEXE} PROPERTY EXCLUDE_FROM_ALL TRUE)
	  add_dependencies(${DAEXE} ${DALIB}_x)
	  target_link_libraries(${DAEXE} "-L${CMAKE_LIBRARY_OUTPUT_DIRECTORY}" "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${DALIB}.a" ${EXTRAROOTLIB} "${ROOTLIBDIR}/libRoot.a" "${ROOTLIBDIR}/libfreetype.a" "${ROOTLIBDIR}/libpcre.a" ${SYSLIBS} ${DAQDALIB} ${MONITORLIBS} ${AMOREDALIBS})
	  add_custom_command(TARGET ${DAEXE} 
		   	PRE_LINK
			COMMAND ${CMAKE_AR} r ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${EXTRAROOTLIB} ${_extraroot})
#	  message("${DALIB} --> ${ALIROOTALIBS}")
	endif(match)
  endforeach(dafile)
endforeach(detector)
