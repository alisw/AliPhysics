# -*- mode: cmake -*-

# AliRoot Build System Module to find and configure IRST ALICE Coding Coventions RuleChecker
#
# Author: Johny Jose (johny.jose@cern.ch)
#         Port of previous Makefile build to cmake

cmake_minimum_required(VERSION 2.8.4 FATAL_ERROR)

find_program(RULECHECKER_SRCML NAMES src2srcml)
message(STATUS "Check for src2srcml: ${RULECHECKER_SRCML}")
find_package(Java)
message(STATUS "Check for Java: ${JAVA_RUNTIME}")

find_file(RULECHECKER_JAR   NAMES NewRuleChecker.jar    PATHS  ${CMAKE_SOURCE_DIR}/RuleChecker NO_DEFAULT_PATH)
find_file(RULECHECKER_RULES NAMES CodingConventions.xml PATHS  ${CMAKE_SOURCE_DIR}/RuleChecker NO_DEFAULT_PATH)
find_file(FACTEXTRACTOR_JAR NAMES FactExtractor.jar     PATHS  ${CMAKE_SOURCE_DIR}/RuleChecker NO_DEFAULT_PATH)
find_file(SMELLDETECTOR_JAR NAMES SmellDetector.jar     PATHS  ${CMAKE_SOURCE_DIR}/RuleChecker NO_DEFAULT_PATH)
if(RULECHECKER_JAR AND RULECHECKER_RULES AND RULECHECKER_SRCML AND JAVA_RUNTIME AND SMELLDETECTOR_JAR)
  set(RULECHECKER_FOUND TRUE)
  message(STATUS "RuleChecker jar :           ${RULECHECKER_JAR}")
  message(STATUS "RuleChecker rules :         ${RULECHECKER_RULES}")
  message(STATUS "RuleChecker factextractor : ${FACTEXTRACTOR_JAR}")
  message(STATUS "RuleChecker smelldetector : ${SMELLDETECTOR_JAR}")
  message(STATUS "RuleChecker found on the system")
  
  if(NOT EXISTS ${CMAKE_BINARY_DIR}/check-hxml-touchfile)
    file(WRITE ${CMAKE_BINARY_DIR}/check-hxml-touchfile "Dummy dependency for factfile")
  endif(NOT EXISTS ${CMAKE_BINARY_DIR}/check-hxml-touchfile)
  set(FACTFILE ${CMAKE_BINARY_DIR}/factFile.xml)
  set(_factfile_deps)
  
  file(GLOB_RECURSE _root_headers  ${ROOTSYS}/include/*.h)
  foreach(_root_header ${_root_headers})
    if(NOT _root_header MATCHES ".*G__ci.h")
      string (REGEX REPLACE "${ROOTSYS}/include/" "" _rel_root_header ${_root_header})
      string (REGEX REPLACE "h$" "h.xml" _rel_root_hxml ${_rel_root_header})
      get_filename_component(_rel_root_header_path ${_rel_root_hxml} PATH)
      set(_root_hxml ${CMAKE_BINARY_DIR}/roothxml/${_rel_root_hxml})
      if(NOT EXISTS ${CMAKE_BINARY_DIR}/roothxml/${_rel_root_header_path})
	file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/roothxml/${_rel_root_header_path})
      endif(NOT EXISTS ${CMAKE_BINARY_DIR}/roothxml/${_rel_root_header_path})
      list(APPEND _factfile_deps ${_root_hxml})
      add_custom_command(OUTPUT  ${_root_hxml}
        COMMAND ${RULECHECKER_SRCML} ${_root_header} ${_root_hxml}
        DEPENDS ${_root_header})
    endif(NOT _root_header MATCHES ".*G__ci.h")
  endforeach(_root_header ${_root_headers})
else()
  message(STATUS "RuleChecker not found on this system")
endif(RULECHECKER_JAR AND RULECHECKER_RULES AND RULECHECKER_SRCML AND JAVA_RUNTIME AND SMELLDETECTOR_JAR)

macro(ALICE_CheckModule)
  if(RULECHECKER_FOUND)
    set(CHECKDIR ${CMAKE_BINARY_DIR}/${MODULE}/check)

    file(GLOB_RECURSE _headers  ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
    file(GLOB_RECURSE _sources_tmp  ${CMAKE_CURRENT_SOURCE_DIR}/*.cxx)
    list(APPEND _sources_tmp ${_headers})
    foreach(_srcfile ${_sources_tmp})
      string(REPLACE ".h"   "" _srcfile_tmp ${_srcfile})
      string(REPLACE ".cxx" "" _srcfile ${_srcfile_tmp})
      list(APPEND _sources ${_srcfile})
    endforeach(_srcfile ${_sources_tmp})
    list(REMOVE_DUPLICATES _sources)
    list(SORT _sources)

    set(_violfiles)
    set(_smellfiles)
    set(_module_factfile_deps)
    foreach(_srcfile ${_sources})
      if(NOT _srcfile MATCHES "^.*LinkDef$" AND NOT _srcfile MATCHES ".*PYTHIA8/pythia8.*")
	set(_violdeps)
	string(REPLACE "${CMAKE_CURRENT_SOURCE_DIR}/" "" _srcfile_short ${_srcfile})
	set(_viol ${CHECKDIR}/${_srcfile_short}.viol)
	set(_smell ${CHECKDIR}/${_srcfile_short}.smell)
	get_filename_component(_viol_path ${_viol} PATH)
	list(APPEND _violfiles ${_viol})
	list(APPEND _smellfiles ${_smell})
	if(EXISTS ${_srcfile}.h)
	  add_custom_command(OUTPUT ${CHECKDIR}/${_srcfile_short}.h.xml
	                     COMMAND ${CMAKE_COMMAND} -E make_directory ${_viol_path}
                             COMMAND ${RULECHECKER_SRCML} ${_srcfile}.h ${CHECKDIR}/${_srcfile_short}.h.xml
			     COMMAND ${CMAKE_COMMAND} -E touch ${CMAKE_BINARY_DIR}/check-hxml-touchfile
			     DEPENDS ${_srcfile}.h)
	  list(APPEND _violdeps ${CHECKDIR}/${_srcfile_short}.h.xml)
	  list(APPEND _module_factfile_deps ${CHECKDIR}/${_srcfile_short}.h.xml)
	endif(EXISTS ${_srcfile}.h)
	if(EXISTS ${_srcfile}.cxx)
	  add_custom_command(OUTPUT ${CHECKDIR}/${_srcfile_short}.cxx.xml
	                     COMMAND ${CMAKE_COMMAND} -E make_directory ${_viol_path}
                             COMMAND ${RULECHECKER_SRCML} ${_srcfile}.cxx ${CHECKDIR}/${_srcfile_short}.cxx.xml
			     DEPENDS ${_srcfile}.cxx)
	  list(APPEND _violdeps ${CHECKDIR}/${_srcfile_short}.cxx.xml)
	endif(EXISTS ${_srcfile}.cxx)
	add_custom_command( OUTPUT ${_viol}
                            COMMAND ${JAVA_RUNTIME} -Xmx1024M -jar ${RULECHECKER_JAR} ${CHECKDIR}/${_srcfile_short}.cxx.xml ${CHECKDIR}/${_srcfile_short}.h.xml ${FACTFILE} ${RULECHECKER_RULES} > ${_viol}
                            DEPENDS ${_violdeps}
                            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
	add_custom_command( OUTPUT ${_smell}
                            COMMAND ${JAVA_RUNTIME} -Xmx1024M -jar ${SMELLDETECTOR_JAR} ${CHECKDIR}/${_srcfile_short}.cxx.xml ${CHECKDIR}/${_srcfile_short}.h.xml > ${_smell}
                            DEPENDS ${_violdeps}
                            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
      endif(NOT _srcfile MATCHES "^.*LinkDef$" AND NOT _srcfile MATCHES ".*PYTHIA8/pythia8.*")
    endforeach(_srcfile ${_sources})

    if(_violfiles)
      add_custom_target(${MODULE}-check DEPENDS ${_violfiles})
      add_dependencies(${MODULE}-check factfile)
      add_dependencies(check-all ${MODULE}-check)

      add_custom_target(${MODULE}-smell DEPENDS ${_smellfiles})
      add_dependencies(smell-all ${MODULE}-smell)

      if(_module_factfile_deps)
	add_custom_target(${MODULE}-hxml DEPENDS ${_module_factfile_deps})
	add_dependencies(check-hxml ${MODULE}-hxml)
      endif(_module_factfile_deps)
    endif(_violfiles)


  endif(RULECHECKER_FOUND)
endmacro(ALICE_CheckModule)

