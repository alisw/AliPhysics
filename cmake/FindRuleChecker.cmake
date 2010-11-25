# AliRoot Build System Module to find and configure IRST ALICE Coding Coventions RuleChecker
#
# Author: Johny Jose (johny.jose@cern.ch)
#         Port of previous Makefile build to cmake

cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

find_program(RULECHECKER_SRCML NAMES src2srcml)
message(STATUS "Check for src2srcml: ${RULECHECKER_SRCML}")
find_package(Java)
message(STATUS "Check for Java: ${JAVA_RUNTIME}")

set(IRST_INSTALLDIR $ENV{IRST_INSTALLDIR})
if(NOT IRST_INSTALLDIR)
  if(ALICE)
    message(STATUS "Setting IRST_INSTALLDIR to ${ALICE}/local/ALICENewRuleChecker")
    set(IRST_INSTALLDIR ${ALICE_ROOT}/cmakelocal/ALICENewRuleChecker)
  endif(ALICE)
endif(NOT IRST_INSTALLDIR)

if(IRST_INSTALLDIR)
  find_file(RULECHECKER_JAR NAMES NewRuleChecker.jar PATHS ${IRST_INSTALLDIR}/NewRuleChecker ${ALICE_ROOT}/cmake/RuleChecker)
  find_file(RULECHECKER_RULES NAMES CodingConventions.xml PATHS ${IRST_INSTALLDIR}/NewRuleChecker/config ${ALICE_ROOT}/cmake/RuleChecker)
  find_file(FACTEXTRACTOR_JAR NAME FactExtractor.jar PATHS ${IRST_INSTALLDIR}/FactExtractor ${ALICE_ROOT}/cmake/RuleChecker)
  if(RULECHECKER_JAR AND RULECHECKER_RULES AND RULECHECKER_SRCML AND JAVA_RUNTIME)
    set(RULECHECKER_FOUND TRUE)
    message(STATUS "RuleChecker jar : ${RULECHECKER_JAR}")
    message(STATUS "RuleChecker rules : ${RULECHECKER_RULES}")
    message(STATUS "RuleChecker factextractor : ${FACTEXTRACTOR_JAR}")
    message(STATUS "RuleChecker found on the system")
  else()
    message(STATUS "RuleChecker not found on this system")
  endif(RULECHECKER_JAR AND RULECHECKER_RULES AND RULECHECKER_SRCML AND JAVA_RUNTIME)
else()
  message(STATUS "RuleChecker not found on this system")
endif(IRST_INSTALLDIR)

macro(ALICE_CheckModule)
  if(RULECHECKER_FOUND)
    set(CHECKDIR ${CMAKE_BINARY_DIR}/${MODULE}/check)
    set(violFiles)
    add_custom_command( OUTPUT ${MODULE}-${PACKAGE}-factFile
                        COMMAND ${CMAKE_COMMAND} -E make_directory ${CHECKDIR}
                        COMMAND ${CMAKE_COMMAND} -E make_directory ${CHECKDIR}/viols
                        COMMAND ${JAVA_RUNTIME} -jar ${FACTEXTRACTOR_JAR} ${CHECKDIR} ${CHECKDIR} 
                        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

    foreach(_srcfile ${SRCS})
      string (REGEX REPLACE "cxx$" "h" _header ${_srcfile})
      get_filename_component(_srcname ${_srcfile} NAME)
      string (REGEX REPLACE "cxx$" "viol" _viol ${_srcname})
      string (REGEX REPLACE "cxx$" "cxx.xml" _srcxml ${_srcname})
      string (REGEX REPLACE "cxx$" "h.xml" _hxml ${_srcname})
      string (REGEX REPLACE ".cxx$" "" _class ${_srcname})
      set(_depends ${_srcfile} ${MODULE}-${PACKAGE}-factFile)

      if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${_header})
        list(APPEND depends ${_header})
        add_custom_command( OUTPUT ${_viol}
                          COMMAND ${RULECHECKER_SRCML} ${_srcfile} ${CHECKDIR}/${_srcxml}
                          COMMAND ${RULECHECKER_SRCML} ${_header} ${CHECKDIR}/${_hxml}
                          COMMAND ${JAVA_RUNTIME} -jar ${RULECHECKER_JAR} ${CHECKDIR}/${_srcxml} ${CHECKDIR}/${_hxml} ${CHECKDIR}/factFile.xml ${RULECHECKER_RULES} > ${CHECKDIR}/viols/${_viol}
                          DEPENDS ${_depends}
                          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
        list(APPEND violFiles ${_viol})
      else()
        add_custom_command( OUTPUT ${_viol}
                          COMMAND ${RULECHECKER_SRCML} ${_srcfile} ${CHECKDIR}/${_srcxml}
                          COMMAND ${JAVA_RUNTIME} -jar ${RULECHECKER_JAR} ${CHECKDIR}/${_srcxml} ${CHECKDIR}/${_hxml} ${CHECKDIR}/factFile.xml ${RULECHECKER_RULES} > ${CHECKDIR}/viols/${_viol}
                          DEPENDS ${_depends}
                          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
        list(APPEND violFiles ${_viol})
      endif(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${_header})
      if(CLASSCHECK STREQUAL "YES")
        add_custom_target(${MODULE}-${_class}-check DEPENDS ${_viol})
      endif(CLASSCHECK STREQUAL "YES")
    endforeach(_srcfile ${SRCS})
    if(violFiles)
      add_custom_target(${PACKAGE}-check DEPENDS ${violFiles})
      add_dependencies(${PACKAGE}-check ${violFiles})
      add_dependencies(${MODULE}-check-all ${PACKAGE}-check)
    endif(violFiles)
    add_custom_command(TARGET clean
                       COMMAND ${CMAKE_COMMAND} -E remove_directory ${CHECKDIR})

  endif(RULECHECKER_FOUND)
endmacro(ALICE_CheckModule)

