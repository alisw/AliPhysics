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
    message(STATUS "Setting IRST_INSTALLDIR to ${ALICE}/local/IRST")
    set(IRST_INSTALLDIR ${ALICE}/local/ALICENewRuleChecker)
  endif(ALICE)
endif(NOT IRST_INSTALLDIR)

if(IRST_INSTALLDIR)
  find_file(RULECHECKER_JAR NAMES NewRuleChecker.jar PATHS ${IRST_INSTALLDIR})
  find_file(RULECHECKER_RULES NAMES AliceCodingConventions.xml PATHS ${IRST_INSTALLDIR})
  if(RULECHECKER_JAR AND RULECHECKER_RULES AND RULECHECKER_SRCML AND JAVA_RUNTIME)
    set(RULECHECKER_FOUND TRUE)
    message(STATUS "RuleChecker found on the system")
  else()
    message(STATUS "RuleChecker not found on this system")
  endif(RULECHECKER_JAR AND RULECHECKER_RULES AND RULECHECKER_SRCML AND JAVA_RUNTIME)
else()
  message(STATUS "RuleChecker not found on this system")
endif(IRST_INSTALLDIR)

