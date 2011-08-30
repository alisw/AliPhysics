# AliRoot Build System Module to find and configure DATE
#
# Author: Johny Jose (johny.jose@cern.ch)
#         Port of previous Makefile build to cmake

cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

find_program(DATE_CONFIG date-config)
message(STATUS "Check for Date : ${DATE_CONFIG}")
if(DATE_CONFIG)
  set(DATE_FOUND TRUE)
  set(DATE_ROOT $ENV{DATE_ROOT})
  set(DATE_COMMON_DEFS $ENV{DATE_COMMON_DEFS})
  set(DATE_MONITOR_DIR $ENV{DATE_MONITOR_DIR})
  execute_process(COMMAND date-config --cflags OUTPUT_VARIABLE DATEFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
  string(REGEX REPLACE "\"" "'" DATEFLAGS ${DATEFLAGS})
  set(DATEFLAGS "-DALI_DATE ${DATEFLAGS}")
  execute_process(COMMAND date-config --monitorlibs=dyn OUTPUT_VARIABLE DMONLIBS OUTPUT_STRIP_TRAILING_WHITESPACE)
else()
  set(DATEFLAGS "-D${CMAKE_SYSTEM_NAME} -DDATE_SYS=${CMAKE_SYSTEM_NAME} -Dlong32='int' -Dlong64='long long' -DdatePointer='long'")
endif(DATE_CONFIG)

