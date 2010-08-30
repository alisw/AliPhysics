# -*- mode: cmake -*-

# - Find the DATE system
# Finds if the RuleChecker is installed and sets the following variables:
#
# DATE_FOUND = Boolean defining if DATE is installed
#

# Check if DATE is installed and in the path

Find_program( DATE_PROGRAM date-config )

If(DATE_PROGRAM)
  Set(DATE_FOUND YES)
  Message(STATUS "DATE is installed on this system")
  Execute_process(COMMAND date-config --cflags OUTPUT_VARIABLE DATEFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
  string(REGEX REPLACE "\"" "'" DATEFLAGS ${DATEFLAGS})
  #separate_arguments(DATEFLAGS)
  execute_process(COMMAND date-config --monitorlibs=dyn OUTPUT_VARIABLE DMONLIBS OUTPUT_STRIP_TRAILING_WHITESPACE)
  separate_arguments(DMONLIBS)
  list(APPEND DMONLIBS "-L$ENV{DIMDIR}/$ENV{ODIR} -ldim")
  set(DATE_COMMON_DEFS $ENV{DATE_COMMON_DEFS})
  set(DATE_MONITOR_DIR $ENV{DATE_MONITOR_DIR})
Else(DATE_PROGRAM)
  Set(DATE_FOUND NO)
  Message(STATUS "DATE is not installed on this system")
  If(UNIX)
    Execute_process(
      COMMAND uname 
      OUTPUT_VARIABLE _uname
      OUTPUT_STRIP_TRAILING_WHITESPACE)
  Else(UNIX)
    Set(_uname "Unknown")
  Endif(UNIX)
  Set(DATEFLAGS "-D${_uname} -DDATE_SYS=${_uname} -Dlong32='int' -Dlong64='long long' -DdatePointer='long'")
Endif(DATE_PROGRAM)
set(__DATEFLAGS ${DATEFLAGS})
separate_arguments(__DATEFLAGS)
Set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${DATEFLAGS})


