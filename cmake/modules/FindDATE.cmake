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
  Execute_process(COMMAND date-config --cflags OUTPUT_VARIABLE DATEFLAGS)
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
  Set(DATEFLAGS "-D${_uname} -DDATE_SYS=${_uname} -Dlong32=\"int\" -Dlong64=\"long long\" -DdatePointer=\"long\"")
Endif(DATE_PROGRAM)

Set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${DATEFLAGS})


