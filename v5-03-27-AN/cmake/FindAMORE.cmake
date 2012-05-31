# AliRoot Build System Module to find and configure AMORE
#
# Author: Johny Jose (johny.jose@cern.ch)
#         Port of previous Makefile build to cmake

cmake_minimum_required(VERSION 2.8.4 FATAL_ERROR)

find_program(AMORE_CONFIG NAMES amore-config)
if(AMORE_CONFIG)
  set(AMORE_FOUND TRUE)
  set(AMOREDEFINITIONS "-DALI_AMORE")
  execute_process(COMMAND ${AMORE_CONFIG} --cflags --includes OUTPUT_VARIABLE AMOREFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND ${AMORE_CONFIG} --ldflags --ldflags-da-static --libs OUTPUT_VARIABLE _lddaflags OUTPUT_STRIP_TRAILING_WHITESPACE)
  ALICE_CleanOutput(_lddaflags "${_lddaflags}")
  set(AMOREFLAGS "-DALI_AMORE ${AMOREFLAGS}")
  if(DATE_FOUND)
    execute_process(COMMAND ${DATE_CONFIG} --rcproxylibs OUTPUT_VARIABLE _rcproxylib OUTPUT_STRIP_TRAILING_WHITESPACE)
  else()
    message(STATUS "AMORE requires DATE")
  endif(DATE_FOUND)
  set(AMOREDALIBS "-static ${_lddaflags} ${_rcproxylib}")
  else()
  message(STATUS "AMORE not found")
endif(AMORE_CONFIG)
