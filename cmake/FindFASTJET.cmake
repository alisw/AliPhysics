# AliRoot Build System Module to find and configure FASTJET
#
# Author: Marco van Leeuwen
#
# Simple first version; checks env variables and fastjet-config; does not verify that includes and libs are really present

cmake_minimum_required(VERSION 2.8.4 FATAL_ERROR)

set(FASTJET_ROOT $ENV{FASTJET_ROOT})

if (NOT FASTJET_ROOT)
  set(FASTJET_ROOT $ENV{FASTJET})
endif(NOT FASTJET_ROOT)

if (NOT FASTJET_ROOT)
  execute_process (COMMAND fastjet-config --prefix OUTPUT_VARIABLE FASTJET_ROOT)
  string (STRIP "${FASTJET_ROOT}" FASTJET_ROOT)
endif(NOT FASTJET_ROOT)

if (FASTJET_ROOT) 
  message(STATUS "FASTJET found in ${FASTJET_ROOT}")
  set (FASTJET_FOUND true)
  set (FASTJET_INCLUDE_DIR ${FASTJET_ROOT}/include)
  message(STATUS "FASTJET include directory: ${FASTJET_INCLUDE_DIR}")
  set (FASTJET_DEFINITIONS -DHAVE_FASTJET)
else()
  message(STATUS "FASTJET not found; make sure to set the FASTJET_ROOT environment variable or have fastjet-config in your path to use fastjet")
endif(FASTJET_ROOT)

