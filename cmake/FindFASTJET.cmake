# AliRoot Build System Module to find and configure FASTJET
#
# Author: Marco van Leeuwen
#
# Simple first version; checks env variables and fastjet-config; verifies that fastjet/PseudoJet.hh exists

cmake_minimum_required(VERSION 2.8.4 FATAL_ERROR)

# FASTJET_ROOT is a legacy env variable; could be removed
if (NOT FASTJET)
  set(FASTJET $ENV{FASTJET_ROOT})
endif(NOT FASTJET)

if (NOT FASTJET)
  set(FASTJET $ENV{FASTJET})
endif(NOT FASTJET)

if (NOT FASTJET)
  execute_process (COMMAND fastjet-config --prefix OUTPUT_VARIABLE FASTJET)
  string (STRIP "${FASTJET}" FASTJET)
endif(NOT FASTJET)

# Check for one of the header files
if (FASTJET)
  find_path(FASTJET include/fastjet/PseudoJet.hh PATHS ${FASTJET})
  if (FASTJET-NOTFOUND)
    message(STATUS "Header file fastjet/PseudoJet.hh not found in ${FASTJET}/include")
  endif (FASTJET-NOTFOUND)
endif (FASTJET)

if (FASTJET AND NOT FASTJET-NOTFOUND) 
  message(STATUS "FASTJET found in ${FASTJET}")
  set (FASTJET_FOUND true)
  set (FASTJET_INCLUDE_DIR ${FASTJET}/include)
  message(STATUS "FASTJET include directory: ${FASTJET_INCLUDE_DIR}")
  set (FASTJET_DEFINITIONS -DHAVE_FASTJET)
else()
  message(STATUS "FASTJET not found; make sure to set the FASTJET_ROOT environment variable or have fastjet-config in your path to use fastjet")
  set (FASTJET)
endif(FASTJET AND NOT FASTJET-NOTFOUND)

