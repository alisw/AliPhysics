# AliRoot Build System CMakeLists for JETAN
#
# Author: Johny Jose m(johny.jose@cern.ch)
#         Port of previous Makefile build to cmake

cmake_minimum_required(VERSION 2.8.4 FATAL_ERROR)

if(NOT FASTJET_FOUND)
  set(PACKAGES CMakelibJETAN.pkg)
  message(STATUS "skipping FASTJETAN")
else()
  set(PACKAGES CMakelibFASTJETAN.pkg CMakelibJETAN.pkg)
endif(NOT FASTJET_FOUND)


ALICE_BuildModule()
