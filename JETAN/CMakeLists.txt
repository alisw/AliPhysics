# AliRoot Build System CMakeLists for JETAN
#
# Author: Johny Jose m(johny.jose@cern.ch)
#         Port of previous Makefile build to cmake

cmake_minimum_required(VERSION 2.8.4 FATAL_ERROR)

# file(GLOB PACKAGES CMake*.pkg)
set(FASTJET_ROOT $ENV{FASTJET_ROOT})

if(NOT FASTJET_ROOT)
  set(PACKAGES CMakelibJETAN.pkg)
  message("JETAN: Fastjet not found; skipping FASTJETAN")
else()
  set(PACKAGES CMakelibFASTJETAN.pkg CMakelibJETAN.pkg)
endif(NOT FASTJET_ROOT)


ALICE_BuildModule()
