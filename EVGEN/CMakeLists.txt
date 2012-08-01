# AliRoot Build System CMakeLists for EVGEN
#
# Author: Johny Jose m(johny.jose@cern.ch)
#         Port of previous Makefile build to cmake

cmake_minimum_required(VERSION 2.8.4 FATAL_ERROR)

file(GLOB PACKAGES CMake*.pkg)

ALICE_BuildModule()
