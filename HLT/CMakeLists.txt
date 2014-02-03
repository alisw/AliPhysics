# AliRoot Build System CMakeLists for HLT
#
# Author: Johny Jose m(johny.jose@cern.ch)
#         Port of previous Makefile build to cmake

cmake_minimum_required(VERSION 2.8.4 FATAL_ERROR)

option(HLT_CUDA   "Build HLT GPU tracker using CUDA"   OFF)
option(HLT_OPENCL "Build HLT GPU tracker using OpenCL" OFF)

include("CMakehltdoc.conf")

file(GLOB PACKAGES CMake*.pkg)

if(HLT_CUDA)
  list(APPEND PACKAGES TPCLib/tracking-ca/cagpu/CMakelibAliHLTTPCCAGPU.pkg)
endif(HLT_CUDA)

if(HLT_OPENCL)
  list(APPEND PACKAGES TPCLib/tracking-ca/cagpu/CMakelibAliHLTTPCCAGPUOpenCL.pkg)
endif()

ALICE_BuildModule()
