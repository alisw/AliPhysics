//**************************************************************************\
//* This file is property of and copyright by the ALICE Project            *\
//* ALICE Experiment at CERN, All rights reserved.                         *\
//*                                                                        *\
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *\
//*                  for The ALICE HLT Project.                            *\
//*                                                                        *\
//* Permission to use, copy, modify and distribute this software and its   *\
//* documentation strictly for non-commercial purposes is hereby granted   *\
//* without fee, provided that the above copyright notice appears in all   *\
//* copies and that both the copyright notice and this permission notice   *\
//* appear in the supporting documentation. The authors make no claims     *\
//* about the suitability of this software for any purpose. It is          *\
//* provided "as is" without express or implied warranty.                  *\
//**************************************************************************

/// \file GPUCommonDef.h
/// \author David Rohr

#ifndef GPUCOMMONDEF_H
#define GPUCOMMONDEF_H

// clang-format off

//Some GPU configuration settings, must be included first
#include "GPUCommonDefSettings.h"

#if !(defined(__CINT__) || defined(__ROOTCINT__) || defined(__CLING__) || defined(__ROOTCLING__) || defined(G__ROOT)) //No GPU code for ROOT
  #if defined(__CUDACC__) || defined(__OPENCL__) || defined(__HIPCC__)
    #define GPUCA_GPUCODE //Compiled by GPU compiler
  #endif

  #if defined(__CUDA_ARCH__) || defined(__OPENCL__) || defined(__HIP_DEVICE_COMPILE__)
    #define GPUCA_GPUCODE_DEVICE //Executed on device
  #endif
#endif

//Definitions for C++11 features not supported by CINT / OpenCL
#if ((defined(__CINT__) || defined(__ROOTCINT__)) && !defined(__CLING__)) || (defined(__OPENCL__) && !defined(__OPENCLCPP__))
  #define CON_DELETE
  #define CON_DEFAULT
  #define CONSTEXPR const
#else
  #define CON_DELETE = delete
  #define CON_DEFAULT = default
  #define CONSTEXPR constexpr
#endif

//Set AliRoot / O2 namespace
#if defined(GPUCA_STANDALONE) || defined(GPUCA_O2_LIB) || defined(GPUCA_GPULIBRARY)
  #define GPUCA_ALIGPUCODE
#endif
#ifdef GPUCA_ALIROOT_LIB
  #define GPUCA_NAMESPACE AliGPU
#else
  #define GPUCA_NAMESPACE o2
#endif

//API Definitions for GPU Compilation
#include "GPUCommonDefAPI.h"

//Definitions steering enabling of GPU processing components
#if (!defined(__OPENCL__) || defined(__OPENCLCPP__)) && !defined(GPUCA_ALIROOT_LIB)
  #define GPUCA_BUILD_MERGER
  #define GPUCA_BUILD_DEDX
  #if defined(HAVE_O2HEADERS)
    #define GPUCA_BUILD_TRD
    #define GPUCA_BUILD_ITS
  #endif
#endif

// clang-format on

#endif
