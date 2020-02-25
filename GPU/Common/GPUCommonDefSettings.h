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

/// \file GPUCommonDefSettings.h
/// \author David Rohr

// clang
#ifndef GPUCOMMONDEFSETTINGS_H
#define GPUCOMMONDEFSETTINGS_H

// clang-format off

#ifndef GPUCOMMONDEF_H
  #error Please include GPUCommonDef.h!
#endif

//#define GPUCA_OPENCL_CPP_CLANG_C11_ATOMICS     // Use C11 atomic instead of old style atomics for OpenCL C++ in clang (OpenCL 2.2 C++ will use C++11 atomics irrespectively)

//#define GPUCA_CUDA_NO_CONSTANT_MEMORY          // Do not use constant memory for CUDA
//#define GPUCA_HIP_NO_CONSTANT_MEMORY           // Do not use constant memory for HIP - MANDATORY for now since all AMD GPUs have insufficient constant memory with HIP
//#define GPUCA_OPENCL_NO_CONSTANT_MEMORY        // Do not use constant memory for OpenCL 1.2
#define GPUCA_OPENCLCPP_NO_CONSTANT_MEMORY       // Do not use constant memory for OpenCL C++ - MANDATORY as OpenCL cannot cast between __constant and __generic yet!

// clang-format on

#endif // GPUCOMMONDEFSETTINGS_H
