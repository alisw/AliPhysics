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

/// \file GPUCommonLogger.h
/// \author David Rohr

#ifndef GPUCOMMONFAIRLOGGER_H
#define GPUCOMMONFAIRLOGGER_H

#if defined(__OPENCL__)
#define LOG(...)
#define LOGF(...)
#elif defined(GPUCA_STANDALONE) ||                    \
  defined(GPUCA_ALIROOT_LIB) ||                       \
  defined(GPUCA_GPUCODE_DEVICE) ||                    \
  (!defined(__cplusplus) || __cplusplus < 201703L) || \
  (defined(__HIPCC__) && (!defined(_GLIBCXX_USE_CXX11_ABI) || _GLIBCXX_USE_CXX11_ABI == 0))

#include <iostream>
#include <cstdio>
#define LOG(type) std::cout
#define LOGF(type, string, ...)         \
  {                                     \
    printf(string "\n", ##__VA_ARGS__); \
  }

#else

#include <Framework/Logger.h>

#endif

#endif
