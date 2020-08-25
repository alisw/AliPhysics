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

/// \file GPULogging.h
/// \author David Rohr

#ifndef GPULOGGING_H
#define GPULOGGING_H

#include "GPUCommonDef.h"
// clang-format off
#if !defined(GPUCA_NOCOMPAT)
  // Cannot do anything for ROOT5CINT / OpenCL1, so just disable
  #define GPUInfo(...)
  #define GPUImportant(...)
  #define GPUWarning(...)
  #define GPUError(...)
  #define GPUFatal(...)
#elif defined(GPUCA_STANDALONE) && !defined(GPUCA_GPUCODE_DEVICE) && !defined(GPUCA_NO_FMT) && !defined(__HIPCC__)
  #include <fmt/printf.h>
  #define GPUInfo(string, ...)                 \
    {                                          \
      fmt::printf(string "\n", ##__VA_ARGS__); \
    }
  #define GPUImportant(...) GPUInfo(__VA_ARGS__)
  #define GPUWarning(string, ...)                       \
    {                                                   \
      fmt::fprintf(stderr, string "\n", ##__VA_ARGS__); \
    }
  #define GPUError(...) GPUWarning(__VA_ARGS__)
  #define GPUFatal(string, ...)                         \
    {                                                   \
      fmt::fprintf(stderr, string "\n", ##__VA_ARGS__); \
      throw std::exception();                           \
    }
#elif defined(GPUCA_STANDALONE) || defined(GPUCA_GPUCODE_DEVICE) || (defined(GPUCA_ALIROOT_LIB) && defined(GPUCA_GPUCODE) && defined(__cplusplus) && __cplusplus < 201703L) || defined(__HIPCC__)
  // For standalone / CUDA / HIP, we just use printf, which should be available
  // Temporarily, we also have to handle CUDA on AliRoot with O2 defaults due to ROOT / CUDA incompatibilities
  #define GPUInfo(string, ...)            \
    {                                     \
      printf(string "\n", ##__VA_ARGS__); \
    }
  #define GPUImportant(...) GPUInfo(__VA_ARGS__)
  #if defined(GPUCA_GPUCODE_DEVICE) || defined(__HIPCC__)
    #define GPUWarning(...) GPUInfo(__VA_ARGS__)
    #define GPUError(...) GPUInfo(__VA_ARGS__)
    #define GPUFatal(...) GPUInfo(__VA_ARGS__)
  #else
    #define GPUWarning(string, ...)                  \
      {                                              \
        fprintf(stderr, string "\n", ##__VA_ARGS__); \
      }
    #define GPUError(...) GPUWarning(__VA_ARGS__)
    #ifdef GPUCA_NOCOMPAT
      #define GPUFatal(string, ...)                    \
        {                                              \
          fprintf(stderr, string "\n", ##__VA_ARGS__); \
          throw std::exception();                      \
        }
    #else
      #define GPUFatal(string, ...)                  \
        {                                            \
          fprintf(stderr, string "\n", __VA_ARGS__); \
          exit(1);                                   \
        }
    #endif
  #endif
#elif defined(GPUCA_ALIROOT_LIB)
  // Forward to HLT Logging functions for AliRoot
  #include "AliHLTLogging.h"
  #define GPUInfo(...) HLTInfo(__VA_ARGS__)
  #define GPUImportant(...) HLTImportant(__VA_ARGS__)
  #define GPUWarning(...) HLTWarning(__VA_ARGS__)
  #define GPUError(...) HLTError(__VA_ARGS__)
  #define GPUFatal(...) HLTFatal(__VA_ARGS__)
  // Workaround for static functions / classes not deriving from AliHLTLogging
  namespace AliGPU
  {
  namespace gpu
  {
  // We pollute the AliGPU::gpu namespace with some anonymous functions that catch the HLT...() magic
  namespace
  {
  AliHLTLogging gAliGPULog; // This creates a couple of bogus instances, but there are plenty anyway
  template <typename... Args>
  void LoggingVarargs(Args... args)
  {
    gAliGPULog.LoggingVarargs(args...);
  }
  template <typename... Args>
  bool CheckFilter(Args... args)
  {
    return gAliGPULog.CheckFilter(args...);
  }
  const char* Class_Name() { return "GPU"; };
  } // namespace
  } // namespace gpu
  } // namespace AliGPU
#elif defined(GPUCA_O2_LIB) || defined(GPUCA_O2_INTERFACE)
  // Forward to O2 LOGF logginf for O2
  #include "GPUCommonLogger.h"
  #define GPUInfo(...) LOGF(info, __VA_ARGS__)
  #define GPUImportant(...) LOGF(info, __VA_ARGS__)
  #define GPUWarning(...) LOGF(warning, __VA_ARGS__)
  #define GPUError(...) LOGF(error, __VA_ARGS__)
  #define GPUFatal(...) LOGF(fatal, __VA_ARGS__)
#endif

// clang-format on

#endif // GPULOGGING_H
