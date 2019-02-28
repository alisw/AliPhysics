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

#define GPUError(...)    \
  {                      \
    printf(__VA_ARGS__); \
    printf("\n");        \
  }
#define GPUWarning(...)  \
  {                      \
    printf(__VA_ARGS__); \
    printf("\n");        \
  }
#define GPUInfo(...)     \
  {                      \
    printf(__VA_ARGS__); \
    printf("\n");        \
  }
#define GPUImportant(...) \
  {                       \
    printf(__VA_ARGS__);  \
    printf("\n");         \
  }
#define GPUFatal(...)    \
  {                      \
    printf(__VA_ARGS__); \
    printf("\n");        \
    exit(1);             \
  }

#endif // GPULOGGING_H
