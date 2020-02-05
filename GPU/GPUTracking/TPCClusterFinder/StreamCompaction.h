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

/// \file StreamCompaction.h
/// \author Felix Weiglhofer

#ifndef O2_GPU_STREAM_COMPACTION_H
#define O2_GPU_STREAM_COMPACTION_H

#include "clusterFinderDefs.h"
#include "GPUTPCClusterFinderKernels.h"

namespace GPUCA_NAMESPACE
{
namespace gpu
{

class StreamCompaction
{

 public:
  static GPUd() void nativeScanUpStartImpl(int, int, int, int, GPUTPCClusterFinderKernels::GPUTPCSharedMemory&,
                                           const uchar*, int*, int*,
                                           int);

  static GPUd() void nativeScanUpImpl(int, int, int, int, GPUTPCClusterFinderKernels::GPUTPCSharedMemory&,
                                      int*, int*, int);

  static GPUd() void nativeScanTopImpl(int, int, int, int, GPUTPCClusterFinderKernels::GPUTPCSharedMemory&,
                                       int*, int);

  static GPUd() void nativeScanDownImpl(int, int, int, int, GPUTPCClusterFinderKernels::GPUTPCSharedMemory&,
                                        int*, const int*, unsigned int, int);

  static GPUd() void compactDigitImpl(int, int, int, int, GPUTPCClusterFinderKernels::GPUTPCSharedMemory&,
                                      const deprecated::Digit*, deprecated::Digit*,
                                      const uchar*, int*, const int*,
                                      int);
};

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
