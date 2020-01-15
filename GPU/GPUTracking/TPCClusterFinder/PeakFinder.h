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

/// \file PeakFinder.h
/// \author Felix Weiglhofer

#ifndef O2_GPU_PEAK_FINDER_H
#define O2_GPU_PEAK_FINDER_H

#include "clusterFinderDefs.h"
#include "GPUTPCClusterFinderKernels.h"

namespace GPUCA_NAMESPACE
{
namespace gpu
{

class PeakFinder
{

 public:
  static GPUd() void findPeaksImpl(int, int, int, int, GPUTPCClusterFinderKernels::GPUTPCSharedMemory&, GPUglobalref() const PackedCharge*, GPUglobalref() const deprecated::Digit*, uint, GPUglobalref() uchar*, GPUglobalref() uchar*);

 private:
  static GPUd() bool isPeakScratchPad(GPUTPCClusterFinderKernels::GPUTPCSharedMemory&, const deprecated::Digit*, ushort, GPUglobalref() const PackedCharge*, GPUsharedref() ChargePos*, GPUsharedref() PackedCharge*);

  static GPUd() bool isPeak(const deprecated::Digit*, GPUglobalref() const PackedCharge*);
};

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
