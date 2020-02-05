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

/// \file Deconvolution.h
/// \author Felix Weiglhofer

#ifndef O2_GPU_DECONVOLUTION_H
#define O2_GPU_DECONVOLUTION_H

#include "clusterFinderDefs.h"
#include "GPUTPCClusterFinderKernels.h"
#include "Array2D.h"
#include "PackedCharge.h"

namespace GPUCA_NAMESPACE
{
namespace gpu
{

class Deconvolution
{

 public:
  static GPUd() void countPeaksImpl(int, int, int, int, GPUTPCClusterFinderKernels::GPUTPCSharedMemory&, const Array2D<uchar>&, Array2D<PackedCharge>&, const deprecated::Digit*, const uint);

 private:
  static GPUd() char countPeaksAroundDigit(const ChargePos&, const Array2D<uchar>&);
  static GPUd() char countPeaksScratchpadInner(ushort, const uchar*, uchar*);
  static GPUd() char countPeaksScratchpadOuter(ushort, ushort, uchar, const uchar*);
};

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
