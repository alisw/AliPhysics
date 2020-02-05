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

/// \file ChargeMapFiller.h
/// \author Felix Weiglhofer

#ifndef O2_GPU_CHARGE_MAP_FILLER_H
#define O2_GPU_CHARGE_MAP_FILLER_H

#include "clusterFinderDefs.h"
#include "GPUTPCClusterFinderKernels.h"
#include "Array2D.h"

namespace GPUCA_NAMESPACE
{
namespace gpu
{

class ChargeMapFiller
{

 public:
  static GPUd() void fillChargeMapImpl(int, int, int, int, GPUTPCClusterFinderKernels::GPUTPCSharedMemory&, const deprecated::Digit*, Array2D<PackedCharge>&, size_t);

  static GPUd() void resetMapsImpl(int, int, int, int, GPUTPCClusterFinderKernels::GPUTPCSharedMemory&, const deprecated::Digit*, Array2D<PackedCharge>&, Array2D<uchar>&);
};

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
