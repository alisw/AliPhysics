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

/// \file Clusterizer.h
/// \author Felix Weiglhofer

#ifndef O2_GPU_CLUSTERIZER_H
#define O2_GPU_CLUSTERIZER_H

#include "clusterFinderDefs.h"
#include "GPUTPCClusterFinderKernels.h"
#include "PackedCharge.h"

namespace GPUCA_NAMESPACE
{
namespace gpu
{

class ClusterAccumulator;

namespace deprecated
{
class ClusterNavite;
}

class Clusterizer
{

 public:
  static GPUd() void computeClustersImpl(int, int, int, int, GPUTPCClusterFinderKernels::GPUTPCSharedMemory&, GPUglobalref() const PackedCharge*, GPUglobalref() const deprecated::Digit*, uint, uint, GPUglobalref() uint*, GPUglobalref() deprecated::ClusterNative*);

 private:
  static GPUd() void addOuterCharge(GPUglobalref() const PackedCharge*, ClusterAccumulator*, GlobalPad, Timestamp, Delta, Delta);

  static GPUd() Charge addInnerCharge(GPUglobalref() const PackedCharge*, ClusterAccumulator*, GlobalPad, Timestamp, Delta, Delta);

  static GPUd() void addCorner(GPUglobalref() const PackedCharge*, ClusterAccumulator*, GlobalPad, Timestamp, Delta, Delta);

  static GPUd() void addLine(GPUglobalref() const PackedCharge*, ClusterAccumulator*, GlobalPad, Timestamp, Delta, Delta);

  static GPUd() void updateClusterScratchpadInner(ushort, ushort, GPUsharedref() const PackedCharge*, ClusterAccumulator*, GPUsharedref() uchar*);

  static GPUd() void updateClusterScratchpadOuter(ushort, ushort, ushort, ushort, GPUsharedref() const PackedCharge*, ClusterAccumulator*);

  static GPUd() void buildClusterScratchPad(GPUglobalref() const PackedCharge*, ChargePos, GPUsharedref() ChargePos*, GPUsharedref() PackedCharge*, GPUsharedref() uchar*, ClusterAccumulator*);

  static GPUd() void buildClusterNaive(GPUglobalref() const PackedCharge*, ClusterAccumulator*, GlobalPad, Timestamp);

  static GPUd() void sortIntoBuckets(const deprecated::ClusterNative*, const uint, const uint, GPUglobalref() uint*, GPUglobalref() deprecated::ClusterNative*);
};

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
