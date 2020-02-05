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
#include "Array2D.h"
#include "PackedCharge.h"

namespace GPUCA_NAMESPACE
{

namespace tpc
{
struct ClusterNative;
}

namespace gpu
{

class ClusterAccumulator;

class Clusterizer
{

 public:
  static GPUd() void computeClustersImpl(int, int, int, int, GPUTPCClusterFinderKernels::GPUTPCSharedMemory&, const Array2D<PackedCharge>&, const deprecated::Digit*, uint, uint, uint*, tpc::ClusterNative*);

 private:
  static GPUd() void addOuterCharge(const Array2D<PackedCharge>&, ClusterAccumulator*, const ChargePos&, Delta2);

  static GPUd() Charge addInnerCharge(const Array2D<PackedCharge>&, ClusterAccumulator*, const ChargePos&, Delta2);

  static GPUd() void addCorner(const Array2D<PackedCharge>&, ClusterAccumulator*, const ChargePos&, Delta2);

  static GPUd() void addLine(const Array2D<PackedCharge>&, ClusterAccumulator*, const ChargePos&, Delta2);

  static GPUd() void updateClusterScratchpadInner(ushort, ushort, const PackedCharge*, ClusterAccumulator*, uchar*);

  static GPUd() void updateClusterScratchpadOuter(ushort, ushort, ushort, ushort, const PackedCharge*, ClusterAccumulator*);

  static GPUd() void buildClusterScratchPad(const Array2D<PackedCharge>&, ChargePos, ChargePos*, PackedCharge*, uchar*, ClusterAccumulator*);

  static GPUd() void buildClusterNaive(const Array2D<PackedCharge>&, ClusterAccumulator*, const ChargePos&);

  static GPUd() void sortIntoBuckets(const tpc::ClusterNative&, const uint, const uint, uint*, tpc::ClusterNative*);
};

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
