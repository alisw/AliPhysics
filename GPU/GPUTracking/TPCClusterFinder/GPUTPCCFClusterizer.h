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
#include "GPUGeneralKernels.h"
#include "GPUConstantMem.h"
#include "GPUTPCClusterFinder.h"
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
class MCLabelAccumulator;

class GPUTPCCFClusterizer : public GPUKernelTemplate
{

 public:
  struct GPUSharedMemory {
    ChargePos posBcast[SCRATCH_PAD_WORK_GROUP_SIZE];
    PackedCharge buf[SCRATCH_PAD_WORK_GROUP_SIZE * SCRATCH_PAD_BUILD_N];
    uchar innerAboveThreshold[SCRATCH_PAD_WORK_GROUP_SIZE];
  };

  enum K : int {
    computeClusters = 0,
  };

  static GPUd() void computeClustersImpl(int, int, int, int, GPUSharedMemory&, const Array2D<PackedCharge>&, const ChargePos*, uint, MCLabelAccumulator*, uint, uint*, tpc::ClusterNative*);

#ifdef HAVE_O2HEADERS
  typedef GPUTPCClusterFinder processorType;
  GPUhdi() static processorType* Processor(GPUConstantMem& processors)
  {
    return processors.tpcClusterer;
  }
#endif

  GPUhdi() CONSTEXPR static GPUDataTypes::RecoStep GetRecoStep()
  {
    return GPUDataTypes::RecoStep::TPCClusterFinding;
  }

  template <int iKernel = defaultKernel, typename... Args>
  GPUd() static void Thread(int nBlocks, int nThreads, int iBlock, int iThread, GPUSharedMemory& smem, processorType& clusterer, Args... args);

  static GPUd() void computeClustersImpl(int, int, int, int, GPUSharedMemory&, const Array2D<PackedCharge>&, const ChargePos*, MCLabelAccumulator*, uint, uint, uint*, tpc::ClusterNative*);

 private:
  static GPUd() void addOuterCharge(const Array2D<PackedCharge>&, ClusterAccumulator*, const ChargePos&, Delta2);

  static GPUd() Charge addInnerCharge(const Array2D<PackedCharge>&, ClusterAccumulator*, const ChargePos&, Delta2);

  static GPUd() void addCorner(const Array2D<PackedCharge>&, ClusterAccumulator*, const ChargePos&, Delta2);

  static GPUd() void addLine(const Array2D<PackedCharge>&, ClusterAccumulator*, const ChargePos&, Delta2);

  static GPUd() void updateClusterScratchpadInner(ushort, ushort, const PackedCharge*, const ChargePos&, ClusterAccumulator*, MCLabelAccumulator*, uchar*);

  static GPUd() void updateClusterScratchpadOuter(ushort, ushort, ushort, ushort, const PackedCharge*, const ChargePos&, ClusterAccumulator*, MCLabelAccumulator*);

  static GPUd() void buildClusterScratchPad(const Array2D<PackedCharge>&, ChargePos, ChargePos*, PackedCharge*, uchar*, ClusterAccumulator*, MCLabelAccumulator*);

  static GPUd() void buildClusterNaive(const Array2D<PackedCharge>&, ClusterAccumulator*, const ChargePos&);

  static GPUd() uint sortIntoBuckets(const tpc::ClusterNative&, uint, uint, uint*, tpc::ClusterNative*);
};

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
