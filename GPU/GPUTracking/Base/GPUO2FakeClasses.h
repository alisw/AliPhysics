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

/// \file GPUO2FakeClasses.h
/// \author David Rohr

#ifndef O2_GPU_GPUO2FAKECLASSES_H
#define O2_GPU_GPUO2FAKECLASSES_H

#include "GPUCommonDef.h"
#include "GPUDataTypes.h"

// These are some dummies of O2 classes needed by AliGPU, to be used when O2 header unavailable

namespace o2
{
namespace gpu
{
} // namespace gpu
} // namespace o2

namespace o2
{
namespace tpc
{
struct ClusterNative {
  GPUd() static float getTime() { return 0.f; }
  GPUd() static float getPad() { return 0.f; }
  GPUd() static int getFlags() { return 0; }
  GPUd() static void setTimeFlags(float t, int f) {}
  GPUd() static void setPad(float p) {}
  GPUd() static void setSigmaTime(float s) {}
  GPUd() static void setSigmaPad(float s) {}

  unsigned char qTot, qMax;
};
struct ClusterNativeAccess {
  const ClusterNative* clustersLinear;
  const ClusterNative* clusters[GPUCA_NSLICES][GPUCA_ROW_COUNT];
  unsigned int nClusters[GPUCA_NSLICES][GPUCA_ROW_COUNT];
  unsigned int nClustersSector[GPUCA_NSLICES];
  unsigned int clusterOffset[GPUCA_NSLICES][GPUCA_ROW_COUNT];
  unsigned int nClustersTotal;
  void setOffsetPtrs() {}
};
} // namespace tpc
namespace base
{
struct MatBudget {
};
class MatLayerCylSet
{
};
} // namespace base
namespace trd
{
class TRDGeometryFlat
{
};
} // namespace trd
} // namespace o2

namespace GPUCA_NAMESPACE
{
namespace gpu
{
namespace deprecated
{
struct PackedDigit {
};
} // namespace deprecated
class GPUFakeEmpty
{
};
class GPUITSFitter
{
};
class GPUTPCConvert
{
};
class GPUTPCCompression
{
 public:
  GPUFakeEmpty mOutput;
};
class GPUTPCClusterFinder
{
};
#ifndef __OPENCL__
struct GPUParam;
class GPUTPCClusterStatistics
{
 public:
  void Finish() {}
  void RunStatistics(const o2::tpc::ClusterNativeAccess* clustersNative, const GPUFakeEmpty* clustersCompressed, const GPUParam& param) {}
};
#endif
} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
