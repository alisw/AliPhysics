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

/// \file GPUO2Interface.h
/// \author David Rohr

#ifndef GPUO2INTERFACE_H
#define GPUO2INTERFACE_H

// Some defines denoting that we are compiling for O2
#ifndef HAVE_O2HEADERS
#define HAVE_O2HEADERS
#endif
#ifndef GPUCA_TPC_GEOMETRY_O2
#define GPUCA_TPC_GEOMETRY_O2
#endif
#ifndef GPUCA_O2_INTERFACE
#define GPUCA_O2_INTERFACE
#endif

#include <memory>
#include "GPUCommonDef.h"
#include "GPUDataTypes.h"
namespace o2::tpc
{
struct ClusterNativeAccess;
struct ClusterNative;
} // namespace o2::tpc

namespace o2::gpu
{
class GPUReconstruction;
class GPUChainTracking;
struct GPUO2InterfaceConfiguration;
struct GPUInterfaceOutputs;
struct GPUOutputControl;

class GPUTPCO2Interface
{
 public:
  GPUTPCO2Interface();
  ~GPUTPCO2Interface();

  int Initialize(const GPUO2InterfaceConfiguration& config);
  void Deinitialize();

  int RunTracking(GPUTrackingInOutPointers* data, GPUInterfaceOutputs* outputs = nullptr);
  void Clear(bool clearOutputs);

  bool GetParamContinuous() { return (mContinuous); }
  void GetClusterErrors2(int row, float z, float sinPhi, float DzDs, short clusterState, float& ErrY2, float& ErrZ2) const;

  int registerMemoryForGPU(const void* ptr, size_t size);
  int unregisterMemoryForGPU(const void* ptr);

  const GPUO2InterfaceConfiguration& getConfig() const { return *mConfig; }

 private:
  GPUTPCO2Interface(const GPUTPCO2Interface&);
  GPUTPCO2Interface& operator=(const GPUTPCO2Interface&);

  bool mInitialized = false;
  bool mContinuous = false;

  std::unique_ptr<GPUReconstruction> mRec;
  GPUChainTracking* mChain = nullptr;
  std::unique_ptr<GPUO2InterfaceConfiguration> mConfig;
  std::unique_ptr<GPUOutputControl> mOutputCompressedClusters;
};
} // namespace o2::gpu

#endif
