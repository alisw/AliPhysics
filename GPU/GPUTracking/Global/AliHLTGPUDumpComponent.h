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

/// \file AliHLTGPUDumpComponent.h
/// \author David Rohr

#ifndef ALIHLTGPUDUMPCOMPONENT_H
#define ALIHLTGPUDUMPCOMPONENT_H

#include "GPUCommonDef.h"
#include "AliHLTProcessor.h"

class AliTPCcalibDB;
class AliTPCRecoParam;
#include "AliRecoParam.h"
class AliTPCTransform;
namespace GPUCA_NAMESPACE
{
namespace gpu
{
class TPCFastTransform;
class TPCFastTransformManager;
class GPUReconstruction;
class GPUChainTracking;
class GPUTPCClusterData;
} // namespace gpu
} // namespace GPUCA_NAMESPACE

class AliHLTGPUDumpComponent : public AliHLTProcessor
{
 public:
  static const unsigned int NSLICES = 36;
  static const unsigned int NPATCHES = 6;

  AliHLTGPUDumpComponent();

  AliHLTGPUDumpComponent(const AliHLTGPUDumpComponent&) CON_DELETE;
  AliHLTGPUDumpComponent& operator=(const AliHLTGPUDumpComponent&) CON_DELETE;

  virtual ~AliHLTGPUDumpComponent();

  const char* GetComponentID();
  void GetInputDataTypes(vector<AliHLTComponentDataType>& list);
  AliHLTComponentDataType GetOutputDataType();
  virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);
  AliHLTComponent* Spawn();

 protected:
  int DoInit(int argc, const char** argv);
  int DoDeinit();
  int Reconfigure(const char* cdbEntry, const char* chainId);
  int DoEvent(const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks);

 private:
  float fSolenoidBz;
  GPUCA_NAMESPACE::gpu::GPUReconstruction* fRec;
  GPUCA_NAMESPACE::gpu::GPUChainTracking* fChain;
  GPUCA_NAMESPACE::gpu::TPCFastTransformManager* fFastTransformManager;
  AliTPCcalibDB* fCalib;
  AliTPCRecoParam* fRecParam;
  AliRecoParam fOfflineRecoParam;
  AliTPCTransform* fOrigTransform;
  bool fIsMC;
};

#endif
