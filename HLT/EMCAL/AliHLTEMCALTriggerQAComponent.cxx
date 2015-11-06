/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Salvatore Aiola                                       *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
#include <TObjArray.h>
#include <THashList.h>
#include <TH1.h>

#include "AliHLTCaloTriggerPatchDataStruct.h"

#include "AliHLTEMCALDefinitions.h"
#include "AliEMCALTriggerQA.h"
#include "AliEMCALTriggerPatchInfo.h"
#include "AliEMCALTriggerBitConfig.h"
#include "AliHLTEMCALGeometry.h"
#include "AliEMCALGeometry.h"

#include "AliHLTEMCALTriggerQAComponent.h"

ClassImp(AliHLTEMCALTriggerQAComponent)

AliHLTEMCALTriggerQAComponent::AliHLTEMCALTriggerQAComponent() :
  AliHLTCaloProcessor(),
  fTriggerBitConfig(NULL),
  fHistoResetOnPush(kTRUE),
  fLocalEventCount(0),
  fGeometry(NULL),
  fTriggerQAPtr(NULL)
{
}

AliHLTEMCALTriggerQAComponent::~AliHLTEMCALTriggerQAComponent()
{
  if (fTriggerQAPtr) delete fTriggerQAPtr;
  if (fGeometry) delete fGeometry;
}

int AliHLTEMCALTriggerQAComponent::DoEvent(const AliHLTComponentEventData& /*evtData*/, const AliHLTComponentBlockData* /*blocks*/,
					   AliHLTComponentTriggerData& /*trigData*/, AliHLTUInt8_t* /*outputPtr*/, AliHLTUInt32_t& /*size*/,
					   std::vector<AliHLTComponentBlockData>& /*outputBlocks*/)
{
  //patch in order to skip calib events
  if (!IsDataEvent()) return 0;

  // This component subscribes to a single input component
  const AliHLTComponentBlockData* block = GetFirstInputBlock(AliHLTEMCALDefinitions::fgkTriggerPatchDataType);
  if (!block) return 0;
  
  AliHLTCaloTriggerPatchDataStruct* hltpatchPtr = reinterpret_cast<AliHLTCaloTriggerPatchDataStruct*>(block->fPtr);
  if (!hltpatchPtr) return 0;

  UInt_t nPatches = block->fSize/sizeof(AliHLTCaloTriggerPatchDataStruct);
  UInt_t specification = block->fSpecification;
    
  AliEMCALTriggerPatchInfo patch;
  patch.SetTriggerBitConfig(fTriggerBitConfig);
  
  for(UInt_t ipatch = 0; ipatch < nPatches; ipatch++) {
    HLTPatch2Patch(hltpatchPtr[ipatch], patch);
    fTriggerQAPtr->ProcessPatch(&patch);
  }
  
  fTriggerQAPtr->EventCompleted();
  
  fLocalEventCount++;

  TIter next(fTriggerQAPtr->GetListOfHistograms());

  TH1* histo = 0;

  while ((histo = static_cast<TH1*>(next()))) {
    if (histo->GetEntries() > 0) {
      Int_t nbytes = PushBack(histo, kAliHLTDataTypeHistogram | kAliHLTDataOriginEMCAL , specification);
      if (fHistoResetOnPush && nbytes > 0) {
	histo->Reset();
      }
    }
  }

  return 0;
}

void AliHLTEMCALTriggerQAComponent::HLTPatch2Patch(const AliHLTCaloTriggerPatchDataStruct& htlpatch, AliEMCALTriggerPatchInfo& patch) const
{
  const TVector3 vect(0,0,0);
  
  patch.Initialize(htlpatch.fCol, htlpatch.fRow, htlpatch.fSize, htlpatch.fADC, 0, 0., htlpatch.fBitMask, vect, fGeometry->GetGeometryPtr());
}

const char* AliHLTEMCALTriggerQAComponent::GetComponentID()
{
  //See headerfile for documentation
  return "EmcalTriggerQA";
}

void AliHLTEMCALTriggerQAComponent::GetInputDataTypes(std::vector<AliHLTComponentDataType>& list)
{
  list.clear();
  list.push_back(AliHLTEMCALDefinitions::fgkTriggerPatchDataType);
}

AliHLTComponentDataType AliHLTEMCALTriggerQAComponent::GetOutputDataType()
{
  return kAliHLTDataTypeHistogram | kAliHLTDataOriginEMCAL;
}

void AliHLTEMCALTriggerQAComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
{
  //see header file for documentation
  constBase = 250000000;
  // to be reviewed later
  inputMultiplier = 0;
}

AliHLTComponent* AliHLTEMCALTriggerQAComponent::Spawn()
{
  return new AliHLTEMCALTriggerQAComponent();
}


int AliHLTEMCALTriggerQAComponent::DoInit(int argc, const char** argv)
{
  InitialiseGeometry();
  fTriggerQAPtr = new AliEMCALTriggerQA;
  fTriggerQAPtr->Init();

  for (int i = 0; i < argc; i++) {
    TString option(argv[i]);
    if (option.Contains("-newTriggerBitConfig")) fTriggerBitConfig = new AliEMCALTriggerBitConfigNew;
    if (option.Contains("-oldTriggerBitConfig")) fTriggerBitConfig = new AliEMCALTriggerBitConfigOld;
    if (option.Contains("-noHistoReset")) fHistoResetOnPush = kFALSE;
  }

  // if not specified use new trigger bit config
  if (!fTriggerBitConfig) fTriggerBitConfig = new AliEMCALTriggerBitConfigNew;
  
  return 0;
}

int AliHLTEMCALTriggerQAComponent::Deinit()
{
  if (fTriggerQAPtr) {
    delete fTriggerQAPtr;
    fTriggerQAPtr = 0;
  }

  if (fGeometry) {
    delete fGeometry;
    fGeometry = 0;
  }
  
  return 0;
}

void AliHLTEMCALTriggerQAComponent::InitialiseGeometry()
{
  fGeometry = new AliHLTEMCALGeometry(GetRunNo());
}
