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
#include <vector>

#include <TObjArray.h>
#include <THashList.h>
#include <TH1.h>
#include <TString.h>

#include "AliEMCALTriggerConstants.h"
#include "AliEMCALTriggerFastOR.h"
#include "AliEMCALTriggerQA.h"
#include "AliEMCALTriggerPatchInfo.h"
#include "AliEMCALTriggerBitConfig.h"
#include "AliEMCALGeometry.h"

#include "AliHLTCTPData.h"

#include "AliHLTCaloTriggerPatchDataStruct.h"
#include "AliHLTCaloTriggerHeaderStruct.h"
#include "AliHLTCaloTriggerDataStruct.h"

#include "AliHLTEMCALDefinitions.h"
#include "AliHLTEMCALGeometry.h"

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

int AliHLTEMCALTriggerQAComponent::DoEvent(const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks,
    AliHLTComponentTriggerData& /*trigData*/, AliHLTUInt8_t* /*outputPtr*/, AliHLTUInt32_t& /*size*/,
    std::vector<AliHLTComponentBlockData>& /*outputBlocks*/)
{
  //patch in order to skip calib events
  if (!IsDataEvent()) return 0;

  TString triggerstring = "";
  const AliHLTCTPData * ctpdata = this->CTPData();
  if (ctpdata) {
    AliHLTComponentTriggerData trgdat;
    AliHLTTriggerMask_t mask = ctpdata->ActiveTriggers(trgdat);
    for (int index=0; index<gkNCTPTriggerClasses; index++) {
      if ((mask&(AliHLTTriggerMask_t(0x1)<<index)) == 0) continue;
      if(triggerstring.Length()) triggerstring += " ";
      triggerstring += ctpdata->Name(index);
    }
  }

  if (!blocks) return 0;

  const AliHLTComponentBlockData* iter = 0;

  for (Int_t ndx = 0; ndx < evtData.fBlockCnt; ndx++) {
    iter = blocks + ndx;

    if(!CheckInputDataType(iter->fDataType)) continue;

    if (iter->fDataType == AliHLTEMCALDefinitions::fgkTriggerPatchDataType) {
      ProcessTriggerPatches(iter);
    }
    else if (iter->fDataType == (kAliHLTDataTypeCaloTrigger | kAliHLTDataOriginEMCAL)) {
      ProcessTriggerFastors(iter);
    }
  }

  fTriggerQAPtr->EventCompleted();

  fLocalEventCount++;

  TIter next(fTriggerQAPtr->GetListOfHistograms());

  TH1* histo = 0;

  while ((histo = static_cast<TH1*>(next()))) {
    if (histo->GetEntries() > 0) {
      Int_t nbytes = PushBack(histo, kAliHLTDataTypeHistogram | kAliHLTDataOriginEMCAL , 0);
      if (fHistoResetOnPush && nbytes > 0) {
        histo->Reset();
      }
    }
  }

  return 0;
}

void AliHLTEMCALTriggerQAComponent::ProcessTriggerPatches(const AliHLTComponentBlockData* block)
{
  AliHLTCaloTriggerPatchDataStruct* hltpatchPtr = reinterpret_cast<AliHLTCaloTriggerPatchDataStruct*>(block->fPtr);
  if (!hltpatchPtr) return;

  UInt_t nPatches = block->fSize/sizeof(AliHLTCaloTriggerPatchDataStruct);

  AliEMCALTriggerPatchInfo patch;
  patch.SetTriggerBitConfig(fTriggerBitConfig);

  for(UInt_t ipatch = 0; ipatch < nPatches; ipatch++) {
    HLTPatch2Patch(hltpatchPtr[ipatch], patch);
    fTriggerQAPtr->ProcessPatch(&patch);
  }
}

void AliHLTEMCALTriggerQAComponent::ProcessTriggerFastors(const AliHLTComponentBlockData* block)
{
  AliHLTCaloTriggerPatchDataStruct* hltpatchPtr = reinterpret_cast<AliHLTCaloTriggerPatchDataStruct*>(block->fPtr);
  if (!hltpatchPtr) return;

  AliHLTCaloTriggerHeaderStruct* triggerhead = reinterpret_cast<AliHLTCaloTriggerHeaderStruct *>(block->fPtr);
  AliHLTCaloTriggerDataStruct *dataptr = reinterpret_cast<AliHLTCaloTriggerDataStruct *>(reinterpret_cast<AliHLTUInt8_t* >(block->fPtr) + sizeof(AliHLTCaloTriggerHeaderStruct));

  AliEMCALTriggerFastOR fastor;

  HLTDebug("Received %d fastor triggers", triggerhead->fNfastor);
  for(Int_t datacount = 0; datacount < triggerhead->fNfastor; datacount++) {
    HLTFastor2Fastor(dataptr[datacount], fastor);
    fTriggerQAPtr->ProcessFastor(&fastor);
  }
}

bool AliHLTEMCALTriggerQAComponent::CheckInputDataType(const AliHLTComponentDataType &datatype)
{
  vector <AliHLTComponentDataType> validTypes;
  GetInputDataTypes(validTypes);

  for(UInt_t i = 0; i < validTypes.size(); i++) {
    if (datatype == validTypes.at(i)) {
      return true;
    }
  }

  HLTDebug("Invalid Datatype");
  return false;
}

void AliHLTEMCALTriggerQAComponent::HLTPatch2Patch(const AliHLTCaloTriggerPatchDataStruct& htlpatch, AliEMCALTriggerPatchInfo& patch) const
{
  const TVector3 vect(0,0,0);
  patch.Initialize(htlpatch.fCol, htlpatch.fRow, htlpatch.fSize, htlpatch.fADC, htlpatch.fOfflineADC, EMCALTrigger::kEMCL1ADCtoGeV*htlpatch.fADC, htlpatch.fBitMask, vect, fGeometry->GetGeometryPtr());
}

void AliHLTEMCALTriggerQAComponent::HLTFastor2Fastor(const AliHLTCaloTriggerDataStruct& htlfastor, AliEMCALTriggerFastOR& fastor) const
{  
  // TO DO: here I always pass 8 as L0 time, need to check whether L0 times are in fact processed upstream
  fastor.Initialize(htlfastor.fAmplitude, htlfastor.fL1TimeSum, htlfastor.fRow, htlfastor.fCol, 8, fGeometry->GetGeometryPtr());
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
  list.push_back(kAliHLTDataTypeCaloTrigger | kAliHLTDataOriginEMCAL);
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
    if (option == "-newTriggerBitConfig") fTriggerBitConfig = new AliEMCALTriggerBitConfigNew;
    if (option == "-oldTriggerBitConfig") fTriggerBitConfig = new AliEMCALTriggerBitConfigOld;
    if (option == "-noHistoReset") fHistoResetOnPush = kFALSE;
    if (option.BeginsWith("-debugLevel")) {
      option.Remove(0, 11);
      Int_t dl = option.Atoi();
      if (dl >= 0) fTriggerQAPtr->SetDebugLevel(dl);
    }
    if (option.BeginsWith("-disableOffline")) {
      fTriggerQAPtr->EnablePatchType(AliEMCALTriggerQA::kOfflinePatch, kFALSE);
    }
    if (option.BeginsWith("-disableOnline")) {
      fTriggerQAPtr->EnablePatchType(AliEMCALTriggerQA::kOnlinePatch, kFALSE);
    }
    if (option.BeginsWith("-disableRecalc")) {
      fTriggerQAPtr->EnablePatchType(AliEMCALTriggerQA::kRecalcPatch, kFALSE);
    }
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
