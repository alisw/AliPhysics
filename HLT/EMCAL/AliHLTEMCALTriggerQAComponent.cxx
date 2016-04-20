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
#include <TStopwatch.h>

#include "AliEMCALTriggerConstants.h"
#include "AliEMCALTriggerFastOR.h"
#include "AliEMCALTriggerOnlineQAPbPb.h"
#include "AliEMCALTriggerOnlineQAPP.h"
#include "AliEMCALTriggerPatchInfo.h"
#include "AliEMCALTriggerBitConfig.h"
#include "AliEMCALGeometry.h"

#include "AliHLTCTPData.h"

#include "AliHLTCaloTriggerPatchDataStruct.h"
#include "AliHLTCaloTriggerHeaderStruct.h"
#include "AliHLTCaloTriggerDataStruct.h"

#include "AliHLTEMCALDefinitions.h"
#include "AliHLTEMCALGeometry.h"
#include "AliHLTEMCALCaloCells.h"

#include "AliHLTEMCALTriggerQAComponent.h"

ClassImp(AliHLTEMCALTriggerQAComponent)

AliHLTEMCALTriggerQAComponent::AliHLTEMCALTriggerQAComponent() :
  AliHLTCaloProcessor(),
  fTriggerBitConfig(NULL),
  fHistoResetOnPush(kTRUE),
  fFilterTrgClass(""),
  fBeamType(kPP),
  fLocalEventCount(0),
  fGeometry(NULL),
  fTriggerQAPtr(NULL),
  fFiredTriggerClasses()
{
  SetPP2016TriggerClasses();
}

AliHLTEMCALTriggerQAComponent::~AliHLTEMCALTriggerQAComponent()
{
  if (fTriggerQAPtr) delete fTriggerQAPtr;
  if (fGeometry) delete fGeometry;
}

int AliHLTEMCALTriggerQAComponent::RetrieveFiredTriggerClasses()
{
  int nTrgClasses = 0;
  fFiredTriggerClasses.clear();
  const AliHLTCTPData * ctpdata = this->CTPData();
  if (ctpdata) {
    AliHLTComponentTriggerData trgdat;
    AliHLTTriggerMask_t mask = ctpdata->ActiveTriggers(trgdat);
    for (int index=0; index<gkNCTPTriggerClasses; index++) {
      if ((mask&(AliHLTTriggerMask_t(0x1)<<index)) == 0) continue;
      TString trgClass = ctpdata->Name(index);
      if (!fFilterTrgClass.IsNull() && !fFilterTrgClass.Contains(trgClass)) continue;
      fFiredTriggerClasses.push_back(trgClass);
      nTrgClasses++;
    }
  }
  return nTrgClasses;
}

int AliHLTEMCALTriggerQAComponent::DoEvent(const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks,
    AliHLTComponentTriggerData& /*trigData*/, AliHLTUInt8_t* /*outputPtr*/, AliHLTUInt32_t& /*size*/,
    std::vector<AliHLTComponentBlockData>& /*outputBlocks*/)
{
  static AliHLTEMCALCaloCells cells;
  cells.Clear();

  //patch in order to skip calib events
  if (!IsDataEvent()) return 0;

  if (!blocks) return 0;

#ifdef __PROFILE__
  TStopwatch profile;
  profile.Start();
#endif

  if (!RetrieveFiredTriggerClasses()) {
    HLTDebug("No trigger classes received\n");
  }

  const AliHLTComponentBlockData* iter = 0;

  const AliHLTComponentBlockData* patchData = 0;
  const AliHLTComponentBlockData* fastorData = 0;
  const AliHLTComponentBlockData* cellData = 0;

  for (Int_t ndx = 0; ndx < evtData.fBlockCnt; ndx++) {
    iter = blocks + ndx;

    if(!CheckInputDataType(iter->fDataType)) continue;

    if (iter->fDataType == AliHLTEMCALDefinitions::fgkDigitDataType) {
      cellData = iter;
    }
    else if (iter->fDataType == AliHLTEMCALDefinitions::fgkTriggerPatchDataType) {
      patchData = iter;
    }
    else if (iter->fDataType == (kAliHLTDataTypeCaloTrigger | kAliHLTDataOriginEMCAL)) {
      fastorData = iter;
    }
  }

  ProcessCells(cellData, cells);
  ProcessTriggerPatches(patchData);
  ProcessTriggerFastors(fastorData, cells);

  fTriggerQAPtr->EventCompleted();

  fLocalEventCount++;

#ifdef __PROFILE__
  profile.Stop();
  printf("End of trigger maker QA: %f (Wall) / %f (CPU)\n", profile.RealTime(), profile.CpuTime());
#endif

  PushHistograms(fTriggerQAPtr->GetListOfHistograms());

  return 0;
}

void AliHLTEMCALTriggerQAComponent::ProcessCells(const AliHLTComponentBlockData* block, AliHLTEMCALCaloCells& cells)
{
  if (!block) return;
  AliHLTCaloDigitDataStruct* digit = reinterpret_cast<AliHLTCaloDigitDataStruct*>(block->fPtr);
  if (!digit) return;

#ifdef __PROFILE__
  TStopwatch profile;
  profile.Start();
#endif

  Int_t nDigits = block->fSize/sizeof(AliHLTCaloDigitDataStruct);
  cells.CreateContainer(nDigits);

  for (Int_t idigit = 0; idigit < nDigits; idigit++){
    cells.SetCell(idigit, digit->fID, digit->fEnergy, digit->fTime);
    digit++;
  }

#ifdef __PROFILE__
  profile.Stop();
  printf("End of ProcessCells: %f (Wall) / %f (CPU)\n", profile.RealTime(), profile.CpuTime());
#endif
}

void AliHLTEMCALTriggerQAComponent::PushHistograms(TCollection* list)
{
  TIter next(list);

  TH1* histo = 0;

  while ((histo = static_cast<TH1*>(next()))) {
    if (histo->GetEntries() > 0) {
      Int_t nbytes = PushBack(histo, kAliHLTDataTypeHistogram | kAliHLTDataOriginEMCAL , 0);
      if (fHistoResetOnPush && nbytes > 0) {
        histo->Reset();
      }
    }
  }
}

void AliHLTEMCALTriggerQAComponent::ProcessTriggerPatches(const AliHLTComponentBlockData* block)
{
  if (!block) return;
  AliHLTCaloTriggerPatchDataStruct* hltpatchPtr = reinterpret_cast<AliHLTCaloTriggerPatchDataStruct*>(block->fPtr);
  if (!hltpatchPtr) return;

#ifdef __PROFILE__
  TStopwatch profile;
  profile.Start();
#endif

  UInt_t nPatches = block->fSize/sizeof(AliHLTCaloTriggerPatchDataStruct);

  AliEMCALTriggerPatchInfo patch;
  patch.SetTriggerBitConfig(fTriggerBitConfig);
  patch.SetOffSet(fTriggerBitConfig->GetTriggerTypesEnd());

  if (fBeamType == kPbPb) {
    for(UInt_t ipatch = 0; ipatch < nPatches; ipatch++) {
      HLTPatch2Patch(hltpatchPtr[ipatch], patch);
      fTriggerQAPtr->ProcessBkgPatch(&patch);
    }

    fTriggerQAPtr->ComputeBackground();
  }

  HLTDebug("Number of patches: %d", nPatches);
  for(UInt_t ipatch = 0; ipatch < nPatches; ipatch++) {
    HLTPatch2Patch(hltpatchPtr[ipatch], patch);
    fTriggerQAPtr->ProcessPatch(&patch);
  }
#ifdef __PROFILE__
  profile.Stop();
  printf("End of ProcessTriggerPatches: %f (Wall) / %f (CPU)\n", profile.RealTime(), profile.CpuTime());
#endif
}

void AliHLTEMCALTriggerQAComponent::ProcessTriggerFastors(const AliHLTComponentBlockData* block, AliHLTEMCALCaloCells& cells)
{
  if (!block) return;
  AliHLTCaloTriggerPatchDataStruct* hltpatchPtr = reinterpret_cast<AliHLTCaloTriggerPatchDataStruct*>(block->fPtr);
  if (!hltpatchPtr) return;

#ifdef __PROFILE__
  TStopwatch profile;
  profile.Start();
#endif

  AliHLTCaloTriggerHeaderStruct* triggerhead = reinterpret_cast<AliHLTCaloTriggerHeaderStruct *>(block->fPtr);
  AliHLTCaloTriggerDataStruct *dataptr = reinterpret_cast<AliHLTCaloTriggerDataStruct *>(reinterpret_cast<AliHLTUInt8_t* >(block->fPtr) + sizeof(AliHLTCaloTriggerHeaderStruct));

  AliEMCALTriggerFastOR fastor;

  HLTDebug("Received %d fastor triggers", triggerhead->fNfastor);
  for(Int_t datacount = 0; datacount < triggerhead->fNfastor; datacount++) {
    HLTFastor2Fastor(dataptr[datacount], fastor);
    fTriggerQAPtr->ProcessFastor(&fastor, &cells);
  }

#ifdef __PROFILE__
  profile.Stop();
  printf("End of ProcessTriggerFastors: %f (Wall) / %f (CPU)\n", profile.RealTime(), profile.CpuTime());
#endif
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

void AliHLTEMCALTriggerQAComponent::HLTFastor2Fastor(const AliHLTCaloTriggerDataStruct& hltfastor, AliEMCALTriggerFastOR& fastor) const
{  
  Int_t l0time = 0;
  if (hltfastor.fNL0Times > 0) {
    l0time = hltfastor.fL0Times[0];
  }
  fastor.Initialize(hltfastor.fAmplitude, hltfastor.fL1TimeSum, hltfastor.fRow, hltfastor.fCol, l0time, fGeometry->GetGeometryPtr());
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
  list.push_back(AliHLTEMCALDefinitions::fgkDigitDataType);
}

AliHLTComponentDataType AliHLTEMCALTriggerQAComponent::GetOutputDataType()
{
  return kAliHLTDataTypeHistogram | kAliHLTDataOriginEMCAL;
}

void AliHLTEMCALTriggerQAComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
{
  //see header file for documentation
  //constBase = 250000000;
  constBase = 250000;
  // to be reviewed later
  inputMultiplier = 0;
}

AliHLTComponent* AliHLTEMCALTriggerQAComponent::Spawn()
{
  return new AliHLTEMCALTriggerQAComponent();
}


int AliHLTEMCALTriggerQAComponent::DoInit(int argc, const char** argv)
{
  //Init the CTP data
  if (SetupCTPData() == -ENOMEM) {
    HLTError("could not SetupCTPData(); ENOMEM");
    return -ENOMEM;
  }
  InitialiseGeometry();

  Int_t debugLevel = 0;
  Bool_t enabledPatchType[3] = {kTRUE};

  Bool_t isPbPb = GetRunNo() > 244823 && GetRunNo() < 246995; // For the moment quick hack to distinguish PbPb from pp
  fBeamType = isPbPb ? kPbPb : kPP;

  for (int i = 0; i < argc; i++) {
    TString option(argv[i]);
    if (option == "-pp") fBeamType = kPP;
    if (option == "-PbPb") fBeamType = kPbPb;
    if (option == "-newTriggerBitConfig") fTriggerBitConfig = new AliEMCALTriggerBitConfigNew;
    if (option == "-oldTriggerBitConfig") fTriggerBitConfig = new AliEMCALTriggerBitConfigOld;
    if (option == "-noHistoReset") fHistoResetOnPush = kFALSE;
    if (option.BeginsWith("-debugLevel")) {
      option.Remove(0, 11);
      debugLevel = option.Atoi();
      if (debugLevel < 0) debugLevel = 0;
    }
    if (option.BeginsWith("-disableOffline")) {
      enabledPatchType[AliEMCALTriggerQA::kOfflinePatch] = kFALSE;
    }
    if (option.BeginsWith("-disableOnline")) {
      enabledPatchType[AliEMCALTriggerQA::kOnlinePatch] = kFALSE;
    }
    if (option.BeginsWith("-disableRecalc")) {
      enabledPatchType[AliEMCALTriggerQA::kRecalcPatch] = kFALSE;
    }
  }

  switch (fBeamType) {
  case kPP:
    fTriggerQAPtr = new AliEMCALTriggerOnlineQAPP("PPTriggerQA");
    break;
  case kPbPb:
    fTriggerQAPtr = new AliEMCALTriggerOnlineQAPbPb("PbPbTriggerQA");
    break;
  }

  fTriggerQAPtr->SetEMCALGeometry(fGeometry->GetGeometryPtr());
  fTriggerQAPtr->Init();

  fTriggerQAPtr->SetDebugLevel(debugLevel);

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

void AliHLTEMCALTriggerQAComponent::SetPbPb2015TriggerClasses()
{
  fFilterTrgClass = "CINT7-B-NOPF-CENT"
      "CINT7EG1-B-NOPF-CENTNOPMD CINT7EG2-B-NOPF-CENTNOPMD CINT7EJ1-B-NOPF-CENTNOPMD CINT7EJ2-B-NOPF-CENTNOPMD"
      "CINT7DG1-B-NOPF-CENTNOPMD CINT7DG2-B-NOPF-CENTNOPMD CINT7DJ1-B-NOPF-CENTNOPMD CINT7DJ2-B-NOPF-CENTNOPMD";
}

void AliHLTEMCALTriggerQAComponent::SetPP2016TriggerClasses()
{
  fFilterTrgClass = "";
}
