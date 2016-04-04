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
#include "AliEMCALTriggerOnlineQAPbPb.h"
//#include "AliEMCALTriggerOnlineQAPP.h"
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
  fFilterTrgClass(""),
  fLocalEventCount(0),
  fGeometry(NULL),
  fTriggerQAPtr(NULL),
  fTrgClassHistos(0),
  fFiredTriggerClasses()
{
  SetPbPb2015TriggerClasses();

  memset(fEMCalBkg, 0, sizeof(Double_t)*3);
  memset(fDCalBkg, 0, sizeof(Double_t)*3);
}

AliHLTEMCALTriggerQAComponent::~AliHLTEMCALTriggerQAComponent()
{
  if (fTriggerQAPtr) delete fTriggerQAPtr;
  if (fGeometry) delete fGeometry;
}

void AliHLTEMCALTriggerQAComponent::SetPbPb2015TriggerClasses()
{
  fFilterTrgClass = "CINT7-B-NOPF-CENT"
      "CINT7EG1-B-NOPF-CENTNOPMD CINT7EG2-B-NOPF-CENTNOPMD CINT7EJ1-B-NOPF-CENTNOPMD CINT7EJ2-B-NOPF-CENTNOPMD"
      "CINT7DG1-B-NOPF-CENTNOPMD CINT7DG2-B-NOPF-CENTNOPMD CINT7DJ1-B-NOPF-CENTNOPMD CINT7DJ2-B-NOPF-CENTNOPMD";
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
      if (!fFilterTrgClass.Contains(trgClass)) continue;
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
  //patch in order to skip calib events
  if (!IsDataEvent()) return 0;

  if (!blocks) return 0;

  //if (!RetrieveFiredTriggerClasses()) return 0;
  if (!RetrieveFiredTriggerClasses()) {
    HLTDebug("No trigger classes received\n");
  }

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

  PushHistograms(fTriggerQAPtr->GetListOfHistograms());
  PushHistograms(fTrgClassHistos);

  return 0;
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

void AliHLTEMCALTriggerQAComponent::CreateTrgClassHistograms(const TString& trgClass, TString patchTag)
{
  TString hname;

  TString patchTypes[3] = {"Online", "Recalc", "Offline"};

  Int_t max = 2000;
  if (patchTag.Contains("JE")) max = 6000;

  for (Int_t i = 0; i < 3; i++) {
    hname = Form("EMCTRQA_%s_%s_%s_PatchAmp", trgClass.Data(), patchTag.Data(), patchTypes[i].Data());
    TH1* hTrgPatchAmp = new TH1F(hname, hname, max/4, 0, max);
    fTrgClassHistos->Add(hTrgPatchAmp);

    hname = Form("EMCTRQA_%s_%s_%s_PatchAmpSubtracted", trgClass.Data(), patchTag.Data(), patchTypes[i].Data());
    TH1* hTrgPatchAmpSub = new TH1F(hname, hname, max/4, -max/2, max/2);
    fTrgClassHistos->Add(hTrgPatchAmpSub);
  }
}

void AliHLTEMCALTriggerQAComponent::FillTrgClassHistograms(const AliEMCALTriggerPatchInfo& patch, const TString& trgClass, TString patchTag)
{
  TString hname;
  TH1* hist = 0;

  // Should find a more elegant way... for the moment it works!
  hname = Form("EMCTRQA_%s_%s_Offline_PatchAmp", trgClass.Data(), patchTag.Data());
  hist = static_cast<TH1*>(fTrgClassHistos->FindObject(hname));
  if (!hist) CreateTrgClassHistograms(trgClass, patchTag);

  Double_t bkg[3] = {0};
  if (patch.IsDCalPHOS()) {
    bkg[0] = fEMCalBkg[0] * patch.GetPatchSize()*patch.GetPatchSize();
    bkg[1] = fEMCalBkg[1] * patch.GetPatchSize()*patch.GetPatchSize();
    bkg[2] = fEMCalBkg[2] * patch.GetPatchSize()*patch.GetPatchSize();
  }
  else {
    bkg[0] = fDCalBkg[0] * patch.GetPatchSize()*patch.GetPatchSize();
    bkg[1] = fDCalBkg[1] * patch.GetPatchSize()*patch.GetPatchSize();
    bkg[2] = fDCalBkg[2] * patch.GetPatchSize()*patch.GetPatchSize();
  }

  if (patch.IsOnline()) {
    hname = Form("EMCTRQA_%s_%s_Online_PatchAmp", trgClass.Data(), patchTag.Data());
    hist = static_cast<TH1*>(fTrgClassHistos->FindObject(hname));
    hist->Fill(patch.GetADCAmp());

    hname = Form("EMCTRQA_%s_%s_Online_PatchAmpSubtracted", trgClass.Data(), patchTag.Data());
    hist = static_cast<TH1*>(fTrgClassHistos->FindObject(hname));

    hist->Fill(patch.GetADCAmp() - bkg[0]);
  }

  if (patch.IsRecalc()) {
    hname = Form("EMCTRQA_%s_%s_Recalc_PatchAmp", trgClass.Data(), patchTag.Data());
    hist = static_cast<TH1*>(fTrgClassHistos->FindObject(hname));
    hist->Fill(patch.GetADCAmp());

    hname = Form("EMCTRQA_%s_%s_Recalc_PatchAmpSubtracted", trgClass.Data(), patchTag.Data());
    hist = static_cast<TH1*>(fTrgClassHistos->FindObject(hname));

    hist->Fill(patch.GetADCAmp() - bkg[1]);
  }

  if (patch.IsOfflineSimple()) {
    hname = Form("EMCTRQA_%s_%s_Offline_PatchAmp", trgClass.Data(), patchTag.Data());
    hist = static_cast<TH1*>(fTrgClassHistos->FindObject(hname));
    hist->Fill(patch.GetADCOfflineAmp());

    hname = Form("EMCTRQA_%s_%s_Offline_PatchAmpSubtracted", trgClass.Data(), patchTag.Data());
    hist = static_cast<TH1*>(fTrgClassHistos->FindObject(hname));

    hist->Fill(patch.GetADCOfflineAmp() - bkg[2]);
  }
}

void AliHLTEMCALTriggerQAComponent::ProcessTriggerClasses(const AliEMCALTriggerPatchInfo& patch)
{
  TString histname;

  // TO DO: a lot of hard-coded stuff... only relevant for PbPb 2015, will need to change it after heavy ion data taking

  for (unsigned int i = 0; i < fFiredTriggerClasses.size(); i++) {
    TString trgClass = fFiredTriggerClasses[i];
    if (patch.IsJetHigh() || patch.IsJetHighSimple() || patch.IsJetHighRecalc()) {
      if ((trgClass == "CINT7-B-NOPF-CENT") ||
          (trgClass.Contains("DJ1") && patch.IsDCalPHOS()) ||
          (trgClass.Contains("EJ1") && patch.IsEMCal())) {
        FillTrgClassHistograms(patch, trgClass, "JEH");
      }
    }

    if (patch.IsJetLow() || patch.IsJetLowSimple() || patch.IsJetLowRecalc()) {
      if ((trgClass == "CINT7-B-NOPF-CENT") ||
          (trgClass.Contains("DJ1") && patch.IsDCalPHOS()) ||
          (trgClass.Contains("EJ1") && patch.IsEMCal())) {
        FillTrgClassHistograms(patch, trgClass, "JEL");
      }
    }

    if (patch.IsGammaHigh() || patch.IsGammaHighSimple() || patch.IsGammaHighRecalc()) {
      if ((trgClass == "CINT7-B-NOPF-CENT") ||
          (trgClass.Contains("DG1") && patch.IsDCalPHOS()) ||
          (trgClass.Contains("EG1") && patch.IsEMCal())) {
        FillTrgClassHistograms(patch, trgClass, "GAH");
      }
    }

    if (patch.IsGammaLow() || patch.IsGammaLowSimple() || patch.IsGammaLowRecalc()) {
      if ((trgClass == "CINT7-B-NOPF-CENT") ||
          (trgClass.Contains("DG1") && patch.IsDCalPHOS()) ||
          (trgClass.Contains("EG1") && patch.IsEMCal())) {
        FillTrgClassHistograms(patch, trgClass, "GAL");
      }
    }
  }
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
    fTriggerQAPtr->ProcessBkgPatch(&patch);
  }

  fTriggerQAPtr->ComputeBackground();

  fTriggerQAPtr->GetEMCalBkg(fEMCalBkg);
  fTriggerQAPtr->GetEMCalBkg(fDCalBkg);

  for(UInt_t ipatch = 0; ipatch < nPatches; ipatch++) {
    HLTPatch2Patch(hltpatchPtr[ipatch], patch);
    fTriggerQAPtr->ProcessPatch(&patch);
    ProcessTriggerClasses(patch);
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
  //Init the CTP data
  if (SetupCTPData() == -ENOMEM) {
    HLTError("could not SetupCTPData(); ENOMEM");
    return -ENOMEM;
  }
  InitialiseGeometry();

  Int_t debugLevel = 0;
  Bool_t enabledPatchType[3] = {kTRUE};

  enum BeamType { kPP, kPbPb } beam = kPP;

  for (int i = 0; i < argc; i++) {
    TString option(argv[i]);
    if (option == "-pp") beam = kPP;
    if (option == "-PbPb") beam = kPbPb;
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

  switch (beam) {
  case kPP:
    //fTriggerQAPtr = new AliEMCALTriggerOnlineQAPP("PPTriggerQA");
    break;
  case kPbPb:
    fTriggerQAPtr = new AliEMCALTriggerOnlineQAPbPb("PbPbTriggerQA");
    break;
  }

  fTriggerQAPtr->Init();

  fTriggerQAPtr->SetDebugLevel(debugLevel);
  fTriggerQAPtr->EnablePatchType(AliEMCALTriggerQA::kOfflinePatch, enabledPatchType[AliEMCALTriggerQA::kOfflinePatch]);
  fTriggerQAPtr->EnablePatchType(AliEMCALTriggerQA::kOnlinePatch, enabledPatchType[AliEMCALTriggerQA::kOnlinePatch]);
  fTriggerQAPtr->EnablePatchType(AliEMCALTriggerQA::kRecalcPatch, enabledPatchType[AliEMCALTriggerQA::kRecalcPatch]);

  // if not specified use new trigger bit config
  if (!fTriggerBitConfig) fTriggerBitConfig = new AliEMCALTriggerBitConfigNew;

  fTrgClassHistos = new THashList();
  fTrgClassHistos->SetName("trgClassHistos");

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
