/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
#include <algorithm>
#include <iostream>
#include <memory>

#include <TClonesArray.h>
#include <TGrid.h>
#include <THashList.h>
#include <THistManager.h>
#include <TMath.h>
#include <TParameter.h>

#include "AliAnalysisUtils.h"
#include "AliOADBContainer.h"
#include "AliAnalysisTaskEmcalPatchMasking.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALTriggerMapping.h"
#include "AliEMCALTriggerPatchInfo.h"
#include "AliEMCALTriggerPatchADCInfoAP.h"
#include "AliESDEvent.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliVCaloTrigger.h"
#include "AliVEvent.h"

ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalPatchMasking)

using namespace PWGJE::EMCALJetTasks;

AliAnalysisTaskEmcalPatchMasking::AliAnalysisTaskEmcalPatchMasking() :
  AliAnalysisTaskEmcal(),
  fHistos(nullptr),
  fL1ADC(),
  fTriggerBits(0),
  fTriggerPattern(""),
  fNameMaskedFastorOADB(""),
  fMaskedFastorOADB(nullptr),
  fListMaskedFastors()
{
  SetNeedEmcalGeom(true);
  SetCaloTriggerPatchInfoName("EmcalTriggers");
}

AliAnalysisTaskEmcalPatchMasking::AliAnalysisTaskEmcalPatchMasking(const char *name) :
  AliAnalysisTaskEmcal(name, true),
  fHistos(nullptr),
  fL1ADC(),
  fTriggerBits(0),
  fTriggerPattern(""),
  fNameMaskedFastorOADB(""),
  fMaskedFastorOADB(nullptr),
  fListMaskedFastors()
{
  SetNeedEmcalGeom(true);
  SetCaloTriggerPatchInfoName("EmcalTriggers");
}


AliAnalysisTaskEmcalPatchMasking::~AliAnalysisTaskEmcalPatchMasking() {
  if(fMaskedFastorOADB) delete fMaskedFastorOADB;
  if(fHistos) delete fHistos;
}

void AliAnalysisTaskEmcalPatchMasking::UserCreateOutputObjects(){
  AliAnalysisTaskEmcal::UserCreateOutputObjects();

  fAliAnalysisUtils = new AliAnalysisUtils;

  std::vector<TString> patchtypes = {"EG1", "EG2", "EJ1", "EJ2", "REGA", "REJE"};
  fHistos = new THistManager(Form("Task%s", GetName()));
  for(const auto &pt : patchtypes){
    fHistos->CreateTH1(Form("hFracMaxFastorTotalAll%s", pt.Data()), Form("Contribution of max fastor to patch ADC for all %s patches", pt.Data()), 100, 0., 1.);
    fHistos->CreateTH1(Form("hFracMaxFastorMeanAll%s", pt.Data()), Form("Ratio max ADC to the mean for all %s patches", pt.Data()), 500, 0., 50.);
    fHistos->CreateTH2(Form("hPatchADCvsADCgoodAll%s", pt.Data()), Form("Patch ADC vs patch good ADC for all %s patches", pt.Data()), 500, 0., 2500., 500, 0., 2500.);
    fHistos->CreateTH1(Form("hFracMaxMaskedFastorAll%s", pt.Data()), Form("Contribution of the max masked fastor to the patch ADC for all %s patches", pt.Data()), 100, 0., 1.);
    fHistos->CreateTH1(Form("hFracMaxGoodFastorAll%s", pt.Data()), Form("Contribution of the max good fastor to the patch ADC for all %s patches", pt.Data()), 100, 0., 1.);
    fHistos->CreateTH1(Form("hFracGoodADCAll%s", pt.Data()), Form("Contribution of the good adc sum to the patch ADC for all %s patches", pt.Data()), 100, 0., 1.);
    fHistos->CreateTH1(Form("hFracMaskedADCAll%s", pt.Data()), Form("Contribution of the masked adc sum to the patch ADC for all %s patches", pt.Data()), 100, 0., 1.);
    fHistos->CreateTH2(Form("hFracMaxMaskedvsNmaskedcontribAll%s", pt.Data()), Form("Number of non-zero masked fastors vs. fraction of max. masked for all %s patches", pt.Data()), 1000, 0., 1000, 100., 0., 1.);
    fHistos->CreateTH1(Form("hFracMaxFastorTotalBad%s", pt.Data()), Form("Contribution of max fastor to patch ADC for bad %s patches", pt.Data()), 100, 0., 1.);
    fHistos->CreateTH1(Form("hFracMaxFastorMeanBad%s", pt.Data()), Form("Ratio max ADC to the mean for bad %s patches", pt.Data()), 500, 0., 50.);
    fHistos->CreateTH2(Form("hPatchADCvsADCgoodBad%s", pt.Data()), Form("Patch ADC vs patch good ADC for bad %s patches", pt.Data()), 500, 0., 2500., 500, 0., 2500.);
    fHistos->CreateTH1(Form("hFracMaxMaskedFastorBad%s", pt.Data()), Form("Contribution of the max masked fastor to the patch ADC for bad %s patches", pt.Data()), 100, 0., 1.);
    fHistos->CreateTH1(Form("hFracMaxGoodFastorBad%s", pt.Data()), Form("Contribution of the max good fastor to the patch ADC for bad %s patches", pt.Data()), 100, 0., 1.);
    fHistos->CreateTH1(Form("hFracGoodADCBad%s", pt.Data()), Form("Contribution of the good adc sum to the patch ADC for bad %s patches", pt.Data()), 100, 0., 1.);
    fHistos->CreateTH1(Form("hFracMaskedADCBad%s", pt.Data()), Form("Contribution of the masked adc sum to the patch ADC bad all %s patches", pt.Data()), 100, 0., 1.);
    fHistos->CreateTH2(Form("hFracMaxMaskedvsNmaskedcontribBad%s", pt.Data()), Form("Number of non-zero masked fastors vs. fraction of max. masked for bad %s patches", pt.Data()), 1000, 0., 1000, 100., 0., 1.);
    fHistos->CreateTH1(Form("hFracMaxMaskedSumMaskedBad%s", pt.Data()), Form("Contibution of the max. masked fastor to the ADC sum of masked fastors for bad %s patches", pt.Data()), 100, 0., 1.);
    fHistos->CreateTH1(Form("hFracMaxFastorTotalGood%s", pt.Data()), Form("Contribution of max fastor to patch ADC for good %s patches", pt.Data()), 100, 0., 1.);
    fHistos->CreateTH1(Form("hFracMaxFastorMeanGood%s", pt.Data()), Form("Ratio max ADC to the mean for good %s patches", pt.Data()), 500, 0., 10.);
  }

  std::vector<TString> maxtypes = {"Online", "Recalc"};
  for(const auto &mt : maxtypes){
    fHistos->CreateTH1(Form("hFracMaxFastorTotalAllMax%s", mt.Data()), Form("Contribution of max fastor to patch ADC for all max %s patches", mt.Data()), 100, 0., 1.);
    fHistos->CreateTH1(Form("hFracMaxFastorMeanAllMax%s", mt.Data()), Form("Ratio max ADC to the mean for all max %s patches", mt.Data()), 500, 0., 50.);
    fHistos->CreateTH2(Form("hPatchADCvsADCgoodAllMax%s", mt.Data()), Form("Patch ADC vs patch good ADC for all max %s patches", mt.Data()), 500, 0., 2500., 500, 0., 2500.);
    fHistos->CreateTH1(Form("hFracMaxMaskedFastorAllMax%s", mt.Data()), Form("Contribution of the max masked fastor to the patch ADC for all max %s patches", mt.Data()), 100, 0., 1.);
    fHistos->CreateTH1(Form("hFracMaxGoodFastorAllMax%s", mt.Data()), Form("Contribution of the max good fastor to the patch ADC for all max %s patches", mt.Data()), 100, 0., 1.);
    fHistos->CreateTH1(Form("hFracGoodADCAllMax%s", mt.Data()), Form("Contribution of the good adc sum to the patch ADC for all max %s patches", mt.Data()), 100, 0., 1.);
    fHistos->CreateTH1(Form("hFracMaskedADCAllMax%s", mt.Data()), Form("Contribution of the masked adc sum to the patch ADC for all max %s patches", mt.Data()), 100, 0., 1.);
    fHistos->CreateTH2(Form("hFracMaxMaskedvsNmaskedcontribAllMax%s", mt.Data()), Form("Number of non-zero masked fastors vs. fraction of max. masked for all max %s patches", mt.Data()), 1000, 0., 1000, 100., 0., 1.);
    fHistos->CreateTH1(Form("hFracMaxFastorTotalBadMax%s", mt.Data()), Form("Contribution of max fastor to patch ADC for bad max %s patches", mt.Data()), 100, 0., 1.);
    fHistos->CreateTH1(Form("hFracMaxFastorMeanBadMax%s", mt.Data()), Form("Ratio max ADC to the mean for bad max %s patches", mt.Data()), 500, 0., 50.);
    fHistos->CreateTH2(Form("hPatchADCvsADCgoodBadMax%s", mt.Data()), Form("Patch ADC vs patch good ADC for bad max %s patches", mt.Data()), 500, 0., 2500., 500, 0., 2500.);
    fHistos->CreateTH1(Form("hFracMaxMaskedFastorBadMax%s", mt.Data()), Form("Contribution of the max masked fastor to the patch ADC for bad max %s patches", mt.Data()), 100, 0., 1.);
    fHistos->CreateTH1(Form("hFracMaxGoodFastorBadMax%s", mt.Data()), Form("Contribution of the max good fastor to the patch ADC for bad max %s patches", mt.Data()), 100, 0., 1.);
    fHistos->CreateTH1(Form("hFracGoodADCBadMax%s", mt.Data()), Form("Contribution of the good adc sum to the patch ADC for bad max %s patches", mt.Data()), 100, 0., 1.);
    fHistos->CreateTH1(Form("hFracMaskedADCBadMax%s", mt.Data()), Form("Contribution of the masked adc sum to the patch ADC bad all max %s patches", mt.Data()), 100, 0., 1.);
    fHistos->CreateTH2(Form("hFracMaxMaskedvsNmaskedcontribBadMax%s", mt.Data()), Form("Number of non-zero masked fastors vs. fraction of max. masked for bad max %s patches", mt.Data()), 1000, 0., 1000, 100., 0., 1.);
    fHistos->CreateTH1(Form("hFracMaxMaskedSumMaskedBadMax%s", mt.Data()), Form("Contibution of the max. masked fastor to the ADC sum of masked fastors for bad max %s patches", mt.Data()), 100, 0., 1.);
    fHistos->CreateTH1(Form("hFracMaxFastorTotalGoodMax%s", mt.Data()), Form("Contribution of max fastor to patch ADC for good max %s patches", mt.Data()), 100, 0., 1.);
    fHistos->CreateTH1(Form("hFracMaxFastorMeanGoodMax%s", mt.Data()), Form("Ratio max ADC to the mean for good max %s patches", mt.Data()), 500, 0., 50.);
  }

  for(auto h : *(fHistos->GetListOfHistograms())) fOutput->Add(h);

  PostData(1, fOutput);
}

bool AliAnalysisTaskEmcalPatchMasking::IsEventSelected(){
  if(fTriggerBits){
    if(!(fInputHandler->IsEventSelected() & fTriggerBits)) return false;
  }
  if(fTriggerPattern.Length()){
    if(!fInputEvent->GetFiredTriggerClasses().Contains(fTriggerPattern)) return false;
  }

  if(fAliAnalysisUtils){
    if(fInputEvent->IsA() == AliESDEvent::Class() && ! fAliAnalysisUtils->IsFirstEventInChunk(fInputEvent)) return false;
    if(fAliAnalysisUtils->IsPileUpEvent(fInputEvent)) return false;
    if(fAliAnalysisUtils->IsVertexSelected2013pA(fInputEvent)) return false;
  }

  AliDebugStream(1) << GetName() << ": Event is selected" << std::endl;

  return true;
}

bool AliAnalysisTaskEmcalPatchMasking::Run(){
  PrepareL1FastorADC();

  // Loop over trigger patches
  AliEMCALTriggerPatchInfo *recpatch(nullptr), *maxOnline(nullptr), *maxRecalc(nullptr);
  for(auto p : *(fTriggerPatchInfo)){
    recpatch = static_cast<AliEMCALTriggerPatchInfo *>(p);

    // select only online and recalc patches
    if(!(recpatch->IsGammaHigh() || recpatch->IsGammaLow() ||
        recpatch->IsJetHigh() || recpatch->IsJetLow() ||
        recpatch->IsGammaLowRecalc() || recpatch->IsJetLowRecalc())) continue;

    // Select only firing patches (patches - unmasked - above threshold for the given trigger)
    // find the highest energetic online and recalc patch at the same step
    if(fTriggerBits & AliVEvent::kEMCEGA){
      if(recpatch->IsJetLow() || recpatch->IsJetHigh() || recpatch->IsJetLowRecalc()) continue;
      if(fTriggerPattern.Contains("EG1")){
        if(!(recpatch->IsGammaHigh() || recpatch->GetADCAmp() > 140)) continue;
        if(recpatch->IsGammaHigh()){
          if(!maxOnline) maxOnline = recpatch;
          else {
            if(recpatch->GetADCAmp() > maxOnline->GetADCAmp()) maxOnline = recpatch;
          }
        }
        if(recpatch->IsGammaLowRecalc()){
          if(!maxRecalc) maxRecalc = recpatch;
          else {
            if(recpatch->GetADCAmp() > maxRecalc->GetADCAmp()) maxRecalc = recpatch;
          }
        }
      }
      if(fTriggerPattern.Contains("EG2")){
        if(!(recpatch->IsGammaLow() || recpatch->GetADCAmp() > 89)) continue;
        if(recpatch->IsGammaLow()){
          if(!maxOnline) maxOnline = recpatch;
          else {
            if(recpatch->GetADCAmp() > maxOnline->GetADCAmp()) maxOnline = recpatch;
          }
        }
        if(recpatch->IsGammaLowRecalc()){
          if(!maxRecalc) maxRecalc = recpatch;
          else {
            if(recpatch->GetADCAmp() > maxRecalc->GetADCAmp()) maxRecalc = recpatch;
          }
        }
      }
    }
    if(fTriggerBits & AliVEvent::kEMCEJE){
      if(recpatch->IsGammaLow() || recpatch->IsGammaHigh() || recpatch->IsGammaLowRecalc()) continue;
      if(fTriggerPattern.Contains("EJ1")){
        if(!(recpatch->IsJetHigh() || recpatch->GetADCAmp() > 260)) continue;
        if(recpatch->IsJetHigh()){
          if(!maxOnline) maxOnline = recpatch;
          else {
            if(recpatch->GetADCAmp() > maxOnline->GetADCAmp()) maxOnline = recpatch;
          }
        }
        if(recpatch->IsJetLowRecalc()){
          if(!maxRecalc) maxRecalc = recpatch;
          else {
            if(recpatch->GetADCAmp() > maxRecalc->GetADCAmp()) maxRecalc = recpatch;
          }
        }
      }
      if(fTriggerPattern.Contains("EJ2")){
        if(!(recpatch->IsJetLow() || recpatch->GetADCAmp() > 127)) continue;
        if(recpatch->IsJetLow()){
          if(!maxOnline) maxOnline = recpatch;
          else {
            if(recpatch->GetADCAmp() > maxOnline->GetADCAmp()) maxOnline = recpatch;
          }
        }
        if(recpatch->IsJetLowRecalc()){
          if(!maxRecalc) maxRecalc = recpatch;
          else {
            if(recpatch->GetADCAmp() > maxRecalc->GetADCAmp()) maxRecalc = recpatch;
          }
        }
      }
    }

    AliDebugStream(2) << GetName() << ": Found firing patch" << std::endl;
    ProcessPatch(*recpatch);
  }

  if(maxOnline) ProcessMaxPatch(*maxOnline, "Online");
  if(maxRecalc) ProcessMaxPatch(*maxRecalc, "Recalc");

  return true;
}

void AliAnalysisTaskEmcalPatchMasking::ExecOnce(){
  AliAnalysisTaskEmcal::ExecOnce();
  if(!fLocalInitialized) return;

  if(fNameMaskedFastorOADB.Length()){
    AliInfoStream() << "Initializing masked fastors from OADB container " << fNameMaskedFastorOADB.Data() << std::endl;
    if(fNameMaskedFastorOADB.Contains("alien://") && !gGrid) TGrid::Connect("alien://");
    fMaskedFastorOADB = new AliOADBContainer("AliEmcalMaskedFastors");
    fMaskedFastorOADB->InitFromFile(fNameMaskedFastorOADB.Data(), "AliEmcalMaskedFastors");
  }

  UChar_t nrow(64);
  if(fGeom->GetTriggerMappingVersion() == 2) nrow = 104;
  fL1ADC.Allocate(48, nrow);
}

void AliAnalysisTaskEmcalPatchMasking::RunChanged(int newrun){
  if(fMaskedFastorOADB){
    fListMaskedFastors.clear();
    TObjArray *badchannelmap = static_cast<TObjArray *>(fMaskedFastorOADB->GetObject(newrun));
    if(!badchannelmap || !badchannelmap->GetEntries()) return;
    for(TIter citer = TIter(badchannelmap).Begin(); citer != TIter::End(); ++citer){
      TParameter<int> *channelID = static_cast<TParameter<int> *>(*citer);
      AliDebugStream(3) << GetName() << ": Found masked fastor channel " << channelID->GetVal() << std::endl;
      fListMaskedFastors.push_back(channelID->GetVal());
    }
  }
}

void AliAnalysisTaskEmcalPatchMasking::PrepareL1FastorADC(){
  fL1ADC.Reset();
  AliVCaloTrigger *emctrigger = fInputEvent->GetCaloTrigger("EMCAL");
  emctrigger->Reset();

  Int_t globCol=-1, globRow=-1, adcAmp=-1;
  while(emctrigger->Next()){
    // get position in global 2x2 tower coordinates
    // A0 left bottom (0,0)
    emctrigger->GetPosition(globCol, globRow);
    emctrigger->GetL1TimeSum(adcAmp);
    AliDebugStream(1) << GetName() << "Fastor at (" << globCol << "," << globRow << ") with ADC " << adcAmp << std::endl;
    if (adcAmp < 0) adcAmp = 0;

    try {
      (fL1ADC)(globCol,globRow) = adcAmp;
    }
    catch (AliEMCALTriggerDataGrid<Int_t>::OutOfBoundsException &e) {
      std::string dirstring = e.GetDirection() == AliEMCALTriggerDataGrid<Int_t>::OutOfBoundsException::kColDir ? "Col" : "Row";
      AliErrorStream() << "Trigger maker task - filling trigger bit grid - index out-of-bounds in " << dirstring << ": " << e.GetIndex() << std::endl;
    }
  }
}

AliEMCALTriggerPatchADCInfoAP *AliAnalysisTaskEmcalPatchMasking::MakeFastorADCValuesForPatch(const AliEMCALTriggerPatchInfo &patch ) const {
  AliEMCALTriggerPatchADCInfoAP *adcpatch = new AliEMCALTriggerPatchADCInfoAP(patch.GetPatchSize());
  for(unsigned char icol = 0; icol < patch.GetPatchSize(); icol++){
    for(unsigned char irow = 0; irow < patch.GetPatchSize(); irow++){
      Int_t adc = 0;
      try{
        adc = fL1ADC(icol + patch.GetColStart(), irow + patch.GetRowStart());
      } catch (AliEMCALTriggerDataGrid<Int_t>::OutOfBoundsException &e){
        adc = 0;
      }
      adcpatch->SetADC(adc, icol, irow);
    }
  }
  return adcpatch;
}

void AliAnalysisTaskEmcalPatchMasking::ProcessPatch(const AliEMCALTriggerPatchInfo &patch){
  std::vector<TString> patchtypes;
  if(patch.IsGammaHigh()) patchtypes.push_back("EG1");
  if(patch.IsGammaLow()) patchtypes.push_back("EG2");
  if(patch.IsGammaLowRecalc()) patchtypes.push_back("REGA");
  if(patch.IsJetHigh()) patchtypes.push_back("EJ1");
  if(patch.IsJetLow()) patchtypes.push_back("EJ2");
  if(patch.IsJetLowRecalc()) patchtypes.push_back("REJE");

  std::unique_ptr<AliEMCALTriggerPatchADCInfoAP> patchADC(MakeFastorADCValuesForPatch(patch));

  // Find the max fastor ADC in the patch
  Int_t maxFastorADC = 0, maxFastorADCmasked = 0, maxFastorADCgood = 0, tmp = 0;
  Int_t nFastorNonZero = 0, nfastorMasked = 0, nMaskedFastorsNonZero = 0;
  Int_t sumADCgood = 0, sumADCbad = 0, absFastorIndex = 0;
  std::vector<Int_t> adcs;
  for(UChar_t icol = 0; icol < patchADC->GetPatchSize(); icol++){
    for(UChar_t irow = 0; irow < patchADC->GetPatchSize(); irow++){
      tmp = patchADC->GetADC(icol, irow);
      adcs.push_back(tmp);
      if(tmp > 0) nFastorNonZero++;
      if(tmp > maxFastorADC) maxFastorADC = tmp;

      // Get masked fastors in the patch
      fGeom->GetTriggerMapping()->GetAbsFastORIndexFromPositionInEMCAL(icol + patch.GetColStart(), irow + patch.GetRowStart(), absFastorIndex);
      if(std::find(fListMaskedFastors.begin(), fListMaskedFastors.end(), absFastorIndex) != fListMaskedFastors.end()){
        // Fastor is masked
        nfastorMasked++;
        if(tmp > 0) nMaskedFastorsNonZero++;
        sumADCbad += tmp;
        if(tmp > maxFastorADCmasked) maxFastorADCmasked = tmp;
      } else {
        // Fastor is not masked
        sumADCgood += tmp;
        if(tmp > maxFastorADCgood) maxFastorADCgood = tmp;
      }
    }
  }
  AliDebugStream(2) << GetName() << ": Max Fastor: " <<  maxFastorADC << ", Sum: " << patch.GetADCAmp() << std::endl;
  AliDebugStream(2) << GetName() << ": Mean ADC: " << TMath::Mean(adcs.begin(), adcs.end()) << std::endl;
  AliDebugStream(2) << GetName() << ": Sum good " << sumADCgood << ", bad " << sumADCbad << std::endl;

  Float_t fracMaxFastorPatchADC = static_cast<Float_t>(maxFastorADC)/static_cast<Float_t>(patch.GetADCAmp()),
          fracMaxFastorMaskedPatchADC = static_cast<Float_t>(maxFastorADCmasked)/static_cast<Float_t>(patch.GetADCAmp()),
          fracMaxFastorGoodPatchADC = static_cast<Float_t>(maxFastorADCgood)/static_cast<Float_t>(patch.GetADCAmp()),
          fracMaxFastorMeanADC = static_cast<Float_t>(maxFastorADC)/TMath::Mean(adcs.begin(), adcs.end()),
          fracGoodADC = static_cast<Float_t>(sumADCgood)/static_cast<Float_t>(patch.GetADCAmp()),
          fracMaskedADC = static_cast<Float_t>(sumADCbad)/static_cast<Float_t>(patch.GetADCAmp());

  // Fill histograms
  for(const auto &pt : patchtypes){
    fHistos->FillTH1(Form("hFracMaxFastorTotalAll%s", pt.Data()), fracMaxFastorPatchADC);
    fHistos->FillTH1(Form("hFracMaxFastorMeanAll%s", pt.Data()), fracMaxFastorMeanADC);
    fHistos->FillTH2(Form("hPatchADCvsADCgoodAll%s", pt.Data()), patch.GetADCAmp(), sumADCgood);
    fHistos->FillTH1(Form("hFracMaxMaskedFastorAll%s", pt.Data()), fracMaxFastorMaskedPatchADC);
    fHistos->FillTH1(Form("hFracMaxGoodFastorAll%s", pt.Data()), fracMaxFastorGoodPatchADC);
    fHistos->FillTH1(Form("hFracGoodADCAll%s", pt.Data()), fracGoodADC);
    fHistos->FillTH1(Form("hFracMaskedADCAll%s", pt.Data()), fracMaskedADC);
    fHistos->FillTH2(Form("hFracMaxMaskedvsNmaskedcontribAll%s", pt.Data()), nMaskedFastorsNonZero, fracMaxFastorMaskedPatchADC);

    if(nfastorMasked){
      Float_t fracMaxMaskedSumMasked = 0;
      if(sumADCbad) fracMaxMaskedSumMasked = static_cast<Float_t>(maxFastorADCmasked)/static_cast<Float_t>(sumADCbad);
      fHistos->FillTH1(Form("hFracMaxFastorTotalBad%s", pt.Data()), fracMaxFastorPatchADC);
      fHistos->FillTH1(Form("hFracMaxFastorMeanBad%s", pt.Data()), fracMaxFastorMeanADC);
      fHistos->FillTH2(Form("hPatchADCvsADCgoodBad%s", pt.Data()), patch.GetADCAmp(), sumADCgood);
      fHistos->FillTH1(Form("hFracMaxMaskedFastorBad%s", pt.Data()), fracMaxFastorMaskedPatchADC);
      fHistos->FillTH1(Form("hFracMaxGoodFastorBad%s", pt.Data()), fracMaxFastorGoodPatchADC);
      fHistos->FillTH1(Form("hFracGoodADCBad%s", pt.Data()), fracGoodADC);
      fHistos->FillTH1(Form("hFracMaskedADCBad%s", pt.Data()), fracMaskedADC);
      fHistos->FillTH2(Form("hFracMaxMaskedvsNmaskedcontribBad%s", pt.Data()), nMaskedFastorsNonZero, fracMaxFastorMaskedPatchADC);
      fHistos->FillTH1(Form("hFracMaxMaskedSumMaskedBad%s", pt.Data()), fracMaxMaskedSumMasked);
    } else {
      // For jet patches can almost never happen
      fHistos->FillTH1(Form("hFracMaxFastorTotalGood%s", pt.Data()), fracMaxFastorPatchADC);
      fHistos->FillTH1(Form("hFracMaxFastorMeanGood%s", pt.Data()), fracMaxFastorMeanADC);
    }
  }
}

void AliAnalysisTaskEmcalPatchMasking::ProcessMaxPatch(const AliEMCALTriggerPatchInfo &patch, const TString &patchtype){
  std::unique_ptr<AliEMCALTriggerPatchADCInfoAP> patchADC(MakeFastorADCValuesForPatch(patch));

  // Find the max fastor ADC in the patch
  UShort_t maxFastorADC = 0, maxFastorADCmasked = 0, maxFastorADCgood = 0, tmp = 0;
  UShort_t nFastorNonZero = 0, nfastorMasked = 0, nMaskedFastorsNonZero = 0;
  Int_t sumADCgood = 0, sumADCbad = 0, absFastorIndex = 0;
  std::vector<UShort_t> adcs;
  for(UChar_t icol = 0; icol < patchADC->GetPatchSize(); icol++){
    for(UChar_t irow = 0; irow < patchADC->GetPatchSize(); irow++){
      tmp = patchADC->GetADC(icol, irow);
      adcs.push_back(tmp);
      if(tmp > 0) nFastorNonZero++;
      if(tmp > maxFastorADC) maxFastorADC = tmp;

      // Get masked fastors in the patch
      fGeom->GetTriggerMapping()->GetAbsFastORIndexFromPositionInEMCAL(icol + patch.GetColStart(), irow + patch.GetRowStart(), absFastorIndex);
      if(std::find(fListMaskedFastors.begin(), fListMaskedFastors.end(), absFastorIndex) != fListMaskedFastors.end()){
        // Fastor is masked
        nfastorMasked++;
        if(tmp > 0) nMaskedFastorsNonZero++;
        sumADCbad += tmp;
        if(tmp > maxFastorADCmasked) maxFastorADCmasked = tmp;
      } else {
        // Fastor is not masked
        sumADCgood += tmp;
        if(tmp > maxFastorADCgood) maxFastorADCgood = tmp;
      }
    }
  }
  Float_t fracMaxFastorPatchADC = static_cast<Float_t>(maxFastorADC)/static_cast<Float_t>(patch.GetADCAmp()),
          fracMaxFastorMaskedPatchADC = static_cast<Float_t>(maxFastorADCmasked)/static_cast<Float_t>(patch.GetADCAmp()),
          fracMaxFastorGoodPatchADC = static_cast<Float_t>(maxFastorADCgood)/static_cast<Float_t>(patch.GetADCAmp()),
          fracMaxFastorMeanADC = static_cast<Float_t>(maxFastorADC)/TMath::Mean(adcs.begin(), adcs.end()),
          fracGoodADC = static_cast<Float_t>(sumADCgood)/static_cast<Float_t>(patch.GetADCAmp()),
          fracMaskedADC = static_cast<Float_t>(sumADCbad)/static_cast<Float_t>(patch.GetADCAmp());

  fHistos->FillTH1(Form("hFracMaxFastorTotalAllMax%s", patchtype.Data()), fracMaxFastorPatchADC);
  fHistos->FillTH1(Form("hFracMaxFastorMeanAllMax%s", patchtype.Data()), fracMaxFastorMeanADC);
  fHistos->FillTH2(Form("hPatchADCvsADCgoodAllMax%s", patchtype.Data()), patch.GetADCAmp(), sumADCgood);
  fHistos->FillTH1(Form("hFracMaxMaskedFastorAllMax%s", patchtype.Data()), fracMaxFastorMaskedPatchADC);
  fHistos->FillTH1(Form("hFracMaxGoodFastorAllMax%s", patchtype.Data()), fracMaxFastorGoodPatchADC);
  fHistos->FillTH1(Form("hFracGoodADCAllMax%s", patchtype.Data()), fracGoodADC);
  fHistos->FillTH1(Form("hFracMaskedADCAllMax%s", patchtype.Data()), fracMaskedADC);
  fHistos->FillTH2(Form("hFracMaxMaskedvsNmaskedcontribAllMax%s", patchtype.Data()), nMaskedFastorsNonZero, fracMaxFastorMaskedPatchADC);

  if(nfastorMasked){
    Float_t fracMaxMaskedSumMasked = 0;
    if(sumADCbad) fracMaxMaskedSumMasked = static_cast<Float_t>(maxFastorADCmasked)/static_cast<Float_t>(sumADCbad);
    fHistos->FillTH1(Form("hFracMaxFastorTotalBadMax%s", patchtype.Data()), fracMaxFastorPatchADC);
    fHistos->FillTH1(Form("hFracMaxFastorMeanBadMax%s", patchtype.Data()), fracMaxFastorMeanADC);
    fHistos->FillTH2(Form("hPatchADCvsADCgoodBadMax%s", patchtype.Data()), patch.GetADCAmp(), sumADCgood);
    fHistos->FillTH1(Form("hFracMaxMaskedFastorBadMax%s", patchtype.Data()), fracMaxFastorMaskedPatchADC);
    fHistos->FillTH1(Form("hFracMaxGoodFastorBadMax%s", patchtype.Data()), fracMaxFastorGoodPatchADC);
    fHistos->FillTH1(Form("hFracGoodADCBadMax%s", patchtype.Data()), fracGoodADC);
    fHistos->FillTH1(Form("hFracMaskedADCBadMax%s", patchtype.Data()), fracMaskedADC);
    fHistos->FillTH2(Form("hFracMaxMaskedvsNmaskedcontribBadMax%s", patchtype.Data()), nMaskedFastorsNonZero, fracMaxFastorMaskedPatchADC);
    fHistos->FillTH1(Form("hFracMaxMaskedSumMaskedBadMax%s", patchtype.Data()), fracMaxMaskedSumMasked);
  } else {
    // For jet patches can almost never happen
    fHistos->FillTH1(Form("hFracMaxFastorTotalGoodMax%s", patchtype.Data()), fracMaxFastorPatchADC);
    fHistos->FillTH1(Form("hFracMaxFastorMeanGoodMax%s", patchtype.Data()), fracMaxFastorMeanADC);
  }

}
