/**************************************************************************
 * Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
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
#include <iostream>
#include <map>
#include <vector>

#include <TClonesArray.h>
#include <TGrid.h>
#include <THistManager.h>
#include <THashList.h>
#include <TLinearBinning.h>
#include <TObjArray.h>
#include <TParameter.h>

#include "AliAnalysisUtils.h"
#include "AliESDEvent.h"
#include "AliEMCALTriggerPatchInfo.h"
#include "AliEmcalTriggerOfflineSelection.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliMultSelection.h"
#include "AliMultEstimator.h"
#include "AliOADBContainer.h"

#include "AliAnalysisTaskEmcalPatchesRef.h"

/// \cond CLASSIMP
ClassImp(EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalPatchesRef)
/// \endcond

namespace EMCalTriggerPtAnalysis {

/**
 * Dummy (I/O) onstructor
 */
AliAnalysisTaskEmcalPatchesRef::AliAnalysisTaskEmcalPatchesRef() :
    AliAnalysisTaskSE(),
    fAnalysisUtil(nullptr),
    fTriggerSelection(nullptr),
    fHistos(nullptr),
    fTriggerPatches(nullptr),
    fRequestAnalysisUtil(kTRUE),
    fTriggerStringFromPatches(kFALSE),
    fCentralityRange(-999., 999.),
    fVertexRange(-999., 999.),
    fRequestCentrality(false),
    fNameDownscaleOADB(""),
    fDownscaleOADB(nullptr),
    fDownscaleFactors(nullptr),
    fCurrentRun(-1),
    fInitialized(false)
{
}

/**
 * Named constructor
 * @param name Name of the task
 */
AliAnalysisTaskEmcalPatchesRef::AliAnalysisTaskEmcalPatchesRef(const char *name):
    AliAnalysisTaskSE(name),
    fAnalysisUtil(nullptr),
    fTriggerSelection(nullptr),
    fHistos(nullptr),
    fTriggerPatches(nullptr),
    fRequestAnalysisUtil(kTRUE),
    fTriggerStringFromPatches(kFALSE),
    fCentralityRange(-999., 999.),
    fVertexRange(-999., 999.),
    fRequestCentrality(false),
    fNameDownscaleOADB(""),
    fDownscaleOADB(nullptr),
    fDownscaleFactors(nullptr),
    fCurrentRun(-1),
    fInitialized(false)
{
  DefineOutput(1, TList::Class());
}

/**
 * Destructor
 */
AliAnalysisTaskEmcalPatchesRef::~AliAnalysisTaskEmcalPatchesRef() {
}

/**
 * Creating output histograms:
 * + Patch (calibrated) energy spectrum - separated by patch type - for different trigger classes
 * + Patch eta-phi map - separated by patch type - for different trigger classes and different min. energies
 */
void AliAnalysisTaskEmcalPatchesRef::UserCreateOutputObjects(){
  AliInfoStream() <<  "Creating histograms for task " << GetName() << std::endl;
  fAnalysisUtil = new AliAnalysisUtils;

  EnergyBinning energybinning;
  TLinearBinning etabinning(100, -0.7, 0.7);
  fHistos = new THistManager("Ref");
  TString triggers[21] = {"MB", "EMC7", "DMC7",
      "EJ1", "EJ2", "EG1", "EG2", "DJ1", "DJ2", "DG1", "DG2",
      "EMC7excl", "DMC7excl", "EG2excl", "EJ2excl", "DG2excl", "DJ2excl",
      "EG1excl", "EJ1excl", "DG1excl", "DJ1excl"
  };
  TString patchtype[10] = {"EG1", "EG2", "EJ1", "EJ2", "EMC7", "DG1", "DG2", "DJ1", "DJ2", "DMC7"};
  Double_t encuts[5] = {1., 2., 5., 10., 20.};
  for(TString *trg = triggers; trg < triggers+18; trg++){
    fHistos->CreateTH1(Form("hEventCount%s", trg->Data()), Form("Event count for trigger class %s", trg->Data()), 1, 0.5, 1.5);
    fHistos->CreateTH1(Form("hEventCentrality%s", trg->Data()), Form("Event centrality for trigger class %s", trg->Data()), 103, -2., 101.);
    fHistos->CreateTH1(Form("hVertexZ%s", trg->Data()), Form("z-position of the primary vertex for trigger class %s", trg->Data()), 200, -40., 40.);
    for(int ipatch = 0; ipatch < 10; ipatch++){
      fHistos->CreateTH1(Form("h%sPatchEnergy%s", patchtype[ipatch].Data(), trg->Data()), Form("%s-patch energy for trigger class %s", patchtype[ipatch].Data(), trg->Data()), energybinning);
      fHistos->CreateTH1(Form("h%sPatchET%s", patchtype[ipatch].Data(), trg->Data()), Form("%s-patch transverse energy for trigger class %s", patchtype[ipatch].Data(), trg->Data()), energybinning);
      fHistos->CreateTH2(Form("h%sPatchEnergyEta%s", patchtype[ipatch].Data(), trg->Data()), Form("%s-patch energy for trigger class %s", patchtype[ipatch].Data(), trg->Data()), energybinning, etabinning);
      fHistos->CreateTH2(Form("h%sPatchETEta%s", patchtype[ipatch].Data(), trg->Data()), Form("%s-patch transverse energy for trigger class %s", patchtype[ipatch].Data(), trg->Data()), energybinning, etabinning);
      for(int ien = 0; ien < 5; ien++){
        fHistos->CreateTH2(Form("h%sEtaPhi%dG%s", patchtype[ipatch].Data(), static_cast<int>(encuts[ien]), trg->Data()), Form("%s-patch #eta-#phi map for patches with energy larger than %f GeV/c for trigger class %s", patchtype[ipatch].Data(), encuts[ien], trg->Data()), 100, -0.7, 0.7, 200, 0, TMath::TwoPi());
        fHistos->CreateTH2(Form("h%sColRow%dG%s", patchtype[ipatch].Data(), static_cast<int>(encuts[ien]), trg->Data()), Form("%s-patch col-row map for patches with energy larger than %f GeV/c for trigger class %s", patchtype[ipatch].Data(), encuts[ien], trg->Data()), 48, -0.5, 47.5, 104, -0.5, 103.5);
      }
    }
  }
  PostData(1, fHistos->GetListOfHistograms());
  AliDebugStream(1) << "Histograms done" << std::endl;
}

/**
 * Run event loop
 * @param Not used
 */
void AliAnalysisTaskEmcalPatchesRef::UserExec(Option_t *){
  AliDebugStream(1) << GetName() << ": Start function" << std::endl;
  if(fInitialized){
    ExecOnce();
    fInitialized = kTRUE;
  }
  if(fCurrentRun != InputEvent()->GetRunNumber()){
    RunChanged(InputEvent()->GetRunNumber());
    fCurrentRun = InputEvent()->GetRunNumber();
  }
  TString triggerstring = "";
  if(fTriggerStringFromPatches){
    triggerstring = GetFiredTriggerClassesFromPatches(fTriggerPatches);
  } else {
    triggerstring = fInputEvent->GetFiredTriggerClasses();
  }
  UInt_t selectionstatus = fInputHandler->IsEventSelected();
  Bool_t isMinBias = selectionstatus & AliVEvent::kINT7,
      isEMC7 = (selectionstatus & AliVEvent::kEMC7) && triggerstring.Contains("EMC7"),
      isEJ1 = (selectionstatus & AliVEvent::kEMCEJE) && triggerstring.Contains("EJ1"),
      isEJ2 = (selectionstatus & AliVEvent::kEMCEJE) && triggerstring.Contains("EJ2"),
      isEG1 = (selectionstatus & AliVEvent::kEMCEGA) && triggerstring.Contains("EG1"),
      isEG2 = (selectionstatus & AliVEvent::kEMCEGA) && triggerstring.Contains("EG2"),
      isDMC7 = (selectionstatus & AliVEvent::kEMC7) && triggerstring.Contains("DMC7"),
      isDJ1 = (selectionstatus & AliVEvent::kEMCEJE) && triggerstring.Contains("DJ1"),
      isDJ2 = (selectionstatus & AliVEvent::kEMCEJE) && triggerstring.Contains("DJ2"),
      isDG1 = (selectionstatus & AliVEvent::kEMCEGA) && triggerstring.Contains("DG1"),
      isDG2 = (selectionstatus & AliVEvent::kEMCEGA) && triggerstring.Contains("DG2");
  if(fTriggerPatches && fTriggerSelection){
      isEMC7 &= fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgEL0, fTriggerPatches);
      isEJ1 &=  fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgEJ1, fTriggerPatches);
      isEJ2 &= fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgEJ2, fTriggerPatches);
      isEG1 &= fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgEG1, fTriggerPatches);
      isEG2 &= fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgEG2, fTriggerPatches);
      isDMC7 &= fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgDL0, fTriggerPatches);
      isDJ1 &=  fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgDJ1, fTriggerPatches);
      isDJ2 &= fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgDJ2, fTriggerPatches);
      isDG1 &= fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgDG1, fTriggerPatches);
      isDG2 &= fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgDG2, fTriggerPatches);

  }
  if(!(isMinBias || isEMC7 || isEG1 || isEG2 || isEJ1 || isEJ2 || isDMC7 || isDG1 || isDG2 || isDJ1 || isDJ2)){
    AliDebugStream(1) << GetName() << ": Reject trigger" << std::endl;
    return;
  }
  AliDebugStream(1) << "Event selected" << std::endl;
  AliMultSelection *mult = dynamic_cast<AliMultSelection *>(InputEvent()->FindListObject("MultSelection"));
  // In case a centrality estimator is used, event selection,
  // otherwise ignore event selection from multiplicity task
  if(fRequestCentrality){
    if(mult && !mult->IsEventSelected()) return;
  }
  double centrality =  mult ? mult->GetEstimator("V0M")->GetPercentile() : -1;
  AliDebugStream(1) << GetName()  << ": Centrality " << centrality << std::endl;
  if(!fCentralityRange.IsInRange(centrality)){
    AliDebugStream(1) << GetName() << ": reject centrality: " << centrality << std::endl;
    return;
  } else {
    AliDebugStream(1) << GetName() << ": select centrality " << centrality << std::endl;
  }
  const AliVVertex *vtx = fInputEvent->GetPrimaryVertex();
  if(!vtx) vtx = fInputEvent->GetPrimaryVertexSPD();
  //if(!fInputEvent->IsPileupFromSPD(3, 0.8, 3., 2., 5.)) return;         // reject pileup event
  if(vtx->GetNContributors() < 1){
    AliDebug(1, Form("%s: Reject contributor\n", GetName()));
    return;
  }
  // Fill reference distribution for the primary vertex before any z-cut
  if(fRequestAnalysisUtil){
    AliDebugStream(1) << GetName() << ": Reject analysis util" << std::endl;
    if(fInputEvent->IsA() == AliESDEvent::Class() && fAnalysisUtil->IsFirstEventInChunk(fInputEvent)) return;
    if(!fAnalysisUtil->IsVertexSelected2013pA(fInputEvent)) return;       // Apply new vertex cut
    if(fAnalysisUtil->IsPileUpEvent(fInputEvent)) return;       // Apply new vertex cut
  }
  // Apply vertex z cut
  if(!fVertexRange.IsInRange(vtx->GetZ())){
    AliDebugStream(1) <<  GetName() << ": Reject Z(" << vtx->GetZ() << ")" << std::endl;
    return;
  }
  AliDebugStream(1) << GetName() << ": Event Selected" << std::endl;

  // Fill Event counter and reference vertex distributions for the different trigger classes
  if(isMinBias) FillEventHistograms("MB", centrality, vtx->GetZ());

  // L0 triggers
  if(isEMC7){
    FillEventHistograms("EMC7", centrality, vtx->GetZ());
    // Check for exclusive classes
    if(!(isMinBias || isEMC7)) FillEventHistograms("EMC7excl", centrality, vtx->GetZ());
  }
  if(isDMC7){
    FillEventHistograms("DMC7", centrality, vtx->GetZ());
    // Check for exclusive classes
    if(!(isMinBias  || isDMC7)) FillEventHistograms("DMC7excl", centrality, vtx->GetZ());
  }

  // L1 jet triggers
  if(isEJ1){
    FillEventHistograms("EJ1", centrality, vtx->GetZ());
    // Check for exclusive classes
    if(!(isMinBias || isEMC7 || isEJ2)) FillEventHistograms("EJ1excl", centrality, vtx->GetZ());
  }
  if(isDJ1){
    FillEventHistograms("DJ1", centrality, vtx->GetZ());
    // Check for exclusive classes
    if(!(isMinBias || isDMC7 || isDJ2)) FillEventHistograms("DJ1excl", centrality, vtx->GetZ());
  }

  if(isEJ2){
    FillEventHistograms("EJ2", centrality, vtx->GetZ());
    // Check for exclusive classes
    if(!(isMinBias || isEMC7)) FillEventHistograms("EJ2excl", centrality, vtx->GetZ());
  }
  if(isDJ2){
    FillEventHistograms("DJ2", centrality, vtx->GetZ());
    // Check for exclusive classes
    if(!(isMinBias || isDMC7)) FillEventHistograms("DJ2excl", centrality, vtx->GetZ());
  }

  // L1 gamma triggers
  if(isEG1){
    FillEventHistograms("EG1", centrality, vtx->GetZ());
    // Check for exclusive classes
    if(!(isMinBias || isEMC7 || isEG2)) FillEventHistograms("EG1excl", centrality, vtx->GetZ());
  }
  if(isDG1){
    FillEventHistograms("DG1", centrality, vtx->GetZ());
    // Check for exclusive classes
    if(!(isMinBias || isDMC7 || isDG2)) FillEventHistograms("DG1excl", centrality, vtx->GetZ());
  }
  if(isEG2){
    FillEventHistograms("EG2", centrality, vtx->GetZ());
    // Check for exclusive classes
    if(!(isMinBias || isEMC7)) FillEventHistograms("EG2excl", centrality, vtx->GetZ());
  }
  if(isDG2){
    FillEventHistograms("DG2", centrality, vtx->GetZ());
    // Check for exclusive classes
    if(!(isMinBias || isDMC7)) FillEventHistograms("DG2excl", centrality, vtx->GetZ());
  }

  if(!fTriggerPatches){
    AliErrorStream() << GetName() << ": Trigger patch container not available" << std::endl;
    return;
  }

  AliDebugStream(1) << GetName() << ": Number of trigger patches " << fTriggerPatches->GetEntries() << std::endl;

  Double_t vertexpos[3];
  fInputEvent->GetPrimaryVertex()->GetXYZ(vertexpos);

  Double_t energy, eta, phi, et;
  Int_t col, row;
  for(TIter patchIter = TIter(fTriggerPatches).Begin(); patchIter != TIter::End(); ++patchIter){
    if(!IsOfflineSimplePatch(*patchIter)) continue;
    AliEMCALTriggerPatchInfo *patch = static_cast<AliEMCALTriggerPatchInfo *>(*patchIter);

    bool isDCAL         = SelectDCALPatch(*patchIter),
        isSingleShower  = SelectSingleShowerPatch(*patchIter),
        isJetPatch      = SelectJetPatch(*patchIter);

    std::vector<TString> patchnames;
    if(isJetPatch){
      if(isDCAL){
        patchnames.push_back("DJ1");
        patchnames.push_back("DJ2");
      } else {
        patchnames.push_back("EJ1");
        patchnames.push_back("EJ2");
      }
    }
    if(isSingleShower){
      if(isDCAL){
        patchnames.push_back("DMC7");
        patchnames.push_back("DG1");
        patchnames.push_back("DG2");
      } else {
        patchnames.push_back("EMC7");
        patchnames.push_back("EG1");
        patchnames.push_back("EG2");
      }
    }
    if(!patchnames.size()){
      // Undefined patch type - ignore
      continue;
    }

    TLorentzVector posvec;
    energy = patch->GetPatchE();
    eta = patch->GetEtaGeo();
    phi = patch->GetPhiGeo();
    col = patch->GetColStart();
    row = patch->GetRowStart();
    et = patch->GetLorentzVectorCenterGeo().Et();

    // fill histograms allEta
    for(std::vector<TString>::iterator nameit = patchnames.begin(); nameit != patchnames.end(); ++nameit){
      if(isMinBias){
        FillPatchHistograms("MB", *nameit, energy, et, eta, phi, col, row);
      }

      if(isEMC7){
        FillPatchHistograms("EMC7", *nameit, energy, et, eta, phi, col, row);
        // check for exclusive classes
        if(!isMinBias) FillPatchHistograms("EMC7excl", *nameit, energy, et, eta, phi, col, row);
      }
      if(isDMC7){
        FillPatchHistograms("DMC7", *nameit, energy, et, eta, phi, col, row);
        // check for exclusive classes
        if(!isMinBias) FillPatchHistograms("DMC7excl", *nameit, energy, et, eta, phi, col, row);
      }

      if(isEJ1){
        FillPatchHistograms("EJ1", *nameit, energy, et, eta, phi, col, row);
        // check for exclusive classes
        if(!(isMinBias || isEMC7 || isEJ2))
          FillPatchHistograms("EJ1excl", *nameit, energy, et, eta, phi, col, row);
      }
      if(isDJ1){
        FillPatchHistograms("DJ1", *nameit, energy, et, eta, phi, col, row);
        // check for exclusive classes
        if(!(isMinBias || isDMC7 || isDJ2))
          FillPatchHistograms("DJ1excl", *nameit, energy, et, eta, phi, col, row);
      }

      if(isEJ2){
        FillPatchHistograms("EJ2", *nameit, energy, et, eta, phi, col, row);
        // check for exclusive classes
        if(!(isMinBias || isEMC7)) FillPatchHistograms("EJ2excl", *nameit, energy, et, eta, phi, col, row);
      }
      if(isDJ2){
        FillPatchHistograms("DJ2", *nameit, energy, et, eta, phi, col, row);
        // check for exclusive classes
        if(!(isMinBias || isDMC7)) FillPatchHistograms("DJ2excl", *nameit, energy, et, eta, phi, col, row);
      }

      if(isEG1){
        FillPatchHistograms("EG1", *nameit, energy, et, eta, phi, col, row);
        // check for exclusive classes
        if(!(isMinBias || isEMC7 || isEG2)) FillPatchHistograms("EG1excl", *nameit, energy, et, eta, phi, col, row);
      }
      if(isDG1){
        FillPatchHistograms("DG1", *nameit, energy, et, eta, phi, col, row);
        // check for exclusive classes
        if(!(isMinBias || isEMC7 || isDG2)) FillPatchHistograms("DG1excl", *nameit, energy, et, eta, phi, col, row);
      }

      if(isEG2){
        FillPatchHistograms("EG2", *nameit, energy, et, eta, phi, col, row);
        // check for exclusive classes
        if(!(isMinBias || isEMC7)) FillPatchHistograms("EG2excl", *nameit, energy, et, eta, phi, col, row);
      }
      if(isDG2){
        FillPatchHistograms("DG2", *nameit, energy, et, eta, phi, col, row);
        // check for exclusive classes
        if(!(isMinBias || isDMC7)) FillPatchHistograms("DG2excl", *nameit, energy, et, eta, phi, col, row);
      }
    }
  }
  PostData(1, fHistos->GetListOfHistograms());
}

void AliAnalysisTaskEmcalPatchesRef::ExecOnce(){
  fTriggerPatches = dynamic_cast<TClonesArray *>(fInputEvent->FindListObject("EmcalTriggers"));
  // Handle OADB container with downscaling factors
  if(fNameDownscaleOADB.Length()){
    if(fNameDownscaleOADB.Contains("alien://") && ! gGrid) TGrid::Connect("alien://");
    fDownscaleOADB = new AliOADBContainer("AliEmcalDownscaleFactors");
    fDownscaleOADB->InitFromFile(fNameDownscaleOADB.Data(), "AliEmcalDownscaleFactors");
  }
}

void AliAnalysisTaskEmcalPatchesRef::RunChanged(Int_t runnumber){
 if(fDownscaleOADB){
    fDownscaleFactors = static_cast<TObjArray *>(fDownscaleOADB->GetObject(runnumber));
  }
}

/**
 * Get a trigger class dependent event weight. The weight
 * is defined as 1/downscalefactor. The downscale factor
 * is taken from the OADB. For triggers which are not downscaled
 * the weight is always 1.
 * @param[in] triggerclass Class for which to obtain the trigger.
 * @return Downscale facror for the trigger class (1 if trigger is not downscaled or no OADB container is available)
 */
Double_t AliAnalysisTaskEmcalPatchesRef::GetTriggerWeight(const TString &triggerclass) const {
  if(fDownscaleFactors){
    TParameter<double> *result(nullptr);
    // Downscaling only done on MB, L0 and the low threshold triggers
    if(triggerclass.Contains("MB")) result = static_cast<TParameter<double> *>(fDownscaleFactors->FindObject("MB"));
    else if(triggerclass.Contains("EMC7")) result = static_cast<TParameter<double> *>(fDownscaleFactors->FindObject("EMC7"));
    else if(triggerclass.Contains("EJ2")) result = static_cast<TParameter<double> *>(fDownscaleFactors->FindObject("EJ2"));
    else if(triggerclass.Contains("EG2")) result = static_cast<TParameter<double> *>(fDownscaleFactors->FindObject("EG2"));
    if(result) return 1./result->GetVal();
  }
  return 1.;
}


/**
 * Filling patch related histogram. In case a downscaling correction is
 * available it is applied to the histograms as weight
 * @param[in] triggerclass Name of the trigger class firing the event
 * @param[in] patchname Name of the patchtype
 * @param[in] energy Calibrated energy of the patch
 * @param[in] eta Patch eta at the geometrical center
 * @param[in] phi Patch phi at the geometrical center
 */
void AliAnalysisTaskEmcalPatchesRef::FillPatchHistograms(TString triggerclass, TString patchname, double energy, double transverseenergy, double eta, double phi, int col, int row){
  Double_t weight = GetTriggerWeight(triggerclass);
  fHistos->FillTH1(Form("h%sPatchEnergy%s", patchname.Data(), triggerclass.Data()), energy, weight);
  fHistos->FillTH1(Form("h%sPatchET%s", patchname.Data(), triggerclass.Data()), transverseenergy, weight);
  fHistos->FillTH2(Form("h%sPatchEnergyEta%s", patchname.Data(), triggerclass.Data()), energy, eta, weight);
  fHistos->FillTH2(Form("h%sPatchETEta%s", patchname.Data(), triggerclass.Data()), transverseenergy, eta, weight);
  Double_t encuts[5] = {1., 2., 5., 10., 20.};
  for(int ien = 0; ien < 5; ien++){
    if(energy > encuts[ien]){
      fHistos->FillTH2(Form("h%sEtaPhi%dG%s", patchname.Data(), static_cast<int>(encuts[ien]), triggerclass.Data()), eta, phi, weight);
      fHistos->FillTH2(Form("h%sColRow%dG%s", patchname.Data(), static_cast<int>(encuts[ien]), triggerclass.Data()), col, row, weight);
    }
  }
}

/**
 * Fill event-based histograms. Monitored are
 * - Number of events
 * - Centrality percentile (if available)
 * - z-position of the primary vertex
 * In case a downscaling correction is avaiable it is applied to all
 * histograms as a weight.
 * @param[in] triggerclass Name of the trigger class
 * @param[in] centrality Centrality percentile of the event
 * @param[in] vertexz z-position of the
 */
void AliAnalysisTaskEmcalPatchesRef::FillEventHistograms(const TString &triggerclass, double centrality, double vertexz){
  Double_t weight = GetTriggerWeight(triggerclass);
  fHistos->FillTH1(Form("hEventCount%s", triggerclass.Data()), 1, weight);
  fHistos->FillTH1(Form("hEventCentrality%s", triggerclass.Data()), centrality, weight);
  fHistos->FillTH1(Form("hVertexZ%s", triggerclass.Data()), vertexz, weight);
}

/**
 * Apply trigger selection using offline patches and trigger thresholds based on offline ADC Amplitude
 * @param triggerpatches Trigger patches found by the trigger maker
 * @return String with EMCAL trigger decision
 */
TString AliAnalysisTaskEmcalPatchesRef::GetFiredTriggerClassesFromPatches(const TClonesArray* triggerpatches) const {
  TString triggerstring = "";
  if(!triggerpatches) return triggerstring;
  Int_t nEJ1 = 0, nEJ2 = 0, nEG1 = 0, nEG2 = 0, nDJ1 = 0, nDJ2 = 0, nDG1 = 0, nDG2 = 0;
  double  minADC_J1 = 260.,
          minADC_J2 = 127.,
          minADC_G1 = 140.,
          minADC_G2 = 89.;
  for(TIter patchIter = TIter(triggerpatches).Begin(); patchIter != TIter::End(); ++patchIter){
    AliEMCALTriggerPatchInfo *patch = dynamic_cast<AliEMCALTriggerPatchInfo *>(*patchIter);
    if(!patch->IsOfflineSimple()) continue;
    if(patch->IsJetHighSimple() && patch->GetADCOfflineAmp() > minADC_J1){
      if(patch->IsDCalPHOS()) nDJ1++;
      else nEJ1++;
    }
    if(patch->IsJetLowSimple() && patch->GetADCOfflineAmp() > minADC_J2){
      if(patch->IsDCalPHOS()) nDJ2++;
      else nEJ2++;
    }
    if(patch->IsGammaHighSimple() && patch->GetADCOfflineAmp() > minADC_G1){
      if(patch->IsDCalPHOS()) nDG1++;
      else nEG1++;
    }
    if(patch->IsGammaLowSimple() && patch->GetADCOfflineAmp() > minADC_G2){
      if(patch->IsDCalPHOS()) nDG2++;
      else nEG2++;
    }
  }
  if(nEJ1) triggerstring += "EJ1";
  if(nEJ2){
    if(triggerstring.Length()) triggerstring += ",";
    triggerstring += "EJ2";
  }
  if(nEG1){
    if(triggerstring.Length()) triggerstring += ",";
    triggerstring += "EG1";
  }
  if(nEG2){
    if(triggerstring.Length()) triggerstring += ",";
    triggerstring += "EG2";
  }
  if(nDJ1){
    if(triggerstring.Length()) triggerstring += ",";
    triggerstring += "DJ1";
  }
  if(nEJ2){
    if(triggerstring.Length()) triggerstring += ",";
    triggerstring += "DJ2";
  }
  if(nDG1){
    if(triggerstring.Length()) triggerstring += ",";
    triggerstring += "DG1";
  }
  if(nDG2){
    if(triggerstring.Length()) triggerstring += ",";
    triggerstring += "DG2";
  }
  return triggerstring;
}

void AliAnalysisTaskEmcalPatchesRef::GetPatchBoundaries(TObject *o, Double_t *boundaries) const {
  AliEMCALTriggerPatchInfo *patch = dynamic_cast<AliEMCALTriggerPatchInfo *>(o);
  boundaries[0] = patch->GetEtaMin();
  boundaries[1] = patch->GetEtaMax();
  boundaries[2] = patch->GetPhiMin();
  boundaries[3] = patch->GetPhiMax();
}

bool AliAnalysisTaskEmcalPatchesRef::IsOfflineSimplePatch(TObject *o) const {
  AliEMCALTriggerPatchInfo *patch = dynamic_cast<AliEMCALTriggerPatchInfo *>(o);
  return patch->IsOfflineSimple();
}

bool AliAnalysisTaskEmcalPatchesRef::SelectDCALPatch(TObject *o) const {
  AliEMCALTriggerPatchInfo *patch = dynamic_cast<AliEMCALTriggerPatchInfo *>(o);
  return patch->GetRowStart() >= 64;
}

bool AliAnalysisTaskEmcalPatchesRef::SelectSingleShowerPatch(TObject *o) const{
  AliEMCALTriggerPatchInfo *patch = dynamic_cast<AliEMCALTriggerPatchInfo *>(o);
  if(!patch->IsOfflineSimple()) return false;
  return patch->IsGammaLowSimple();
}

bool AliAnalysisTaskEmcalPatchesRef::SelectJetPatch(TObject *o) const{
  AliEMCALTriggerPatchInfo *patch = dynamic_cast<AliEMCALTriggerPatchInfo *>(o);
  if(!patch->IsOfflineSimple()) return false;
  return patch->IsJetLowSimple();
}

double AliAnalysisTaskEmcalPatchesRef::GetPatchEnergy(TObject *o) const {
  double energy = 0.;
  AliEMCALTriggerPatchInfo *patch = dynamic_cast<AliEMCALTriggerPatchInfo *>(o);
  energy = patch->GetPatchE();
  return energy;
}

AliAnalysisTaskEmcalPatchesRef::EnergyBinning::EnergyBinning():
    TCustomBinning()
{
  this->SetMinimum(0.);
  this->AddStep(1., 0.05);
  this->AddStep(2., 0.1);
  this->AddStep(4, 0.2);
  this->AddStep(7, 0.5);
  this->AddStep(16, 1);
  this->AddStep(32, 2);
  this->AddStep(40, 4);
  this->AddStep(50, 5);
  this->AddStep(100, 10);
  this->AddStep(200, 20);
}


} /* namespace EMCalTriggerPtAnalysis */
