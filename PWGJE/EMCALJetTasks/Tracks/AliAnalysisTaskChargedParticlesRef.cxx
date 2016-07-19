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
#include <array>
#include <iostream>
#include <vector>
#include <map>

#include <TArrayD.h>
#include <TMath.h>
#include <TGrid.h>
#include <THashList.h>
#include <THistManager.h>
#include <TObjArray.h>
#include <TParameter.h>

#include "AliAnalysisUtils.h"
#include "AliAODInputHandler.h"
#include "AliAODTrack.h"
#include "AliEmcalAnalysisFactory.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALRecoUtils.h"
#include "AliEMCALTriggerPatchInfo.h"
#include "AliEmcalTrackSelection.h"
#include "AliEmcalTriggerOfflineSelection.h"
#include "AliESDtrackCuts.h"
#include "AliESDEvent.h"
#include "AliInputEventHandler.h"
#include "AliOADBContainer.h"
#include "AliPIDResponse.h"
#include "AliTOFPIDResponse.h"
#include "AliVVertex.h"

#include "AliEMCalTriggerExtraCuts.h"
#include "AliAnalysisTaskChargedParticlesRef.h"

/// \cond CLASSIMP
ClassImp(EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRef)
/// \endcond

namespace EMCalTriggerPtAnalysis {

/**
 * Dummy constructor
 */
AliAnalysisTaskChargedParticlesRef::AliAnalysisTaskChargedParticlesRef() :
    AliAnalysisTaskSE(),
    fTrackCuts(nullptr),
    fAnalysisUtil(nullptr),
    fTriggerSelection(nullptr),
    fHistos(nullptr),
    fGeometry(nullptr),
    fTriggerPatches(nullptr),
    fTriggerStringFromPatches(kFALSE),
    fYshift(0.465),
    fEtaSign(1),
    fNameDownscaleOADB(""),
    fDownscaleOADB(nullptr),
    fDownscaleFactors(nullptr),
    fCurrentRun(-1),
    fInitialized(false)
{
  // Restrict analysis to the EMCAL acceptance
  fEtaLabCut[0] = -0.6;
  fEtaLabCut[1] = 0.6;
  fEtaCmsCut[0] = -0.13;
  fEtaCmsCut[1] = 0.13;
}

/**
 * Main constructor
 * @param[in] name Name of the task
 */
AliAnalysisTaskChargedParticlesRef::AliAnalysisTaskChargedParticlesRef(const char *name) :
    AliAnalysisTaskSE(name),
    fTrackCuts(nullptr),
    fAnalysisUtil(nullptr),
    fTriggerSelection(nullptr),
    fHistos(nullptr),
    fGeometry(nullptr),
    fTriggerPatches(nullptr),
    fTriggerStringFromPatches(kFALSE),
    fYshift(0.465),
    fEtaSign(1),
    fNameDownscaleOADB(""),
    fDownscaleOADB(nullptr),
    fDownscaleFactors(nullptr),
    fCurrentRun(-1),
    fInitialized(false)
{
  // Restrict analysis to the EMCAL acceptance
  fEtaLabCut[0] = -0.6;
  fEtaLabCut[1] = 0.6;
  fEtaCmsCut[0] = -0.13;
  fEtaCmsCut[1] = 0.13;
  DefineOutput(1, TList::Class());
}

/**
 * Destuctor
 */
AliAnalysisTaskChargedParticlesRef::~AliAnalysisTaskChargedParticlesRef() {
  //if(fTrackCuts) delete fTrackCuts;
  if(fAnalysisUtil) delete fAnalysisUtil;
  if(fTriggerSelection) delete fTriggerSelection;
  if(fHistos) delete fHistos;
  if(fDownscaleOADB) delete fDownscaleOADB;
}

/**
 * Create the output histograms
 */
void AliAnalysisTaskChargedParticlesRef::UserCreateOutputObjects() {
  fAnalysisUtil = new AliAnalysisUtils;

  if(!fTrackCuts) InitializeTrackCuts("standard", fInputHandler->IsA() == AliAODInputHandler::Class());

  OldPtBinning oldbinning;
  NewPtBinning newbinning;

  fHistos = new THistManager("Ref");
  // Exclusive means classes without lower trigger classes (which are downscaled) - in order to make samples statistically independent:
  // MBExcl means MinBias && !EMCAL trigger
  std::array<TString, 11> triggers = {
      "MB", "EMC7",
      "EJ1", "EJ2", "EG1", "EG2",
      "EMC7excl", "EG1excl", "EG2excl", "EJ1excl", "EJ2excl"
  };
  std::array<Double_t, 5> ptcuts = {1., 2., 5., 10., 20.};
  // Binning for the PID histogram
  const int kdimPID = 3;
  const int knbinsPID[kdimPID] = {1000, 200, 300};
  const double kminPID[kdimPID] = {-100., 0.,  0.}, kmaxPID[kdimPID] = {100., 200., 1.5};
  for(auto trg : triggers){
    fHistos->CreateTH1(Form("hEventCount%s", trg.Data()), Form("Event Counter for trigger class %s", trg.Data()), 1, 0.5, 1.5);
    fHistos->CreateTH1(Form("hVertexBefore%s", trg.Data()), Form("Vertex distribution before z-cut for trigger class %s", trg.Data()), 500, -50, 50);
    fHistos->CreateTH1(Form("hVertexAfter%s", trg.Data()), Form("Vertex distribution after z-cut for trigger class %s", trg.Data()), 100, -10, 10);
    fHistos->CreateTH1(Form("hPtEtaAllOldBinning%s", trg.Data()), Form("Charged particle pt distribution all eta old binning trigger %s", trg.Data()), oldbinning);
    fHistos->CreateTH1(Form("hPtEtaCentOldBinning%s", trg.Data()), Form("Charged particle pt distribution central eta old binning trigger %s", trg.Data()), oldbinning);
    fHistos->CreateTH1(Form("hPtEtaAllNewBinning%s", trg.Data()), Form("Charged particle pt distribution all eta new binning trigger %s", trg.Data()), newbinning);
    fHistos->CreateTH1(Form("hPtEtaCentNewBinning%s", trg.Data()), Form("Charged particle pt distribution central eta new binning trigger %s", trg.Data()), newbinning);
    fHistos->CreateTH1(Form("hPtEMCALEtaAllOldBinning%s", trg.Data()), Form("Charged particle in EMCAL pt distribution all eta old binning trigger %s", trg.Data()), oldbinning);
    fHistos->CreateTH1(Form("hPtEMCALEtaCentOldBinning%s", trg.Data()), Form("Charged particle in EMCAL pt distribution central eta old binning trigger %s", trg.Data()), oldbinning);
    fHistos->CreateTH1(Form("hPtEMCALEtaAllNewBinning%s", trg.Data()), Form("Charged particle in EMCAL pt distribution all eta new binning trigger %s", trg.Data()), newbinning);
    fHistos->CreateTH1(Form("hPtEMCALEtaCentNewBinning%s", trg.Data()), Form("Charged particle in EMCAL pt distribution central eta new binning trigger %s", trg.Data()), newbinning);
    fHistos->CreateTH2(Form("hTPCdEdxEMCAL%s", trg.Data()), Form("TPC dE/dx of charged particles in the EMCAL region for trigger %s", trg.Data()), 400, -20., 20., 200, 0., 200.);
    fHistos->CreateTH2(Form("hTOFBetaEMCAL%s", trg.Data()), Form("TOF beta  of charged particles in the EMCAL region for trigger %s", trg.Data()), 400, -20., 20., 150, 0., 1.5);
    fHistos->CreateTHnSparse(Form("hPIDcorrEMCAL%s", trg.Data()), Form("Correlation of PID observables for Trigger %s", trg.Data()), kdimPID, knbinsPID, kminPID, kmaxPID);
    fHistos->CreateTH1(Form("hPtEMCALNoTRDEtaAllOldBinning%s", trg.Data()), Form("Charged particle in EMCAL (no TRD in front) pt distribution all eta old binning trigger %s", trg.Data()), oldbinning);
    fHistos->CreateTH1(Form("hPtEMCALNoTRDEtaCentOldBinning%s", trg.Data()), Form("Charged particle in EMCAL (no TRD in front) pt distribution central eta old binning trigger %s", trg.Data()), oldbinning);
    fHistos->CreateTH1(Form("hPtEMCALNoTRDEtaAllNewBinning%s", trg.Data()), Form("Charged particle in EMCAL (no TRD in front) pt distribution all eta new binning trigger %s", trg.Data()), newbinning);
    fHistos->CreateTH1(Form("hPtEMCALNoTRDEtaCentNewBinning%s", trg.Data()), Form("Charged particle in EMCAL (no TRD in front) pt distribution central eta new binning trigger %s", trg.Data()), newbinning);
    fHistos->CreateTH1(Form("hPtEMCALWithTRDEtaAllOldBinning%s", trg.Data()), Form("Charged particle in EMCAL (with TRD in front) pt distribution all eta old binning trigger %s", trg.Data()), oldbinning);
    fHistos->CreateTH1(Form("hPtEMCALWithTRDEtaCentOldBinning%s", trg.Data()), Form("Charged particle in EMCAL (with TRD in front) pt distribution central eta old binning trigger %s", trg.Data()), oldbinning);
    fHistos->CreateTH1(Form("hPtEMCALWithTRDEtaAllNewBinning%s", trg.Data()), Form("Charged particle in EMCAL (with TRD in front) pt distribution all eta new binning trigger %s", trg.Data()), newbinning);
    fHistos->CreateTH1(Form("hPtEMCALWithTRDEtaCentNewBinning%s", trg.Data()), Form("Charged particle in EMCAL (with TRD in front) pt distribution central eta new binning trigger %s", trg.Data()), newbinning);
    for(auto ptcut : ptcuts) {
      fHistos->CreateTH1(
          Form("hEtaLabDistAllPt%d%s", static_cast<Int_t>(ptcut), trg.Data()),
          Form("Eta (lab) distribution without etacut for tracks with Pt above %.1f GeV/c trigger %s", ptcut, trg.Data()),
          100,
          -1.,
          1.
          );
      fHistos->CreateTH1(
          Form("hEtaLabDistCutPt%d%s", static_cast<Int_t>(ptcut), trg.Data()),
          Form("Eta (lab) distribution with etacut for tracks with Pt above %.1f GeV/c trigger %s", ptcut, trg.Data()),
          100,
          -1.,
          1.
          );
      fHistos->CreateTH1(
          Form("hEtaCentDistAllPt%d%s", static_cast<Int_t>(ptcut), trg.Data()),
          Form("Eta (cent) distribution without etacut for tracks with Pt above %.1f GeV/c trigger %s",
              ptcut, trg.Data()),
              160,
              -1.3,
              1.3
              );
      fHistos->CreateTH1(
          Form("hEtaCentDistCutPt%d%s", static_cast<Int_t>(ptcut), trg.Data()),
          Form("Eta (cent) distribution with etacut for tracks with Pt above %.1f GeV/c trigger %s", ptcut, trg.Data()),
          160,
          -1.3,
          1.3
          );
      fHistos->CreateTH1(
          Form("hEtaLabDistAllEMCALPt%d%s", static_cast<Int_t>(ptcut), trg.Data()),
          Form("Eta (lab) distribution without etacut for tracks in EMCAL with Pt above %.1f GeV/c trigger %s", ptcut, trg.Data()),
          100,
          -1.,
          1.
          );
      fHistos->CreateTH1(
          Form("hEtaLabDistCutEMCALPt%d%s", static_cast<Int_t>(ptcut), trg.Data()),
          Form("Eta (lab) distribution with etacut for tracks in EMCAL with Pt above %.1f GeV/c trigger %s", ptcut, trg.Data()),
          100,
          -1.,
          1.
          );
      fHistos->CreateTH1(
          Form("hEtaCentDistAllEMCALPt%d%s", static_cast<Int_t>(ptcut), trg.Data()),
          Form("Eta (cent) distribution without etacut for tracks in EMCAL with Pt above %.1f GeV/c trigger %s", ptcut, trg.Data()),
          160,
          -1.3,
          1.3
          );
      fHistos->CreateTH1(
          Form("hEtaCentDistCutEMCALPt%d%s", static_cast<Int_t>(ptcut), trg.Data()),
          Form("Eta (cent) distribution with etacut for tracks in EMCAL with Pt above %.1f GeV/c trigger %s", ptcut, trg.Data()),
          160,
          -1.3,
          1.3
          );
      fHistos->CreateTH1(
          Form("hPhiDistAllPt%d%s", static_cast<Int_t>(ptcut), trg.Data()),
          Form("#phi distribution of particles with Pt above %.1f GeV/c trigger %s", ptcut, trg.Data()),
          300,
          0.,
          2*TMath::Pi()
          );
    }
  }
  fHistos->GetListOfHistograms()->Add(fTrackCuts);
  PostData(1, fHistos->GetListOfHistograms());
}

/**
 * Simple unit test framework
 * - Select event using AliAnalysisUtil
 * - Assing trigger type (Request INT7, EJ*, EG*)
 * - Loop over tracks, select particles
 * - Fill distributions
 * @param option Not used
 */
void AliAnalysisTaskChargedParticlesRef::UserExec(Option_t*) {
  if(!fInitialized){
    AliInfoStream() << GetName() << ": Initializing ..." << std::endl;
    ExecOnce();
    fInitialized = kTRUE;
  }
  if(InputEvent()->GetRunNumber() != fCurrentRun){
    AliInfoStream() << GetName() << ": Changing run from " <<  fCurrentRun << " to " << InputEvent()->GetRunNumber() << std::endl;
    RunChanged(InputEvent()->GetRunNumber());
    fCurrentRun = InputEvent()->GetRunNumber();
  }
  // Select event
  TString triggerstring = "";
  if(fTriggerStringFromPatches){
    triggerstring = GetFiredTriggerClassesFromPatches(fTriggerPatches);
  } else {
    triggerstring = fInputEvent->GetFiredTriggerClasses();
  }
  UInt_t selectionstatus = fInputHandler->IsEventSelected();
  Bool_t isMinBias = selectionstatus & AliVEvent::kINT7,
      isEJ1 = (selectionstatus & AliVEvent::kEMCEJE) && triggerstring.Contains("EJ1"),
      isEJ2 = (selectionstatus & AliVEvent::kEMCEJE) && triggerstring.Contains("EJ2"),
      isEG1 = (selectionstatus & AliVEvent::kEMCEGA) && triggerstring.Contains("EG1"),
      isEG2 = (selectionstatus & AliVEvent::kEMCEGA) && triggerstring.Contains("EG2"),
      isEMC7 = (selectionstatus & AliVEvent::kEMC7) && triggerstring.Contains("CEMC7");
  if(fTriggerPatches && fTriggerSelection){
      isEJ1 &= fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgEJ1, fTriggerPatches);
      isEJ2 &= fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgEJ2, fTriggerPatches);
      isEG1 &= fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgEG1, fTriggerPatches);
      isEG2 &= fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgEG2, fTriggerPatches);
      isEMC7 &= fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgEL0, fTriggerPatches);
  }
  if(!(isMinBias || isEMC7 || isEG1 || isEG2 || isEJ1 || isEJ2)) return;
  const AliVVertex *vtx = fInputEvent->GetPrimaryVertex();
  //if(!fInputEvent->IsPileupFromSPD(3, 0.8, 3., 2., 5.)) return;         // reject pileup event
  if(vtx->GetNContributors() < 1) return;
  if(fInputEvent->IsA() == AliESDEvent::Class() && fAnalysisUtil->IsFirstEventInChunk(fInputEvent)) return;
  bool isSelected  = kTRUE;
  if(!fAnalysisUtil->IsVertexSelected2013pA(fInputEvent)) isSelected = kFALSE;       // Apply new vertex cut
  if(fAnalysisUtil->IsPileUpEvent(fInputEvent)) isSelected = kFALSE;       // Apply new vertex cut
  // Apply vertex z cut
  if(vtx->GetZ() < -10. || vtx->GetZ() > 10.) isSelected = kFALSE;

  // Fill Event counter and reference vertex distributions for the different trigger classes
  if(isMinBias){
    FillEventCounterHists("MB", vtx->GetZ(), isSelected);
  }
  if(isEMC7){
    FillEventCounterHists("EMC7", vtx->GetZ(), isSelected);
    if(!isMinBias){
      FillEventCounterHists("EMC7excl", vtx->GetZ(), isSelected);
    }
  }
  if(isEJ2){
    FillEventCounterHists("EJ2", vtx->GetZ(), isSelected);
    // Check for exclusive classes
    if(!(isMinBias || isEMC7)){
      FillEventCounterHists("EJ2excl", vtx->GetZ(), isSelected);
    }
  }
  if(isEJ1){
    FillEventCounterHists("EJ1", vtx->GetZ(), isSelected);
    // Check for exclusive classes
    if(!(isMinBias || isEMC7 || isEJ2)){
      FillEventCounterHists("EJ1excl", vtx->GetZ(), isSelected);
    }
  }
  if(isEG2){
    FillEventCounterHists("EG2", vtx->GetZ(), isSelected);
    // Check for exclusive classes
    if(!(isMinBias || isEMC7)){
      FillEventCounterHists("EG2excl", vtx->GetZ(), isSelected);
    }
  }
  if(isEG1){
    FillEventCounterHists("EG1", vtx->GetZ(), isSelected);
    // Check for exclusive classes
    if(!(isMinBias || isEMC7 || isEG2)){
      FillEventCounterHists("EG1excl", vtx->GetZ(), isSelected);
    }
  }

  if(!isSelected) return;

  // Loop over tracks, fill select particles
  // Histograms
  // - Full eta_{lab} (-0.8, 0.8), new binning
  // - Full eta_{lab} (-0.8, 0.8), old binning
  // - Eta distribution for tracks above 1, 2, 5, 10 GeV/c without eta cut
  // - Central eta_{cms} (-0.3, -0.2), new binning,
  // - Central eta_{cms} (-0.8, -0.2), old binning,
  // - Eta distribution for tracks above 1, 2, 5, 10 GeV/c
  // - Eta distribution for tracks above 1, 2, 5, 10 GeV/c with eta cut
  AliVTrack *checktrack(nullptr);
  int ptmin[5] = {1,2,5,10,20}; // for eta distributions
  Bool_t isEMCAL(kFALSE), hasTRD(kFALSE);
  Double_t etaEMCAL(0.), phiEMCAL(0.);
  for(int itrk = 0; itrk < fInputEvent->GetNumberOfTracks(); ++itrk){
    checktrack = dynamic_cast<AliVTrack *>(fInputEvent->GetTrack(itrk));
    if(!checktrack) continue;
    if((checktrack->Eta() < fEtaLabCut[0]) || (checktrack->Eta() > fEtaLabCut[1])) continue;
    if(TMath::Abs(checktrack->Pt()) < 0.1) continue;
    if(checktrack->IsA() == AliESDtrack::Class()){
      AliESDtrack copytrack(*(static_cast<AliESDtrack *>(checktrack)));
      AliEMCALRecoUtils::ExtrapolateTrackToEMCalSurface(&copytrack);
      etaEMCAL = copytrack.GetTrackEtaOnEMCal();
      phiEMCAL = copytrack.GetTrackPhiOnEMCal();
    } else {
      AliAODTrack copytrack(*(static_cast<AliAODTrack *>(checktrack)));
      AliEMCALRecoUtils::ExtrapolateTrackToEMCalSurface(&copytrack);
      etaEMCAL = copytrack.GetTrackEtaOnEMCal();
      phiEMCAL = copytrack.GetTrackPhiOnEMCal();
    }
    Int_t supermoduleID = -1;
    isEMCAL = fGeometry->SuperModuleNumberFromEtaPhi(etaEMCAL, phiEMCAL, supermoduleID);
    // Exclude supermodules 10 and 11 as they did not participate in the trigger
    isEMCAL = isEMCAL && supermoduleID < 10;
    hasTRD = isEMCAL && supermoduleID >= 4;  // supermodules 4 - 10 have TRD in front in the 2012-2013 ALICE setup

    // Calculate eta in cms frame according
    // EPJC74 (2014) 3054:
    // eta_cms = - eta_lab - |yshift|
    Double_t etacent = -1. * checktrack->Eta() - TMath::Abs(fYshift);
    etacent *= fEtaSign;

    Bool_t etacentcut = etacent > fEtaCmsCut[0] && etacent < fEtaCmsCut[1];

    if(!fTrackCuts->IsTrackAccepted(checktrack)) continue;

    // fill histograms allEta
    if(isMinBias){
      FillTrackHistos("MB", checktrack->Pt(), checktrack->Eta() * fEtaSign, etacent, checktrack->Phi(), etacentcut, isEMCAL, hasTRD);
      if(isEMCAL) FillPIDHistos("MB", *checktrack);
    }
    if(isEMC7){
      FillTrackHistos("EMC7", checktrack->Pt(), checktrack->Eta() * fEtaSign, etacent, checktrack->Phi(), etacentcut, isEMCAL, hasTRD);
      if(isEMCAL) FillPIDHistos("EMC7", *checktrack);
      if(!isMinBias){
        FillTrackHistos("EMC7excl", checktrack->Pt(), checktrack->Eta() * fEtaSign, etacent, checktrack->Phi(), etacentcut, isEMCAL, hasTRD);
        if(isEMCAL) FillPIDHistos("EMC7excl", *checktrack);
      }
    }
    if(isEJ2){
      FillTrackHistos("EJ2", checktrack->Pt(), checktrack->Eta() * fEtaSign, etacent, checktrack->Phi(), etacentcut, isEMCAL, hasTRD);
      if(isEMCAL) FillPIDHistos("EJ2", *checktrack);
      // check for exclusive classes
      if(!(isMinBias || isEMC7)){
        FillTrackHistos("EJ2excl", checktrack->Pt(), checktrack->Eta() * fEtaSign, etacent, checktrack->Phi(), etacentcut, isEMCAL, hasTRD);
        if(isEMCAL) FillPIDHistos("EJ2excl", *checktrack);
      }
    }
    if(isEJ1){
      FillTrackHistos("EJ1", checktrack->Pt(), checktrack->Eta() * fEtaSign, etacent, checktrack->Phi(), etacentcut, isEMCAL, hasTRD);
      if(isEMCAL) FillPIDHistos("EJ1", *checktrack);
      // check for exclusive classes
      if(!(isMinBias ||isEMC7 || isEJ2)){
        FillTrackHistos("EJ1excl", checktrack->Pt(), checktrack->Eta() * fEtaSign, etacent, checktrack->Phi(), etacentcut, isEMCAL, hasTRD);
        if(isEMCAL) FillPIDHistos("EJ1excl", *checktrack);
      }
    }
    if(isEG2){
      FillTrackHistos("EG2", checktrack->Pt(), checktrack->Eta() * fEtaSign, etacent, checktrack->Phi(), etacentcut, isEMCAL, hasTRD);
      if(isEMCAL) FillPIDHistos("EG2", *checktrack);
      // check for exclusive classes
      if(!(isMinBias || isEMC7)){
        FillTrackHistos("EG2excl", checktrack->Pt(), checktrack->Eta() * fEtaSign, etacent, checktrack->Phi(), etacentcut, isEMCAL, hasTRD);
        if(isEMCAL) FillPIDHistos("EG2excl", *checktrack);
      }
    }
    if(isEG1){
      FillTrackHistos("EG1", checktrack->Pt(), checktrack->Eta() * fEtaSign, etacent, checktrack->Phi(), etacentcut, isEMCAL, hasTRD);
      if(isEMCAL) FillPIDHistos("EG1", *checktrack);
      if(!(isMinBias || isEMC7 || isEG2)){
        FillTrackHistos("EG1excl", checktrack->Pt(), checktrack->Eta() * fEtaSign, etacent, checktrack->Phi(), etacentcut, isEMCAL, hasTRD);
        if(isEMCAL) FillPIDHistos("EG1excl", *checktrack);
      }
    }
  }
  PostData(1, fHistos->GetListOfHistograms());
}

/**
 * Fill event counter histogram for a given trigger class
 * @param triggerclass Trigger class firing the trigger
 * @param vtxz z-position of the primary vertex
 * @param isSelected Check whether track is selected
 */
void AliAnalysisTaskChargedParticlesRef::FillEventCounterHists(
    const std::string &triggerclass,
    double vtxz,
    bool isSelected
)
{
  Double_t weight = GetTriggerWeight(triggerclass);
  AliDebugStream(1) << GetName() << ": Using weight " << weight << " for trigger " << triggerclass << " in event histograms." << std::endl;
  // Fill reference distribution for the primary vertex before any z-cut
  fHistos->FillTH1(Form("hVertexBefore%s", triggerclass.c_str()), vtxz, weight);
  if(isSelected){
    // Fill Event counter and reference vertex distributions after event selection
    fHistos->FillTH1(Form("hEventCount%s", triggerclass.c_str()), 1, weight);
    fHistos->FillTH1(Form("hVertexAfter%s", triggerclass.c_str()), vtxz, weight);
  }
}

/**
 * Perform gloabl initializations
 * Used for the moment for
 * - Connect EMCAL geometry
 * - Initialize OADB container with downscaling factors
 */
void AliAnalysisTaskChargedParticlesRef::ExecOnce(){
  // Handle geometry
  if(!fGeometry){
    fGeometry = AliEMCALGeometry::GetInstance();
    if(!fGeometry)
      fGeometry = AliEMCALGeometry::GetInstanceFromRunNumber(InputEvent()->GetRunNumber());
  }
  fTriggerPatches = static_cast<TClonesArray *>(fInputEvent->FindListObject("EmcalTriggers"));
  // Handle OADB container with downscaling factors
  if(fNameDownscaleOADB.Length()){
    if(fNameDownscaleOADB.Contains("alien://") && ! gGrid) TGrid::Connect("alien://");
    fDownscaleOADB = new AliOADBContainer("AliEmcalDownscaleFactors");
    fDownscaleOADB->InitFromFile(fNameDownscaleOADB.Data(), "AliEmcalDownscaleFactors");
  }
}

/**
 * Run change method. Called when the run number of the new event
 * is different compared to the run number of the previous event.
 * Used for loading of the downscale factor for a given
 * run from the downscale OADB.
 * @param[in] runnumber Number of the new run.
 */
void AliAnalysisTaskChargedParticlesRef::RunChanged(Int_t runnumber){
 if(fDownscaleOADB){
    fDownscaleFactors = static_cast<TObjArray *>(fDownscaleOADB->GetObject(runnumber));
  }
}

/**
 * Fill track histograms
 * @param eventclass Trigger class fired
 * @param pt track \f$ p_{t} \f$
 * @param etalab Track \f$ \eta \f$ in lab frame
 * @param etacent Track \f$ \eta \f$ in cms frame
 * @param phi Track \f$ \eta \f$ in lab frame
 * @param etacut Track accepted by \f$ \eta \f$ cut
 * @param inEmcal Track in EMCAL \f$ \phi \f$ acceptance
 */
void AliAnalysisTaskChargedParticlesRef::FillTrackHistos(
    const std::string &eventclass,
    Double_t pt,
    Double_t etalab,
    Double_t etacent,
    Double_t phi,
    Bool_t etacut,
    Bool_t inEmcal,
    Bool_t hasTRD
    )
{
  Double_t weight = GetTriggerWeight(eventclass);
  AliDebugStream(1) << GetName() << ": Using weight " << weight << " for trigger " << eventclass << " in particle histograms." << std::endl;
  fHistos->FillTH1(Form("hPtEtaAllNewBinning%s", eventclass.c_str()), TMath::Abs(pt), weight);
  fHistos->FillTH1(Form("hPtEtaAllOldBinning%s", eventclass.c_str()), TMath::Abs(pt), weight);
  if(inEmcal){
    fHistos->FillTH1(Form("hPtEMCALEtaAllNewBinning%s", eventclass.c_str()), TMath::Abs(pt), weight);
    fHistos->FillTH1(Form("hPtEMCALEtaAllOldBinning%s", eventclass.c_str()), TMath::Abs(pt), weight);
    if(hasTRD){
      fHistos->FillTH1(Form("hPtEMCALWithTRDEtaAllNewBinning%s", eventclass.c_str()), TMath::Abs(pt), weight);
      fHistos->FillTH1(Form("hPtEMCALWithTRDEtaAllOldBinning%s", eventclass.c_str()), TMath::Abs(pt), weight);
    } else {
      fHistos->FillTH1(Form("hPtEMCALNoTRDEtaAllNewBinning%s", eventclass.c_str()), TMath::Abs(pt), weight);
      fHistos->FillTH1(Form("hPtEMCALNoTRDEtaAllOldBinning%s", eventclass.c_str()), TMath::Abs(pt), weight);
    }
  }

  std::array<int, 5> ptmin = {1,2,5,10,20}; // for eta distributions
  for(auto ptmincut : ptmin){
    if(TMath::Abs(pt) > static_cast<double>(ptmincut)){
      fHistos->FillTH1(Form("hPhiDistAllPt%d%s", ptmincut, eventclass.c_str()), phi, weight);
      fHistos->FillTH1(Form("hEtaLabDistAllPt%d%s", ptmincut, eventclass.c_str()), etalab, weight);
      fHistos->FillTH1(Form("hEtaCentDistAllPt%d%s", ptmincut, eventclass.c_str()), etacent, weight);
      if(inEmcal){
        fHistos->FillTH1(Form("hEtaLabDistAllEMCALPt%d%s", ptmincut, eventclass.c_str()), etalab, weight);
        fHistos->FillTH1(Form("hEtaCentDistAllEMCALPt%d%s", ptmincut, eventclass.c_str()), etacent, weight);
      }
    }
  }

  if(etacut){
    fHistos->FillTH1(Form("hPtEtaCentNewBinning%s", eventclass.c_str()), TMath::Abs(pt), weight);
    fHistos->FillTH1(Form("hPtEtaCentOldBinning%s", eventclass.c_str()), TMath::Abs(pt), weight);
    if(inEmcal){
      fHistos->FillTH1(Form("hPtEMCALEtaCentNewBinning%s", eventclass.c_str()), TMath::Abs(pt), weight);
      fHistos->FillTH1(Form("hPtEMCALEtaCentOldBinning%s", eventclass.c_str()), TMath::Abs(pt), weight);
      if(hasTRD){
        fHistos->FillTH1(Form("hPtEMCALWithTRDEtaCentNewBinning%s", eventclass.c_str()), TMath::Abs(pt), weight);
        fHistos->FillTH1(Form("hPtEMCALWithTRDEtaCentOldBinning%s", eventclass.c_str()), TMath::Abs(pt), weight);
      } else {
        fHistos->FillTH1(Form("hPtEMCALNoTRDEtaCentNewBinning%s", eventclass.c_str()), TMath::Abs(pt), weight);
        fHistos->FillTH1(Form("hPtEMCALNoTRDEtaCentOldBinning%s", eventclass.c_str()), TMath::Abs(pt), weight);
      }
    }
    for(auto ptmincut : ptmin){
      if(TMath::Abs(pt) > static_cast<double>(ptmincut)){
        fHistos->FillTH1(Form("hEtaLabDistCutPt%d%s", ptmincut, eventclass.c_str()), etalab, weight);
        fHistos->FillTH1(Form("hEtaCentDistCutPt%d%s", ptmincut, eventclass.c_str()), etacent, weight);
        if(inEmcal){
          fHistos->FillTH1(Form("hEtaLabDistCutEMCALPt%d%s", ptmincut, eventclass.c_str()), etalab, weight);
          fHistos->FillTH1(Form("hEtaCentDistCutEMCALPt%d%s", ptmincut, eventclass.c_str()), etacent, weight);
        }
      }
    }
  }
}

void AliAnalysisTaskChargedParticlesRef::FillPIDHistos(
    const std::string &eventclass,
    const AliVTrack &trk
) {
  Double_t weight = GetTriggerWeight(eventclass);
  AliDebugStream(1) << GetName() << ": Using weight " << weight << " for trigger " << eventclass << " in PID histograms." << std::endl;
  AliPIDResponse *pid = fInputHandler->GetPIDResponse();
  if(TMath::Abs(trk.Eta()) > 0.5) return;
  if(!((trk.GetStatus() & AliVTrack::kTOFout) && (trk.GetStatus() & AliVTrack::kTIME))) return;

  double poverz = TMath::Abs(trk.P())/static_cast<double>(trk.Charge());
  fHistos->FillTH2(Form("hTPCdEdxEMCAL%s", eventclass.c_str()), poverz, trk.GetTPCsignal(), weight);
  // correct for units - TOF in ps, track length in cm
  Double_t trtime = (trk.GetTOFsignal() - pid->GetTOFResponse().GetTimeZero()) * 1e-12;
  Double_t v = trk.GetIntegratedLength()/(100. * trtime);
  Double_t beta =  v / TMath::C();
  fHistos->FillTH2(Form("hTOFBetaEMCAL%s", eventclass.c_str()), poverz, beta, weight);
  double datapoint[3] = {poverz, trk.GetTPCsignal(), beta};
  fHistos->FillTHnSparse(Form("hPIDcorrEMCAL%s", eventclass.c_str()), datapoint, weight);
}

/**
 * Get a trigger class dependent event weight. The weight
 * is defined as 1/downscalefactor. The downscale factor
 * is taken from the OADB. For triggers which are not downscaled
 * the weight is always 1.
 * @param[in] triggerclass Class for which to obtain the trigger.
 * @return Downscale facror for the trigger class (1 if trigger is not downscaled or no OADB container is available)
 */
Double_t AliAnalysisTaskChargedParticlesRef::GetTriggerWeight(const std::string &triggerclass) const {
  if(fDownscaleFactors){
    TParameter<double> *result(nullptr);
    // Downscaling only done on MB, L0 and the low threshold triggers
    if(triggerclass.find("MB") != std::string::npos) result = static_cast<TParameter<double> *>(fDownscaleFactors->FindObject("INT7"));
    else if(triggerclass.find("EMC7") != std::string::npos) result = static_cast<TParameter<double> *>(fDownscaleFactors->FindObject("EMC7"));
    else if(triggerclass.find("EJ2") != std::string::npos) result = static_cast<TParameter<double> *>(fDownscaleFactors->FindObject("EJ2"));
    else if(triggerclass.find("EG2") != std::string::npos) result = static_cast<TParameter<double> *>(fDownscaleFactors->FindObject("EG2"));
    if(result) return 1./result->GetVal();
  }
  return 1.;
}

/**
 * Set the track selection
 * @param cutname Name of the track cuts
 * @param isAOD check whether we run on ESDs or AODs
 */
void AliAnalysisTaskChargedParticlesRef::InitializeTrackCuts(TString cutname, bool isAOD){
  SetTrackSelection(AliEmcalAnalysisFactory::TrackCutsFactory(cutname, isAOD));
}

/**
 * Apply trigger selection using offline patches and trigger thresholds based on offline ADC Amplitude
 * @param triggerpatches Trigger patches found by the trigger maker
 * @return String with EMCAL trigger decision
 */
TString AliAnalysisTaskChargedParticlesRef::GetFiredTriggerClassesFromPatches(const TClonesArray* triggerpatches) const {
  TString triggerstring = "";
  Int_t nEJ1 = 0, nEJ2 = 0, nEG1 = 0, nEG2 = 0;
  double  minADC_EJ1 = 260.,
          minADC_EJ2 = 127.,
          minADC_EG1 = 140.,
          minADC_EG2 = 89.;
  for(TIter patchIter = TIter(triggerpatches).Begin(); patchIter != TIter::End(); ++patchIter){
    AliEMCALTriggerPatchInfo *patch = dynamic_cast<AliEMCALTriggerPatchInfo *>(*patchIter);
    if(!patch->IsOfflineSimple()) continue;
    if(patch->IsJetHighSimple() && patch->GetADCOfflineAmp() > minADC_EJ1) nEJ1++;
    if(patch->IsJetLowSimple() && patch->GetADCOfflineAmp() > minADC_EJ2) nEJ2++;
    if(patch->IsGammaHighSimple() && patch->GetADCOfflineAmp() > minADC_EG1) nEG1++;
    if(patch->IsGammaLowSimple() && patch->GetADCOfflineAmp() > minADC_EG2) nEG2++;
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
  return triggerstring;
}

/**
 * @brief Constructor
 * Create old pt binning
 */
AliAnalysisTaskChargedParticlesRef::OldPtBinning::OldPtBinning() :
  TCustomBinning()
{
  this->SetMinimum(0.);
  this->AddStep(2.5, 0.1);
  this->AddStep(7., 0.25);
  this->AddStep(15., 0.5);
  this->AddStep(25., 1.);
  this->AddStep(40., 2.5);
  this->AddStep(50., 5.);
  this->AddStep(100., 10.);
}

/**
 * @brief Constructor
 * Create new Pt binning
 */
AliAnalysisTaskChargedParticlesRef::NewPtBinning::NewPtBinning():
    TCustomBinning()
{
  this->SetMinimum(0.);
  this->AddStep(1, 0.05);
  this->AddStep(2, 0.1);
  this->AddStep(4, 0.2);
  this->AddStep(7, 0.5);
  this->AddStep(16, 1);
  this->AddStep(36, 2);
  this->AddStep(40, 4);
  this->AddStep(50, 5);
  this->AddStep(100, 10);
  this->AddStep(200, 20);
}

} /* namespace EMCalTriggerPtAnalysis */
