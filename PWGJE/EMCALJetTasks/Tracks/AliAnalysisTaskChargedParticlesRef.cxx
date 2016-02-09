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
#include <vector>
#include <map>

#include <TArrayD.h>
#include <TMath.h>
#include <THashList.h>
#include <TString.h>

#include "AliAnalysisUtils.h"
#include "AliAODTrack.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALRecoUtils.h"
#include "AliESDtrackCuts.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliInputEventHandler.h"
#include "AliVVertex.h"

#include "AliEMCALTriggerPatchInfo.h"
#include "AliEMCalHistoContainer.h"
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
    fTrackCuts(NULL),
    fAnalysisUtil(NULL),
    fHistos(NULL),
    fGeometry(NULL),
    fTriggerStringFromPatches(kFALSE),
    fYshift(0.465),
    fEtaSign(1),
    fSwitchoffSPDcut(kFALSE),
    fSwitchoffITScut(kFALSE)
{
  // Restrict analysis to the EMCAL acceptance
  fEtaLabCut[0] = -0.6;
  fEtaLabCut[1] = 0.6;
  fEtaCmsCut[0] = -0.13;
  fEtaCmsCut[1] = 0.13;
  for(int itrg = 0; itrg < kCPRntrig; itrg++){
    fOfflineEnergyThreshold[itrg] = -1.;
  }
}

/**
 * Main constructor
 * @param name Name of the task
 */
AliAnalysisTaskChargedParticlesRef::AliAnalysisTaskChargedParticlesRef(const char *name) :
    AliAnalysisTaskSE(name),
    fTrackCuts(NULL),
    fAnalysisUtil(NULL),
    fHistos(NULL),
    fGeometry(NULL),
    fTriggerStringFromPatches(kFALSE),
    fYshift(0.465),
    fEtaSign(1),
    fSwitchoffSPDcut(kFALSE),
    fSwitchoffITScut(kFALSE)
{
  // Restrict analysis to the EMCAL acceptance
  fEtaLabCut[0] = -0.6;
  fEtaLabCut[1] = 0.6;
  fEtaCmsCut[0] = -0.13;
  fEtaCmsCut[1] = 0.13;
  for(int itrg = 0; itrg < kCPRntrig; itrg++){
    fOfflineEnergyThreshold[itrg] = -1.;
  }
  DefineOutput(1, TList::Class());
}

/**
 * Destuctor
 */
AliAnalysisTaskChargedParticlesRef::~AliAnalysisTaskChargedParticlesRef() {
  //if(fTrackCuts) delete fTrackCuts;
  if(fAnalysisUtil) delete fAnalysisUtil;
  if(fHistos) delete fHistos;
}

/**
 * Create the output histograms
 */
void AliAnalysisTaskChargedParticlesRef::UserCreateOutputObjects() {
  fAnalysisUtil = new AliAnalysisUtils;

  fTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(true, 1);
  fTrackCuts->DefineHistograms(kRed);
  fTrackCuts->SetName("Standard Track cuts");
  fTrackCuts->SetMinNCrossedRowsTPC(120);
  fTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
  if(fSwitchoffSPDcut || fSwitchoffITScut){
    fTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kOff);
  }
  if(fSwitchoffITScut){
    fTrackCuts->SetRequireITSRefit(kFALSE);           // Allows also TPC-only tracks
    fTrackCuts->SetMinNClustersITS(0);
  }

  TArrayD oldbinning, newbinning;
  CreateOldPtBinning(oldbinning);
  CreateNewPtBinning(newbinning);

  fHistos = new AliEMCalHistoContainer("Ref");
  // Exclusive means classes without higher trigger classes:
  // EG2excl means EG2 && !EG1
  // EJ2excl means EJ2 && !EJ1
  // MBExcl means MinBias && !EMCAL trigger
  // Combined means: gamma and ANY jet class, jet and ANY gamma class
  // Jonly means: No gamma class fired at the same time
  // Gonly means: No jet class fired at the same time
  TString triggers[17] = {
      "MB", "EMC7",
      "EJ1", "EJ2", "EG1", "EG2",
      "EMC7excl", "EG1excl", "EG2excl", "EJ1excl", "EJ2excl",
      "E1combined", "E1Jonly", "E1Gonly", "E2combined", "E2Jonly", "E2Gonly"
  };
  Double_t ptcuts[5] = {1., 2., 5., 10., 20.};
  for(TString *trg = triggers; trg < triggers + sizeof(triggers)/sizeof(TString); trg++){
    fHistos->CreateTH1(Form("hEventCount%s", trg->Data()), Form("Event Counter for trigger class %s", trg->Data()), 1, 0.5, 1.5);
    fHistos->CreateTH1(Form("hVertexBefore%s", trg->Data()), Form("Vertex distribution before z-cut for trigger class %s", trg->Data()), 500, -50, 50);
    fHistos->CreateTH1(Form("hVertexAfter%s", trg->Data()), Form("Vertex distribution after z-cut for trigger class %s", trg->Data()), 100, -10, 10);
    fHistos->CreateTH1(Form("hPtEtaAllOldBinning%s", trg->Data()), Form("Charged particle pt distribution all eta old binning trigger %s", trg->Data()), oldbinning);
    fHistos->CreateTH1(Form("hPtEtaCentOldBinning%s", trg->Data()), Form("Charged particle pt distribution central eta old binning trigger %s", trg->Data()), oldbinning);
    fHistos->CreateTH1(Form("hPtEtaAllNewBinning%s", trg->Data()), Form("Charged particle pt distribution all eta new binning trigger %s", trg->Data()), newbinning);
    fHistos->CreateTH1(Form("hPtEtaCentNewBinning%s", trg->Data()), Form("Charged particle pt distribution central eta new binning trigger %s", trg->Data()), newbinning);
    fHistos->CreateTH1(Form("hPtEMCALEtaAllOldBinning%s", trg->Data()), Form("Charged particle in EMCAL pt distribution all eta old binning trigger %s", trg->Data()), oldbinning);
    fHistos->CreateTH1(Form("hPtEMCALEtaCentOldBinning%s", trg->Data()), Form("Charged particle in EMCAL pt distribution central eta old binning trigger %s", trg->Data()), oldbinning);
    fHistos->CreateTH1(Form("hPtEMCALEtaAllNewBinning%s", trg->Data()), Form("Charged particle in EMCAL pt distribution all eta new binning trigger %s", trg->Data()), newbinning);
    fHistos->CreateTH1(Form("hPtEMCALEtaCentNewBinning%s", trg->Data()), Form("Charged particle in EMCAL pt distribution central eta new binning trigger %s", trg->Data()), newbinning);
    fHistos->CreateTH1(Form("hPtEMCALNoTRDEtaAllOldBinning%s", trg->Data()), Form("Charged particle in EMCAL (no TRD in front) pt distribution all eta old binning trigger %s", trg->Data()), oldbinning);
    fHistos->CreateTH1(Form("hPtEMCALNoTRDEtaCentOldBinning%s", trg->Data()), Form("Charged particle in EMCAL (no TRD in front) pt distribution central eta old binning trigger %s", trg->Data()), oldbinning);
    fHistos->CreateTH1(Form("hPtEMCALNoTRDEtaAllNewBinning%s", trg->Data()), Form("Charged particle in EMCAL (no TRD in front) pt distribution all eta new binning trigger %s", trg->Data()), newbinning);
    fHistos->CreateTH1(Form("hPtEMCALNoTRDEtaCentNewBinning%s", trg->Data()), Form("Charged particle in EMCAL (no TRD in front) pt distribution central eta new binning trigger %s", trg->Data()), newbinning);
    fHistos->CreateTH1(Form("hPtEMCALWithTRDEtaAllOldBinning%s", trg->Data()), Form("Charged particle in EMCAL (with TRD in front) pt distribution all eta old binning trigger %s", trg->Data()), oldbinning);
    fHistos->CreateTH1(Form("hPtEMCALWithTRDEtaCentOldBinning%s", trg->Data()), Form("Charged particle in EMCAL (with TRD in front) pt distribution central eta old binning trigger %s", trg->Data()), oldbinning);
    fHistos->CreateTH1(Form("hPtEMCALWithTRDEtaAllNewBinning%s", trg->Data()), Form("Charged particle in EMCAL (with TRD in front) pt distribution all eta new binning trigger %s", trg->Data()), newbinning);
    fHistos->CreateTH1(Form("hPtEMCALWithTRDEtaCentNewBinning%s", trg->Data()), Form("Charged particle in EMCAL (with TRD in front) pt distribution central eta new binning trigger %s", trg->Data()), newbinning);
    for(int ipt = 0; ipt < 5; ipt++){
      fHistos->CreateTH1(
          Form("hEtaLabDistAllPt%d%s", static_cast<Int_t>(ptcuts[ipt]), trg->Data()),
          Form("Eta (lab) distribution without etacut for tracks with Pt above 1 GeV/c trigger %s", trg->Data()),
          100,
          -1.,
          1.
          );
      fHistos->CreateTH1(
          Form("hEtaLabDistCutPt%d%s", static_cast<Int_t>(ptcuts[ipt]), trg->Data()),
          Form("Eta (lab) distribution with etacut for tracks with Pt above 1 GeV/c trigger %s", trg->Data()),
          100,
          -1.,
          1.
          );
      fHistos->CreateTH1(
          Form("hEtaCentDistAllPt%d%s", static_cast<Int_t>(ptcuts[ipt]), trg->Data()),
          Form("Eta (cent) distribution without etacut for tracks with Pt above 1 GeV/c trigger %s",
              trg->Data()),
              160,
              -1.3,
              1.3
              );
      fHistos->CreateTH1(
          Form("hEtaCentDistCutPt%d%s", static_cast<Int_t>(ptcuts[ipt]), trg->Data()),
          Form("Eta (cent) distribution with etacut for tracks with Pt above 1 GeV/c trigger %s", trg->Data()),
          160,
          -1.3,
          1.3
          );
      fHistos->CreateTH1(
          Form("hEtaLabDistAllEMCALPt%d%s", static_cast<Int_t>(ptcuts[ipt]), trg->Data()),
          Form("Eta (lab) distribution without etacut for tracks in EMCAL with Pt above %.1f GeV/c trigger %s", ptcuts[ipt], trg->Data()),
          100,
          -1.,
          1.
          );
      fHistos->CreateTH1(
          Form("hEtaLabDistCutEMCALPt%d%s", static_cast<Int_t>(ptcuts[ipt]), trg->Data()),
          Form("Eta (lab) distribution with etacut for tracks in EMCAL with Pt above %.1f GeV/c trigger %s", ptcuts[ipt], trg->Data()),
          100,
          -1.,
          1.
          );
      fHistos->CreateTH1(
          Form("hEtaCentDistAllEMCALPt%d%s", static_cast<Int_t>(ptcuts[ipt]), trg->Data()),
          Form("Eta (cent) distribution without etacut for tracks in EMCAL with Pt above %.1f GeV/c trigger %s", ptcuts[ipt], trg->Data()),
          160,
          -1.3,
          1.3
          );
      fHistos->CreateTH1(
          Form("hEtaCentDistCutEMCALPt%d%s", static_cast<Int_t>(ptcuts[ipt]), trg->Data()),
          Form("Eta (cent) distribution with etacut for tracks in EMCAL with Pt above %.1f GeV/c trigger %s", ptcuts[ipt], trg->Data()),
          160,
          -1.3,
          1.3
          );
      fHistos->CreateTH1(
          Form("hPhiDistAllPt%d%s", static_cast<Int_t>(ptcuts[ipt]), trg->Data()),
          Form("#phi distribution of particles with Pt above %.1f GeV/c trigger %s", ptcuts[ipt], trg->Data()),
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
  if(!fGeometry){
    fGeometry = AliEMCALGeometry::GetInstance();
    if(!fGeometry)
      fGeometry = AliEMCALGeometry::GetInstanceFromRunNumber(InputEvent()->GetRunNumber());
  }
  // Select event
  TClonesArray *triggerpatches = static_cast<TClonesArray *>(fInputEvent->FindListObject("EmcalTriggers"));
  TString triggerstring = "";
  if(fTriggerStringFromPatches){
    triggerstring = GetFiredTriggerClassesFromPatches(triggerpatches);
  } else {
    triggerstring = fInputEvent->GetFiredTriggerClasses();
  }
  UInt_t selectionstatus = fInputHandler->IsEventSelected();
  Bool_t isMinBias = selectionstatus & AliVEvent::kINT7,
      isEJ1 = (selectionstatus & AliVEvent::kEMCEJE) && triggerstring.Contains("EJ1") && IsOfflineSelected(kCPREJ1, triggerpatches),
      isEJ2 = (selectionstatus & AliVEvent::kEMCEJE) && triggerstring.Contains("EJ2") && IsOfflineSelected(kCPREJ2, triggerpatches),
      isEG1 = (selectionstatus & AliVEvent::kEMCEGA) && triggerstring.Contains("EG1") && IsOfflineSelected(kCPREG1, triggerpatches),
      isEG2 = (selectionstatus & AliVEvent::kEMCEGA) && triggerstring.Contains("EG2") && IsOfflineSelected(kCPREG2, triggerpatches),
      isEMC7 = (selectionstatus & AliVEvent::kEMC7) && triggerstring.Contains("CEMC7") && IsOfflineSelected(kCPREL0, triggerpatches);
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
    if(!isMinBias){
      FillEventCounterHists("EJ2excl", vtx->GetZ(), isSelected);
    }
    if(isEG1 || isEG2){
      FillEventCounterHists("E2combined", vtx->GetZ(), isSelected);
    } else {
      FillEventCounterHists("E2Jonly", vtx->GetZ(), isSelected);
    }
  }
  if(isEJ1){
    FillEventCounterHists("EJ1", vtx->GetZ(), isSelected);
    // Check for exclusive classes
    if(!(isMinBias || isEJ2)){
      FillEventCounterHists("EJ1excl", vtx->GetZ(), isSelected);
    }
    if(isEG1 || isEG2){
      FillEventCounterHists("E1combined", vtx->GetZ(), isSelected);
    } else {
      FillEventCounterHists("E1Jonly", vtx->GetZ(), isSelected);
    }
  }
  if(isEG2){
    FillEventCounterHists("EG2", vtx->GetZ(), isSelected);
    // Check for exclusive classes
    if(!isMinBias){
      FillEventCounterHists("EG2excl", vtx->GetZ(), isSelected);
    }
    if(!(isEJ1 || isEJ2)){
      FillEventCounterHists("E2Gonly", vtx->GetZ(), isSelected);
    }
  }
  if(isEG1){
    FillEventCounterHists("EG1", vtx->GetZ(), isSelected);
    // Check for exclusive classes
    if(!(isMinBias || isEG2 )){
      FillEventCounterHists("EG1excl", vtx->GetZ(), isSelected);
    }
    if(!(isEJ1 || isEJ2)){
      FillEventCounterHists("E1Gonly", vtx->GetZ(), isSelected);
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
  AliVTrack *checktrack(NULL);
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

    // Distinguish track selection for ESD and AOD tracks
    AliESDtrack *esdtrack(NULL);
    AliAODTrack *aodtrack(NULL);
    if((esdtrack = dynamic_cast<AliESDtrack *>(checktrack))){
      if(!TrackSelectionESD(esdtrack)) continue;
    } else if((aodtrack = dynamic_cast<AliAODTrack *>(checktrack))){
      if(!TrackSelectionAOD(aodtrack)) continue;
    } else {
      continue;
    }

    // fill histograms allEta
    if(isMinBias){
      FillTrackHistos("MB", checktrack->Pt(), checktrack->Eta() * fEtaSign, etacent, checktrack->Phi(), etacentcut, isEMCAL, hasTRD);
    }
    if(isEMC7){
      FillTrackHistos("EMC7", checktrack->Pt(), checktrack->Eta() * fEtaSign, etacent, checktrack->Phi(), etacentcut, isEMCAL, hasTRD);
      if(!isMinBias){
        FillTrackHistos("EMC7excl", checktrack->Pt(), checktrack->Eta() * fEtaSign, etacent, checktrack->Phi(), etacentcut, isEMCAL, hasTRD);
      }
    }
    if(isEJ2){
      FillTrackHistos("EJ2", checktrack->Pt(), checktrack->Eta() * fEtaSign, etacent, checktrack->Phi(), etacentcut, isEMCAL, hasTRD);
      // check for exclusive classes
      if(!isMinBias){
        FillTrackHistos("EJ2excl", checktrack->Pt(), checktrack->Eta() * fEtaSign, etacent, checktrack->Phi(), etacentcut, isEMCAL, hasTRD);
      }
      if(isEG1 || isEG2){
        FillTrackHistos("E2combined", checktrack->Pt(), checktrack->Eta() * fEtaSign, etacent, checktrack->Phi(), etacentcut, isEMCAL, hasTRD);
      } else {
        FillTrackHistos("E2Jonly", checktrack->Pt(), checktrack->Eta() * fEtaSign, etacent, checktrack->Phi(), etacentcut, isEMCAL, hasTRD);
      }
    }
    if(isEJ1){
      FillTrackHistos("EJ1", checktrack->Pt(), checktrack->Eta() * fEtaSign, etacent, checktrack->Phi(), etacentcut, isEMCAL, hasTRD);
      // check for exclusive classes
      if(!(isMinBias || isEJ2)){
        FillTrackHistos("EJ1excl", checktrack->Pt(), checktrack->Eta() * fEtaSign, etacent, checktrack->Phi(), etacentcut, isEMCAL, hasTRD);
      }
      if(isEG1 || isEG2) {
        FillTrackHistos("E1combined", checktrack->Pt(), checktrack->Eta() * fEtaSign, etacent, checktrack->Phi(), etacentcut, isEMCAL, hasTRD);
      } else {
        FillTrackHistos("E1Jonly", checktrack->Pt(), checktrack->Eta() * fEtaSign, etacent, checktrack->Phi(), etacentcut, isEMCAL, hasTRD);
      }
    }
    if(isEG2){
      FillTrackHistos("EG2", checktrack->Pt(), checktrack->Eta() * fEtaSign, etacent, checktrack->Phi(), etacentcut, isEMCAL, hasTRD);
      // check for exclusive classes
      if(!isMinBias){
        FillTrackHistos("EG2excl", checktrack->Pt(), checktrack->Eta() * fEtaSign, etacent, checktrack->Phi(), etacentcut, isEMCAL, hasTRD);
      }
      if(!(isEJ2 || isEJ1)){
        FillTrackHistos("E2Gonly", checktrack->Pt(), checktrack->Eta() * fEtaSign, etacent, checktrack->Phi(), etacentcut, isEMCAL, hasTRD);
      }
    }
    if(isEG1){
      FillTrackHistos("EG1", checktrack->Pt(), checktrack->Eta() * fEtaSign, etacent, checktrack->Phi(), etacentcut, isEMCAL, hasTRD);
      if(!(isMinBias || isEG2)){
        FillTrackHistos("EG1excl", checktrack->Pt(), checktrack->Eta() * fEtaSign, etacent, checktrack->Phi(), etacentcut, isEMCAL, hasTRD);
      }
      if(!(isEJ1 || isEJ2)){
        FillTrackHistos("E1Gonly", checktrack->Pt(), checktrack->Eta() * fEtaSign, etacent, checktrack->Phi(), etacentcut, isEMCAL, hasTRD);
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
    const char *triggerclass,
    double vtxz,
    bool isSelected
)
{
  // Fill reference distribution for the primary vertex before any z-cut
  fHistos->FillTH1(Form("hVertexBefore%s", triggerclass), vtxz);
  if(isSelected){
    // Fill Event counter and reference vertex distributions after event selection
    fHistos->FillTH1(Form("hEventCount%s", triggerclass), 1);
    fHistos->FillTH1(Form("hVertexAfter%s", triggerclass), vtxz);
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
    const char *eventclass,
    Double_t pt,
    Double_t etalab,
    Double_t etacent,
    Double_t phi,
    Bool_t etacut,
    Bool_t inEmcal,
    Bool_t hasTRD
    )
{
  fHistos->FillTH1(Form("hPtEtaAllNewBinning%s", eventclass), TMath::Abs(pt));
  fHistos->FillTH1(Form("hPtEtaAllOldBinning%s", eventclass), TMath::Abs(pt));
  if(inEmcal){
    fHistos->FillTH1(Form("hPtEMCALEtaAllNewBinning%s", eventclass), TMath::Abs(pt));
    fHistos->FillTH1(Form("hPtEMCALEtaAllOldBinning%s", eventclass), TMath::Abs(pt));
    if(hasTRD){
      fHistos->FillTH1(Form("hPtEMCALWithTRDEtaAllNewBinning%s", eventclass), TMath::Abs(pt));
      fHistos->FillTH1(Form("hPtEMCALWithTRDEtaAllOldBinning%s", eventclass), TMath::Abs(pt));
    } else {
      fHistos->FillTH1(Form("hPtEMCALNoTRDEtaAllNewBinning%s", eventclass), TMath::Abs(pt));
      fHistos->FillTH1(Form("hPtEMCALNoTRDEtaAllOldBinning%s", eventclass), TMath::Abs(pt));
    }
  }

  int ptmin[5] = {1,2,5,10,20}; // for eta distributions
  for(int icut = 0; icut < 5; icut++){
    if(TMath::Abs(pt) > static_cast<double>(ptmin[icut])){
      fHistos->FillTH1(Form("hPhiDistAllPt%d%s", ptmin[icut], eventclass), phi);
      fHistos->FillTH1(Form("hEtaLabDistAllPt%d%s", ptmin[icut], eventclass), etalab);
      fHistos->FillTH1(Form("hEtaCentDistAllPt%d%s", ptmin[icut], eventclass), etacent);
      if(inEmcal){
        fHistos->FillTH1(Form("hEtaLabDistAllEMCALPt%d%s", ptmin[icut], eventclass), etalab);
        fHistos->FillTH1(Form("hEtaCentDistAllEMCALPt%d%s", ptmin[icut], eventclass), etacent);
      }
    }
  }

  if(etacut){
    fHistos->FillTH1(Form("hPtEtaCentNewBinning%s", eventclass), TMath::Abs(pt));
    fHistos->FillTH1(Form("hPtEtaCentOldBinning%s", eventclass), TMath::Abs(pt));
    if(inEmcal){
      fHistos->FillTH1(Form("hPtEMCALEtaCentNewBinning%s", eventclass), TMath::Abs(pt));
      fHistos->FillTH1(Form("hPtEMCALEtaCentOldBinning%s", eventclass), TMath::Abs(pt));
      if(hasTRD){
        fHistos->FillTH1(Form("hPtEMCALWithTRDEtaCentNewBinning%s", eventclass), TMath::Abs(pt));
        fHistos->FillTH1(Form("hPtEMCALWithTRDEtaCentOldBinning%s", eventclass), TMath::Abs(pt));
      } else {
        fHistos->FillTH1(Form("hPtEMCALNoTRDEtaCentNewBinning%s", eventclass), TMath::Abs(pt));
        fHistos->FillTH1(Form("hPtEMCALNoTRDEtaCentOldBinning%s", eventclass), TMath::Abs(pt));
      }
    }
    for(int icut = 0; icut < 5; icut++){
      if(TMath::Abs(pt) > static_cast<double>(ptmin[icut])){
        fHistos->FillTH1(Form("hEtaLabDistCutPt%d%s", ptmin[icut], eventclass), etalab);
        fHistos->FillTH1(Form("hEtaCentDistCutPt%d%s", ptmin[icut], eventclass), etacent);
        if(inEmcal){
          fHistos->FillTH1(Form("hEtaLabDistCutEMCALPt%d%s", ptmin[icut], eventclass), etalab);
          fHistos->FillTH1(Form("hEtaCentDistCutEMCALPt%d%s", ptmin[icut], eventclass), etacent);
        }
      }
    }
  }
}


/**
 * Create old pt binning
 * @param binning
 */
void AliAnalysisTaskChargedParticlesRef::CreateOldPtBinning(TArrayD &binning) const{
 std::vector<double> mybinning;
 std::map<double,double> definitions;
 definitions.insert(std::pair<double,double>(2.5, 0.1));
 definitions.insert(std::pair<double,double>(7., 0.25));
 definitions.insert(std::pair<double,double>(15., 0.5));
 definitions.insert(std::pair<double,double>(25., 1.));
 definitions.insert(std::pair<double,double>(40., 2.5));
 definitions.insert(std::pair<double,double>(50., 5.));
 definitions.insert(std::pair<double,double>(100., 10.));
 double currentval = 0;
 for(std::map<double,double>::iterator id = definitions.begin(); id != definitions.end(); ++id){
   double limit = id->first, binwidth = id->second;
   while(currentval < limit){
     currentval += binwidth;
     mybinning.push_back(currentval);
    }
  }
  binning.Set(mybinning.size());
  int ib = 0;
  for(std::vector<double>::iterator it = mybinning.begin(); it != mybinning.end(); ++it)
    binning[ib++] = *it;
}

/**
 * Create new Pt binning
 * @param binning
 */
void AliAnalysisTaskChargedParticlesRef::CreateNewPtBinning(TArrayD& binning) const {
  std::vector<double> mybinning;
  std::map<double,double> definitions;
  definitions.insert(std::pair<double, double>(1, 0.05));
  definitions.insert(std::pair<double, double>(2, 0.1));
  definitions.insert(std::pair<double, double>(4, 0.2));
  definitions.insert(std::pair<double, double>(7, 0.5));
  definitions.insert(std::pair<double, double>(16, 1));
  definitions.insert(std::pair<double, double>(36, 2));
  definitions.insert(std::pair<double, double>(40, 4));
  definitions.insert(std::pair<double, double>(50, 5));
  definitions.insert(std::pair<double, double>(100, 10));
  definitions.insert(std::pair<double, double>(200, 20));
  double currentval = 0.;
  mybinning.push_back(currentval);
  for(std::map<double,double>::iterator id = definitions.begin(); id != definitions.end(); ++id){
    double limit = id->first, binwidth = id->second;
    while(currentval < limit){
      currentval += binwidth;
      mybinning.push_back(currentval);
    }
  }
  binning.Set(mybinning.size());
  int ib = 0;
  for(std::vector<double>::iterator it = mybinning.begin(); it != mybinning.end(); ++it)
    binning[ib++] = *it;
}

/**
 * Run track selection for ESD tracks
 * @param track The track to check
 * @return True if the track is selected, false otherwise
 */
Bool_t AliAnalysisTaskChargedParticlesRef::TrackSelectionESD(AliESDtrack* track) {
  return fTrackCuts->AcceptTrack(track);
}

/**
 * Run track selection for AOD tracks
 * @param track The track to check
 * @return True if the track is selected, false otherwise
 */
Bool_t AliAnalysisTaskChargedParticlesRef::TrackSelectionAOD(AliAODTrack* track) {
  if(!track->TestFilterBit(AliAODTrack::kTrkGlobal)) return false;
  if(track->GetTPCNCrossedRows() < 120) return false;
  return true;
}

/**
 * Apply additional cut requiring at least one offline patch above a given energy (not fake ADC!)
 * Attention: This task groups into single shower triggers (L0, EG1, EG2) and jet triggers (EJ1 and EJ2).
 * Per convention the low threshold patch is selected. No energy cut should be applied in the trigger maker
 * @param trgcls Trigger class for which to apply additional offline patch selection
 * @param triggerpatches Array of trigger patches
 * @return True if at least on patch above threshold is found or no cut is applied
 */
Bool_t AliAnalysisTaskChargedParticlesRef::IsOfflineSelected(EmcalTriggerClass trgcls, const TClonesArray * const triggerpatches) const {
  if(fOfflineEnergyThreshold[trgcls] < 0) return true;
  bool isSingleShower = ((trgcls == kCPREL0) || (trgcls == kCPREG1) || (trgcls == kCPREG2));
  int nfound = 0;
  AliEMCALTriggerPatchInfo *patch = NULL;
  for(TIter patchIter = TIter(triggerpatches).Begin(); patchIter != TIter::End(); ++patchIter){
    patch = static_cast<AliEMCALTriggerPatchInfo *>(*patchIter);
    if(!patch->IsOfflineSimple()) continue;
    if(isSingleShower){
     if(!patch->IsGammaLowSimple()) continue;
    } else {
      if(!patch->IsJetLowSimple()) continue;
    }
    if(patch->GetPatchE() > fOfflineEnergyThreshold[trgcls]) nfound++;
  }
  return nfound > 0;
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

} /* namespace EMCalTriggerPtAnalysis */
