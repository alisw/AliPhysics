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
#include <memory>

#include <TArrayD.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <THashList.h>
#include <TH1.h>
#include <THistManager.h>
#include <TKey.h>
#include <TList.h>
#include <TPDGCode.h>
#include <TProfile.h>
#include <TMath.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>

#include "AliAnalysisManager.h"
#include "AliAnalysisUtils.h"
#include "AliAODMCHeader.h"
#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h"
#include "AliAODTrack.h"
#include "AliEmcalAnalysisFactory.h"
#include "AliEmcalTrackSelection.h"
#include "AliEmcalTriggerOfflineSelection.h"
#include "AliEMCALTriggerPatchInfo.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALRecoUtils.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliGenPythiaEventHeader.h"
#include "AliInputEventHandler.h"
#include "AliMCEvent.h"
#include "AliVVertex.h"

#include "AliAnalysisTaskChargedParticlesRefMC.h"
#include "AliEMCalTriggerWeightHandler.h"

/// \cond CLASSIMP
ClassImp(EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRefMC)
/// \endcond

namespace EMCalTriggerPtAnalysis {

/**
 * Dummy constructor
 */
AliAnalysisTaskChargedParticlesRefMC::AliAnalysisTaskChargedParticlesRefMC():
        AliAnalysisTaskSE(),
        fTrackCuts(NULL),
        fAnalysisUtil(NULL),
        fTriggerSelection(NULL),
        fHistos(NULL),
        fGeometry(NULL),
        fWeightHandler(NULL),
        fUsePythiaHard(kFALSE),
        fPtHard(0),
        fPtHardBin(0),
        fNTrials(0),
        fXsection(0),
        fYshift(0.465),
        fEtaSign(1),
        fEtaLabCut(-0.6, 0.6),
        fEtaCmsCut(-0.13, 0.13),
        fFracPtHard(-1)
{
}

/**
 * Main constructor
 * @param name Name of the task
 */
AliAnalysisTaskChargedParticlesRefMC::AliAnalysisTaskChargedParticlesRefMC(const char* name):
        AliAnalysisTaskSE(name),
        fTrackCuts(NULL),
        fAnalysisUtil(NULL),
        fTriggerSelection(NULL),
        fHistos(NULL),
        fGeometry(NULL),
        fWeightHandler(NULL),
        fUsePythiaHard(kFALSE),
        fPtHard(0),
        fPtHardBin(0),
        fNTrials(0),
        fXsection(0),
        fYshift(0.465),
        fEtaSign(1),
        fEtaLabCut(-0.6, 0.6),
        fEtaCmsCut(-0.13, 0.13),
        fFracPtHard(-1)
{
  DefineOutput(1, TList::Class());
}

/**
 * Destuctor
 */
AliAnalysisTaskChargedParticlesRefMC::~AliAnalysisTaskChargedParticlesRefMC() {
  //if(fTrackCuts) delete fTrackCuts;
  if(fAnalysisUtil) delete fAnalysisUtil;
  if(fTriggerSelection) delete fTriggerSelection;
  if(fHistos) delete fHistos;
}
/**
 * Create the output histograms
 */
void AliAnalysisTaskChargedParticlesRefMC::UserCreateOutputObjects() {
  if(!fAnalysisUtil) fAnalysisUtil = new AliAnalysisUtils;
  fHistos = new THistManager("Ref");

  if(!fTrackCuts) InitializeTrackCuts("standard",fInputHandler->IsA() == AliAODInputHandler::Class());

  PtBinning newbinning;

  fHistos->CreateTH1("hNtrials", "Number of trials", 1, 0.5, 1.5);
  fHistos->CreateTProfile("hCrossSection", "PYTHIA cross section", 1, 0.5, 1.5);
  fHistos->CreateTH1("hNtrialsNoSelect", "Number of trials (without event selection)", 1, 0.5, 1.5);
  fHistos->CreateTH1("hNtrialsEvent", "Number of trials (from header, after event selection)", 1, 0.5, 1.5);
  fHistos->CreateTProfile("hCrossSectionEvent", "PYTHIA cross section (from header, after event selection)", 1, 0.5, 1.5);
  fHistos->CreateTH1("hPtHard", "Pt of the hard interaction", 1000, 0., 500);
  fHistos->CreateTH1("hTriggerJetPtNoCut", "pt of trigger jets wihtout cuts", 1000, 0., 500);
  fHistos->CreateTH1("hRatioPtJetPtHardNoCut", "Ratio of pt jet / pt hard without cut", 1000, 0., 20.);
  fHistos->CreateTH1("hTriggerJetPtWithCut", "pt of trigger jets after cuts", 1000, 0., 500);
  fHistos->CreateTH1("hRatioPtJetPtHardWithCut", "Ratio of pt jet / pt hard with cut on this ratio", 1000, 0., 20.);
  TString triggers[7] = {"True", "MB", "EMC7", "EJ1", "EJ2", "EG1", "EG2"};
  Double_t ptcuts[5] = {1., 2., 5., 10., 20.};
  TString species[6] = {"El", "Mu", "Pi", "Ka", "Pr", "Ot"};
  for(TString *trg = triggers; trg < triggers + sizeof(triggers)/sizeof(TString); trg++){
    fHistos->CreateTH1(Form("hEventCount%s", trg->Data()), Form("Event Counter for trigger class %s", trg->Data()), 1, 0.5, 1.5);
    fHistos->CreateTH1(Form("hVertexBefore%s", trg->Data()), Form("Vertex distribution before z-cut for trigger class %s", trg->Data()), 500, -50, 50);
    fHistos->CreateTH1(Form("hVertexAfter%s", trg->Data()), Form("Vertex distribution after z-cut for trigger class %s", trg->Data()), 100, -10, 10);
    fHistos->CreateTH1(Form("hPtEtaAll%s", trg->Data()), Form("Charged particle pt distribution all eta trigger %s", trg->Data()), newbinning);
    fHistos->CreateTH1(Form("hPtEtaCent%s", trg->Data()), Form("Charged particle pt distribution central eta trigger %s", trg->Data()), newbinning);
    fHistos->CreateTH1(Form("hPtEMCALEtaAll%s", trg->Data()), Form("Charged particle in EMCAL pt distribution all eta trigger %s", trg->Data()), newbinning);
    fHistos->CreateTH1(Form("hPtEMCALEtaCent%s", trg->Data()), Form("Charged particle in EMCAL pt distribution central eta trigger %s", trg->Data()), newbinning);
    fHistos->CreateTH1(Form("hPtEMCALNoTRDEtaAll%s", trg->Data()), Form("Charged particle in EMCAL (no TRD in front) pt distribution all eta trigger %s", trg->Data()), newbinning);
    fHistos->CreateTH1(Form("hPtEMCALNoTRDEtaCent%s", trg->Data()), Form("Charged particle in EMCAL (no TRD in front) pt distribution central eta trigger %s", trg->Data()), newbinning);
    fHistos->CreateTH1(Form("hPtEMCALWithTRDEtaAll%s", trg->Data()), Form("Charged particle in EMCAL (with TRD in front) pt distribution all eta trigger %s", trg->Data()), newbinning);
    fHistos->CreateTH1(Form("hPtEMCALWithTRDEtaCent%s", trg->Data()), Form("Charged particle in EMCAL (with TRD in front) pt distribution central eta trigger %s", trg->Data()), newbinning);
    for(TString *piditer = species; piditer < species + sizeof(species)/sizeof(TString); ++piditer){
      fHistos->CreateTH1(Form("hPtEtaAll%s%s", piditer->Data(), trg->Data()), Form("Charged %s pt distribution all eta trigger %s", piditer->Data(), trg->Data()), newbinning);
      fHistos->CreateTH1(Form("hPtEtaCent%s%s", piditer->Data(), trg->Data()), Form("Charged %s pt distribution central eta trigger %s", piditer->Data(), trg->Data()), newbinning);
      fHistos->CreateTH1(Form("hPtEMCALEtaAll%s%s", piditer->Data(), trg->Data()), Form("Charged %s in EMCAL pt distribution all eta trigger %s", piditer->Data(), trg->Data()), newbinning);
      fHistos->CreateTH1(Form("hPtEMCALEtaCent%s%s", piditer->Data(), trg->Data()), Form("Charged %s in EMCAL pt distribution central eta trigger %s", piditer->Data(), trg->Data()), newbinning);
      fHistos->CreateTH1(Form("hPtEMCALNoTRDEtaAll%s%s", piditer->Data(), trg->Data()), Form("Charged %s in EMCAL (no TRD in front) pt distribution all eta trigger %s", piditer->Data(), trg->Data()), newbinning);
      fHistos->CreateTH1(Form("hPtEMCALNoTRDEtaCent%s%s", piditer->Data(), trg->Data()), Form("Charged %s in EMCAL (no TRD in front) pt distribution central eta trigger %s", piditer->Data(), trg->Data()), newbinning);
      fHistos->CreateTH1(Form("hPtEMCALWithTRDEtaAll%s%s", piditer->Data(), trg->Data()), Form("Charged %s in EMCAL (with TRD in front) pt distribution all eta trigger %s", piditer->Data(), trg->Data()), newbinning);
      fHistos->CreateTH1(Form("hPtEMCALWithTRDEtaCent%s%s", piditer->Data(), trg->Data()), Form("Charged %s in EMCAL (with TRD in front) pt distribution central eta trigger %s", piditer->Data(), trg->Data()), newbinning);
    }
    for(int ipt = 0; ipt < 5; ipt++){
      fHistos->CreateTH1(
          Form("hEtaLabDistAllPt%d%s", static_cast<Int_t>(ptcuts[ipt]), trg->Data()),
          Form("Eta (lab) distribution without etacut for tracks with Pt above %.1f GeV/c trigger %s", ptcuts[ipt], trg->Data()),
          100,
          -1.,
          1.
          );
      fHistos->CreateTH1(
          Form("hEtaLabDistCutPt%d%s", static_cast<Int_t>(ptcuts[ipt]), trg->Data()),
          Form("Eta (lab) distribution with etacut for tracks with Pt above %.1f GeV/c trigger %s", ptcuts[ipt], trg->Data()),
          100,
          -1.,
          1.
          );
      fHistos->CreateTH1(
          Form("hEtaCentDistAllPt%d%s", static_cast<Int_t>(ptcuts[ipt]), trg->Data()),
          Form("Eta (cent) distribution without etacut for tracks with Pt above %.1f GeV/c trigger %s", ptcuts[ipt], trg->Data()),
          160,
          -1.3,
          1.3
          );
      fHistos->CreateTH1(
          Form("hEtaCentDistCutPt%d%s", static_cast<Int_t>(ptcuts[ipt]), trg->Data()),
          Form("Eta (cent) distribution with etacut for tracks with Pt above %.1f GeV/c trigger %s", ptcuts[ipt], trg->Data()),
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
  //fHistos->GetListOfHistograms()->Add(fTrackCuts);
  PostData(1, fHistos->GetListOfHistograms());
}

void AliAnalysisTaskChargedParticlesRefMC::UserExec(Option_t*) {  // Select event
  if(MCEvent()) return;
  if(!fGeometry){
    fGeometry = AliEMCALGeometry::GetInstanceFromRunNumber(InputEvent()->GetRunNumber());
  }
  AliGenPythiaEventHeader *pyheader = GetPythiaHeader();
  if(pyheader){
    FillTriggerJetHistograms(kFALSE, pyheader);
    fHistos->FillTH1("hNtrialsNoSelect",1,pyheader->Trials());
    if(fFracPtHard > 0){
      // Apply outlier cut in case of a positive fraction of pthard
      if(IsOutlierJet(pyheader))
        return;
    }
    FillTriggerJetHistograms(kTRUE, pyheader);

    if(IsOutlierParticle(pyheader)) return;
  }

  Double_t weight = fWeightHandler ? fWeightHandler->GetEventWeight(pyheader) : 1.;

  //TString triggerstring = GetFiredTriggerClasses(fTriggerPatches);
  Bool_t isMinBias = fInputHandler->IsEventSelected() & AliVEvent::kINT7,
      isEMC7 = kFALSE, isEJ1 = kFALSE, isEJ2 = kFALSE, isEG1 = kFALSE, isEG2 = kFALSE;
  TClonesArray *fTriggerPatches = dynamic_cast<TClonesArray *>(fInputEvent->FindListObject("EmcalTriggers"));
  if(fTriggerPatches && fTriggerSelection){
      isEMC7 = isMinBias && fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgEL0, fTriggerPatches); // triggerstring.Contains("EMC7"),
      isEJ1 = isMinBias && fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgEJ1, fTriggerPatches); // triggerstring.Contains("EJ1"),
      isEJ2 = isMinBias && fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgEJ2, fTriggerPatches); // triggerstring.Contains("EJ2"),
      isEG1 = isMinBias && fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgEG1, fTriggerPatches); // triggerstring.Contains("EG1"),
      isEG2 = isMinBias && fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgEG2, fTriggerPatches); // triggerstring.Contains("EG2");               // In simulations triggered events are a subset of min. bias events

  }
  if(!(isMinBias || isEG1 || isEG2 || isEJ1 || isEJ2)) return;
  const AliVVertex *vtx = fInputEvent->GetPrimaryVertex();
  //if(!fInputEvent->IsPileupFromSPD(3, 0.8, 3., 2., 5.)) return;         // reject pileup event
  if(vtx->GetNContributors() < 1) return;
  // Fill reference distribution for the primary vertex before any z-cut
  fHistos->FillTH1("hVertexBeforeTrue", vtx->GetZ(), weight);
  if(isMinBias) fHistos->FillTH1("hVertexBeforeMB", vtx->GetZ(), weight);
  if(isEMC7) fHistos->FillTH1("hVertexBeforeEMC7", vtx->GetZ(), weight);
  if(isEJ1) fHistos->FillTH1("hVertexBeforeEJ1", vtx->GetZ(), weight);
  if(isEJ2) fHistos->FillTH1("hVertexBeforeEJ2", vtx->GetZ(), weight);
  if(isEG1) fHistos->FillTH1("hVertexBeforeEG1", vtx->GetZ(), weight);
  if(isEG2) fHistos->FillTH1("hVertexBeforeEG2", vtx->GetZ(), weight);
  if(!fAnalysisUtil->IsVertexSelected2013pA(fInputEvent)) return;       // Apply new vertex cut
  if(fAnalysisUtil->IsPileUpEvent(fInputEvent)) return;                 // Apply new vertex cut
  // Apply vertex z cut
  if(vtx->GetZ() < -10. || vtx->GetZ() > 10.) return;

  // Fill Event counter and reference vertex distributions for the different trigger classes
  fHistos->FillTH1("hEventCountTrue", 1, weight);
  fHistos->FillTH1("hVertexAfterTrue", vtx->GetZ(), weight);
  if(isMinBias){
    fHistos->FillTH1("hEventCountMB", 1, weight);
    fHistos->FillTH1("hVertexAfterMB", vtx->GetZ(), weight);
  }
  if(isEMC7){
    fHistos->FillTH1("hEventCountEMC7", 1, weight);
    fHistos->FillTH1("hVertexAfterEMC7", vtx->GetZ(), weight);
  }
  if(isEJ1){
    fHistos->FillTH1("hEventCountEJ1", 1, weight);
    fHistos->FillTH1("hVertexAfterEJ1", vtx->GetZ(), weight);
  }
  if(isEJ2){
    fHistos->FillTH1("hEventCountEJ2", 1, weight);
    fHistos->FillTH1("hVertexAfterEJ2", vtx->GetZ(), weight);
  }
  if(isEG1){
    fHistos->FillTH1("hEventCountEG1", 1, weight);
    fHistos->FillTH1("hVertexAfterEG1", vtx->GetZ(), weight);
  }
  if(isEG2){
    fHistos->FillTH1("hEventCountEG2", 1, weight);
    fHistos->FillTH1("hVertexAfterEG2", vtx->GetZ(), weight);
  }

  // Fill PYTHIA histograms from event header
  if(pyheader){
   fHistos->FillTH1("hNtrialsEvent", 1., pyheader->Trials());
   fHistos->FillProfile("hCrossSectionEvent", 1., pyheader->GetXsection());
   fHistos->FillTH1("hPtHard", pyheader->GetPtHard());
  }

  // MonteCarlo Loop
  // Histograms
  // - Full eta (-0.8, 0.8), new binning
  // - Eta distribution for tracks above 1, 2, 5, 10 GeV/c without eta cut
  // - Central eta_{cms} (-0.3, 0.3), new binning,
  // - Eta distribution for tracks above 1, 2, 5, 10 GeV/c
  // - Eta distribution for tracks above 1, 2, 5, 10 GeV/c with eta cut
  AliVParticle *truepart = NULL;
  Bool_t isEMCAL(kFALSE);
  for(int ipart = 0; ipart < fMCEvent->GetNumberOfTracks(); ipart++){
    truepart = fMCEvent->GetTrack(ipart);

    // Select only particles within ALICE acceptance
    if(!fEtaLabCut.IsInRange(truepart->Eta())) continue;
    if(TMath::Abs(truepart->Pt()) < 0.1) continue;
    if(!truepart->Charge()) continue;

    if(!IsPhysicalPrimary(truepart, fMCEvent)) continue;
    isEMCAL = (truepart->Phi() > 1.5 && truepart->Phi() < 3.1) ? kTRUE : kFALSE;

    // Calculate eta in cms frame according
    // EPJC74 (2014) 3054:
    // eta_cms = - eta_lab - |yshift|
    Double_t etacent = -1. * truepart->Eta() - TMath::Abs(fYshift);
    etacent *= fEtaSign;

    Bool_t etacentcut = fEtaCmsCut.IsInRange(etacent);

    // Get PID
    TString pid = "";
    switch(TMath::Abs(truepart->PdgCode())){
    case kPiPlus: pid = "Pi"; break;
    case kMuonMinus: pid = "Mu"; break;
    case kElectron: pid = "El"; break;
    case kKPlus: pid = "Ka"; break;
    case kProton: pid = "Pr"; break;
    default: pid = "Ot"; break;
    };

    // Particle selected (do not filter TRD sectors for MC truth)
    FillTrackHistos("True", weight, truepart->Pt(), truepart->Eta() * fEtaSign, etacent, truepart->Phi(), etacentcut, isEMCAL, kFALSE, pid);
  }

  // Loop over tracks, fill select particles
  // Histograms
  // - Full eta (-0.8, 0.8), new binning
  // - Eta distribution for tracks above 1, 2, 5, 10 GeV/c without eta cut
  // - Central eta (-0.8, -0.2), new binning,
  // - Eta distribution for tracks above 1, 2, 5, 10 GeV/c
  // - Eta distribution for tracks above 1, 2, 5, 10 GeV/c with eta cut
  AliVTrack *checktrack(NULL);
  AliVParticle *assocMC(NULL);
  double ptparticle(-1.), etaparticle(-100.), etaEMCAL(0.), phiEMCAL(0.);
  Bool_t hasTRD = kFALSE;
  for(int itrk = 0; itrk < fInputEvent->GetNumberOfTracks(); ++itrk){
    checktrack = dynamic_cast<AliVTrack *>(fInputEvent->GetTrack(itrk));
    if(!checktrack) continue;
    // Find associated particle
    assocMC = fMCEvent->GetTrack(TMath::Abs(checktrack->GetLabel()));
    if(!assocMC) continue;        // Fake track
    if(!IsPhysicalPrimary(assocMC, fMCEvent)) continue;

    // Select only particles within ALICE acceptance
    if(!fEtaLabCut.IsInRange(checktrack->Eta())) continue;
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

    if(!fTrackCuts->IsTrackAccepted(checktrack)) continue;

    ptparticle = TMath::Abs(assocMC->Pt());
    etaparticle = assocMC->Eta();

    // Calculate eta in cms frame according
    // EPJC74 (2014) 3054:
    // eta_cms = - eta_lab - |yshift|
    Double_t etacent = -1. * checktrack->Eta() - TMath::Abs(fYshift);
    etacent *= fEtaSign;

    Bool_t etacentcut = fEtaCmsCut.IsInRange(etacent);

    // Get PID
    TString assocpid = "";
    switch(TMath::Abs(assocMC->PdgCode())){
    case kPiPlus: assocpid = "Pi"; break;
    case kMuonMinus: assocpid = "Mu"; break;
    case kElectron: assocpid = "El"; break;
    case kKPlus: assocpid = "Ka"; break;
    case kProton: assocpid = "Pr"; break;
    default: assocpid = "Ot"; break;
    };

    // Go through the trigger classes and fill histograms
    if(isMinBias){
      FillTrackHistos("MB", weight,  ptparticle, checktrack->Eta() * fEtaSign, etacent, checktrack->Phi(), etacentcut, isEMCAL, hasTRD, assocpid);
    }
    if(isEMC7){
      FillTrackHistos("EMC7", weight, ptparticle, checktrack->Eta() * fEtaSign, etacent, checktrack->Phi(), etacentcut, isEMCAL, hasTRD, assocpid);
    }
    if(isEJ1){
      FillTrackHistos("EJ1", weight, ptparticle, checktrack->Eta() * fEtaSign, etacent, checktrack->Phi(), etacentcut, isEMCAL, hasTRD, assocpid);
    }
    if(isEJ2){
      FillTrackHistos("EJ2", weight, ptparticle, checktrack->Eta() * fEtaSign, etacent, checktrack->Phi(), etacentcut, isEMCAL, hasTRD, assocpid);
    }
    if(isEG1){
      FillTrackHistos("EG1", weight, ptparticle, checktrack->Eta() * fEtaSign, etacent, checktrack->Phi(), etacentcut, isEMCAL, hasTRD, assocpid);
    }
    if(isEG2){
      FillTrackHistos("EG2", weight, ptparticle, checktrack->Eta() * fEtaSign, etacent, checktrack->Phi(), etacentcut, isEMCAL, hasTRD, assocpid);
    }
  }
  PostData(1, fHistos->GetListOfHistograms());
}

/**
 * Fill track histograms
 * @param eventclass Trigger class fired
 * @param weight \f$ p_{t} \f$-hard dependent weight
 * @param pt track \f$ p_{t} \f$
 * @param etalab Track \f$ \eta \f$ in lab frame
 * @param etacent Track \f$ \eta \f$ in cms frame
 * @param phi Track \f$ \eta \f$ in lab frame
 * @param etacut Track accepted by \f$ \eta \f$ cut
 * @param inEmcal Track in EMCAL \f$ \phi \f$ acceptance
 */
void AliAnalysisTaskChargedParticlesRefMC::FillTrackHistos(
    const char *eventclass,
    Double_t weight,
    Double_t pt,
    Double_t etalab,
    Double_t etacent,
    Double_t phi,
    Bool_t etacut,
    Bool_t inEmcal,
    Bool_t hasTRD,
    const char *pid
    )
{
  fHistos->FillTH1(Form("hPtEtaAll%s", eventclass), TMath::Abs(pt), weight);
  fHistos->FillTH1(Form("hPtEtaAll%s%s", pid, eventclass), TMath::Abs(pt), weight);
  if(inEmcal){
    fHistos->FillTH1(Form("hPtEMCALEtaAll%s", eventclass), TMath::Abs(pt), weight);
    fHistos->FillTH1(Form("hPtEMCALEtaAll%s%s", pid, eventclass), TMath::Abs(pt), weight);
    if(hasTRD){
      fHistos->FillTH1(Form("hPtEMCALWithTRDEtaAll%s", eventclass), TMath::Abs(pt), weight);
      fHistos->FillTH1(Form("hPtEMCALWithTRDEtaAll%s%s", pid, eventclass), TMath::Abs(pt), weight);
    } else {
      fHistos->FillTH1(Form("hPtEMCALNoTRDEtaAll%s", eventclass), TMath::Abs(pt), weight);
      fHistos->FillTH1(Form("hPtEMCALNoTRDEtaAll%s%s", pid, eventclass), TMath::Abs(pt), weight);
    }
  }

  int ptmin[5] = {1,2,5,10,20}; // for eta distributions
  for(int icut = 0; icut < 5; icut++){
    if(TMath::Abs(pt) > static_cast<double>(ptmin[icut])){
      fHistos->FillTH1(Form("hPhiDistAllPt%d%s", ptmin[icut], eventclass), phi, weight);
      fHistos->FillTH1(Form("hEtaLabDistAllPt%d%s", ptmin[icut], eventclass), etalab, weight);
      fHistos->FillTH1(Form("hEtaCentDistAllPt%d%s", ptmin[icut], eventclass), etacent, weight);
      if(inEmcal){
        fHistos->FillTH1(Form("hEtaLabDistAllEMCALPt%d%s", ptmin[icut], eventclass), etalab, weight);
        fHistos->FillTH1(Form("hEtaCentDistAllEMCALPt%d%s", ptmin[icut], eventclass), etacent, weight);
      }
    }
  }

  if(etacut){
    fHistos->FillTH1(Form("hPtEtaCent%s", eventclass), TMath::Abs(pt), weight);
    fHistos->FillTH1(Form("hPtEtaCent%s%s", pid, eventclass), TMath::Abs(pt), weight);
    if(inEmcal){
      fHistos->FillTH1(Form("hPtEMCALEtaCent%s", eventclass), TMath::Abs(pt), weight);
      fHistos->FillTH1(Form("hPtEMCALEtaCent%s%s", pid, eventclass), TMath::Abs(pt), weight);
      if(hasTRD){
        fHistos->FillTH1(Form("hPtEMCALWithTRDEtaCent%s", eventclass), TMath::Abs(pt), weight);
        fHistos->FillTH1(Form("hPtEMCALWithTRDEtaCent%s%s", pid, eventclass), TMath::Abs(pt), weight);
      } else {
        fHistos->FillTH1(Form("hPtEMCALNoTRDEtaCent%s", eventclass), TMath::Abs(pt), weight);
        fHistos->FillTH1(Form("hPtEMCALNoTRDEtaCent%s%s", pid, eventclass), TMath::Abs(pt), weight);
      }
    }
    for(int icut = 0; icut < 5; icut++){
      if(TMath::Abs(pt) > static_cast<double>(ptmin[icut])){
        fHistos->FillTH1(Form("hEtaLabDistCutPt%d%s", ptmin[icut], eventclass), etalab, weight);
        fHistos->FillTH1(Form("hEtaCentDistCutPt%d%s", ptmin[icut], eventclass), etacent, weight);
        if(inEmcal){
          fHistos->FillTH1(Form("hEtaLabDistCutEMCALPt%d%s", ptmin[icut], eventclass), etalab, weight);
          fHistos->FillTH1(Form("hEtaCentDistCutEMCALPt%d%s", ptmin[icut], eventclass), etacent, weight);
        }
      }
    }
  }
}

/**
 * Fill histogram with pt of all triggering jets
 * @param histname Name of the histogram to fill
 * @param header PYTHIA event header with jet information
 */
void AliAnalysisTaskChargedParticlesRefMC::FillTriggerJetHistograms(Bool_t aftercut, AliGenPythiaEventHeader *const header){
  TString ending = aftercut ? "WithCut" : "NoCut";
  Float_t pbuf[4];
  TLorentzVector jetvec;
  TString hnameSpec = "hTriggerJetPt" + ending,
          hnameptratio = "hRatioPtJetPtHard" + ending;
  for(int ijet = 0; ijet < header->NTriggerJets(); ijet++){
    memset(pbuf, 0, sizeof(Float_t) * 4);
    header->TriggerJet(ijet, pbuf);
    jetvec.SetPxPyPzE(pbuf[0], pbuf[1], pbuf[2], pbuf[3]);
    fHistos->FillTH1(hnameSpec.Data(), TMath::Abs(jetvec.Pt()));
    fHistos->FillTH1(hnameptratio.Data(), TMath::Abs(jetvec.Pt()/header->GetPtHard()));
  }
}

/**
 * Perform actions when giles change
 * @return
 */
Bool_t AliAnalysisTaskChargedParticlesRefMC::UserNotify()
{
  // Called when file changes.
  if(!fUsePythiaHard) return kTRUE;

  TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();
  if (!tree) {
    AliErrorStream() << GetName() << " - UserNotify: No current tree!" << std::endl;
    return kFALSE;
  }

  Float_t xsection    = 0;
  Float_t trials      = 0;
  Int_t   pthardbin   = 0;

  TFile *curfile = tree->GetCurrentFile();
  if (!curfile) {
    AliErrorStream() << GetName() << " - UserNotify: No current file!" << std::endl;
    return kFALSE;
  }

  TChain *chain = dynamic_cast<TChain*>(tree);
  if (chain) tree = chain->GetTree();

  Int_t nevents = tree->GetEntriesFast();

  PythiaInfoFromFile(curfile->GetName(), xsection, trials, pthardbin);

  fHistos->FillTH1("hNtrials", 1., trials);
  fHistos->FillProfile("hCrossSection", 1., xsection);
  //fHistos->FillTH1(pthardbin, nevents);

  return kTRUE;
}

/**
 * Get the cross section and the trails either from pyxsec.root or from pysec_hists.root
 * Get the pt hard bin from the file path
 *  This is to called in Notify and should provide the path to the AOD/ESD file
 * (Partially copied from AliAnalysisHelperJetTasks)
 * From AliAnalysisTaskEmcal
 * @param currFile File name with PYTHIA hard cross section
 * @param fXsec Output storage for the cross section
 * @param fTrials Output storage for the number of trials
 * @param pthard Output storage of the pthardbin
 * @return True if reading was successful, false in case of errors
 */
Bool_t AliAnalysisTaskChargedParticlesRefMC::PythiaInfoFromFile(const char* currFile, Float_t &fXsec, Float_t &fTrials, Int_t &pthard) const {

  TString file(currFile);
  fXsec = 0;
  fTrials = 1;

  if (file.Contains(".zip#")) {
    Ssiz_t pos1 = file.Index("root_archive",12,0,TString::kExact);
    Ssiz_t pos = file.Index("#",1,pos1,TString::kExact);
    Ssiz_t pos2 = file.Index(".root",5,TString::kExact);
    file.Replace(pos+1,pos2-pos1,"");
  } else {
    // not an archive take the basename....
    file.ReplaceAll(gSystem->BaseName(file.Data()),"");
  }
  AliDebugStream(1) << "File name: " << file << std::endl;

  // Get the pt hard bin
  TString strPthard(file);

  strPthard.Remove(strPthard.Last('/'));
  strPthard.Remove(strPthard.Last('/'));
  if (strPthard.Contains("AOD")) strPthard.Remove(strPthard.Last('/'));
  strPthard.Remove(0,strPthard.Last('/')+1);
  if (strPthard.IsDec())
    pthard = strPthard.Atoi();
  else
    AliWarningStream() << "Could not extract file number from path " << strPthard << std::endl;

  // problem that we cannot really test the existance of a file in a archive so we have to live with open error message from root
  std::unique_ptr<TFile> fxsec(TFile::Open(Form("%s%s",file.Data(),"pyxsec.root")));

  if (!fxsec) {
    // next trial fetch the histgram file
    fxsec = std::unique_ptr<TFile>(TFile::Open(Form("%s%s",file.Data(),"pyxsec_hists.root")));
    if (!fxsec) return kFALSE;// not a severe condition but inciate that we have no information
    else {
      // find the tlist we want to be independtent of the name so use the Tkey
      TKey* key = (TKey*)fxsec->GetListOfKeys()->At(0);
      if (!key) return kFALSE;
      TList *list = dynamic_cast<TList*>(key->ReadObj());
      if (!list) return kFALSE;
      fXsec = ((TProfile*)list->FindObject("h1Xsec"))->GetBinContent(1);
      fTrials  = ((TH1F*)list->FindObject("h1Trials"))->GetBinContent(1);
    }
  } else { // no tree pyxsec.root
    TTree *xtree = (TTree*)fxsec->Get("Xsection");
    if (!xtree) return kFALSE;
    UInt_t   ntrials  = 0;
    Double_t  xsection  = 0;
    xtree->SetBranchAddress("xsection",&xsection);
    xtree->SetBranchAddress("ntrials",&ntrials);
    xtree->GetEntry(0);
    fTrials = ntrials;
    fXsec = xsection;
  }
  return kTRUE;
}

/**
 * Access PYTHIA event header
 * @return pythia event header (if existing)
 */
AliGenPythiaEventHeader *AliAnalysisTaskChargedParticlesRefMC::GetPythiaHeader() const {
  AliGenPythiaEventHeader *pythiaHeader = dynamic_cast<AliGenPythiaEventHeader*>(MCEvent()->GenEventHeader());
  if (!pythiaHeader) {
    // Check if AOD
    AliAODMCHeader* aodMCH = dynamic_cast<AliAODMCHeader*>(InputEvent()->FindListObject(AliAODMCHeader::StdBranchName()));
    if (aodMCH) {
      for (UInt_t i = 0;i<aodMCH->GetNCocktailHeaders();i++) {
        pythiaHeader = dynamic_cast<AliGenPythiaEventHeader*>(aodMCH->GetCocktailHeader(i));
        if (pythiaHeader) break;
      }
    }
  }
  return pythiaHeader;
}

/**
 * Set the track selection
 * @param cutname Name of the track cuts
 * @param isAOD check whether we run on ESDs or AODs
 */
void AliAnalysisTaskChargedParticlesRefMC::InitializeTrackCuts(TString cutname, bool isAOD){
  SetTrackSelection(AliEmcalAnalysisFactory::TrackCutsFactory(cutname, isAOD));
}



/**
 * Apply trigger selection using offline patches and trigger thresholds based on offline ADC Amplitude
 * @param triggerpatches Trigger patches found by the trigger maker
 * @return String with EMCAL trigger decision
 */
TString AliAnalysisTaskChargedParticlesRefMC::GetFiredTriggerClasses(const TClonesArray* triggerpatches) {
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
 * Check in a transparent way for ESDs and AODs whether the particle is physical primary or not
 * -# AOD: Information stored in the AliAODMCParticle
 * -# ESD: Information needs to be retrieved from the stack via the label of the MC particle
 * @param part The particle to check
 * @param mcevent The MC event containing the stack (ESD only)
 * @return True if particle is a physical primary particle, false otherwise
 */
Bool_t AliAnalysisTaskChargedParticlesRefMC::IsPhysicalPrimary(const AliVParticle* const part, AliMCEvent* const mcevent) {
  Bool_t physprim = false;
  const AliAODMCParticle *aodmc = dynamic_cast<const AliAODMCParticle *>(part);
  if(aodmc){
    physprim = aodmc->IsPhysicalPrimary();
  } else {
    physprim = mcevent->IsPhysicalPrimary(part->GetLabel());
  }
  return physprim;
}

/**
 * Find outlier jets compared to the pt hard
 * @param header PYTHIA header with trigger jets and pt hard
 * @return True if event has at least one outlier, false otherwise
 */
Bool_t AliAnalysisTaskChargedParticlesRefMC::IsOutlierJet(AliGenPythiaEventHeader * const header) const {
  Bool_t hasOutlier = kFALSE;
  Float_t pbuf[4];
  TLorentzVector jetvec;
  for(int ijet = 0; ijet < header->NTriggerJets(); ijet++){
    memset(pbuf, 0, sizeof(Float_t) * 4);
    header->TriggerJet(ijet, pbuf);
    jetvec.SetPxPyPzE(pbuf[0], pbuf[1], pbuf[2], pbuf[3]);
    if(TMath::Abs(jetvec.Pt()) >= this->fFracPtHard * header->GetPtHard()){
      hasOutlier = true;
      break;
    }
  }
  return hasOutlier;
}

/**
 * Find outlier events in second category: Events that contain a particle
 * with a higher pt than the pthard
 * @return True if event has at least one outlier, false otherwise
 */
Bool_t AliAnalysisTaskChargedParticlesRefMC::IsOutlierParticle(AliGenPythiaEventHeader *const header) const {
  Bool_t hasOutlier = false;
  const double pthard = header->GetPtHard();
  const AliVParticle  *part(nullptr);
  for(int ipart = 0; ipart < MCEvent()->GetNumberOfTracks(); ipart++){
    part = MCEvent()->GetTrack(ipart);
    if(TMath::Abs(part->Pt() > pthard)){
      hasOutlier = true;
      break;
    }
  }
  return hasOutlier;
}

/**
 * Create \f$ p_{t} \f$ binning
 */
AliAnalysisTaskChargedParticlesRefMC::PtBinning::PtBinning() :
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
