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
        AliAnalysisTaskEmcal(),
        fTrackCuts(nullptr),
        fTriggerSelection(nullptr),
        fHistos(nullptr),
        fWeightHandler(nullptr),
        fEventTriggers(),
        fEventWeight(1.),
        fYshift(0.465),
        fEtaSign(1),
        fEtaLabCut(-0.6, 0.6),
        fEtaCmsCut(-0.13, 0.13),
        fFracPtHard(-1)
{
  SetCaloTriggerPatchInfoName("EmcalTriggers");
  SetNeedEmcalGeom(true);
}

/**
 * Main constructor
 * @param name Name of the task
 */
AliAnalysisTaskChargedParticlesRefMC::AliAnalysisTaskChargedParticlesRefMC(const char* name):
        AliAnalysisTaskEmcal(name, true),
        fTrackCuts(nullptr),
        fTriggerSelection(nullptr),
        fHistos(nullptr),
        fWeightHandler(nullptr),
        fEventTriggers(),
        fEventWeight(1.),
        fYshift(0.465),
        fEtaSign(1),
        fEtaLabCut(-0.6, 0.6),
        fEtaCmsCut(-0.13, 0.13),
        fFracPtHard(-1)
{
  SetCaloTriggerPatchInfoName("EmcalTriggers");
  SetNeedEmcalGeom(true);
}

/**
 * Destuctor
 */
AliAnalysisTaskChargedParticlesRefMC::~AliAnalysisTaskChargedParticlesRefMC() {
  //if(fTrackCuts) delete fTrackCuts;
  if(fTriggerSelection) delete fTriggerSelection;
  if(fHistos) delete fHistos;
}
/**
 * Create the output histograms
 */
void AliAnalysisTaskChargedParticlesRefMC::UserCreateOutputObjects() {
  if(!fAliAnalysisUtils) fAliAnalysisUtils = new AliAnalysisUtils;
  fHistos = new THistManager("Ref");

  if(!fTrackCuts) InitializeTrackCuts("standard",fInputHandler->IsA() == AliAODInputHandler::Class());

  PtBinning newbinning;

  fHistos->CreateTH1("hPtHard", "Pt of the hard interaction", 1000, 0., 500);
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
  for(auto hist : *(fHistos->GetListOfHistograms())){
    fOutput->Add(hist);
  }
}

bool AliAnalysisTaskChargedParticlesRefMC::IsEventSelected(){
  AliDebugStream(1) << GetName() << ": Using custom event selection" << std::endl;
  if(!MCEvent()) return false;
  if(!fTriggerPatchInfo) return false;
  fEventWeight = fWeightHandler ? fWeightHandler->GetEventWeight(fPythiaHeader) : 1.;

  // Do MC outlier cut
  if(fIsPythia){
    if(!CheckMCOutliers()) return false;
  }


  // select trigger
  fEventTriggers.clear();
  bool isMinBias;
  if((isMinBias = fInputHandler->IsEventSelected() & AliVEvent::kINT7)) fEventTriggers.push_back("MB");
  // In simulations triggered events are a subset of min. bias events
  if(fTriggerPatchInfo && fTriggerSelection){
    if(isMinBias && fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgEL0, fTriggerPatchInfo))
      fEventTriggers.push_back("EMC7"); // triggerstring.Contains("EMC7"),
    if(isMinBias && fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgEJ1, fTriggerPatchInfo))
      fEventTriggers.push_back("EJ1"); // triggerstring.Contains("EJ1"),
    if(isMinBias && fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgEJ2, fTriggerPatchInfo))
      fEventTriggers.push_back("EJ2"); // triggerstring.Contains("EJ2"),
    if(isMinBias && fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgEG1, fTriggerPatchInfo))
      fEventTriggers.push_back("EG1"); // triggerstring.Contains("EG1"),
    if(isMinBias && fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgEG2, fTriggerPatchInfo))
      fEventTriggers.push_back("EG2"); // triggerstring.Contains("EG2");
  }
  if(!fEventTriggers.size()){
    AliDebugStream(1) << GetName() << ": No trigger selected" << std::endl;
  }

  const AliVVertex *vtx = fInputEvent->GetPrimaryVertex();
  if(vtx->GetNContributors() < 1) return false;
  if(!fAliAnalysisUtils->IsVertexSelected2013pA(fInputEvent)) return false;       // Apply new vertex cut
  if(fAliAnalysisUtils->IsPileUpEvent(fInputEvent)) return false;                 // Apply new vertex cut
  // Fill reference distribution for the primary vertex before any z-cut
  fHistos->FillTH1("hVertexBeforeTrue", vtx->GetZ(), fEventWeight);
  for(const auto &trg : fEventTriggers) fHistos->FillTH1(Form("hVertexBefore%s", trg.c_str()), vtx->GetZ(), fEventWeight);
  // Apply vertex z cut
  if(vtx->GetZ() < -10. || vtx->GetZ() > 10.) return false;


  // Fill Event counter and reference vertex distributions for the different trigger classes
  fHistos->FillTH1("hEventCountTrue", 1, fEventWeight);
  fHistos->FillTH1("hVertexAfterTrue", vtx->GetZ(), fEventWeight);
  for(const auto &trg : fEventTriggers){
      fHistos->FillTH1(Form("hEventCount%s", trg.c_str()), 1, fEventWeight);
      fHistos->FillTH1(Form("hVertexAfter%s", trg.c_str()), vtx->GetZ(), fEventWeight);
  }

  return true;
}

bool AliAnalysisTaskChargedParticlesRefMC::Run() {  // Select event

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
    FillTrackHistos("True", fEventWeight, truepart->Pt(), truepart->Eta() * fEtaSign, etacent, truepart->Phi(), etacentcut, isEMCAL, kFALSE, pid);
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
    isEMCAL = fGeom->SuperModuleNumberFromEtaPhi(etaEMCAL, phiEMCAL, supermoduleID);
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
    for(const auto &trg : fEventTriggers)
      FillTrackHistos(trg.c_str(), fEventWeight,  ptparticle, checktrack->Eta() * fEtaSign, etacent, checktrack->Phi(), etacentcut, isEMCAL, hasTRD, assocpid);
  }
  return true;
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
