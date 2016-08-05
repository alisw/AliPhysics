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
#include <memory>
#include <string>

#include <TClonesArray.h>
#include <TFile.h>
#include <THistManager.h>
#include <TMath.h>
#include <TPDGCode.h>
#include <TTree.h>

#include "AliAnalysisUtils.h"
#include "AliAODInputHandler.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODTrack.h"
#include "AliEmcalAnalysisFactory.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALRecoUtils.h"
#include "AliEMCALTriggerPatchInfo.h"
#include "AliEMCalTriggerWeightHandler.h"
#include "AliESDtrack.h"
#include "AliGenPythiaEventHeader.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliVEvent.h"
#include "AliVParticle.h"

#include "AliAnalysisTaskChargedParticlesMCTriggerMimic.h"

/// \cond CLASSIMP
ClassImp(EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesMCTriggerMimic)
/// \endcond

namespace EMCalTriggerPtAnalysis {

/**
 * Dummy constructor
 */
AliAnalysisTaskChargedParticlesMCTriggerMimic::AliAnalysisTaskChargedParticlesMCTriggerMimic():
            AliAnalysisTaskEmcal(),
            fTrackCuts(NULL),
            fHistos(NULL),
            fWeightHandler(NULL),
            fYshift(0.465),
            fEtaSign(1),
            fEtaLabCut(-0.6, 0.6),
            fEtaCmsCut(-0.13, 0.13),
            fPatchType(kUndef),
            fEnergyThreshold(0.)
{
  SetCaloTriggerPatchInfoName("EmcalTriggers");
  SetNeedEmcalGeom(true);
}

/**
 * Main constructor
 * @param name Name of the task
 */
AliAnalysisTaskChargedParticlesMCTriggerMimic::AliAnalysisTaskChargedParticlesMCTriggerMimic(const char *name):
            AliAnalysisTaskEmcal(name, true),
            fTrackCuts(nullptr),
            fHistos(nullptr),
            fWeightHandler(nullptr),
            fYshift(0.465),
            fEtaSign(1),
            fPatchType(kUndef),
            fEnergyThreshold(0.)
{
  SetCaloTriggerPatchInfoName("EmcalTriggers");
  SetNeedEmcalGeom(true);
}

/**
 * Destructor - cleaning up
 */
AliAnalysisTaskChargedParticlesMCTriggerMimic::~AliAnalysisTaskChargedParticlesMCTriggerMimic() {
  //if(fTrackCuts) delete fTrackCuts;
  if(fHistos) delete fHistos;
}

void AliAnalysisTaskChargedParticlesMCTriggerMimic::UserCreateOutputObjects(){
  AliAnalysisTaskEmcal::UserCreateOutputObjects();

  if(!fAliAnalysisUtils) fAliAnalysisUtils = new AliAnalysisUtils;
  if(!fTrackCuts) InitializeTrackCuts("standard",fInputHandler->IsA() == AliAODInputHandler::Class());

  PtBinning newbinning;
  fHistos = new THistManager("Ref");

  // The first two histograms should lead to a re-implementation inside AliAnalysisTaskEMCAL:
  // Users should be able to set event weights.
  fHistos->CreateTH1("hUserEventCount", "User event counter", 1, 0.5, 1.5);
  fHistos->CreateTH1("hUserVertexZ", "User vertex distribution after z-cut", 100, -10, 10);
  const std::array<std::string, 2> kInputs = {"True", "Accept"};
  const std::array<std::string, 6> kSpecies = {"El", "Mu", "Pi", "Ka", "Pr", "Ot"};
  const std::array<double, 5> kPtCuts = {1., 2., 5., 10., 20.};
  for(const auto &input : kInputs){
    fHistos->CreateTH1(Form("hPtEtaAll%s", input.c_str()), Form("Charged particle p_{t} distribution all #eta %s", input.c_str()), newbinning, "s");
    fHistos->CreateTH1(Form("hPtEtaCent%s", input.c_str()), Form("Charged particle p_{t} distribution central #eta  %s", input.c_str()), newbinning, "s");
    fHistos->CreateTH1(Form("hPtEMCALEtaAll%s", input.c_str()), Form("Charged particle in EMCAL p_{t} distribution all #eta trigger %s", input.c_str()), newbinning);
    fHistos->CreateTH1(Form("hPtEMCALEtaCent%s", input.c_str()), Form("Charged particle in EMCAL p_{t} distribution central eta trigger %s", input.c_str()), newbinning);
    for(const auto &piditer : kSpecies){
      fHistos->CreateTH1(Form("hPtEtaAll%s%s", piditer.c_str(), input.c_str()), Form("Charged %s p_{t} distribution all #eta %s", piditer.c_str(), input.c_str()), newbinning);
      fHistos->CreateTH1(Form("hPtEtaCent%s%s", piditer.c_str(), input.c_str()), Form("Charged %s p_{t} distribution central #eta %s", piditer.c_str(), input.c_str()), newbinning);
      fHistos->CreateTH1(Form("hPtEMCALEtaAll%s%s", piditer.c_str(), input.c_str()), Form("Charged %s in EMCAL p_{t} distribution all #eta %s", piditer.c_str(), input.c_str()), newbinning);
      fHistos->CreateTH1(Form("hPtEMCALEtaCent%s%s", piditer.c_str(), input.c_str()), Form("Charged %s in EMCAL p_{t} distribution central #eta %s", piditer.c_str(), input.c_str()), newbinning);
    }
    for(const auto &ptcut : kPtCuts){
      fHistos->CreateTH1(
        Form("hEtaLabDistAllPt%d%s", static_cast<Int_t>(ptcut), input.c_str()),
        Form("#eta_{lab} distribution without #eta-cut for tracks with p_{t} above %.1f GeV/c %s", ptcut, input.c_str()),
        100,
        -1.,
        1.
        );
      fHistos->CreateTH1(
        Form("hEtaLabDistCutPt%d%s", static_cast<Int_t>(ptcut), input.c_str()),
        Form("#eta_{lab} distribution with #eta-cut for tracks with p_{t} above %.1f GeV/c %s", ptcut, input.c_str()),
        100,
        -1.,
        1.
        );
      fHistos->CreateTH1(
        Form("hEtaCentDistAllPt%d%s", static_cast<Int_t>(ptcut), input.c_str()),
        Form("#eta_{cent} distribution without #eta-cut for tracks with p_{t} above %.1f GeV/c %s", ptcut, input.c_str()),
        160,
        -1.3,
        1.3
        );
    fHistos->CreateTH1(
        Form("hEtaCentDistCutPt%d%s", static_cast<Int_t>(ptcut), input.c_str()),
        Form("#eta_{cent} distribution with #eta-cut for tracks with p_{t} above %.1f GeV/c %s", ptcut, input.c_str()),
        160,
        -1.3,
        1.3
        );
    fHistos->CreateTH1(
        Form("hEtaLabDistAllEMCALPt%d%s", static_cast<Int_t>(ptcut), input.c_str()),
        Form("#eta_{lab} distribution without #eta-cut for tracks in EMCAL with p_{t} above %.1f GeV/c %s", ptcut, input.c_str()),
        100,
        -1.,
        1.
        );
    fHistos->CreateTH1(
        Form("hEtaLabDistCutEMCALPt%d%s", static_cast<Int_t>(ptcut), input.c_str()),
        Form("#eta_{lab} distribution with #eta-cut for tracks in EMCAL with p_{t} above %.1f GeV/c %s", ptcut, input.c_str()),
        100,
        -1.,
        1.
        );
    fHistos->CreateTH1(
        Form("hEtaCentDistAllEMCALPt%d%s", static_cast<Int_t>(ptcut), input.c_str()),
        Form("#eta_{cent} distribution without #eta-cut for tracks in EMCAL with p_{t} above %.1f GeV/c %s", ptcut, input.c_str()),
        160,
        -1.3,
        1.3
        );
    fHistos->CreateTH1(
        Form("hEtaCentDistCutEMCALPt%d%s", static_cast<Int_t>(ptcut), input.c_str()),
        Form("Eta (cent) distribution with #eta-cut for tracks in EMCAL with p_{t} above %.1f GeV/c %s", ptcut, input.c_str()),
        160,
        -1.3,
        1.3
        );
    fHistos->CreateTH1(
        Form("hPhiDistAllPt%d%s", static_cast<Int_t>(ptcut), input.c_str()),
        Form("#phi distribution of particles with p_{t} above %.1f GeV/c trigger %s", ptcut, input.c_str()),
        300,
        0.,
        2*TMath::Pi()
        );
  }
  }


  for(auto hist : *(fHistos->GetListOfHistograms())) fOutput->Add(hist);
}

/**
 * Perform event selection: This function overwrites the
 * event selection of the AliAnalysisTaskEmcal
 * @return Result of the event selection (true if event is selected)
 */
Bool_t AliAnalysisTaskChargedParticlesMCTriggerMimic::IsEventSelected(){
  AliDebugStream(1) << GetName() << "Using custom event selection method" << std::endl;
  if(!fTriggerPatchInfo){
    AliErrorStream() << GetName() << "Trigger patch container not found but required" << std::endl;
    return false;
  }
  if(!(fInputHandler->IsEventSelected() & AliVEvent::kINT7)) return false;

  // MC outlier cut
  if(fIsPythia){
    if(!CheckMCOutliers()) return false;
  }

  // Generall event quality cuts
  // The vertex cut also contains cuts on the number
  // of contributors and the position in z
  if(fAliAnalysisUtils){
    if(!fAliAnalysisUtils->IsVertexSelected2013pA(InputEvent())) return false;
    if(!fAliAnalysisUtils->IsPileUpEvent(InputEvent())) return false;
  }

  if(fPatchType != kUndef){
    if(!SelectEmcalTrigger(fTriggerPatchInfo)) return false;
  }
  return true;
}

/**
 * Run function
 * @return Always true
 */
Bool_t AliAnalysisTaskChargedParticlesMCTriggerMimic::Run(){
  Double_t weight = 1.;
  if(fIsPythia){
    weight = fWeightHandler ? fWeightHandler->GetEventWeight(fPythiaHeader) : 1.;
  }

  AliVParticle *truepart(nullptr);
  Bool_t isEMCAL(kFALSE);
  if(MCEvent()){
    for(int ipart = 0; ipart < fMCEvent->GetNumberOfTracks(); ipart++){
      truepart = MCEvent()->GetTrack(ipart);

      // Select only particles within ALICE acceptance
      if(!fEtaLabCut.IsInRange(truepart->Eta())) continue;
      if(TMath::Abs(truepart->Pt()) < 0.1) continue;
      if(!truepart->Charge()) continue;

      if(!IsPhysicalPrimary(truepart, fMCEvent)) continue;
      AliAODMCParticle *aodmc = static_cast<AliAODMCParticle *>(truepart);
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
      FillTrackHistos("True", weight, TMath::Abs(truepart->Pt()), truepart->Eta() * fEtaSign, etacent, truepart->Phi(), etacentcut, isEMCAL, pid);
    }
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
  for(int itrk = 0; itrk < fInputEvent->GetNumberOfTracks(); ++itrk){
    checktrack = dynamic_cast<AliVTrack *>(fInputEvent->GetTrack(itrk));
    if(!checktrack) continue;
    // Find associated particle
    assocMC = MCEvent()->GetTrack(TMath::Abs(checktrack->GetLabel()));
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

    FillTrackHistos("Accepted", weight,  ptparticle, checktrack->Eta() * fEtaSign, etacent, checktrack->Phi(), etacentcut, isEMCAL, assocpid);
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
void AliAnalysisTaskChargedParticlesMCTriggerMimic::FillTrackHistos(
    const char *eventclass,
    Double_t weight,
    Double_t pt,
    Double_t etalab,
    Double_t etacent,
    Double_t phi,
    Bool_t etacut,
    Bool_t inEmcal,
    const char *pid
    )
{
  fHistos->FillTH1(Form("hPtEtaAll%s", eventclass), TMath::Abs(pt), weight);
  fHistos->FillTH1(Form("hPtEtaAll%s%s", pid, eventclass), TMath::Abs(pt), weight);
  if(inEmcal){
    fHistos->FillTH1(Form("hPtEMCALEtaAll%s", eventclass), TMath::Abs(pt), weight);
    fHistos->FillTH1(Form("hPtEMCALEtaAll%s%s", pid, eventclass), TMath::Abs(pt), weight);
  }

  const std::array<int, 5> kPtMin = {1,2,5,10,20}; // for eta distributions
  for(const auto &ptmin : kPtMin){
    if(TMath::Abs(pt) > static_cast<double>(ptmin)){
      fHistos->FillTH1(Form("hPhiDistAllPt%d%s", ptmin, eventclass), phi, weight);
      fHistos->FillTH1(Form("hEtaLabDistAllPt%d%s", ptmin, eventclass), etalab, weight);
      fHistos->FillTH1(Form("hEtaCentDistAllPt%d%s", ptmin, eventclass), etacent, weight);
      if(inEmcal){
        fHistos->FillTH1(Form("hEtaLabDistAllEMCALPt%d%s", ptmin, eventclass), etalab, weight);
        fHistos->FillTH1(Form("hEtaCentDistAllEMCALPt%d%s", ptmin, eventclass), etacent, weight);
      }
    }
  }

  if(etacut){
    fHistos->FillTH1(Form("hPtEtaCent%s", eventclass), TMath::Abs(pt), weight);
    fHistos->FillTH1(Form("hPtEtaCent%s%s", pid, eventclass), TMath::Abs(pt), weight);
    if(inEmcal){
      fHistos->FillTH1(Form("hPtEMCALEtaCent%s", eventclass), TMath::Abs(pt), weight);
      fHistos->FillTH1(Form("hPtEMCALEtaCent%s%s", pid, eventclass), TMath::Abs(pt), weight);
    }
    for(const auto &ptmin : kPtMin){
      if(TMath::Abs(pt) > static_cast<double>(ptmin)){
        fHistos->FillTH1(Form("hEtaLabDistCutPt%d%s", ptmin, eventclass), etalab, weight);
        fHistos->FillTH1(Form("hEtaCentDistCutPt%d%s", ptmin, eventclass), etacent, weight);
        if(inEmcal){
          fHistos->FillTH1(Form("hEtaLabDistCutEMCALPt%d%s", ptmin, eventclass), etalab, weight);
          fHistos->FillTH1(Form("hEtaCentDistCutEMCALPt%d%s", ptmin, eventclass), etacent, weight);
        }
      }
    }
  }
}


/**
 * Check in a transparent way for ESDs and AODs whether the particle is physical primary or not
 * -# AOD: Information stored in the AliAODMCParticle
 * -# ESD: Information needs to be retrieved from the stack via the label of the MC particle
 * @param part The particle to check
 * @param mcevent The MC event containing the stack (ESD only)
 * @return True if particle is a physical primary particle, false otherwise
 */
Bool_t AliAnalysisTaskChargedParticlesMCTriggerMimic::IsPhysicalPrimary(const AliVParticle* const part, AliMCEvent* const mcevent) {
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
 * Set the track selection
 * @param cutname Name of the track cuts
 * @param isAOD check whether we run on ESDs or AODs
 */
void AliAnalysisTaskChargedParticlesMCTriggerMimic::InitializeTrackCuts(TString cutname, bool isAOD){
  SetTrackSelection(AliEmcalAnalysisFactory::TrackCutsFactory(cutname, isAOD));
}

/**
 * Select EMCAL-triggered event based on the presence of a trigger patch above
 * (offline) energy threshold. The threshold is settable, as well as at the patch
 * type.
 * @param[in] triggerpatches Array of trigger patches used for the trigger decision
 * @return True if the event is selected as triggered, false otherwise
 */
bool AliAnalysisTaskChargedParticlesMCTriggerMimic::SelectEmcalTrigger(const TClonesArray *triggerpatches){
  bool selected = false;
  AliEMCALTriggerPatchInfo *patch(nullptr);
  for(TIter patchiter = TIter(triggerpatches).Begin(); patchiter != TIter::End(); ++patchiter){
    patch = static_cast<AliEMCALTriggerPatchInfo *>(*patchiter);
    if(!patch->IsOfflineSimple()) continue;
    if(fPatchType == kEMCEJE){
      if(!patch->IsJetHighSimple()) continue;
    } else if(fPatchType == kEMCEGA){
      if(!patch->IsGammaHighSimple()) continue;
    }
    if(patch->GetPatchE() > fEnergyThreshold){
      // firing trigger patch found
      selected = true;
      break;
    }
  }
  return true;
}

/**
 * Create \f$ p_{t} \f$ binning
 */
AliAnalysisTaskChargedParticlesMCTriggerMimic::PtBinning::PtBinning() :
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
