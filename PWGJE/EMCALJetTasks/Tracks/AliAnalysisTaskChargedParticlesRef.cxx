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
#include <algorithm>
#include <array>
#include <functional>
#include <iostream>
#include <map>

#include <TMath.h>
#include <THistManager.h>
#include <TLinearBinning.h>

#include "AliAnalysisUtils.h"
#include "AliAODInputHandler.h"
#include "AliAODTrack.h"
#include "AliEmcalAnalysisFactory.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALRecoUtils.h"
#include "AliEMCALTriggerPatchInfo.h"
#include "AliEmcalTrackSelection.h"
#include "AliESDtrackCuts.h"
#include "AliESDEvent.h"
#include "AliInputEventHandler.h"
#include "AliPIDResponse.h"
#include "AliTOFPIDResponse.h"
#include "AliVVertex.h"

#include "AliEMCalTriggerExtraCuts.h"
#include "AliAnalysisTaskChargedParticlesRef.h"

/// \cond CLASSIMP
ClassImp(EMCalTriggerPtAnalysis::AliAnalysisTaskChargedParticlesRef)
/// \endcond

namespace EMCalTriggerPtAnalysis {

AliAnalysisTaskChargedParticlesRef::AliAnalysisTaskChargedParticlesRef() :
    AliAnalysisTaskEmcalTriggerBase(),
    fTrackCuts(nullptr),
    fYshift(0.465),
    fEtaSign(1),
    fEtaLabCut(-0.5, 0.5),
    fEtaCmsCut(-2., 2.),
    fPhiCut(0., TMath::TwoPi()),
    fStudyPID(false)
{
}

AliAnalysisTaskChargedParticlesRef::AliAnalysisTaskChargedParticlesRef(const char *name) :
    AliAnalysisTaskEmcalTriggerBase(name),
    fTrackCuts(nullptr),
    fYshift(0.465),
    fEtaSign(1),
    fEtaLabCut(-0.5, 0.5),
    fEtaCmsCut(-2., 2.),
    fPhiCut(0., TMath::TwoPi()),
    fStudyPID(false)
{
  SetNeedEmcalGeom(true);
  SetCaloTriggerPatchInfoName("EmcalTriggers");
}

AliAnalysisTaskChargedParticlesRef::~AliAnalysisTaskChargedParticlesRef() {
  //if(fTrackCuts) delete fTrackCuts;
}

void AliAnalysisTaskChargedParticlesRef::CreateUserObjects(){
  if(!fTrackCuts) InitializeTrackCuts("standard", fInputHandler->IsA() == AliAODInputHandler::Class());
}

void AliAnalysisTaskChargedParticlesRef::CreateUserHistos() {

  PtBinning newbinning;

  // Binning for the PID histogram
  const int kdimPID = 3;
  const int knbinsPID[kdimPID] = {1000, 200, 300};
  const double kminPID[kdimPID] = {-100., 0.,  0.}, kmaxPID[kdimPID] = {100., 200., 1.5};
  for(auto trg : GetSupportedTriggers()){
    fHistos->CreateTH1(Form("hEventCount%s", trg.Data()), Form("Event Counter for trigger class %s", trg.Data()), 1, 0.5, 1.5);
    fHistos->CreateTH1(Form("hVertexBefore%s", trg.Data()), Form("Vertex distribution before z-cut for trigger class %s", trg.Data()), 500, -50, 50);
    fHistos->CreateTH1(Form("hVertexAfter%s", trg.Data()), Form("Vertex distribution after z-cut for trigger class %s", trg.Data()), 100, -10, 10);

    fHistos->CreateTH3(Form("hPtEtaPhiAll%s", trg.Data()), Form("p_{t}-#eta-#phi distribution of all accepted tracks for trigger %s; p_{t} (GeV/c); #eta; #phi", trg.Data()), newbinning, TLinearBinning(64, -0.8, 0.8), TLinearBinning(100, 0., 2*TMath::Pi()));
    fHistos->CreateTH3(Form("hPtEtaPhiEMCALAll%s", trg.Data()), Form("p_{t}-#eta-#phi distribution of all accepted tracks pointing to the EMCAL for trigger %s; p_{t} (GeV/c); #eta; #phi", trg.Data()), newbinning, TLinearBinning(64, -0.8, 0.8), TLinearBinning(100, 0., 2*TMath::Pi()));
    fHistos->CreateTH3(Form("hPtEtaPhiCent%s", trg.Data()), Form("p_{t}-#eta-#phi distribution of all accepted tracks for trigger %s; p_{t} (GeV/c); #eta; #phi", trg.Data()), newbinning, TLinearBinning(64, -0.8, 0.8), TLinearBinning(100, 0., 2*TMath::Pi()));
    fHistos->CreateTH3(Form("hPtEtaPhiEMCALCent%s", trg.Data()), Form("p_{t}-#eta-#phi distribution of all accepted tracks pointing to the EMCAL for trigger %s; p_{t} (GeV/c); #eta; #phi", trg.Data()), newbinning, TLinearBinning(64, -0.8, 0.8), TLinearBinning(100, 0., 2*TMath::Pi()));

    if(fStudyPID){
      fHistos->CreateTH2(Form("hTPCdEdxEMCAL%s", trg.Data()), Form("TPC dE/dx of charged particles in the EMCAL region for trigger %s", trg.Data()), 400, -20., 20., 200, 0., 200.);
      fHistos->CreateTH2(Form("hTOFBetaEMCAL%s", trg.Data()), Form("TOF beta  of charged particles in the EMCAL region for trigger %s", trg.Data()), 400, -20., 20., 150, 0., 1.5);
      fHistos->CreateTHnSparse(Form("hPIDcorrEMCAL%s", trg.Data()), Form("Correlation of PID observables for Trigger %s", trg.Data()), kdimPID, knbinsPID, kminPID, kmaxPID);
    }
  }
}

Bool_t AliAnalysisTaskChargedParticlesRef::Run() {
  Bool_t hasPIDresponse = fInputHandler->GetPIDResponse() != nullptr;
  if(fStudyPID && !hasPIDresponse) AliErrorStream() << "PID requested but PID response not available" << std::endl;


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
  Bool_t isEMCAL(kFALSE);
  Double_t etaEMCAL(0.), phiEMCAL(0.);
  for(int itrk = 0; itrk < fInputEvent->GetNumberOfTracks(); ++itrk){
    checktrack = dynamic_cast<AliVTrack *>(fInputEvent->GetTrack(itrk));
    if(!checktrack) continue;
    if(!fEtaLabCut.IsInRange(checktrack->Eta())) continue;
    if(!fPhiCut.IsInRange(checktrack->Phi())) continue;
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

    // Calculate eta in cms frame according
    // EPJC74 (2014) 3054:
    // eta_cms = - eta_lab - |yshift|
    Double_t etacent = -1. * checktrack->Eta() - TMath::Abs(fYshift);
    etacent *= fEtaSign;

    if(!fEtaCmsCut.IsInRange(etacent)) continue;    // Apply eta-cent cut
    if(!fTrackCuts->IsTrackAccepted(checktrack)) continue;

    for(const auto &t : fSelectedTriggers){
      FillTrackHistos(t, checktrack->Pt(), checktrack->Eta() * fEtaSign, etacent, checktrack->Phi(), isEMCAL);
      if(fStudyPID && hasPIDresponse)
        if(isEMCAL) FillPIDHistos(t, *checktrack);
    }
  }
  return true;
}

void AliAnalysisTaskChargedParticlesRef::UserFillHistosBeforeEventSelection(){
  // Apply vertex z cut
  for(const auto &t : fSelectedTriggers){
    Double_t weight = GetTriggerWeight(t);
    fHistos->FillTH1(Form("hVertexBefore%s", t.Data()), fVertex[2], weight);
  }
}

void AliAnalysisTaskChargedParticlesRef::UserFillHistosAfterEventSelection(){
  for(const auto &t : fSelectedTriggers) {
    Double_t weight = GetTriggerWeight(t);
    // Fill Event counter and reference vertex distributions after event selection
    fHistos->FillTH1(Form("hEventCount%s", t.Data()), 1, weight);
    fHistos->FillTH1(Form("hVertexAfter%s", t.Data()), fVertex[2], weight);
  }
}

void AliAnalysisTaskChargedParticlesRef::FillTrackHistos(
    const TString &eventclass,
    Double_t pt,
    Double_t etalab,
    Double_t etacent,
    Double_t phi,
    Bool_t inEmcal
    )
{
  Double_t weight = GetTriggerWeight(eventclass);
  AliDebugStream(1) << GetName() << ": Using weight " << weight << " for trigger " << eventclass << " in particle histograms." << std::endl;
  double kinepointall[3] = {TMath::Abs(pt), etalab, phi}, kinepointcent[3] = {TMath::Abs(pt), etacent, phi};
  fHistos->FillTH3(Form("hPtEtaPhiAll%s", eventclass.Data()), kinepointall, weight);
  fHistos->FillTH3(Form("hPtEtaPhiCent%s", eventclass.Data()), kinepointcent, weight);
  if(inEmcal){
    fHistos->FillTH3(Form("hPtEtaPhiEMCALAll%s", eventclass.Data()), kinepointall, weight);
    fHistos->FillTH3(Form("hPtEtaPhiEMCALCent%s", eventclass.Data()), kinepointall, weight);
  }
}

void AliAnalysisTaskChargedParticlesRef::FillPIDHistos(
    const TString &eventclass,
    const AliVTrack &trk
) {
  Double_t weight = GetTriggerWeight(eventclass);
  AliDebugStream(1) << GetName() << ": Using weight " << weight << " for trigger " << eventclass << " in PID histograms." << std::endl;
  AliPIDResponse *pid = fInputHandler->GetPIDResponse();
  if(TMath::Abs(trk.Eta()) > 0.5) return;
  if(!((trk.GetStatus() & AliVTrack::kTOFout) && (trk.GetStatus() & AliVTrack::kTIME))) return;

  double poverz = TMath::Abs(trk.P())/static_cast<double>(trk.Charge());
  fHistos->FillTH2(Form("hTPCdEdxEMCAL%s", eventclass.Data()), poverz, trk.GetTPCsignal(), weight);
  // correct for units - TOF in ps, track length in cm
  Double_t trtime = (trk.GetTOFsignal() - pid->GetTOFResponse().GetTimeZero()) * 1e-12;
  Double_t v = trk.GetIntegratedLength()/(100. * trtime);
  Double_t beta =  v / TMath::C();
  fHistos->FillTH2(Form("hTOFBetaEMCAL%s", eventclass.Data()), poverz, beta, weight);
  double datapoint[3] = {poverz, trk.GetTPCsignal(), beta};
  fHistos->FillTHnSparse(Form("hPIDcorrEMCAL%s", eventclass.Data()), datapoint, weight);
}

void AliAnalysisTaskChargedParticlesRef::InitializeTrackCuts(TString cutname, bool isAOD){
  SetEMCALTrackSelection(AliEmcalAnalysisFactory::TrackCutsFactory(cutname, isAOD));
}

AliAnalysisTaskChargedParticlesRef::PtBinning::PtBinning():
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
