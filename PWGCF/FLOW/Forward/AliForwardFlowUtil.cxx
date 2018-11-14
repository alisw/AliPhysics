#include "TString.h"
#include "TMath.h"
#include "AliForwardFlowUtil.h"
#include "TFile.h"

#include <iostream>
#include <TROOT.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <TList.h>
#include <THn.h>

#include "AliLog.h"
#include "AliForwardFlowRun2Task.h"
#include "AliForwardQCumulantRun2.h"
#include "AliForwardGenericFramework.h"

#include "AliAODForwardMult.h"
#include "AliAODCentralMult.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"

#include "AliForwardFlowUtil.h"

#include "AliVVZERO.h"
#include "AliAODVertex.h"
#include "AliCentrality.h"

#include "AliESDEvent.h"
#include "AliVTrack.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliAODTracklets.h"

#include "AliAnalysisFilter.h"
#include "AliMultSelection.h"
#include "AliMultiplicity.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"

#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
//________________________________________________________________________
AliForwardFlowUtil::AliForwardFlowUtil():
fevent(),
fAODevent(),
fMCevent()
{
}

Bool_t AliForwardFlowUtil::ExtraEventCutFMD(TH2D& forwarddNdedp, double cent, Bool_t mc){
  Bool_t useEvent = true;
  Int_t nBadBins = 0;
  Int_t phibins = forwarddNdedp.GetNbinsY();
  Double_t totalFMDpar = 0;

  for (Int_t etaBin = 1; etaBin <= forwarddNdedp.GetNbinsX(); etaBin++) {

    Double_t eta = forwarddNdedp.GetXaxis()->GetBinCenter(etaBin);
    Double_t runAvg = 0;
    Double_t avgSqr = 0;
    Double_t max = 0;
    Int_t nInAvg = 0;

    for (Int_t phiBin = 0; phiBin <= phibins; phiBin++) {
      // if (!mc){
      //   if ( fabs(eta) > 1.7) {
      //     if (phiBin == 0 && forwarddNdedp.GetBinContent(etaBin, 0) == 0) break;
      //   }
      // }
      Double_t weight = forwarddNdedp.GetBinContent(etaBin, phiBin);
      if (!weight){
        weight = 0;
      }
      totalFMDpar += weight;

      // We calculate the average Nch per. bin
      avgSqr += weight*weight;
      runAvg += weight;
      nInAvg++;
      if (weight == 0) continue;
      if (weight > max) {
        max = weight;
      }
    } // End of phi loop

    // Outlier cut calculations
    double fSigmaCut = 4.0;
    if (nInAvg > 0) {
      runAvg /= nInAvg;
      avgSqr /= nInAvg;
      Double_t stdev = (nInAvg > 1 ? TMath::Sqrt(nInAvg/(nInAvg-1))*TMath::Sqrt(avgSqr - runAvg*runAvg) : 0);
      Double_t nSigma = (stdev == 0 ? 0 : (max-runAvg)/stdev);
      if (fSigmaCut > 0. && nSigma >= fSigmaCut && cent < 60) nBadBins++;
      else nBadBins = 0;
      // We still finish the loop, for fOutliers to make sense,
      // but we do no keep the event for analysis
      if (nBadBins > 3) useEvent = false;
     //if (nBadBins > 3) std::cout << "NUMBER OF BAD BINS > 3" << std::endl;
    }
  } // End of eta bin
  if (totalFMDpar < 10) useEvent = false;

  return useEvent;
}


Double_t AliForwardFlowUtil::GetZ(){
  return fevent->GetPrimaryVertex()->GetZ();
}


Double_t AliForwardFlowUtil::GetCentrality(TString centrality_estimator){
  // Get MultSelection
  AliMultSelection *MultSelection;

  //if (fSettings.esd) MultSelection = (AliMultSelection*)dynamic_cast<AliMCEvent*>(InputEvent())->FindListObject("MultSelection");
  //else
  MultSelection  = (AliMultSelection*)fevent->FindListObject("MultSelection");

  return MultSelection->GetMultiplicityPercentile(centrality_estimator);
}


void AliForwardFlowUtil::FillFromTrackrefs(TH2D*& cen, TH2D*& fwd) const
{

  Int_t nTracks = fMCevent->Stack()->GetNtrack();

  for (Int_t iTr = 0; iTr < nTracks; iTr++) {
      AliMCParticle* p = static_cast< AliMCParticle* >(fMCevent->GetTrack(iTr));
    if (p->Charge() == 0) continue;

    for (Int_t iTrRef = 0; iTrRef < p->GetNumberOfTrackReferences(); iTrRef++) {
      AliTrackReference* ref = p->GetTrackReference(iTrRef);
      // Check hit on FMD
      if (!ref) continue;
      if (AliTrackReference::kTPC != ref->DetectorId()){
        Double_t x      = ref->X() - fevent->GetPrimaryVertex()->GetX();
        Double_t y      = ref->Y() - fevent->GetPrimaryVertex()->GetY();
        Double_t z      = ref->Z() - fevent->GetPrimaryVertex()->GetZ();
        Double_t rr     = TMath::Sqrt(x * x + y * y);
        Double_t thetaR = TMath::ATan2(rr, z);
        Double_t phiR   = TMath::ATan2(y,x);

        if (phiR < 0) {
          phiR += 2*TMath::Pi();
        }
        if (thetaR < 0) {
          thetaR += 2*TMath::Pi();
        }
        cen->Fill(-TMath::Log(TMath::Tan(thetaR / 2)),phiR);
      }
      else if (AliTrackReference::kFMD != ref->DetectorId()) {
        Double_t x      = ref->X() - fevent->GetPrimaryVertex()->GetX();
        Double_t y      = ref->Y() - fevent->GetPrimaryVertex()->GetY();
        Double_t z      = ref->Z() - fevent->GetPrimaryVertex()->GetZ();
        Double_t rr     = TMath::Sqrt(x * x + y * y);
        Double_t thetaR = TMath::ATan2(rr, z);
        Double_t phiR   = TMath::ATan2(y,x);

        if (phiR < 0) {
          phiR += 2*TMath::Pi();
        }
        if (thetaR < 0) {
          thetaR += 2*TMath::Pi();
        }
        fwd->Fill(-TMath::Log(TMath::Tan(thetaR / 2)),phiR);
      }
    }
  }
}


void AliForwardFlowUtil::FillFromPrimaries(TH2D*& cen, TH2D*& fwd) const
{
  Int_t nTracksMC = fMCevent->GetNumberOfTracks();

  for (Int_t iTr = 0; iTr < nTracksMC; iTr++) {
    AliMCParticle* p = static_cast< AliMCParticle* >(fMCevent->GetTrack(iTr));
    if (!p->IsPhysicalPrimary()) continue;
    if (p->Charge() == 0) continue;

    Double_t eta = p->Eta();
    if (TMath::Abs(eta) < 1.1) {
      if (p->Pt()>=0.2 && p->Pt()<=5)
        cen->Fill(eta,p->Phi(),1);
    }
    if (eta < 5 /*fwd->GetXaxis()-GetXmax()*/ && eta > -3.5 /*fwd->GetXaxis()-GetXmin()*/) {
      if (TMath::Abs(eta) >= 1.7)
        fwd->Fill(eta,p->Phi(),1);
    }
  }
}


void AliForwardFlowUtil::FillFromTracklets(TH2D*& cen) const {
  AliAODTracklets* aodTracklets = fAODevent->GetTracklets();

  for (Int_t i = 0; i < aodTracklets->GetNumberOfTracklets(); i++) {
    cen->Fill(aodTracklets->GetEta(i),aodTracklets->GetPhi(i), 1);
  }
}


void AliForwardFlowUtil::FillFromTracks(TH2D*& cen, Int_t tracktype) const {
  Int_t  iTracks(fevent->GetNumberOfTracks());
  for(Int_t i(0); i < iTracks; i++) {

  // loop  over  all  the  tracks
    AliAODTrack* track = static_cast<AliAODTrack *>(fAODevent->GetTrack(i));
    if (track->TestFilterBit(tracktype)){
      if (track->Pt() >= 0.2 && track->Pt() <= 5){
        cen->Fill(track->Eta(),track->Phi(), 1);
      }
    }
  }
}


AliTrackReference* AliForwardFlowUtil::IsHitFMD(AliMCParticle* p) {
  //std::cout << "p->GetNumberOfTrackReferences() = " << p->GetNumberOfTrackReferences() << std::endl;
  for (Int_t iTrRef = 0; iTrRef < p->GetNumberOfTrackReferences(); iTrRef++) {
    AliTrackReference* ref = p->GetTrackReference(iTrRef);
    // Check hit on FMD
    if (!ref || AliTrackReference::kFMD != ref->DetectorId()) {
      continue;
    }
    else {
      return ref;
    }
  }
  return 0x0;
}

AliTrackReference* AliForwardFlowUtil::IsHitTPC(AliMCParticle* p) {
  for (Int_t iTrRef = 0; iTrRef < p->GetNumberOfTrackReferences(); iTrRef++) {
    AliTrackReference* ref = p->GetTrackReference(iTrRef);
    // Check hit on FMD
    if (!ref || AliTrackReference::kTPC != ref->DetectorId()) {
      continue;
    }
    else {
      return ref;
    }
  }
  return 0x0;
}
