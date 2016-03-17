//
// Calculation of rho for flow bias studies
//
// Author: S.Aiola

#include "AliAnalysisTaskRhoFlow.h"

#include <TClonesArray.h>
#include <TMath.h>
#include <TH2F.h>

#include "AliVTrack.h"
#include "AliLog.h"
#include "AliRhoParameter.h"
#include "AliParticleContainer.h"

ClassImp(AliAnalysisTaskRhoFlow)

//________________________________________________________________________
AliAnalysisTaskRhoFlow::AliAnalysisTaskRhoFlow() : 
  AliAnalysisTaskRho("AliAnalysisTaskRhoFlow",kTRUE),
  fRhoNearSide(0),
  fRhoAwaySide(0),
  fRhoPerpSide1(0),
  fRhoPerpSide2(0),
  fHistRhoNearVsCent(0),
  fHistDeltaRhoNearVsCent(0),
  fHistRhoAwayVsCent(0),
  fHistDeltaRhoAwayVsCent(0),
  fHistRhoPerp1VsCent(0),
  fHistDeltaRhoPerp1VsCent(0),
  fHistRhoPerp2VsCent(0),
  fHistDeltaRhoPerp2VsCent(0)
{
  // Constructor.
  SetAttachToEvent(kFALSE);
}

//________________________________________________________________________
AliAnalysisTaskRhoFlow::AliAnalysisTaskRhoFlow(const char *name) :
  AliAnalysisTaskRho(name, kTRUE),
  fRhoNearSide(0),
  fRhoAwaySide(0),
  fRhoPerpSide1(0),
  fRhoPerpSide2(0),
  fHistRhoNearVsCent(0),
  fHistDeltaRhoNearVsCent(0),
  fHistRhoAwayVsCent(0),
  fHistDeltaRhoAwayVsCent(0),
  fHistRhoPerp1VsCent(0),
  fHistDeltaRhoPerp1VsCent(0),
  fHistRhoPerp2VsCent(0),
  fHistDeltaRhoPerp2VsCent(0)
{
  // Constructor.
  SetAttachToEvent(kFALSE);
}

//________________________________________________________________________
void AliAnalysisTaskRhoFlow::UserCreateOutputObjects()
{
  // User create output objects, called at the beginning of the analysis.

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  fHistRhoNearVsCent = new TH2F("RhoNearVsCent", "RhoNearVsCent", 101, -1,  100, fNbins, fMinBinPt, fMaxBinPt*2);
  fOutput->Add(fHistRhoNearVsCent);

  fHistRhoAwayVsCent = new TH2F("RhoAwayVsCent", "RhoAwayVsCent", 101, -1,  100, fNbins, fMinBinPt, fMaxBinPt*2);
  fOutput->Add(fHistRhoAwayVsCent);

  fHistRhoPerp1VsCent = new TH2F("RhoPerp1VsCent", "RhoPerp1VsCent", 101, -1,  100, fNbins, fMinBinPt, fMaxBinPt*2);
  fOutput->Add(fHistRhoPerp1VsCent);

  fHistRhoPerp2VsCent = new TH2F("RhoPerp2VsCent", "RhoPerp2VsCent", 101, -1,  100, fNbins, fMinBinPt, fMaxBinPt*2);
  fOutput->Add(fHistRhoPerp2VsCent);

  if (!fCompareRhoName.IsNull()) {
    fHistDeltaRhoNearVsCent = new TH2F("DeltaRhoNearVsCent", "DeltaRhoNearVsCent", 101, -1, 100, fNbins, -fMaxBinPt, fMaxBinPt);
    fOutput->Add(fHistDeltaRhoNearVsCent);

    fHistDeltaRhoAwayVsCent = new TH2F("DeltaRhoAwayVsCent", "DeltaRhoAwayVsCent", 101, -1, 100, fNbins, -fMaxBinPt, fMaxBinPt);
    fOutput->Add(fHistDeltaRhoAwayVsCent);
    
    fHistDeltaRhoPerp1VsCent = new TH2F("DeltaRhoPerp1VsCent", "DeltaRhoPerp1VsCent", 101, -1, 100, fNbins, -fMaxBinPt, fMaxBinPt);
    fOutput->Add(fHistDeltaRhoPerp1VsCent);

    fHistDeltaRhoPerp2VsCent = new TH2F("DeltaRhoPerp2VsCent", "DeltaRhoPerp2VsCent", 101, -1, 100, fNbins, -fMaxBinPt, fMaxBinPt);
    fOutput->Add(fHistDeltaRhoPerp2VsCent);
  }
}

//________________________________________________________________________
Bool_t AliAnalysisTaskRhoFlow::Run() 
{
  AliParticleContainer* tracks = GetParticleContainer(0);
  if (!tracks) return kFALSE;
  
  Double_t jetRadius = GetJetRadius();
  Double_t maxTrackPhi = -1;
  Double_t maxTrackPt  = 0;

  AliVParticle *track = 0;
  tracks->ResetCurrentID();
  while ((track = tracks->GetNextAcceptParticle())) {

    if (track->Pt() > maxTrackPt) {
      maxTrackPt  = track->Pt();
      maxTrackPhi = track->Phi();
    }
  }

  Double_t minPhi = -1;
  Double_t maxPhi = -1;

  // away side
  minPhi = maxTrackPhi + TMath::Pi() - TMath::Pi()/4 + jetRadius;
  maxPhi = maxTrackPhi + TMath::Pi() + TMath::Pi()/4 + jetRadius;
  if (maxPhi > TMath::Pi() * 2) {
    minPhi -= TMath::Pi() * 2;
    maxPhi -= TMath::Pi() * 2;
  }
  SetJetPhiLimits(minPhi, maxPhi);  
  fNExclLeadJets = 1;
  AliAnalysisTaskRho::Run();
  fRhoAwaySide = fOutRho->GetVal();

  // perp 1 side
  minPhi = maxTrackPhi + TMath::Pi()/2 - TMath::Pi()/4 + jetRadius;
  maxPhi = maxTrackPhi + TMath::Pi()/2 + TMath::Pi()/4 + jetRadius;
  if (maxPhi > TMath::Pi() * 2) {
    minPhi -= TMath::Pi() * 2;
    maxPhi -= TMath::Pi() * 2;
  }
  SetJetPhiLimits(minPhi, maxPhi);  
  fNExclLeadJets = 0;
  AliAnalysisTaskRho::Run();
  fRhoPerpSide1 = fOutRho->GetVal();

  // perp 2 side
  minPhi = maxTrackPhi - TMath::Pi()/2 - TMath::Pi()/4 + jetRadius;
  maxPhi = maxTrackPhi - TMath::Pi()/2 + TMath::Pi()/4 + jetRadius;
  if (maxPhi > TMath::Pi() * 2) {
    minPhi -= TMath::Pi() * 2;
    maxPhi -= TMath::Pi() * 2;
  }
  SetJetPhiLimits(minPhi, maxPhi);
  fNExclLeadJets = 0;
  AliAnalysisTaskRho::Run();
  fRhoPerpSide2 = fOutRho->GetVal();

  // near side
  minPhi = maxTrackPhi - TMath::Pi()/4 + jetRadius;
  maxPhi = maxTrackPhi + TMath::Pi()/4 + jetRadius;
  if (maxPhi > TMath::Pi() * 2) {
    minPhi -= TMath::Pi() * 2;
    maxPhi -= TMath::Pi() * 2;
  }
  SetJetPhiLimits(minPhi, maxPhi);
  fNExclLeadJets = 1;
  AliAnalysisTaskRho::Run();
  fRhoNearSide = fOutRho->GetVal();
  
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskRhoFlow::FillHistograms() 
{
  // Fill histograms.

  fHistRhoNearVsCent->Fill(fCent, fRhoNearSide);  
  fHistRhoAwayVsCent->Fill(fCent, fRhoAwaySide);
  fHistRhoPerp1VsCent->Fill(fCent, fRhoPerpSide1);
  fHistRhoPerp2VsCent->Fill(fCent, fRhoPerpSide2);

  if (fCompareRho) {
    fHistDeltaRhoNearVsCent->Fill(fCent, fRhoNearSide - fCompareRho->GetVal());
    fHistDeltaRhoAwayVsCent->Fill(fCent, fRhoAwaySide - fCompareRho->GetVal());
    fHistDeltaRhoPerp1VsCent->Fill(fCent, fRhoPerpSide1 - fCompareRho->GetVal());
    fHistDeltaRhoPerp2VsCent->Fill(fCent, fRhoPerpSide2 - fCompareRho->GetVal());
  }

  return kTRUE;
}      
