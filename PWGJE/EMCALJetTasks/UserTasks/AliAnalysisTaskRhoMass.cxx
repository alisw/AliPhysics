// $Id$
//
// Calculation of rho mass from a collection of jets.
// If scale function is given the scaled rho will be exported
// with the name as "fOutRhoMassName".Apppend("_Scaled").
//
// Authors: M. Verweij

#include "AliAnalysisTaskRhoMass.h"

#include <TClonesArray.h>
#include <TMath.h>

#include "AliAnalysisManager.h"
#include "AliEmcalJet.h"
#include "AliLog.h"
#include "AliRhoParameter.h"

ClassImp(AliAnalysisTaskRhoMass)

//________________________________________________________________________
AliAnalysisTaskRhoMass::AliAnalysisTaskRhoMass() : 
  AliAnalysisTaskRhoMassBase("AliAnalysisTaskRhoMass"),
  fNExclLeadJets(0),
  fJetRhoMassType(kMd),
  fHistMdAreavsCent(0)
{
  // Constructor.
}

//________________________________________________________________________
AliAnalysisTaskRhoMass::AliAnalysisTaskRhoMass(const char *name, Bool_t histo) :
  AliAnalysisTaskRhoMassBase(name, histo),
  fNExclLeadJets(0),
  fJetRhoMassType(kMd),
  fHistMdAreavsCent(0)
{
  // Constructor.
}

//________________________________________________________________________
void AliAnalysisTaskRhoMass::UserCreateOutputObjects()
{
  // User create output objects, called at the beginning of the analysis.

  if (!fCreateHisto)
    return;
  
  AliAnalysisTaskRhoMassBase::UserCreateOutputObjects();

  fHistMdAreavsCent = new TH2F("fHistMdAreavsCent", "fHistMdAreavsCent", 101, -1,  100, fNbins, fMinBinPt, fMaxBinPt/2.);
  fHistMdAreavsCent->GetXaxis()->SetTitle("Centrality (%)");
  fHistMdAreavsCent->GetYaxis()->SetTitle("#rho_{m} (GeV/c * rad^{-1})");
  fOutput->Add(fHistMdAreavsCent);
}


//________________________________________________________________________
Bool_t AliAnalysisTaskRhoMass::Run() 
{
  // Run the analysis.

  fOutRhoMass->SetVal(0);
  if (fOutRhoMassScaled)
    fOutRhoMassScaled->SetVal(0);

  if (!fJets)
    return kFALSE;

  const Int_t Njets   = fJets->GetEntries();

  Int_t maxJetIds[]   = {-1, -1};
  Float_t maxJetPts[] = { 0,  0};

  if (fNExclLeadJets > 0) {
    for (Int_t ij = 0; ij < Njets; ++ij) {
      AliEmcalJet *jet = static_cast<AliEmcalJet*>(fJets->At(ij));
      if (!jet) {
	AliError(Form("%s: Could not receive jet %d", GetName(), ij));
	continue;
      } 

      if (!AcceptJet(jet))
        continue;

      if (jet->Pt() > maxJetPts[0]) {
	maxJetPts[1] = maxJetPts[0];
	maxJetIds[1] = maxJetIds[0];
	maxJetPts[0] = jet->Pt();
	maxJetIds[0] = ij;
      } else if (jet->Pt() > maxJetPts[1]) {
	maxJetPts[1] = jet->Pt();
	maxJetIds[1] = ij;
      }
    }
    if (fNExclLeadJets < 2) {
      maxJetIds[1] = -1;
      maxJetPts[1] = 0;
    }
  }

  static Double_t rhomvec[999];
  static Double_t Evec[999];
  static Double_t Mvec[999];
  Int_t NjetAcc = 0;

  // push all jets within selected acceptance into stack
  for (Int_t iJets = 0; iJets < Njets; ++iJets) {

    // exlcuding lead jets
    if (iJets == maxJetIds[0] || iJets == maxJetIds[1])
      continue;

    AliEmcalJet *jet = static_cast<AliEmcalJet*>(fJets->At(iJets));
    if (!jet) {
      AliError(Form("%s: Could not receive jet %d", GetName(), iJets));
      continue;
    } 

    if (!AcceptJet(jet))
      continue;

    // Double_t sumM = GetSumMConstituents(jet);
    // Double_t sumPt = GetSumPtConstituents(jet);
    if(jet->Area()>0.) {// && (jet->M()*jet->M() + jet->Pt()*jet->Pt())>0.) {
      //rhomvec[NjetAcc] = (TMath::Sqrt(sumM*sumM + sumPt*sumPt) - sumPt ) / jet->Area();
      // rhomvec[NjetAcc] = (TMath::Sqrt(jet->M()*jet->M() + jet->Pt()*jet->Pt()) - jet->Pt() ) / jet->Area();
      rhomvec[NjetAcc] = GetMd(jet) / jet->Area();
      fHistMdAreavsCent->Fill(fCent,rhomvec[NjetAcc]);
      Evec[NjetAcc] = jet->E();
      Mvec[NjetAcc] = jet->M();
      ++NjetAcc;
    }
  }

  if (NjetAcc > 0) {
    //find median value
    Double_t rhom = TMath::Median(NjetAcc, rhomvec);
    fOutRhoMass->SetVal(rhom);

    Int_t Ntracks = fTracks->GetEntries();
    Double_t meanM = TMath::Mean(NjetAcc, Mvec);
    Double_t meanE = TMath::Mean(NjetAcc, Evec);
    Double_t gamma = 0.;
    if(meanM>0.) gamma = meanE/meanM;
    fHistGammaVsNtrack->Fill(Ntracks,gamma);

    if (fOutRhoMassScaled) {
      Double_t rhomScaled = rhom * GetScaleFactor(fCent);
      fOutRhoMassScaled->SetVal(rhomScaled);
    }
  }

  return kTRUE;
} 

//________________________________________________________________________
Double_t AliAnalysisTaskRhoMass::GetSumMConstituents(AliEmcalJet *jet) {
  
  Double_t sum = 0.;
  
  AliVParticle *vp;
  for(Int_t icc=0; icc<jet->GetNumberOfTracks(); icc++) {
    vp = static_cast<AliVParticle*>(jet->TrackAt(icc, fTracks));
    if(!vp) continue;
    sum+=vp->M();
  }
  return sum;
}

//________________________________________________________________________
Double_t AliAnalysisTaskRhoMass::GetSumPtConstituents(AliEmcalJet *jet) {
  
  Double_t sum = 0.;
  
  AliVParticle *vp;
  for(Int_t icc=0; icc<jet->GetNumberOfTracks(); icc++) {
    vp = static_cast<AliVParticle*>(jet->TrackAt(icc, fTracks));
    if(!vp) continue;
    sum+=vp->Pt();
  }
  return sum;
}

//________________________________________________________________________
Double_t AliAnalysisTaskRhoMass::GetMd(AliEmcalJet *jet) {
  //get md as defined in http://arxiv.org/pdf/1211.2811.pdf
  Double_t sum = 0.;
  Double_t px = 0.;
  Double_t py = 0.;
  Double_t pz = 0.;
  Double_t E = 0.;

  if (fTracks) {
    AliVParticle *vp;
    for(Int_t icc=0; icc<jet->GetNumberOfTracks(); icc++) {
      vp = static_cast<AliVParticle*>(jet->TrackAt(icc, fTracks));
      if(!vp) continue;
      if(fJetRhoMassType==kMd) sum += TMath::Sqrt(vp->M()*vp->M() + vp->Pt()*vp->Pt()) - vp->Pt(); //sqrt(E^2-P^2+pt^2)=sqrt(E^2-pz^2)
      else if(fJetRhoMassType==kMdP) sum += TMath::Sqrt(vp->M()*vp->M() + vp->P()*vp->P()) - vp->P();
      else if(fJetRhoMassType==kMd4) {
	px+=vp->Px();
	py+=vp->Py();
	pz+=vp->Pz();
	E+=vp->E();
      }
    }
  }

  if (fCaloClusters) {
    AliVCluster *vp;
    for(Int_t icc=0; icc<jet->GetNumberOfClusters(); icc++) {
      vp = static_cast<AliVCluster*>(jet->ClusterAt(icc, fCaloClusters));
      if(!vp) continue;
      TLorentzVector nPart;
      vp->GetMomentum(nPart, fVertex);
      
      if(fJetRhoMassType==kMd) sum += TMath::Sqrt(nPart.M()*nPart.M() + nPart.Pt()*nPart.Pt()) - nPart.Pt();
      else if(fJetRhoMassType==kMdP) sum += TMath::Sqrt(nPart.M()*nPart.M() + nPart.P()*nPart.P()) - nPart.P();
      else if(fJetRhoMassType==kMd4) {
	px+=nPart.Px();
	py+=nPart.Py();
	pz+=nPart.Pz();
	E+=nPart.E();
      }
    }
  }

  if(fJetRhoMassType==kMd4) {
    Double_t pt = TMath::Sqrt(px*px + py*py);
    Double_t m2 = E*E - pt*pt - pz*pz;
    sum = TMath::Sqrt(m2 + pt*pt) - pt;
  }
  return sum;
}
