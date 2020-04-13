/************************************************************************************
 * Copyright (C) 2014, Copyright Holders of the ALICE Collaboration                 *
 * All rights reserved.                                                             *
 *                                                                                  *
 * Redistribution and use in source and binary forms, with or without               *
 * modification, are permitted provided that the following conditions are met:      *
 *     * Redistributions of source code must retain the above copyright             *
 *       notice, this list of conditions and the following disclaimer.              *
 *     * Redistributions in binary form must reproduce the above copyright          *
 *       notice, this list of conditions and the following disclaimer in the        *
 *       documentation and/or other materials provided with the distribution.       *
 *     * Neither the name of the <organization> nor the                             *
 *       names of its contributors may be used to endorse or promote products       *
 *       derived from this software without specific prior written permission.      *
 *                                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND  *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED    *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE           *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY              *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES       *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND      *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT       *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                     *
 ************************************************************************************/
#include "AliAnalysisTaskRhoMassSparse.h"

#include <TClonesArray.h>
#include <TMath.h>

#include "AliAnalysisManager.h"
#include "AliEmcalJet.h"
#include "AliLog.h"
#include "AliRhoParameter.h"
#include "AliJetContainer.h"

ClassImp(AliAnalysisTaskRhoMassSparse)

AliAnalysisTaskRhoMassSparse::AliAnalysisTaskRhoMassSparse() : 
  AliAnalysisTaskRhoMassBase("AliAnalysisTaskRhoMassSparse"),
  fNExclLeadJets(0),
  fJetRhoMassType(kMd),
  fPionMassClusters(kFALSE),
  fHistMdAreavsCent(0)
{
}

AliAnalysisTaskRhoMassSparse::AliAnalysisTaskRhoMassSparse(const char *name, Bool_t histo) :
  AliAnalysisTaskRhoMassBase(name, histo),
  fNExclLeadJets(0),
  fJetRhoMassType(kMd),
  fPionMassClusters(kFALSE),
  fHistMdAreavsCent(0)
{
}

void AliAnalysisTaskRhoMassSparse::UserCreateOutputObjects()
{
  if (!fCreateHisto)
    return;
  
  AliAnalysisTaskRhoMassBase::UserCreateOutputObjects();

  fHistMdAreavsCent = new TH2F("fHistMdAreavsCent", "fHistMdAreavsCent", 101, -1,  100, fNbins, fMinBinPt, fMaxBinPt/2.);
  fHistMdAreavsCent->GetXaxis()->SetTitle("Centrality (%)");
  fHistMdAreavsCent->GetYaxis()->SetTitle("#rho_{m} (GeV/c * rad^{-1})");
  fOutput->Add(fHistMdAreavsCent);
  
  fHistOccCorrvsCent = new TH2F("OccCorrvsCent", "OccCorrvsCent", 101, -1, 100, 2000, 0 , 2);
  fOutput->Add(fHistOccCorrvsCent);

}

Bool_t AliAnalysisTaskRhoMassSparse::IsJetOverlapping(AliEmcalJet* jet1, AliEmcalJet* jet2)
{
  for (Int_t i = 0; i < jet1->GetNumberOfTracks(); ++i)
  {
    Int_t jet1Track = jet1->TrackAt(i);
    for (Int_t j = 0; j < jet2->GetNumberOfTracks(); ++j)
    {
      Int_t jet2Track = jet2->TrackAt(j);
      if (jet1Track == jet2Track)
        return kTRUE;
    }
  }
  return kFALSE;
}

Bool_t AliAnalysisTaskRhoMassSparse::IsJetSignal(AliEmcalJet* jet)
{
  if(jet->Pt()>5){
      return kTRUE;
  }else{
    return kFALSE;
  }
}


Bool_t AliAnalysisTaskRhoMassSparse::Run() 
{
  fOutRhoMass->SetVal(0);
  if (fOutRhoMassScaled)
    fOutRhoMassScaled->SetVal(0);

  if (!fJets)
    return kFALSE;

  const Int_t Njets   = fJets->GetEntries();
  
  AliJetContainer *sigjets = static_cast<AliJetContainer*>(fJetCollArray.At(1));
  
  Int_t NjetsSig = 0;
  if (sigjets) NjetsSig = sigjets->GetNJets();

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
  Double_t TotaljetArea=0;
  Double_t TotaljetAreaPhys=0;

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

    TotaljetArea+=jet->Area();

    if(jet->Pt()>0.1){
      TotaljetAreaPhys+=jet->Area();
    }

    if (!AcceptJet(jet))
      continue;

      // Search for overlap with signal jets
    Bool_t isOverlapping = kFALSE;
    if (sigjets) {
      for(Int_t j=0;j<NjetsSig;j++)
	{
	  AliEmcalJet* signalJet = sigjets->GetAcceptJet(j);
	  if(!signalJet)
	    continue;
	  if(!IsJetSignal(signalJet))     
	    continue;
	  
	  if(IsJetOverlapping(signalJet, jet))
	    {
	      isOverlapping = kTRUE;
	      break;
	    }
	}
    }

    if(isOverlapping) 
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

  Double_t OccCorr=0.0;
  if(TotaljetArea>0) OccCorr=TotaljetAreaPhys/TotaljetArea;
 
  if (fCreateHisto)
    fHistOccCorrvsCent->Fill(fCent, OccCorr);


  if (NjetAcc > 0) {
    //find median value
    Double_t rhom = TMath::Median(NjetAcc, rhomvec);
    if(fRhoCMS){
      rhom = rhom * OccCorr;
    }

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

Double_t AliAnalysisTaskRhoMassSparse::GetSumMConstituents(AliEmcalJet *jet) {
  
  Double_t sum = 0.;
  
  AliVParticle *vp;
  for(Int_t icc=0; icc<jet->GetNumberOfTracks(); icc++) {
    vp = static_cast<AliVParticle*>(jet->TrackAt(icc, fTracks));
    if(!vp) continue;
    sum+=vp->M();
  }
  return sum;
}

Double_t AliAnalysisTaskRhoMassSparse::GetSumPtConstituents(AliEmcalJet *jet) {
  
  Double_t sum = 0.;
  
  AliVParticle *vp;
  for(Int_t icc=0; icc<jet->GetNumberOfTracks(); icc++) {
    vp = static_cast<AliVParticle*>(jet->TrackAt(icc, fTracks));
    if(!vp) continue;
    sum+=vp->Pt();
  }
  return sum;
}

Double_t AliAnalysisTaskRhoMassSparse::GetMd(AliEmcalJet *jet) {
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
      Double_t m = 0.;
      if(fPionMassClusters) m = 0.13957;
      if(fJetRhoMassType==kMd) sum += TMath::Sqrt(m*m + nPart.Pt()*nPart.Pt()) - nPart.Pt();
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
