/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */
//-------------------------------------------------------------------------
//     Analysis Task that uses the Soft K-Means Algorithm to find clusters in
//     the eta-phi space of Minimum Bias. No pt information is used for the clustering.
//     
//
//     Author: Andreas Morsch (CERN)
//     andreas.morsch@cern.ch
//-------------------------------------------------------------------------



#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"

#include "TList.h"
#include "TParticle.h"
#include "TParticlePDG.h"
#include "TProfile.h"
#include "TMath.h"
#include "TRandom.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliStack.h"
#include "AliESDVertex.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliMultiplicity.h"

#include "AliMCParticle.h"
#include "AliMCEvent.h"
#include "AliAnalysisTaskKMeans.h"
#include "AliTrackReference.h"
#include "AliStack.h"
#include "AliHeader.h"
#include "AliKMeansClustering.h"



ClassImp(AliAnalysisTaskKMeans)

AliAnalysisTaskKMeans::AliAnalysisTaskKMeans() 
    : AliAnalysisTaskSE() 
    ,fHists(0)
    ,fH1CEta(0)
    ,fH1CPhi(0)
    ,fH1CEtaR(0)
    ,fH2N1N2(0)
    ,fH1Pt(0)
    ,fH1PtC(0)
    ,fH1PtC1(0)
    ,fH1PtC2(0)
    ,fH1SPt(0)
    ,fH1SPtC(0)
    ,fH1DPhi(0)
    ,fH1DR(0)
    ,fH1DRR(0)
    ,fH2DPhiEta(0)
    ,fH2DPhiEtaR(0)
    ,fH2DPhiEtaL(0)
    ,fH2DPhiEtaC(0)
    ,fH2DPhiEtaCR(0)
    ,fCuts(0)
{
  //
  // Constructor
  //
}

//________________________________________________________________________
AliAnalysisTaskKMeans::AliAnalysisTaskKMeans(const char *name) 
    : AliAnalysisTaskSE(name) 
      ,fHists(0)
      ,fH1CEta(0)
      ,fH1CPhi(0)
      ,fH1CEtaR(0)
      ,fH2N1N2(0)
      ,fH1Pt(0)
      ,fH1PtC(0)
      ,fH1PtC1(0)
      ,fH1PtC2(0)
      ,fH1SPt(0)
      ,fH1SPtC(0)
      ,fH1DPhi(0)
      ,fH1DR(0)
      ,fH1DRR(0)
      ,fH2DPhiEta(0)
      ,fH2DPhiEtaR(0)
      ,fH2DPhiEtaL(0)
      ,fH2DPhiEtaC(0)
      ,fH2DPhiEtaCR(0)
      ,fCuts(0)
{
  //
  // Constructor
  //
  DefineOutput(1,  TList::Class());
}


//________________________________________________________________________
void AliAnalysisTaskKMeans::UserCreateOutputObjects()
{
  // Create histograms
  // Called once
    fHists   = new TList();
    fH1CEta  = new TH1F("fH1CEta",   "eta distribution of clusters",        90, -1.5, 1.5);
    fH1CEtaR = new TH1F("fH1CEtaR",  "eta distribution of clusters",        90, -1.5, 1.5);
    fH1CPhi  = new TH1F("fH1CPhi",   "phi distribution of clusters",       157,  0.0, 2. * TMath::Pi());
    fH2N1N2  = new TH2F("fH2N1N2",   "multiplicity distribution",          50, 0., 50., 50, 0., 50.);
    
    fH1Pt    = new TH1F("fH1Pt",     "pt distribution",50, 0., 10.);
    fH1PtC   = new TH1F("fH1PtC",    "pt distribution",50, 0., 10.);
    fH1PtC1  = new TH1F("fH1PtC1",   "pt distribution",50, 0., 10.);
    fH1PtC2  = new TH1F("fH1PtC2",   "pt distribution",50, 0., 10.);

    fH1SPt    = new TH1F("fH1SPt",     "sum pt distribution",50, 0., 10.);
    fH1SPtC   = new TH1F("fH1SPtC",    "sum pt distribution",50, 0., 10.);

    fH1DR        = new TH1F("fH1DR",     "dR distribution", 50, 0., 5.);
    fH1DPhi      = new TH1F("fH1DPhi",   "dPhi distribution", 31, 0., TMath::Pi());
    fH2DPhiEta   = new TH2F("fH2DPhiEta","eta phi distribution", 31, 0., TMath::Pi(), 20, 0., 2.);
    fH2DPhiEtaR  = new TH2F("fH2DPhiEtaR","eta phi distribution", 31, 0., TMath::Pi(), 20, 0., 2.);
    fH2DPhiEtaL  = new TH2F("fH2DPhiEtaL","eta phi distribution", 31, 0., TMath::Pi(), 20, 0., 2.);
    fH2DPhiEtaC  = new TH2F("fH2DPhiEtaC","eta phi distribution", 31, 0., TMath::Pi(), 20, 0., 2.);
    fH2DPhiEtaCR = new TH2F("fH2DPhiEtaCR","eta phi distribution", 31, 0., TMath::Pi(), 20, 0., 2.);
    fH1DRR       = new TH1F("fH1DRR",    "dR distribution", 50, 0., 5.);
    
    fHists->SetOwner();

    fHists->Add(fH1CEta);
    fHists->Add(fH1CEtaR);
    fHists->Add(fH1CPhi);
    fHists->Add(fH2N1N2);
    fHists->Add(fH1Pt);
    fHists->Add(fH1PtC);
    fHists->Add(fH1PtC1);
    fHists->Add(fH1PtC2);
    fHists->Add(fH1DR);
    fHists->Add(fH1SPtC);
    fHists->Add(fH1SPt);
    fHists->Add(fH1DPhi);
    fHists->Add(fH2DPhiEta);
    fHists->Add(fH2DPhiEtaR);
    fHists->Add(fH2DPhiEtaL);
    fHists->Add(fH2DPhiEtaC);
    fHists->Add(fH2DPhiEtaCR);
    fHists->Add(fH1DRR);
    
    //
    AliKMeansClustering::SetBeta(20.);
    
}

//________________________________________________________________________
void AliAnalysisTaskKMeans::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event
    Double_t   phi [500];
    Double_t   phiR[500];
    Double_t   eta[500];
    Double_t   etaR[500];
    Double_t   pt [500];
    Double_t   mPhi[20];
    Double_t   mEta[20];
    Double_t     rk[20];
    Int_t       ind[20];

  if (!fInputEvent) {
    Printf("ERROR: fESD not available");
    return;
  }
  
  AliESDEvent* esdE = (AliESDEvent*) fInputEvent;

  Int_t ic = 0;
  //
  // Fill eta-phi positions
  //
  /*
  const AliMultiplicity *spdMult = esdE->GetMultiplicity();
  for (Int_t i = 0; i < spdMult->GetNumberOfTracklets(); ++i) {
	  phi[ic] = spdMult->GetPhi(i);
	  eta[ic] = spdMult->GetEta(i);
	  ic++;
  }
  */
  Double_t ptmax = 0.;
  Int_t    icmax = -1;
  
  for (Int_t iTracks = 0; iTracks < esdE->GetNumberOfTracks(); iTracks++) {
      AliESDtrack* track = esdE->GetTrack(iTracks);
      if ((fCuts->AcceptTrack(track))) 
      {
	  phi[ic]  = track->Phi();
	  phiR[ic] = 2. * TMath::Pi() * gRandom->Rndm();
	  eta[ic]  = track->Eta();
	  etaR[ic] = 2. * gRandom->Rndm() - 1.;
	  pt[ic]   = track->Pt();
	  if (pt[ic] > ptmax) {
	      ptmax = pt[ic];
	      icmax = ic;
	  }
	  ic++;
      }
  } //track loop 

  //
  // Cluster
  if (ic < 10) {
      PostData(1, fHists);
      return;
  }
  //
  AliKMeansClustering::SoftKMeans(2, ic, phi, eta, mPhi, mEta, rk);
  //
  // Sort
  TMath::Sort(2, rk, ind);
  //
  // Analyse
  //
  // Cluster Multiplicity
  Int_t mult[2] = {0, 0};
  //
  Double_t dphic = DeltaPhi(mPhi[0], mPhi[1]);
  Double_t detac = TMath::Abs(mEta[0] - mEta[1]);
  fH2DPhiEtaC->Fill(dphic, detac);  

  //
  Double_t sumPt  = 0;
  Double_t sumPtC = 0;
  for (Int_t i = 0; i < 1; i++) {
      fH1CEta->Fill(mEta[ind[i]]);
      fH1CPhi->Fill(mPhi[ind[i]]);      
      for (Int_t j = 0; j < ic; j++) {
	  Double_t r    = DeltaR(mPhi[ind[i]], mEta[ind[i]], phi[j], eta[j]);
	  Double_t dphi = DeltaPhi(mPhi[ind[i]], phi[j]);
	  Double_t deta = mEta[ind[i]] - eta[j];
	  
	  fH1DR->Fill(r);
	  fH1DPhi->Fill(dphi);
	  fH2DPhiEta->Fill(dphi, TMath::Abs(deta));
	  if (j == icmax) fH2DPhiEtaL->Fill(dphi, TMath::Abs(deta));

	  if (r < 0.2) {
	      fH1PtC2->Fill(pt[j]);
	  }
	  if (r < 0.3) 
	  {
	      fH1PtC1->Fill(pt[j]);
	  }
	  if (r < 0.4) 
	  {
	    sumPtC += pt[j];
	    mult[i]++;
	    fH1PtC->Fill(pt[j]);
	  } 
	  if (r > 0.7) {
	      fH1Pt->Fill(pt[j]);
	  }

	  if (r > 0.7 && r < 1.) {
	    sumPt += pt[j];
	  }
      }
  }

  fH2N1N2->Fill(Float_t(mult[0]), Float_t(mult[1]));
  fH1SPt ->Fill(sumPt);
  fH1SPtC->Fill(sumPtC);
  

  // Randomized phi
  //
  AliKMeansClustering::SoftKMeans(2, ic, phiR, etaR, mPhi, mEta, rk);
  //
  // Sort
  TMath::Sort(2, rk, ind);
  //
  // Analyse
  //
  // Cluster Multiplicity
  for (Int_t i = 0; i < 1; i++) {
      fH1CEtaR->Fill(mEta[ind[i]]);
  }
  for (Int_t i = 0; i < 2; i++) {
      for (Int_t j = 0; j < ic; j++) {
	  Double_t dphi = DeltaPhi(mPhi[ind[i]], phi[j]);
	  Double_t deta = mEta[ind[i]] - eta[j];
	  Double_t r = DeltaR(mPhi[ind[i]], mEta[ind[i]], phi[j], eta[j]);
	  fH1DRR->Fill(r);
	  fH2DPhiEtaR->Fill(dphi, TMath::Abs(deta));
      }
  }
  dphic = DeltaPhi(mPhi[0], mPhi[1]);
  detac = TMath::Abs(mEta[0] - mEta[1]);
  fH2DPhiEtaCR->Fill(dphic, detac);  

  
  //
  // Post output data.
  PostData(1, fHists);
}      

Double_t AliAnalysisTaskKMeans::DeltaPhi(Double_t phi1, Double_t phi2)
{
    Double_t dphi = TMath::Abs(phi1 - phi2);
    if (dphi > TMath::Pi()) dphi = 2. * TMath::Pi() - dphi;
    return dphi;
}

Double_t AliAnalysisTaskKMeans::DeltaR(Double_t phi1, Double_t eta1, Double_t phi2, Double_t eta2)
{
    Double_t dphi = DeltaPhi(phi1, phi2);
    Double_t deta = eta1 - eta2;
    return (TMath::Sqrt(dphi * dphi + deta * deta));
    
}

//________________________________________________________________________
void AliAnalysisTaskKMeans::Terminate(Option_t *) 
{
}  

