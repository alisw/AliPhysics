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
//---------------------------------------------------------------------------------------
//     Analysis Task that uses the Soft K-Means Algorithm to find clusters in
//     the eta-phi space of Minimum Bias. No pt information is used for the clustering.
//     
//
//     Author: Andreas Morsch (CERN)
//     andreas.morsch@cern.ch
//---------------------------------------------------------------------------------------



#include "TChain.h"
#include "TFile.h"
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
#include "TObjArray.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliExternalTrackParam.h"
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
    ,fK(0)
    ,fNMin(0)
    ,fHists(0)
    ,fH1CEta(0)
    ,fH1CPhi(0)
    ,fH1CEtaR(0)
    ,fH2N1N2(0)
    ,fH1Pt(0)
    ,fH1PtC(0)
    ,fH1PtC1(0)
    ,fH1PtC2(0)
    ,fH1PtAS(0)
    ,fH1PtR(0)
    ,fH1SPt(0)
    ,fH1SPtC(0)
    ,fH1DPhi(0)
    ,fH1DR(0)
    ,fH1DRR(0)
    ,fH2DPhiEta(0)
    ,fH2DPhiEtaR(0)
    ,fH2DPhiEtaL(0)
    ,fH2DPhiEtaLR(0)
    ,fH2DPhiEtaC(0)
    ,fH2DPhiEtaCR(0)
    ,fH1Resp(0)
    ,fH1RespR(0)
    ,fH2Sigma(0)
    ,fH2SigmaR(0)
    ,fHDensity(0)
    ,fHCSize(0)
    ,fHNCluster(0)
    ,fHPtDensity(0)
    ,fHDPhi(0)
    ,fH2EtaPhi(0)
    ,fH2REtaPhi(0)
    ,fCuts(0)

{
  //

  // Constructor
  //
    for (Int_t i=0; i< 10; i++) {
	fA[i] = 0;
	fB[i] = 0;
    }
}

//________________________________________________________________________
AliAnalysisTaskKMeans::AliAnalysisTaskKMeans(const char *name) 
    : AliAnalysisTaskSE(name) 
      ,fK(0)
      ,fNMin(0)
      ,fHists(0)
      ,fH1CEta(0)
      ,fH1CPhi(0)
      ,fH1CEtaR(0)
      ,fH2N1N2(0)
      ,fH1Pt(0)
      ,fH1PtC(0)
      ,fH1PtC1(0)
      ,fH1PtC2(0)
      ,fH1PtAS(0)
      ,fH1PtR(0)
      ,fH1SPt(0)
      ,fH1SPtC(0)
      ,fH1DPhi(0)
      ,fH1DR(0)
      ,fH1DRR(0)
      ,fH2DPhiEta(0)
      ,fH2DPhiEtaR(0)
      ,fH2DPhiEtaL(0)
      ,fH2DPhiEtaLR(0)
      ,fH2DPhiEtaC(0)
      ,fH2DPhiEtaCR(0)
      ,fH1Resp(0)
      ,fH1RespR(0)
      ,fH2Sigma(0)
      ,fH2SigmaR(0)
      ,fHDensity(0)
      ,fHCSize(0)
      ,fHNCluster(0)
      ,fHPtDensity(0)
      ,fHDPhi(0)
      ,fH2EtaPhi(0)
      ,fH2REtaPhi(0)
      ,fCuts(0)

{
  //
  // Constructor
  //
    for (Int_t i=0; i< 10; i++) {
	fA[i] = 0;
	fB[i] = 0;
    }
  DefineOutput(1,  TList::Class());
}

AliAnalysisTaskKMeans::AliAnalysisTaskKMeans(const AliAnalysisTaskKMeans &res)
: AliAnalysisTaskSE(res) 
      ,fK(0)
      ,fNMin(0)
      ,fHists(0)
      ,fH1CEta(0)
      ,fH1CPhi(0)
      ,fH1CEtaR(0)
      ,fH2N1N2(0)
      ,fH1Pt(0)
      ,fH1PtC(0)
      ,fH1PtC1(0)
      ,fH1PtC2(0)
      ,fH1PtAS(0)
      ,fH1PtR(0)
      ,fH1SPt(0)
      ,fH1SPtC(0)
      ,fH1DPhi(0)
      ,fH1DR(0)
      ,fH1DRR(0)
      ,fH2DPhiEta(0)
      ,fH2DPhiEtaR(0)
      ,fH2DPhiEtaL(0)
      ,fH2DPhiEtaLR(0)
      ,fH2DPhiEtaC(0)
      ,fH2DPhiEtaCR(0)
      ,fH1Resp(0)
      ,fH1RespR(0)
      ,fH2Sigma(0)
      ,fH2SigmaR(0)
      ,fHDensity(0)
      ,fHCSize(0)  
      ,fHNCluster(0)
      ,fHPtDensity(0)
      ,fHDPhi(0)
      ,fH2EtaPhi(0)
      ,fH2REtaPhi(0)
      ,fCuts(0)
{
    // Dummy copy constructor
    for (Int_t i=0; i< 10; i++) {
	fA[i] = 0;
	fB[i] = 0;
    }
}

AliAnalysisTaskKMeans& AliAnalysisTaskKMeans::operator=(const AliAnalysisTaskKMeans& /*trclass*/)
{
    // Dummy assignment operator
    return *this;
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
    fH1PtAS  = new TH1F("fH1PtAS",   "pt distribution",50, 0., 10.);
    fH1PtR   = new TH1F("fH1PtR",    "pt distribution",50, 0., 10.);

    fH1SPt    = new TH1F("fH1SPt",     "sum pt distribution",50, 0., 10.);
    fH1SPtC   = new TH1F("fH1SPtC",    "sum pt distribution",50, 0., 10.);

    fH1DR         = new TH1F("fH1DR",     "dR distribution", 50, 0., 5.);
    fH1DPhi       = new TH1F("fH1DPhi",   "dPhi distribution", 31, 0., TMath::Pi());
    fH2DPhiEta    = new TH2F("fH2DPhiEta","eta phi distribution", 31, 0., TMath::Pi(), 20, 0., 2.);
    fH2DPhiEtaR   = new TH2F("fH2DPhiEtaR","eta phi distribution", 31, 0., TMath::Pi(), 20, 0., 2.);
    fH2DPhiEtaL   = new TH2F("fH2DPhiEtaL","eta phi distribution", 31, 0., TMath::Pi(), 20, 0., 2.);
    fH2DPhiEtaLR  = new TH2F("fH2DPhiEtaLR","eta phi distribution", 31, 0., TMath::Pi(), 20, 0., 2.);
    fH2DPhiEtaC   = new TH2F("fH2DPhiEtaC","eta phi distribution", 31, 0., TMath::Pi(), 20, 0., 2.);
    fH2DPhiEtaCR  = new TH2F("fH2DPhiEtaCR","eta phi distribution", 31, 0., TMath::Pi(), 20, 0., 2.);
    fH1DRR        = new TH1F("fH1DRR",    "dR distribution", 50, 0., 5.);
    fH1Resp       = new TH1F("fH1Resp",   "Responsibility", 50, 0., 1.);
    fH1RespR      = new TH1F("fH1RespR",  "Responsibility", 50, 0., 1.);
    fH2Sigma      = new TH2F("fH2Sigma",  "Sigma", 31, 0., TMath::Pi(), 20, 0., 2.);
    fH2SigmaR     = new TH2F("fH2SigmaR", "Sigma", 31, 0., TMath::Pi(), 20, 0., 2.);


    fHDensity    = new TH1F("fHDensity",    "density distribution", 100, 0., 20.);
    fHCSize      = new TH1F("fHCSize",      "cluster size", 20, -0.5, 19.5);

    fHNCluster   = new TH1F("fHNCluster",   "Number of clusters", 11, -0.5, 10.5);
    fHPtDensity  = new TH2F("fHPtDensity", "Pt vs density", 100, 0., 20., 50, 0., 10.);
    fHDPhi       = new TH1F("fHDPhi",   "phi correlation", 100, 0., 0.5);
    fH2EtaPhi    = new TH2F("fH2EtaPhi", "Eta Phi", 200, -1., 1., 628, 0., 2. * TMath::Pi());  
    fHists->SetOwner();

    fHists->Add(fH1CEta);
    fHists->Add(fH1CEtaR);
    fHists->Add(fH1CPhi);
    fHists->Add(fH2N1N2);
    fHists->Add(fH1Pt);
    fHists->Add(fH1PtR);
    fHists->Add(fH1PtC);
    fHists->Add(fH1PtC1);
    fHists->Add(fH1PtC2);
    fHists->Add(fH1PtAS);
    fHists->Add(fH1DR);
    fHists->Add(fH1SPtC);
    fHists->Add(fH1SPt);
    fHists->Add(fH1DPhi);
    fHists->Add(fH2DPhiEta);
    fHists->Add(fH2DPhiEtaR);
    fHists->Add(fH2DPhiEtaL);
    fHists->Add(fH2DPhiEtaLR);
    fHists->Add(fH2DPhiEtaC);
    fHists->Add(fH2DPhiEtaCR);
    fHists->Add(fH1DRR);
    fHists->Add(fH1RespR);
    fHists->Add(fH1Resp);    
    fHists->Add(fH2Sigma);    
    fHists->Add(fH2SigmaR);
    fHists->Add(fHCSize);    
    fHists->Add(fHDensity);
    fHists->Add(fHNCluster);
    fHists->Add(fHPtDensity);
    fHists->Add(fHDPhi);
    fHists->Add(fH2EtaPhi);
    //
    for (Int_t i = 0; i < 10; i++) {
      fA[i] = new AliKMeansResult(i+1);
      fB[i] = new AliKMeansResult(i+1);
    }
}

//________________________________________________________________________
void AliAnalysisTaskKMeans::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event
    Double_t   phi [500] = {0};
    Double_t   phiR[500] = {0};
    Double_t    eta[500] = {0};
    Double_t   etaR[500] = {0};
    Double_t    pt [500] = {0};

  if (!fInputEvent) {
    Printf("ERROR: fESD not available");
    return;
  }
  

  Int_t ic = 0;
  
  //
  // Fill eta-phi positions
  //

  Double_t ptmax = 0.;
  Int_t    icmax = -1;
  
  for (Int_t iTracks = 0; iTracks < fInputEvent->GetNumberOfTracks(); iTracks++) {
      AliVParticle* track = fInputEvent->GetTrack(iTracks);
      if ((fCuts->AcceptTrack((AliESDtrack*)track))) 
      {
	const AliExternalTrackParam * tpcT = ((AliESDtrack*) track)->GetTPCInnerParam();
	if (!tpcT) continue;
	if (TMath::Abs(tpcT->Eta()) > 0.9) continue;

	  phi [ic] = tpcT->Phi();
	  eta [ic] = tpcT->Eta();
	  pt  [ic] = tpcT->Pt();

	  if (fH2REtaPhi) {
	    fH2REtaPhi->GetRandom2(etaR[ic], phiR[ic]);
	  } else {
	    phiR[ic] = 2. * TMath::Pi() * gRandom->Rndm();
	    etaR[ic] = 1.8 * gRandom->Rndm() - 0.9;
	  }
	  

	  if (pt[ic] > ptmax) {
	      ptmax = pt[ic];
	      icmax = ic;
	  }

	  fH2EtaPhi->Fill(eta[ic], phi[ic]);
	  ic++;
      }
  } //track loop 

  for (Int_t i = 0; i < ic; i++) {
    for (Int_t j = i+1; j < ic; j++) {
      Double_t dphi = TMath::Abs(phi[i] - phi[j]);
      fHDPhi->Fill(dphi);
    }
  }

  //
  // Cluster
  if (ic < fNMin) {
      PostData(1, fHists);
      return;
  }

  //
  Double_t rk0[10];
  AliKMeansResult* res = 0;
  AliKMeansResult best(10);
  Float_t   rmaxG = -1.;

  for (Int_t k = 0; k < 20; k++) {
    Float_t   rmax   = -1.;
    Int_t     imax   = 0;
    for (Int_t i = 0; i < fK; i++) {
      res = fA[i];
      AliKMeansClustering::SoftKMeans2(i+1, ic, phi, eta, res->GetMx(),res->GetMy(), res->GetSigma2(), res->GetRk());
      res->Sort(ic, phi, eta);
      Int_t j = (res->GetInd())[0];
      rk0[i]  = (res->GetTarget())[j];
      if (rk0[i] > rmax)  {
	rmax = rk0[i];
	imax = i;
      }
    }
    if (rmax > rmaxG) {
      rmaxG = rmax;
      best.CopyResults(fA[imax]);
    }
  }
  Double_t* mPhi     = best.GetMx();
  Double_t* mEta     = best.GetMy();
  Double_t* sigma2   = best.GetSigma2();
  Int_t     nk       = best.GetK();
  Int_t     im       = (best.GetInd())[0];
  Double_t  etaC     = mEta[im];

  fHDensity->Fill(rmaxG / TMath::Pi());
  fHCSize->Fill(2.28 * rmaxG * sigma2[im]);
  fHNCluster->Fill(Float_t(nk));

  Double_t dphic, detac;

  if (rmaxG > 0. && TMath::Abs(etaC) < 0.4) {
    // Analysis
    //
    // Cluster Multiplicity
    Int_t mult[2] = {0, 0};
    //
    if (nk > 1) {
      dphic = DeltaPhi(mPhi[0], mPhi[1]);
      detac = TMath::Abs(mEta[0] - mEta[1]);
      fH2DPhiEtaC->Fill(dphic, detac);  
    }
    //
    // Random cluster position
    Int_t ir = Int_t(Float_t(ic) * gRandom->Rndm());

    Double_t crPhi = phi[ir];
    Double_t crEta = eta[ir];
    //
    Double_t sumPt  = 0;
    Double_t sumPtC = 0;

    for (Int_t i = 0; i < 1; i++) {
      fH1CEta->Fill(mEta[im]);
      fH1CPhi->Fill(mPhi[im]);      
      for (Int_t j = 0; j < ic; j++) {
	Double_t r    = DeltaR(mPhi[im], mEta[im], phi[j], eta[j]);
	Double_t dphi = DeltaPhi(mPhi[im], phi[j]);
	Double_t deta = mEta[im] - eta[j];
	Double_t rr   = DeltaR(crPhi, crEta, phi[j], eta[j]);

	fH1DR->Fill(r);
	fH1DPhi->Fill(dphi);
	fH2DPhiEta->Fill(dphi, TMath::Abs(deta));
	if (j == icmax) fH2DPhiEtaL->Fill(dphi, TMath::Abs(deta));
	
	if (r < 0.2) {
	  fH1PtC2->Fill(pt[j]);
	}
	if (r < 0.3) {
	  fH1PtC1->Fill(pt[j]);
	  fHPtDensity->Fill(rmaxG/TMath::Pi(), pt[j]);
	}
	if (rr < 0.3) {
	    fH1PtR->Fill(pt[j]);
	}

	if (r < 0.4)  {
	  sumPtC += pt[j];
	  mult[i]++;
	  fH1PtC->Fill(pt[j]);
	} 
	if (r > 0.7 && dphi < (TMath::Pi() - 0.3)) {
	  fH1Pt->Fill(pt[j]);
	}
	
	if (r > 0.7 && r < 1.) {
	  sumPt += pt[j];
	}
	
	if (dphi > (TMath::Pi() - 0.3)) {
	  fH1PtAS->Fill(pt[j], 1.);
	}
      }
    }

    fH2N1N2->Fill(Float_t(mult[0]), Float_t(mult[1]));
    fH1SPt ->Fill(sumPt);
    fH1SPtC->Fill(sumPtC);
  }
    //
    // Randomized phi
    //
  rmaxG = -1.;
  for (Int_t k = 0; k < 20; k++) {
    Float_t rmax   = -1.;
    Int_t   imax   =  0;
    for (Int_t i = 0; i < fK; i++) {
      res = fB[i];
      AliKMeansClustering::SoftKMeans2(i+1, ic, phiR, etaR, res->GetMx(),res->GetMy(), res->GetSigma2(), res->GetRk());
      res->Sort(ic, phiR, etaR);
      Int_t j = (res->GetInd())[0];
      rk0[i]  = (res->GetTarget())[j];
      if (rk0[i] > rmax) {
	rmax = rk0[i];
	imax = i;
      }
    }
    if (rmax > rmaxG) {
      rmaxG = rmax;
      best.CopyResults(fB[imax]);
    }    
  }
  
    mPhi    = best.GetMx();
    mEta    = best.GetMy();
    nk      = best.GetK();
    im      = (best.GetInd())[0];
    etaC    = mEta[im];    

    //
    // Cluster Multiplicity
    if (rmaxG > 0. && TMath::Abs(etaC) < 0.4) {
      for (Int_t i = 0; i < 1; i++) {
	im = (best.GetInd())[i];
	fH1CEtaR->Fill(mEta[im]);
      }
      
      for (Int_t i = 0; i < 1; i++) {
	im = (best.GetInd())[i];
	for (Int_t j = 0; j < ic; j++) {
	  Double_t dphi = DeltaPhi(mPhi[im], phiR[j]);
	  Double_t deta = mEta[im] - etaR[j];
	  Double_t r = DeltaR(mPhi[im], mEta[im], phiR[j], etaR[j]);
	  fH1DRR->Fill(r);
	  fH2DPhiEtaR->Fill(dphi, TMath::Abs(deta));
	  if (j == icmax) fH2DPhiEtaLR->Fill(dphi, TMath::Abs(deta));
	}
      }
      if (nk > 1) {
	dphic = DeltaPhi(mPhi[0], mPhi[1]);
	detac = TMath::Abs(mEta[0] - mEta[1]);
	fH2DPhiEtaCR->Fill(dphic, detac);  
      }
    }
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

