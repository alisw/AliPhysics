#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TClonesArray.h"
#include "TCanvas.h"
#include "TList.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"

#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODJet.h"
#include "AliAnalysisTaskJetsTM.h"

ClassImp(AliAnalysisTaskJetsTM)
//
// Thrust Major (TM) analysis of reconstructed jets.
// TM is the thrust in the plane perpendicular to the jet axis
// The present amalysis performs the following steps:
// (a) Construct to orthogonal unit vectors (e1, e2) in the plane perpendicular to the jet axis
// (b) Calculate the components of all particles with jT > 1 GeV with respect to e1, e2
// (c) Construct the sphericity matrix
// (d) Find the two orthogonal eigenvectors of the spericity matrix
// (e) Caluclate the components of all particles with jT < 1 GeV in the reference frame spanned by the eigenvectors
// (f) Calculate the azimuthal angle in this frame
//
//
// Author: andreas.morsch@cern.ch


AliAnalysisTaskJetsTM::AliAnalysisTaskJetsTM(const char *name) 
    : AliAnalysisTaskSE(name) 
      ,fHists(0)
      ,fPtH(0)
      ,fPtTH(0)
      ,fPhiM(0)
      ,fPhiMPt(0)
      ,fPhiMPtJ(0)
      ,fPtSum(0)
{
  //
  // Constructor
  //
  DefineOutput(1,  TList::Class());
}


//________________________________________________________________________
void AliAnalysisTaskJetsTM::UserCreateOutputObjects()
{
  // Create histograms
  // Called once
    fHists    = new TList();
    fPtH      = new TH1F("fPtH" ,     " p_{T}       distribution of jets"  ,            50, 0., 100.0);
    fPtSum    = new TH2F("fPtSum" ,   " p_{T} sum in cone"                 ,            50, 0., 100.0, 50, 0, 100.);
    fPtTH     = new TH1F("fPtTH",     " p_{T}       distribution of tracks",            50, 0., 100.0);
    fPhiM     = new TH1F("fPhiM",     " phi^{major} distribution of tracks",            64, 0.,   3.2);
    fPhiMPt   = new TH2F("fPhiMPt",   " phi^{major} distribution of tracks vs pt",      64, 0.,   3.2, 20, 0., 2.);
    fPhiMPtJ  = new TH2F("fPhiMPtJ",  " phi^{major} distribution of tracks vs pt Jet",  64, 0.,   3.2, 40, 0., 100.);

    fHists->SetOwner();
    fHists->Add(fPtH);
    fHists->Add(fPtSum);
    fHists->Add(fPtTH);
    fHists->Add(fPhiM);
    fHists->Add(fPhiMPt);
    fHists->Add(fPhiMPtJ);
}


//________________________________________________________________________
void AliAnalysisTaskJetsTM::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event

  if (!fInputEvent) {
    Printf("ERROR: AOD not available");
    return;
  }

  AliAODEvent* aodE  = dynamic_cast<AliAODEvent*>  (fInputEvent);

  if (!aodE) {
    Printf("ERROR: AOD not available");
    return;
  }

  TClonesArray* jets = dynamic_cast<TClonesArray*> (aodE->FindListObject("jetsAOD_FASTKT04"));
  if (!jets) {
    Printf("ERROR: Jet branch not available");
    return;
  }

  Int_t nJ = jets->GetEntries();

  Float_t ptmax = 0.;
  Int_t   imax  = -1;
  
//
// Find highest pT jet with pt > 20 GeV
//
  for (Int_t i = 0; i < nJ; i++) {
      AliAODJet* jet = dynamic_cast<AliAODJet*> (jets->At(i));
      if (!jet) continue;
      Float_t ptJ  = jet->Pt();
      Float_t etaJ = TMath::Abs(jet->Eta());
      if ((ptJ > 20.) && (ptJ  > ptmax) && etaJ < 0.5) {
	  ptmax = ptJ;
	  imax = i;
      }
  }

  if (imax == -1) return;
  

  AliAODJet* jet = dynamic_cast<AliAODJet*> (jets->At(imax));
  Float_t ptJ  = jet->Pt();
  Float_t phiJ = jet->Phi();
  Float_t etaJ = jet->Eta();
  
  fPtH->Fill(ptJ);

//
// The transverse plane
//
// 2 normalized vectors perpendicular to the jet direction
//
      Float_t pxJ = jet->Px();
      Float_t pyJ = jet->Py();
      Float_t pzJ = jet->Pz();
      
      TVector3  ppJ1(pxJ, pyJ, pzJ);
      TVector3  ppJ3(- pxJ * pzJ, - pyJ * pzJ, pxJ * pxJ + pyJ * pyJ);
      ppJ3.SetMag(1.);
      TVector3  ppJ2(-pyJ, pxJ, 0);
      ppJ2.SetMag(1.);

//
// 1st track loop to determine the sphericity matrix
//   
      Int_t nT = aodE->GetNumberOfTracks();
      // Initialize Shericity matrix
      Float_t mxx    = 0.;
      Float_t myy    = 0.;
      Float_t mxy    = 0.;
      Int_t   nc     = 0;
      Float_t sump2  = 0.;
      Float_t ptMax  = 0.;
      Float_t etaMax = 0.;
      Float_t phiMax = 0.;
      Int_t   iMax   = -1;
      Float_t ptsum  = 0.;
      
      for (Int_t i = 0; i < nT; i++) {
	  AliAODTrack* track = dynamic_cast<AliAODTrack*>(aodE->GetTrack(i));
	  if(!track) AliFatal("Not a standard AOD");
//
//    Track quality cuts
	  if (!track->TestFilterBit(1<<4)) continue;
	  fPtTH->Fill(track->Pt());

//
//
	  TVector3 pp(track->Px(), track->Py(), track->Pz());
	      
	  Float_t phi = track->Phi();
	  Float_t eta = track->Eta();
	  Float_t pt  = track->Pt();
	  Float_t jT  = pp.Perp(ppJ1);
	  
	  
	  // pT > 2.
	  
	    
	  Float_t deta = eta - etaJ;
	  Float_t dphi = phi - phiJ;

	  if (dphi >   TMath::Pi()) dphi =   2. * TMath::Pi() - dphi;
	  if (dphi < - TMath::Pi()) dphi = - 2. * TMath::Pi() - dphi;
	  Float_t r = TMath::Sqrt(dphi * dphi + deta * deta);

	  if (r < 0.4) ptsum += pt;
	  if (r < 0.5 && jT > 1.) {
	      TVector3 pLong = pp.Dot(ppJ1) / ppJ1.Mag2() * ppJ1;
	      TVector3 pPerp = pp - pLong;

	      Float_t ppjX = pPerp.Dot(ppJ2);
	      Float_t ppjY = pPerp.Dot(ppJ3);
	      Float_t ppjT = TMath::Sqrt(ppjX * ppjX + ppjY * ppjY);

	      mxx += (ppjX * ppjX / ppjT);
	      myy += (ppjY * ppjY / ppjT);
	      mxy += (ppjX * ppjY / ppjT);
	      nc++;
	      sump2 += ppjT;
	      // max pt
	      if (pt > ptMax) {
		  ptMax  = pt;
		  iMax   = i;
		  etaMax = deta;
		  phiMax = dphi;
	      }
	  } // R < 0.4
      } // 1st Track Loop
      fPtSum->Fill(ptJ, ptsum);
      
//
// At this point we have mxx, myy, mxy
      if (nc == 0) return;      
// Shericity Matrix	
      const Double_t ele[4] = {mxx / sump2, mxy / sump2, mxy / sump2, myy / sump2};	
      TMatrixDSym m0(2, ele);
// Find eigenvectors
      TMatrixDSymEigen m(m0);
      TVectorD eval(2);
      TMatrixD evecm = m.GetEigenVectors();
      eval  = m.GetEigenValues();
// Largest eigenvector
      Int_t jev = 0;
      if (eval[0] < eval[1]) jev = 1;
      TVectorD evec0(2);
// Principle axis
      evec0 = TMatrixDColumn(evecm, jev);
      TVector2 evec(evec0[0], evec0[1]); 
// Principle axis from leading partice
      Float_t phiM = TMath::ATan2(phiMax, etaMax);
      TVector2 evecM(TMath::Cos(phiM), TMath::Sin(phiM)); 
      Float_t phistM = evecM.DeltaPhi(evec);
      if (TMath::Abs(phistM) > TMath::Pi()/2.) evec*=(-1.);

//
// 2nd correlation with principal axis
//   
      for (Int_t i = 0; i < nT; i++) {
	  AliAODTrack* track = dynamic_cast<AliAODTrack*>(aodE->GetTrack(i));
	  if(!track) AliFatal("Not a standard AOD");
//
//    Track quality cuts
	  if (!track->TestFilterBit(1<<4)) continue;
//
//
	  TVector3 pp(track->Px(), track->Py(), track->Pz());
	  Float_t phi = track->Phi();
	  Float_t eta = track->Eta();
	  Float_t pt  = track->Pt();
	  Float_t jT  = pp.Perp(ppJ1);
	  // pT < 2.
	  if (jT > 1.) continue;	    
	  Float_t deta = eta - etaJ;
	  Float_t dphi = phi - phiJ;
	  if (dphi >   TMath::Pi()) dphi =   2. * TMath::Pi() - dphi;
	  if (dphi < - TMath::Pi()) dphi = - 2. * TMath::Pi() - dphi;

	  Float_t r = TMath::Sqrt(dphi * dphi + deta * deta);

	  if (r > 0.7) continue;

	  TVector3 pLong = pp.Dot(ppJ1) / ppJ1.Mag2() * ppJ1;
	  TVector3 pPerp = pp - pLong;
	  
	  Float_t ppjX = pPerp.Dot(ppJ2);
	  Float_t ppjY = pPerp.Dot(ppJ3);
	    
	  TVector2 vr(ppjX, ppjY) ;
	  Float_t phistr = evec.DeltaPhi(vr);
	  fPhiM->Fill(phistr);
	  fPhiMPt ->Fill(phistr, pt);
	  fPhiMPtJ->Fill(phistr, ptJ);
	  
      } // 2nd Track loop
      
  
  // Post output data.
  PostData(1, fHists);
}      

//________________________________________________________________________
void AliAnalysisTaskJetsTM::Terminate(Option_t *) 
{
// Terminate
}  
