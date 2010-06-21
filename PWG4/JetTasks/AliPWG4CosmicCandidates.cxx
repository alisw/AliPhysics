/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

//-----------------------------------------------------------------------
// Analysis task looking for cosmic candidates
// Task checks if particles are back-to-back in eta and phi
//             
//
// Author : Marta Verweij - UU - marta.verweij@cern.ch
//-----------------------------------------------------------------------

#include "TVector3.h"
#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TList.h"
#include "TChain.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"

#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliExternalTrackParam.h"

#include "AliLog.h"

#include "AliPWG4CosmicCandidates.h"

//using namespace std; //required for resolving the 'cout' symbol
using namespace std;

ClassImp(AliPWG4CosmicCandidates)

//________________________________________________________________________
AliPWG4CosmicCandidates::AliPWG4CosmicCandidates()
: AliAnalysisTaskSE(),
  fTrackCuts(0), 
  fPtMin(5.),
  fMaxCosmicAngle(0.002),
  fNEventAll(0),
  fNEventSel(0),
  fPtSignedCosmicCandidates(0),
  fDeltaPtCosmicCandidates(0),
  fDeltaPhiSumEta(0),
  fDCAZCosmicCandidates(0),
  fDCARCosmicCandidates(0),
  fTheta(0),
  fThetaZoom(0),
  fThetaPt1Pt2(0),
  fThetaPt1Pt2Signed(0),
  fDeltaPhiSumEtaPt1(0),
  fDeltaPhiSumEtaPt2(0),
  fThetaDCAZ1DCAZ2(0),
  fRisol(0),
  fRisolTheta(0),
  fHistListCosmics(0)
{
  //
  // Default constructor
  //
}

//________________________________________________________________________
AliPWG4CosmicCandidates::AliPWG4CosmicCandidates(const char *name)
  : AliAnalysisTaskSE(name),
    fTrackCuts(0), 
    fPtMin(5.),
    fMaxCosmicAngle(0.002),
    fNEventAll(0),
    fNEventSel(0),
    fPtSignedCosmicCandidates(0),
    fDeltaPtCosmicCandidates(0),
    fDeltaPhiSumEta(0),
    fDCAZCosmicCandidates(0),
    fDCARCosmicCandidates(0),
    fTheta(0),
    fThetaZoom(0),
    fThetaPt1Pt2(0),
    fThetaPt1Pt2Signed(0),
    fDeltaPhiSumEtaPt1(0),
    fDeltaPhiSumEtaPt2(0),
    fThetaDCAZ1DCAZ2(0),
    fRisol(0),
    fRisolTheta(0),
    fHistListCosmics(0)
{
  // Constructor. Initialization of Inputs and Outputs

  // Define input and output slots here
  // Input slot #0 works with a TChain
  //  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TList
  DefineOutput(1, TList::Class());
  // Output slot #2 writes into a AliESDtrackCuts
  DefineOutput(2, AliESDtrackCuts::Class());
}
//________________________________________________________________________
AliPWG4CosmicCandidates::AliPWG4CosmicCandidates(const AliPWG4CosmicCandidates &res)
  : AliAnalysisTaskSE(res),
    fTrackCuts(0), 
    fPtMin(5.),
    fMaxCosmicAngle(0.002),
    fNEventAll(0),
    fNEventSel(0),
    fPtSignedCosmicCandidates(0),
    fDeltaPtCosmicCandidates(0),
    fDeltaPhiSumEta(0),
    fDCAZCosmicCandidates(0),
    fDCARCosmicCandidates(0),
    fTheta(0),
    fThetaZoom(0),
    fThetaPt1Pt2(0),
    fThetaPt1Pt2Signed(0),
    fDeltaPhiSumEtaPt1(0),
    fDeltaPhiSumEtaPt2(0),
    fThetaDCAZ1DCAZ2(0),
    fRisol(0),
    fRisolTheta(0),
    fHistListCosmics(0)
{
    // Dummy copy constructor
}

//________________________________________________________________________
AliPWG4CosmicCandidates& AliPWG4CosmicCandidates::operator=(const AliPWG4CosmicCandidates& /*trclass*/)
{
    // Dummy assignment operator
    return *this;
}

//________________________________________________________________________
void AliPWG4CosmicCandidates::LocalInit()
{
  //
  // Only called once at beginning
  //
  PostData(2,fTrackCuts);
}

//________________________________________________________________________
void AliPWG4CosmicCandidates::UserCreateOutputObjects()
{
  //Create output objects
  // Called once
  AliDebug(2,Form(">> AliPWG4CosmicCandidates::UserCreateOutputObjects \n")); 

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE); 
  
  OpenFile(1);
  fHistListCosmics = new TList();

  fNEventAll = new TH1F("fNEventAll","NEventAll",1,-0.5,0.5);
  fHistListCosmics->Add(fNEventAll);
  fNEventSel = new TH1F("fNEventSel","NEvent Selected for analysis",1,-0.5,0.5);
  fHistListCosmics->Add(fNEventSel);

  Float_t fgkPtMin=0.;
  Float_t fgkPtMax=100.;
  Int_t fgkNPtBins= (int)(fgkPtMax-fgkPtMin);

  Int_t fgkNPhiBins=18;
  Float_t kMinPhi = -0.5*TMath::Pi();
  Float_t kMaxPhi = 3./2.*TMath::Pi();

  Int_t fgkNThetaBins=fgkNPhiBins*8;
  Float_t kMinTheta = -0.5*TMath::Pi();
  Float_t kMaxTheta = 3./2.*TMath::Pi();

  Int_t fgkNDCARBins=80;
  Float_t fgkDCARMin = -0.2;
  Float_t fgkDCARMax = 0.2;
  Double_t *binsDCAR=new Double_t[fgkNDCARBins+1];
  for(Int_t i=0; i<=fgkNDCARBins; i++) binsDCAR[i]=(Double_t)fgkDCARMin + (fgkDCARMax-fgkDCARMin)/fgkNDCARBins*(Double_t)i ;

  Int_t fgkNDCAZBins=80;
  Float_t fgkDCAZMin = -2.;
  Float_t fgkDCAZMax = 2.;
  Double_t *binsDCAZ=new Double_t[fgkNDCAZBins+1];
  for(Int_t i=0; i<=fgkNDCAZBins; i++) binsDCAZ[i]=(Double_t)fgkDCAZMin + (fgkDCAZMax-fgkDCAZMin)/fgkNDCAZBins*(Double_t)i ;

  fPtSignedCosmicCandidates = new TH1F("fPtSignedCosmicCandidates","fPtSignedCosmicCandidates",2*(int)(fgkPtMax-fgkPtMin), -1.*fgkPtMax, fgkPtMax);
  fHistListCosmics->Add(fPtSignedCosmicCandidates);  

  fDeltaPtCosmicCandidates = new TH1F("fDeltaPtCosmicCandidates","fDeltaPtCosmicCandidates",fgkNPtBins, -50., 50.);
  fHistListCosmics->Add(fDeltaPtCosmicCandidates);  

  fDeltaPhiSumEta = new TH2F("fDeltaPhiSumEta","fDeltaPhiSumEta",fgkNPhiBins*4,kMinPhi,kMaxPhi,80, -2.,2.);
  fHistListCosmics->Add(fDeltaPhiSumEta);  

  fDCAZCosmicCandidates = new TH2F("fDCAZCosmicCandidates","fDCAZCosmicCandidates",fgkNDCAZBins,binsDCAZ,fgkNDCAZBins,binsDCAZ);
  fHistListCosmics->Add(fDCAZCosmicCandidates);

  fDCARCosmicCandidates = new TH2F("fDCARCosmicCandidates","fDCARCosmicCandidates",fgkNDCARBins,binsDCAR,fgkNDCARBins,binsDCAR);
  fHistListCosmics->Add(fDCARCosmicCandidates);

  fTheta = new TH1F("fTheta","fTheta",fgkNThetaBins,kMinTheta,kMaxTheta);
  fHistListCosmics->Add(fTheta);

  fThetaZoom = new TH1F("fThetaZoom","fThetaZoom",100,TMath::Pi()-1.,TMath::Pi()+1.);
  fHistListCosmics->Add(fThetaZoom);

  fThetaPt1Pt2 = new TH3F("fThetaPt1Pt2","fThetaPt1Pt2",fgkNThetaBins,kMinTheta,kMaxTheta,(int)(fgkPtMax-fgkPtMin),fgkPtMin,fgkPtMax,(int)(fgkPtMax-fgkPtMin),fgkPtMin,fgkPtMax);
  fHistListCosmics->Add(fThetaPt1Pt2);

  fThetaPt1Pt2Signed = new TH3F("fThetaPt1Pt2Signed","fThetaPt1Pt2Signed",fgkNThetaBins,kMinTheta,kMaxTheta,4*(int)(fgkPtMax-fgkPtMin),-1.*fgkPtMax,fgkPtMax,4*(int)(fgkPtMax-fgkPtMin),-1.*fgkPtMax,fgkPtMax);
  fHistListCosmics->Add(fThetaPt1Pt2Signed);

  fDeltaPhiSumEtaPt1 = new TH3F("fDeltaPhiSumEtaPt1","fDeltaPhiSumEtaPt1",fgkNThetaBins,kMinTheta,kMaxTheta,80, -2.,2.,(int)(fgkPtMax-fgkPtMin),fgkPtMin,fgkPtMax);
  fHistListCosmics->Add(fDeltaPhiSumEtaPt1);

  fDeltaPhiSumEtaPt2 = new TH3F("fDeltaPhiSumEtaPt2","fDeltaPhiSumEtaPt2",fgkNThetaBins,kMinTheta,kMaxTheta,80, -2.,2.,(int)(fgkPtMax-fgkPtMin),fgkPtMin,fgkPtMax);
  fHistListCosmics->Add(fDeltaPhiSumEtaPt2);

  fThetaDCAZ1DCAZ2 = new TH3F("fThetaDCAZ1DCAZ2","fThetaDCAZ1DCAZ2",fgkNThetaBins,kMinTheta,kMaxTheta,fgkNDCAZBins,-2.,2.,fgkNDCAZBins,-2.,2.);
  fHistListCosmics->Add(fThetaDCAZ1DCAZ2);

  fRisol = new TH1F("fRisol","fRisol",100,0.,10.);
  fHistListCosmics->Add(fRisol);

  fRisolTheta = new TH2F("fRisolTheta","fRisolTheta",100,0.,10.,fgkNThetaBins,kMinTheta,kMaxTheta);
  fHistListCosmics->Add(fRisolTheta);

  TH1::AddDirectory(oldStatus);   

  if(binsDCAR) delete [] binsDCAR;
  if(binsDCAZ) delete [] binsDCAZ;

}

//________________________________________________________________________
void AliPWG4CosmicCandidates::UserExec(Option_t *) 
{
//   // Main loop
//   // Called for each event

  // All events without selection
  fNEventAll->Fill(0.);
   
  if (!fInputEvent) {
    AliDebug(2,Form("ERROR: fESD not available"));
    cout << "ERROR: fESD not available" << endl;
    PostData(1, fHistListCosmics);
    return;
  }

  const AliVVertex *vtx = fInputEvent->GetPrimaryVertex();
  // Need vertex cut
  if (!vtx || vtx->GetNContributors() < 2) {
    // Post output data
    PostData(1, fHistListCosmics);
    return;
  }

  //  AliDebug(2,Form("Vertex title %s, status %d, nCont %d\n",vtx->GetTitle(), vtx->GetStatus(), vtx->GetNContributors()));
  double primVtx[3];
  vtx->GetXYZ(primVtx);
  if(TMath::Sqrt(primVtx[0]*primVtx[0] + primVtx[1]*primVtx[1])>1. || TMath::Abs(primVtx[2]>10.)){
    // Post output data
    PostData(1, fHistListCosmics);
    return;
  }
  if(!fInputEvent->GetNumberOfTracks() || fInputEvent->GetNumberOfTracks()<2){ 
    // Post output data
    PostData(1, fHistListCosmics);
    return;
  }
  Int_t nTracks = fInputEvent->GetNumberOfTracks();

  if(!fTrackCuts) {
   // Post output data
    PostData(1, fHistListCosmics);
    return;
  }

 fNEventSel->Fill(0.);

  Float_t dcaR[2] = {0.,0.};
  Float_t dcaZ[2] = {0.,0.};

  // Track loop to fill a pT spectrum
  for (Int_t iTrack1 = 0; iTrack1 < nTracks; iTrack1++) {

    AliESDtrack* track1 = (AliESDtrack*)fInputEvent->GetTrack(iTrack1);
    if (!track1)  continue;
    if(!(fTrackCuts->AcceptTrack(track1))) { continue; }
    if(track1->Pt()<fPtMin) continue;
    //Start 2nd track loop to look for correlations
    for (Int_t iTrack2 = iTrack1+1; iTrack2 < nTracks; iTrack2++) {
      AliESDtrack *track2 = (AliESDtrack*)fInputEvent->GetTrack(iTrack2);
      if(!track2) continue;
      if(!(fTrackCuts->AcceptTrack(track2))) { continue; }
         
      //Check if back-to-back
      Double_t mom1[3],mom2[3];
      track1->GetPxPyPz(mom1);
      track2->GetPxPyPz(mom2);
      //     Double_t cosTheta = (mom1[0]*mom2[0]+mom1[1]*mom2[1]+mom1[2]*mom2[2])/( TMath::Sqrt(mom1[0]*mom1[0]+mom1[1]*mom1[1]+mom1[2]*mom1[2])*TMath::Sqrt(mom2[0]*mom2[0]+mom2[1]*mom2[1]+mom2[2]*mom2[2]) );
      TVector3 momv1(mom1[0],mom1[1],mom1[2]);
      TVector3 momv2(mom2[0],mom2[1],mom2[2]);
      //Double_t theta = momv1.Angle(momv2);
      Double_t theta = momv1.Phi()-momv2.Phi();
      if(theta<-0.5*TMath::Pi()) theta+=2.*TMath::Pi();
    
      fDeltaPtCosmicCandidates->Fill(track1->Pt()-track2->Pt());
      Float_t deltaPhi = track1->Phi()-track2->Phi();
      if(deltaPhi<-0.5*TMath::Pi()) deltaPhi+=2.*TMath::Pi();
      fDeltaPhiSumEta->Fill(deltaPhi,track1->Eta()+track2->Eta());

      track1->GetImpactParameters(dcaR[0],dcaZ[0]);
      track2->GetImpactParameters(dcaR[1],dcaZ[1]);

      if(track2->Pt()<0.5) continue;
      Double_t rIsol = TMath::Sqrt( deltaPhi*deltaPhi+(track1->Eta()-track2->Eta())*(track1->Eta()-track2->Eta()) );
      fRisol->Fill(rIsol); //Fill R histogram
      if(track2->Pt()<fPtMin) continue;

	fTheta->Fill(theta);
	fThetaZoom->Fill(theta);
	fThetaPt1Pt2->Fill(theta,track1->Pt(),track2->Pt());
	fThetaPt1Pt2Signed->Fill(theta,track1->GetSign()*track1->Pt(),track2->GetSign()*track2->Pt());
	fDeltaPhiSumEtaPt1->Fill(deltaPhi,track1->Eta()+track2->Eta(),track1->Pt());
	fDeltaPhiSumEtaPt2->Fill(deltaPhi,track1->Eta()+track2->Eta(),track2->Pt());
	fThetaDCAZ1DCAZ2->Fill(theta,dcaZ[0],dcaZ[1]);
	fRisolTheta->Fill(rIsol,theta);
	if(TMath::Abs(TMath::Pi()-theta)<fMaxCosmicAngle) {
	  fDCAZCosmicCandidates->Fill(dcaZ[0],dcaZ[1]);
	  fDCARCosmicCandidates->Fill(dcaR[0],dcaR[1]);
	  fPtSignedCosmicCandidates->Fill(track1->GetSign()*track1->Pt());
	}

    } // track2 loop
    
  } // track1 loop 

  PostData(1,fHistListCosmics);
}      

//________________________________________________________________________
void AliPWG4CosmicCandidates::Terminate(Option_t *) 
{
  // Called once at the end of the query

  
}
