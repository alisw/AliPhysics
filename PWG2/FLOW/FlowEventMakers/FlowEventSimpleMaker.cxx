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
/* $Id */

#include "Riostream.h"
#include "FlowEventSimpleMaker.h"
#include "AliFlowCommon/AliFlowEventSimple.h"
#include "AliFlowCommon/AliFlowTrackSimple.h"
#include "TTree.h"
#include "TParticle.h"
#include "AliFlowCommon/AliFlowTrackSimpleCuts.h"


// FlowEventSimpleMaker:
// Class to fill the AliFlowEventSimple
// with AliFlowTrackSimple objects
// ouside the AliRoot Framework
// Has fill methods for TTree, 

ClassImp(FlowEventSimpleMaker)
//----------------------------------------------------------------------- 
FlowEventSimpleMaker::FlowEventSimpleMaker():
  fMCReactionPlaneAngle(0.),
  fCount(0)
{
  //constructor
}

//-----------------------------------------------------------------------   
FlowEventSimpleMaker::~FlowEventSimpleMaker()
{
  //destructor
}

//-----------------------------------------------------------------------   
AliFlowEventSimple* FlowEventSimpleMaker::FillTracks(TTree* anInput, AliFlowTrackSimpleCuts* rpCuts, AliFlowTrackSimpleCuts* poiCuts)
{
  //fills the event from a TTree of kinematic.root files
  
  // number of times to use the same particle (trick to introduce nonflow)
  Int_t iLoops = 1;
  
  //flags for particles passing int. and diff. flow cuts
  Bool_t bPassedRPFlowCuts  = kFALSE;
  Bool_t bPassedPOIFlowCuts = kFALSE;
  
  //track cut values
  Double_t dPtMaxRP  = rpCuts->GetPtMax();
  Double_t dPtMinRP  = rpCuts->GetPtMin();
  Double_t dEtaMaxRP = rpCuts->GetEtaMax();
  Double_t dEtaMinRP = rpCuts->GetEtaMin();
  Double_t dPhiMaxRP = rpCuts->GetPhiMax();
  Double_t dPhiMinRP = rpCuts->GetPhiMin();
  Int_t iPIDRP       = rpCuts->GetPID();
  
  Double_t dPtMaxPOI  = poiCuts->GetPtMax();
  Double_t dPtMinPOI  = poiCuts->GetPtMin();
  Double_t dEtaMaxPOI = poiCuts->GetEtaMax();
  Double_t dEtaMinPOI = poiCuts->GetEtaMin();
  Double_t dPhiMaxPOI = poiCuts->GetPhiMax();
  Double_t dPhiMinPOI = poiCuts->GetPhiMin();
  Int_t iPIDPOI       = poiCuts->GetPID();
  
  Int_t iNumberOfInputTracks = anInput->GetEntries() ;

  TParticle* pParticle = new TParticle();
  anInput->SetBranchAddress("Particles",&pParticle);  
  //  AliFlowEventSimple* pEvent = new AliFlowEventSimple(iNumberOfInputTracks);
  AliFlowEventSimple* pEvent = new AliFlowEventSimple(10);
  
  Int_t iN = iNumberOfInputTracks; // additional variable to artificially fix the number of tracks
  //  Int_t iN = 576; //multiplicity for chi=1.5
  //  Int_t iN = 256; //multiplicity for chi=1
  //  Int_t iN = 164; //multiplicity for chi=0.8
  
  Int_t iGoodTracks = 0;
  Int_t itrkN = 0;
  Int_t iSelParticlesRP = 0;
  Int_t iSelParticlesPOI = 0;
  
  
  while (itrkN < iNumberOfInputTracks) {
    anInput->GetEntry(itrkN);   //get input particle
    if (pParticle->IsPrimary()) {
      //checking the cuts for int. and diff. flow
      if (pParticle->Pt() > dPtMinRP && pParticle->Pt() < dPtMaxRP &&
	  pParticle->Eta() > dEtaMinRP && pParticle->Eta() < dEtaMaxRP &&
	  pParticle->Phi() > dPhiMinRP && pParticle->Phi() < dPhiMaxRP &&
	  TMath::Abs(pParticle->GetPdgCode()) == iPIDRP) { 
	bPassedRPFlowCuts = kTRUE;
      } 
    
      if (pParticle->Pt() > dPtMinPOI && pParticle->Pt() < dPtMaxPOI &&
	  pParticle->Eta() > dEtaMinPOI && pParticle->Eta() < dEtaMaxPOI &&
	  pParticle->Phi() > dPhiMinPOI && pParticle->Phi() < dPhiMaxPOI &&
	  TMath::Abs(pParticle->GetPdgCode()) == iPIDPOI){ 
	bPassedPOIFlowCuts = kTRUE; 
      }
    }
    if (bPassedRPFlowCuts || bPassedPOIFlowCuts) {
      for(Int_t d=0;d<iLoops;d++) {
	AliFlowTrackSimple* pTrack = new AliFlowTrackSimple();
	pTrack->SetPt(pParticle->Pt());
	pTrack->SetEta(pParticle->Eta());
	pTrack->SetPhi(pParticle->Phi());
	
	//marking the particles used for int. flow:
	if(bPassedRPFlowCuts && iSelParticlesRP < iN*iLoops) {  
	  pTrack->SetForRPSelection(kTRUE);
	  iSelParticlesRP++;
	}
	//marking the particles used for diff. flow:
	if(bPassedPOIFlowCuts) {
	  pTrack->SetForPOISelection(kTRUE);
	  iSelParticlesPOI++;
	}
	//adding a particles which were used either for int. or diff. flow to the list
	pEvent->AddTrack(pTrack);
	iGoodTracks++;
      }//end of for(Int_t d=0;d<iLoops;d++)
    }//end of if(bPassedIntFlowCuts || bPassedDiffFlowCuts) 
    itrkN++;  
    bPassedRPFlowCuts  = kFALSE;
    bPassedPOIFlowCuts = kFALSE;
  }//end of while (itrkN < iNumberOfInputTracks)
  
  pEvent->SetEventNSelTracksRP(iSelParticlesRP);  
  pEvent->SetMCReactionPlaneAngle(fMCReactionPlaneAngle);

  if (!fMCReactionPlaneAngle == 0) cout<<" MC Reaction Plane Angle = "<<  fMCReactionPlaneAngle << endl;
  else cout<<" MC Reaction Plane Angle = unknown "<< endl;

  cout<<" iGoodTracks = "<< iGoodTracks << endl;
  cout<<" # of RP selected tracks = "<<iSelParticlesRP<<endl;
  cout<<" # of POI selected tracks = "<<iSelParticlesPOI<<endl;  
  cout << "# " << ++fCount << " events processed" << endl;

  delete pParticle;
  return pEvent;
}




 
