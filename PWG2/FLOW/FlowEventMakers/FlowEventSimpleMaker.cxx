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
FlowEventSimpleMaker::FlowEventSimpleMaker()
{
  //constructor
}

//-----------------------------------------------------------------------   
FlowEventSimpleMaker::~FlowEventSimpleMaker()
{
  //destructor
}

//-----------------------------------------------------------------------   
AliFlowEventSimple* FlowEventSimpleMaker::FillTracks(TTree* anInput, AliFlowTrackSimpleCuts* RPCuts, AliFlowTrackSimpleCuts* POICuts)
{
  //fills the event from a TTree of kinematic.root files
  
  // number of times to use the same particle (trick to introduce nonflow)
  Int_t iLoops = 1;
  
  //flags for particles passing RP and POI cuts
  Bool_t bPassedRPCuts  = kFALSE;
  Bool_t bPassedPOICuts = kFALSE;
  
  //track cut values
  Double_t dPtMaxRP  = RPCuts->GetPtMax();
  Double_t dPtMinRP  = RPCuts->GetPtMin();
  Double_t dEtaMaxRP = RPCuts->GetEtaMax();
  Double_t dEtaMinRP = RPCuts->GetEtaMin();
  Double_t dPhiMaxRP = RPCuts->GetPhiMax();
  Double_t dPhiMinRP = RPCuts->GetPhiMin();
  Int_t iPIDRP       = RPCuts->GetPID();
  
  Double_t dPtMaxPOI  = POICuts->GetPtMax();
  Double_t dPtMinPOI  = POICuts->GetPtMin();
  Double_t dEtaMaxPOI = POICuts->GetEtaMax();
  Double_t dEtaMinPOI = POICuts->GetEtaMin();
  Double_t dPhiMaxPOI = POICuts->GetPhiMax();
  Double_t dPhiMinPOI = POICuts->GetPhiMin();
  Int_t iPIDPOI       = POICuts->GetPID();
  
  Int_t iNumberOfInputTracks = anInput->GetEntries() ;
  //cerr<<"iNumberOfInputTracks = "<<iNumberOfInputTracks<<endl;
  TParticle* pParticle = new TParticle();
  anInput->SetBranchAddress("Particles",&pParticle);  
  //  AliFlowEventSimple* pEvent = new AliFlowEventSimple(iNumberOfInputTracks);
  AliFlowEventSimple* pEvent = new AliFlowEventSimple(10);
  //cerr<<pEvent<<" pEvent "<<endl;
  
  Int_t iN = iNumberOfInputTracks; // additional variable to artificially fix the number of tracks
  //  Int_t iN = 576; //multiplicity for chi=1.5
  //  Int_t iN = 256; //multiplicity for chi=1
  //  Int_t iN = 164; //multiplicity for chi=0.8
  
  Int_t iGoodTracks = 0;
  Int_t itrkN = 0;
  Int_t iSelParticlesPOI = 0;
  Int_t iSelParticlesRP  = 0;
  
  while (itrkN < iNumberOfInputTracks) {
    anInput->GetEntry(itrkN);   //get input particle
    //checking the cuts for int. and diff. flow
    if (pParticle->Pt() > dPtMinRP && pParticle->Pt() < dPtMaxRP &&
	pParticle->Eta() > dEtaMinRP && pParticle->Eta() < dEtaMaxRP &&
	pParticle->Phi() > dPhiMinRP && pParticle->Phi() < dPhiMaxRP &&
	TMath::Abs(pParticle->GetPdgCode()) == iPIDRP) { 
      bPassedRPCuts = kTRUE; 
    } 
    
    if (pParticle->Pt() > dPtMinPOI && pParticle->Pt() < dPtMaxPOI &&
	pParticle->Eta() > dEtaMinPOI && pParticle->Eta() < dEtaMaxPOI &&
	pParticle->Phi() > dPhiMinPOI && pParticle->Phi() < dPhiMaxPOI &&
	TMath::Abs(pParticle->GetPdgCode()) == iPIDPOI){ 
      bPassedPOICuts = kTRUE; 
    }
    
    if (bPassedRPCuts || bPassedPOICuts) {
      for(Int_t d=0;d<iLoops;d++) {
	AliFlowTrackSimple* pTrack = new AliFlowTrackSimple();
	pTrack->SetPt(pParticle->Pt());
	pTrack->SetEta(pParticle->Eta());
	pTrack->SetPhi(pParticle->Phi());
	
	//marking the particles used for int. flow:
	if(bPassedRPCuts && iSelParticlesRP < iN*iLoops) {  
	  //pTrack->SetForIntegratedFlow(kTRUE);
	  pTrack->SetForPRSelection(kTRUE);
	  iSelParticlesRP++;
	}
	//marking the particles used for diff. flow:
	if(bPassedPOICuts) {
	  pTrack->SetForPOISelection(kTRUE);
	  iSelParticlesPOI++;
	}
	//adding a particles which were used either for int. or diff. flow to the list
	pEvent->TrackCollection()->Add(pTrack);
	iGoodTracks++;
      }//end of for(Int_t d=0;d<iLoops;d++)
    }//end of if(bPassedIntFlowCuts || bPassedDiffFlowCuts) 
    itrkN++;  
    bPassedRPCuts  = kFALSE;
    bPassedPOICuts = kFALSE;
  }//end of while (itrkN < iNumberOfInputTracks)
  
  pEvent->SetEventNSelTracksRP(iSelParticlesRP);  
  pEvent->SetNumberOfTracks(iGoodTracks);//tracks used either for int. or for diff. flow

  cout<<" iGoodTracks = "<<iGoodTracks<<endl;
  cout<<" # of selected tracks for RP  = "<<iSelParticlesRP<<endl;
  cout<<" # of selected tracks for POI = "<<iSelParticlesPOI<<endl;  

  delete pParticle;
  return pEvent;
}




 
