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
#include "AliFlowEventSimple.h"
#include "AliFlowTrackSimple.h"
#include "TTree.h"
#include "TParticle.h"
#include "AliFlowTrackSimpleCuts.h"


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
AliFlowEventSimple* FlowEventSimpleMaker::FillTracks(TTree* anInput, AliFlowTrackSimpleCuts* intCuts, AliFlowTrackSimpleCuts* diffCuts)
{
  //fills the event from a TTree of kinematic.root files
  
  // number of times to use the same particle (trick to introduce nonflow)
  Int_t iLoops = 1;
  
  //flags for particles passing int. and diff. flow cuts
  Bool_t bPassedIntFlowCuts  = kFALSE;
  Bool_t bPassedDiffFlowCuts = kFALSE;
  
  //track cut values
  Double_t dPtMaxInt  = intCuts->GetPtMax();
  Double_t dPtMinInt  = intCuts->GetPtMin();
  Double_t dEtaMaxInt = intCuts->GetEtaMax();
  Double_t dEtaMinInt = intCuts->GetEtaMin();
  Double_t dPhiMaxInt = intCuts->GetPhiMax();
  Double_t dPhiMinInt = intCuts->GetPhiMin();
  Int_t iPIDInt       = intCuts->GetPID();
  
  Double_t dPtMaxDiff  = diffCuts->GetPtMax();
  Double_t dPtMinDiff  = diffCuts->GetPtMin();
  Double_t dEtaMaxDiff = diffCuts->GetEtaMax();
  Double_t dEtaMinDiff = diffCuts->GetEtaMin();
  Double_t dPhiMaxDiff = diffCuts->GetPhiMax();
  Double_t dPhiMinDiff = diffCuts->GetPhiMin();
  Int_t iPIDDiff       = diffCuts->GetPID();
  
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
  Int_t iSelParticlesDiff = 0;
  Int_t iSelParticlesInt = 0;
  
  while (itrkN < iNumberOfInputTracks) {
    anInput->GetEntry(itrkN);   //get input particle
    //checking the cuts for int. and diff. flow
    if (pParticle->Pt() > dPtMinInt && pParticle->Pt() < dPtMaxInt &&
	pParticle->Eta() > dEtaMinInt && pParticle->Eta() < dEtaMaxInt &&
	pParticle->Phi() > dPhiMinInt && pParticle->Phi() < dPhiMaxInt &&
	TMath::Abs(pParticle->GetPdgCode()) == iPIDInt) { 
      bPassedIntFlowCuts = kTRUE; 
    } 
    
    if (pParticle->Pt() > dPtMinDiff && pParticle->Pt() < dPtMaxDiff &&
	pParticle->Eta() > dEtaMinDiff && pParticle->Eta() < dEtaMaxDiff &&
	pParticle->Phi() > dPhiMinDiff && pParticle->Phi() < dPhiMaxDiff &&
	TMath::Abs(pParticle->GetPdgCode()) == iPIDDiff){ 
      bPassedDiffFlowCuts = kTRUE; 
    }
    
    if (bPassedIntFlowCuts || bPassedDiffFlowCuts) {
      for(Int_t d=0;d<iLoops;d++) {
	AliFlowTrackSimple* pTrack = new AliFlowTrackSimple();
	pTrack->SetPt(pParticle->Pt());
	pTrack->SetEta(pParticle->Eta());
	pTrack->SetPhi(pParticle->Phi());
	
	//marking the particles used for int. flow:
	if(bPassedIntFlowCuts && iSelParticlesInt < iN*iLoops) {  
	  pTrack->SetForIntegratedFlow(kTRUE);
	  iSelParticlesInt++;
	}
	//marking the particles used for diff. flow:
	if(bPassedDiffFlowCuts) {
	  pTrack->SetForDifferentialFlow(kTRUE);
	  iSelParticlesDiff++;
	}
	//adding a particles which were used either for int. or diff. flow to the list
	pEvent->TrackCollection()->Add(pTrack);
	iGoodTracks++;
      }//end of for(Int_t d=0;d<iLoops;d++)
    }//end of if(bPassedIntFlowCuts || bPassedDiffFlowCuts) 
    itrkN++;  
    bPassedIntFlowCuts  = kFALSE;
    bPassedDiffFlowCuts = kFALSE;
  }//end of while (itrkN < iNumberOfInputTracks)
  
  pEvent->SetEventNSelTracksIntFlow(iSelParticlesInt);  
  pEvent->SetNumberOfTracks(iGoodTracks);//tracks used either for int. or for diff. flow

  cout<<" iGoodTracks = "<<iGoodTracks<<endl;
  cout<<" # of selected tracks for int. flow  = "<<iSelParticlesInt<<endl;
  cout<<" # of selected tracks for diff. flow = "<<iSelParticlesDiff<<endl;  

  delete pParticle;
  return pEvent;
}




 
