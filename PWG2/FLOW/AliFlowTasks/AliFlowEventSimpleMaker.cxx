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

/////////////////////////////////////////////////////////////////////////
// AliFlowEventSimpleMaker:
//
// Class to fill the AliFlowEventSimple
// with AliFlowTrackSimple objects
// Has fill methods for TTree, AliMCEvent, AliESDEvent and AliAODEvent
// author: N. van der Kolk (kolk@nikhef.nl)
/////////////////////////////////////////////////////////////////////////


#include "Riostream.h"
#include "AliFlowEventSimpleMaker.h"
#include "AliFlowEventSimple.h"
#include "AliFlowTrackSimple.h"
#include "TTree.h"
#include "TParticle.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliLog.h"
#include "AliCFManager.h"
#include "AliFlowTrackSimpleCuts.h"


ClassImp(AliFlowEventSimpleMaker)
//----------------------------------------------------------------------- 
AliFlowEventSimpleMaker::AliFlowEventSimpleMaker() :
  fMCReactionPlaneAngle(0.),
  fCount(0),
  fNoOfLoops(1),
  fEllipticFlowValue(0.),
  fMultiplicityOfEvent(1000000000)
{
  //constructor
}

//-----------------------------------------------------------------------   
AliFlowEventSimpleMaker::~AliFlowEventSimpleMaker()
{
  //destructor
}

//-----------------------------------------------------------------------   
AliFlowEventSimple* AliFlowEventSimpleMaker::FillTracks(TTree* anInput, AliFlowTrackSimpleCuts* rpCuts, AliFlowTrackSimpleCuts* poiCuts)
{
  //fills the event from a TTree of kinematic.root files
  
  // number of times to use the same particle (trick to introduce nonflow)
  //  Int_t iLoops = 1;
  
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
  
  //  Int_t fMultiplicityOfEvent = 576; //multiplicity for chi=1.5
  //  Int_t fMultiplicityOfEvent = 256; //multiplicity for chi=1
  //  Int_t fMultiplicityOfEvent = 164; //multiplicity for chi=0.8
  
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
      for(Int_t d=0;d<fNoOfLoops;d++) {
	AliFlowTrackSimple* pTrack = new AliFlowTrackSimple();
	pTrack->SetPt(pParticle->Pt());
	pTrack->SetEta(pParticle->Eta());
	pTrack->SetPhi(pParticle->Phi()-fEllipticFlowValue*TMath::Sin(2*(pParticle->Phi()-fMCReactionPlaneAngle)));
	
	//marking the particles used for int. flow:
	if(bPassedRPFlowCuts && iSelParticlesRP < fMultiplicityOfEvent) {  
	  pTrack->SetForRPSelection(kTRUE);
	  iSelParticlesRP++;
	}
	//marking the particles used for diff. flow:
	if(bPassedPOIFlowCuts && iGoodTracks%fNoOfLoops==0) {
	  pTrack->SetForPOISelection(kTRUE);
	  iSelParticlesPOI++;
	}
	//adding a particles which were used either for int. or diff. flow to the list
	pEvent->TrackCollection()->Add(pTrack);
	iGoodTracks++;
      }//end of for(Int_t d=0;d<iLoops;d++)
    }//end of if(bPassedIntFlowCuts || bPassedDiffFlowCuts) 
    itrkN++;  
    bPassedRPFlowCuts  = kFALSE;
    bPassedPOIFlowCuts = kFALSE;
  }//end of while (itrkN < iNumberOfInputTracks)
  
  pEvent->SetEventNSelTracksRP(iSelParticlesRP);  
  pEvent->SetNumberOfTracks(iGoodTracks);//tracks used either for RP or for POI selection
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

//-----------------------------------------------------------------------   
AliFlowEventSimple* AliFlowEventSimpleMaker::FillTracks(AliMCEvent* anInput)
{
  //Fills the event from the MC kinematic information
  
  Int_t iNumberOfInputTracks = anInput->GetNumberOfTracks() ;
 
  AliFlowEventSimple* pEvent = new AliFlowEventSimple(10);
    
  //Int_t iN = 256; //multiplicity for chi=1
  Int_t iN = iNumberOfInputTracks;
  Int_t iGoodTracks = 0;
  Int_t itrkN = 0;
  Int_t iSelParticlesDiff = 0;
  Int_t iSelParticlesInt = 0;

  //normal loop
  while (iGoodTracks < iN && itrkN < iNumberOfInputTracks) {
    AliMCParticle* pParticle = anInput->GetTrack(itrkN);   //get input particle
    //cut on tracks
    if (TMath::Abs(pParticle->Eta()) < 0.9)
      {
	if(
	   TMath::Abs(pParticle->Particle()->GetPdgCode()) == 211
	   //	      TMath::Abs(pParticle->Particle()->GetPdgCode()) == 211 ||
	   //	      TMath::Abs(pParticle->Particle()->GetPdgCode()) == 321 ||
	   //	      TMath::Abs(pParticle->Particle()->GetPdgCode()) == 2212
	   )
	  {
	    AliFlowTrackSimple* pTrack = new AliFlowTrackSimple();
	    pTrack->SetPt(pParticle->Pt() );
	    pTrack->SetEta(pParticle->Eta() );
	    pTrack->SetPhi(pParticle->Phi() );
	    pTrack->SetForRPSelection(kTRUE);
	    pTrack->SetForPOISelection(kTRUE);

	    if (pTrack->InRPSelection())
	      { iSelParticlesInt++; }
	    if (pTrack->InPOISelection())
	      { iSelParticlesDiff++; }
	    iGoodTracks++;
	    pEvent->TrackCollection()->Add(pTrack) ;  	     
	  }
	/*	  else if(
		  TMath::Abs(pParticle->Particle()->GetPdgCode()) == 211
		  )
	    {
	      AliFlowTrackSimple* pTrack = new AliFlowTrackSimple();
	      pTrack->SetPt(pParticle->Pt() );
	      pTrack->SetEta(pParticle->Eta() );
	      pTrack->SetPhi(pParticle->Phi() );
	      pTrack->SetForRPSelection(kFALSE);
	      pTrack->SetForPOISelection(kTRUE);

	      if (pTrack->InRPSelection())
		{ iSelParticlesInt++; }
	      if (pTrack->InPOISelection())
		{ iSelParticlesDiff++; }
	      iGoodTracks++;
	      pEvent->TrackCollection()->Add(pTrack);  	     
	    }
	  */
      }
      
    itrkN++; 
  }
  
  pEvent-> SetEventNSelTracksRP(iSelParticlesInt);  
  pEvent->SetNumberOfTracks(iGoodTracks);

  pEvent->SetMCReactionPlaneAngle(fMCReactionPlaneAngle);
  if (!fMCReactionPlaneAngle == 0) cout<<" MC Reaction Plane Angle = "<<  fMCReactionPlaneAngle << endl;
  else cout<<" MC Reaction Plane Angle = unknown "<< endl;

  cout<<" iGoodTracks = "<<iGoodTracks<<endl;
  cout<<" # of selected tracks for int. flow  = "<<iSelParticlesInt<<endl;
  cout<<" # of selected tracks for diff. flow = "<<iSelParticlesDiff<<endl;  
  cout << "# " << ++fCount << " events processed" << endl;

  return pEvent;

}


//-----------------------------------------------------------------------   
AliFlowEventSimple* AliFlowEventSimpleMaker::FillTracks(AliMCEvent* anInput, AliCFManager* intCFManager, AliCFManager* diffCFManager)
{
  //Fills the event from the MC kinematic information
  
  Int_t iNumberOfInputTracks = anInput->GetNumberOfTracks() ;

  if (iNumberOfInputTracks==-1) {
    cout<<"Skipping Event -- No MC information available for this event"<<endl;
    return 0;
  }

  AliFlowEventSimple* pEvent = new AliFlowEventSimple(10);
    
  Int_t iN = iNumberOfInputTracks; //maximum number of tracks in AliFlowEventSimple
  Int_t iGoodTracks = 0;           //number of good tracks
  Int_t itrkN = 0;                 //track counter
  Int_t iSelParticlesDiff = 0;     //number of tracks selected for Diff
  Int_t iSelParticlesInt = 0;      //number of tracks selected for Int

   
  //normal loop
  while (iGoodTracks < iN && itrkN < iNumberOfInputTracks) {
    //get input particle
    AliMCParticle* pParticle = anInput->GetTrack(itrkN);   
    //make new AliFlowTrackSimple
    AliFlowTrackSimple* pTrack = new AliFlowTrackSimple();
    pTrack->SetPt(pParticle->Pt() );
    pTrack->SetEta(pParticle->Eta() );
    pTrack->SetPhi(pParticle->Phi() );

    //check if pParticle passes the cuts
    if (intCFManager->CheckParticleCuts(AliCFManager::kPartGenCuts,pParticle)) {
      pTrack->SetForRPSelection(kTRUE);
      //cout<<"integrated selection. PID = "<<pParticle->Particle()->GetPdgCode()<<endl; 
    }
    if (diffCFManager->CheckParticleCuts(AliCFManager::kPartGenCuts,pParticle)) {
      pTrack->SetForPOISelection(kTRUE);
      //cout<<"differential selection. PID = "<<pParticle->Particle()->GetPdgCode()<<endl; 
    }

    //check if any bits are set
    TBits bFlowBits = pTrack->GetFlowBits();
    if (bFlowBits.CountBits() ==0) {
      delete pTrack; } //track will not be used anymore
    else {
      pEvent->TrackCollection()->Add(pTrack) ; 
      iGoodTracks++;

      if (pTrack->InRPSelection())
	{ iSelParticlesInt++; }
      if (pTrack->InPOISelection())
	{ iSelParticlesDiff++; }
    }
    
    itrkN++; 
  }
  
  pEvent-> SetEventNSelTracksRP(iSelParticlesInt);  
  pEvent-> SetNumberOfTracks(iGoodTracks);

  pEvent->SetMCReactionPlaneAngle(fMCReactionPlaneAngle);
  if (!fMCReactionPlaneAngle == 0) cout<<" MC Reaction Plane Angle = "<<  fMCReactionPlaneAngle << endl;
  else cout<<" MC Reaction Plane Angle = unknown "<< endl;

  cout<<" iGoodTracks = "<<iGoodTracks<<endl;
  cout<<" # of selected tracks for int. flow  = "<<iSelParticlesInt<<endl;
  cout<<" # of selected tracks for diff. flow = "<<iSelParticlesDiff<<endl;  
  cout << "# " << ++fCount << " events processed" << endl;

  return pEvent;

}


//-----------------------------------------------------------------------   
AliFlowEventSimple* AliFlowEventSimpleMaker::FillTracks(AliESDEvent* anInput)
{
  //Fills the event from the ESD
  
  Int_t iNumberOfInputTracks = anInput->GetNumberOfTracks() ;
  
  AliFlowEventSimple* pEvent = new AliFlowEventSimple(10);
    
  //Int_t iN = 256; //multiplicity for chi=1
  Int_t iN = iNumberOfInputTracks;
  Int_t iGoodTracks = 0;
  Int_t itrkN = 0;
  Int_t iSelParticlesDiff = 0;
  Int_t iSelParticlesInt = 0;

 
  
  //normal loop
  while (iGoodTracks < iN && itrkN < iNumberOfInputTracks) {
    AliESDtrack* pParticle = anInput->GetTrack(itrkN);   //get input particle
    //cut on tracks
    if (TMath::Abs(pParticle->Eta()) < 0.9)
    
   
    
    {
	AliFlowTrackSimple* pTrack = new AliFlowTrackSimple();
	pTrack->SetPt(pParticle->Pt() );
	pTrack->SetEta(pParticle->Eta() );
	pTrack->SetPhi(pParticle->Phi() );
	pTrack->SetForRPSelection(kTRUE);
	pTrack->SetForPOISelection(kTRUE);

	if (pTrack->InRPSelection())
	  { iSelParticlesInt++; }
	if (pTrack->InPOISelection())
	  { iSelParticlesDiff++; }
	iGoodTracks++;
	pEvent->TrackCollection()->Add(pTrack) ;  	     
      }
      
    itrkN++; 
  }
  
  pEvent-> SetEventNSelTracksRP(iSelParticlesInt);  
  pEvent->SetNumberOfTracks(iGoodTracks);

  pEvent->SetMCReactionPlaneAngle(fMCReactionPlaneAngle);
  if (!fMCReactionPlaneAngle == 0) cout<<" MC Reaction Plane Angle = "<<  fMCReactionPlaneAngle << endl;
  else cout<<" MC Reaction Plane Angle = unknown "<< endl;

  cout<<" iGoodTracks = "<<iGoodTracks<<endl;
  cout<<" # of selected tracks for int. flow  = "<<iSelParticlesInt<<endl;
  cout<<" # of selected tracks for diff. flow = "<<iSelParticlesDiff<<endl;  
  cout << "# " << ++fCount << " events processed" << endl;

  return pEvent;
}


//-----------------------------------------------------------------------   
AliFlowEventSimple* AliFlowEventSimpleMaker::FillTracks(AliESDEvent* anInput, AliCFManager* intCFManager, AliCFManager* diffCFManager)
{
  //Fills the event from the ESD
  
  Int_t iNumberOfInputTracks = anInput->GetNumberOfTracks() ;
  
  AliFlowEventSimple* pEvent = new AliFlowEventSimple(10);
    
  Int_t iN = iNumberOfInputTracks; //maximum number of tracks in AliFlowEventSimple
  Int_t iGoodTracks = 0;           //number of good tracks
  Int_t itrkN = 0;                 //track counter
  Int_t iSelParticlesDiff = 0;     //number of tracks selected for Diff
  Int_t iSelParticlesInt = 0;      //number of tracks selected for Int

  //normal loop
  while (iGoodTracks < iN && itrkN < iNumberOfInputTracks) {
    AliESDtrack* pParticle = anInput->GetTrack(itrkN);   //get input particle
    //make new AliFLowTrackSimple
    AliFlowTrackSimple* pTrack = new AliFlowTrackSimple();
    pTrack->SetPt(pParticle->Pt() );
    pTrack->SetEta(pParticle->Eta() );
    pTrack->SetPhi(pParticle->Phi() );
    //check if pParticle passes the cuts

    if (intCFManager->CheckParticleCuts(AliCFManager::kPartRecCuts,pParticle) && 
	intCFManager->CheckParticleCuts(AliCFManager::kPartSelCuts,pParticle)) {
      pTrack->SetForRPSelection(kTRUE);
    }
    if (diffCFManager->CheckParticleCuts(AliCFManager::kPartRecCuts,pParticle) && 
	diffCFManager->CheckParticleCuts(AliCFManager::kPartSelCuts,pParticle)) {
      pTrack->SetForPOISelection(kTRUE);}

    //check if any bits are set
    TBits bFlowBits = pTrack->GetFlowBits();
    if (bFlowBits.CountBits() ==0) {
      delete pTrack; } //track will not be used anymore
    else {
      pEvent->TrackCollection()->Add(pTrack) ;  
      iGoodTracks++;

      if (pTrack->InRPSelection())
	{ iSelParticlesInt++; }
      if (pTrack->InPOISelection())
	{ iSelParticlesDiff++; }
      
    }
    
    itrkN++; 
  }
  
  pEvent-> SetEventNSelTracksRP(iSelParticlesInt);  
  pEvent->SetNumberOfTracks(iGoodTracks);

  pEvent->SetMCReactionPlaneAngle(fMCReactionPlaneAngle);
  if (!fMCReactionPlaneAngle == 0) cout<<" MC Reaction Plane Angle = "<<  fMCReactionPlaneAngle << endl;
  else cout<<" MC Reaction Plane Angle = unknown "<< endl;

  cout<<" iGoodTracks = "<<iGoodTracks<<endl;
  cout<<" # of selected tracks for int. flow  = "<<iSelParticlesInt<<endl;
  cout<<" # of selected tracks for diff. flow = "<<iSelParticlesDiff<<endl;  
  cout << "# " << ++fCount << " events processed" << endl;

  return pEvent;
}

//-----------------------------------------------------------------------   
AliFlowEventSimple* AliFlowEventSimpleMaker::FillTracks(AliAODEvent* anInput)
{
  //Fills the event from the AOD
  
  Int_t iNumberOfInputTracks = anInput->GetNumberOfTracks() ;
  
  AliFlowEventSimple* pEvent = new AliFlowEventSimple(10);
    
  //Int_t iN = 256; //multiplicity for chi=1
  Int_t iN = iNumberOfInputTracks;
  Int_t iGoodTracks = 0;
  Int_t itrkN = 0;
  Int_t iSelParticlesDiff = 0;
  Int_t iSelParticlesInt = 0;
  
  //normal loop
  while (iGoodTracks < iN && itrkN < iNumberOfInputTracks) {
    AliAODTrack* pParticle = anInput->GetTrack(itrkN);   //get input particle
    //cut on tracks
    if (TMath::Abs(pParticle->Eta()) < 0.9)
      {
	AliFlowTrackSimple* pTrack = new AliFlowTrackSimple();
	pTrack->SetPt(pParticle->Pt() );
	pTrack->SetEta(pParticle->Eta() );
	pTrack->SetPhi(pParticle->Phi() );
	pTrack->SetForRPSelection(kTRUE);
	pTrack->SetForPOISelection(kTRUE);

	if (pTrack->InRPSelection())
	  { iSelParticlesInt++; }
	if (pTrack->InPOISelection())
	  { iSelParticlesDiff++; }
	iGoodTracks++;
	pEvent->TrackCollection()->Add(pTrack) ;  	     
      }
      
    itrkN++; 
  }
  
  pEvent-> SetEventNSelTracksRP(iSelParticlesInt);  
  pEvent->SetNumberOfTracks(iGoodTracks);

  pEvent->SetMCReactionPlaneAngle(fMCReactionPlaneAngle);
  if (!fMCReactionPlaneAngle == 0) cout<<" MC Reaction Plane Angle = "<<  fMCReactionPlaneAngle << endl;
  else cout<<" MC Reaction Plane Angle = unknown "<< endl;

  cout<<" iGoodTracks = "<<iGoodTracks<<endl;
  cout<<" # of selected tracks for int. flow  = "<<iSelParticlesInt<<endl;
  cout<<" # of selected tracks for diff. flow = "<<iSelParticlesDiff<<endl;  
  cout << "# " << ++fCount << " events processed" << endl;

  return pEvent;
}


//-----------------------------------------------------------------------   
AliFlowEventSimple* AliFlowEventSimpleMaker::FillTracks(AliAODEvent* anInput,  AliCFManager* intCFManager, AliCFManager* diffCFManager)
{
  //Fills the event from the AOD
  
  Int_t iNumberOfInputTracks = anInput->GetNumberOfTracks() ;
  
  AliFlowEventSimple* pEvent = new AliFlowEventSimple(10);
  
  Int_t iN = iNumberOfInputTracks; //maximum number of tracks in AliFlowEventSimple
  Int_t iGoodTracks = 0;           //number of good tracks
  Int_t itrkN = 0;                 //track counter
  Int_t iSelParticlesDiff = 0;     //number of tracks selected for Diff
  Int_t iSelParticlesInt = 0;      //number of tracks selected for Int

  //normal loop
  while (iGoodTracks < iN && itrkN < iNumberOfInputTracks) {
    AliAODTrack* pParticle = anInput->GetTrack(itrkN);   //get input particle
    //make new AliFlowTrackSimple
    AliFlowTrackSimple* pTrack = new AliFlowTrackSimple();
    pTrack->SetPt(pParticle->Pt() );
    pTrack->SetEta(pParticle->Eta() );
    pTrack->SetPhi(pParticle->Phi() );

    //check if pParticle passes the cuts
    if (intCFManager->CheckParticleCuts(AliCFManager::kPartRecCuts,pParticle) &&
	intCFManager->CheckParticleCuts(AliCFManager::kPartSelCuts,pParticle)) {          
      pTrack->SetForRPSelection(kTRUE); }
    if (diffCFManager->CheckParticleCuts(AliCFManager::kPartRecCuts,pParticle) &&
	diffCFManager->CheckParticleCuts(AliCFManager::kPartSelCuts,pParticle)) {
      pTrack->SetForPOISelection(kTRUE);}	

    //check if any bits are set
    TBits bFlowBits = pTrack->GetFlowBits();
    if (bFlowBits.CountBits() ==0) {
      delete pTrack; } //track will not be used anymore
    else {
      pEvent->TrackCollection()->Add(pTrack) ; 
      iGoodTracks++;

      if (pTrack->InRPSelection())
	{ iSelParticlesInt++; }
      if (pTrack->InPOISelection())
	{ iSelParticlesDiff++; }
      	     
    }
       
    itrkN++; 
  }
  
  pEvent-> SetEventNSelTracksRP(iSelParticlesInt);  
  pEvent->SetNumberOfTracks(iGoodTracks);

  pEvent->SetMCReactionPlaneAngle(fMCReactionPlaneAngle);
  if (!fMCReactionPlaneAngle == 0) cout<<" MC Reaction Plane Angle = "<<  fMCReactionPlaneAngle << endl;
  else cout<<" MC Reaction Plane Angle = unknown "<< endl;

  cout<<" iGoodTracks = "<<iGoodTracks<<endl;
  cout<<" # of selected tracks for int. flow  = "<<iSelParticlesInt<<endl;
  cout<<" # of selected tracks for diff. flow = "<<iSelParticlesDiff<<endl;  
  cout << "# " << ++fCount << " events processed" << endl;

  return pEvent;
}

//-----------------------------------------------------------------------   
AliFlowEventSimple*  AliFlowEventSimpleMaker::FillTracks(AliESDEvent* anInput, AliMCEvent* anInputMc, Int_t anOption)
{
  //fills the event with tracks from the ESD and kinematics from the MC info via the track label

  if (!(anOption ==0 || anOption ==1)) {
    cout<<"WRONG OPTION IN AliFlowEventSimpleMaker::FillTracks(AliESDEvent* anInput, AliMCEvent* anInputMc, Int_t anOption)"<<endl;
    exit(1);
  }

  Int_t iNumberOfInputTracks = anInput->GetNumberOfTracks() ;
  
  AliFlowEventSimple* pEvent = new AliFlowEventSimple(10);
    
  //Int_t iN = 256; //multiplicity for chi=1
  Int_t iN = iNumberOfInputTracks;
  Int_t iGoodTracks = 0;
  Int_t itrkN = 0;
  Int_t iSelParticlesDiff = 0;
  Int_t iSelParticlesInt = 0;

  //normal loop
  while (iGoodTracks < iN && itrkN < iNumberOfInputTracks) {
    AliESDtrack* pParticle = anInput->GetTrack(itrkN);   //get input particle
    //get Label
    Int_t iLabel = pParticle->GetLabel();
    //match to mc particle
    AliMCParticle* pMcParticle = anInputMc->GetTrack(TMath::Abs(iLabel));
    
    //check
    if (TMath::Abs(pParticle->GetLabel())!=pMcParticle->Label()) cout<<"pParticle->GetLabel()!=pMcParticle->Label() "<<pParticle->GetLabel()<<"  "<<pMcParticle->Label()<<endl;
    
    //cut on tracks
    if (TMath::Abs(pParticle->Eta()) < 0.2)
      {
	if(
	   TMath::Abs(pMcParticle->Particle()->GetPdgCode()) == 211 //pions
	   //	      TMath::Abs(pMcParticle->Particle()->GetPdgCode()) == 211 ||
	   //	      TMath::Abs(pMcParticle->Particle()->GetPdgCode()) == 321 ||
	   //	      TMath::Abs(pMcParticle->Particle()->GetPdgCode()) == 2212
	   )
	  {
	    AliFlowTrackSimple* pTrack = new AliFlowTrackSimple();
	    if(anOption == 0) { //take the PID from the MC & the kinematics from the ESD
	      pTrack->SetPt(pParticle->Pt() );
	      pTrack->SetEta(pParticle->Eta() );
	      pTrack->SetPhi(pParticle->Phi() );
	      pTrack->SetForRPSelection(kTRUE);
	      pTrack->SetForPOISelection(kTRUE);
	    }
	    else if (anOption == 1) { //take the PID and kinematics from the MC
	      pTrack->SetPt(pMcParticle->Pt() );
	      pTrack->SetEta(pMcParticle->Eta() );
	      pTrack->SetPhi(pMcParticle->Phi() );
	      pTrack->SetForRPSelection(kTRUE);
	      pTrack->SetForPOISelection(kTRUE);
	    }
	    else { cout<<"Not a valid option"<<endl; }
	    if (pTrack->InRPSelection())
	      { iSelParticlesInt++; }
	    if (pTrack->InPOISelection())
	      { iSelParticlesDiff++; }
	    iGoodTracks++;
	    pEvent->TrackCollection()->Add(pTrack) ;  	     
	  }
      }
    itrkN++; 
  }
  
  pEvent-> SetEventNSelTracksRP(iSelParticlesInt);  
  pEvent->SetNumberOfTracks(iGoodTracks);

  pEvent->SetMCReactionPlaneAngle(fMCReactionPlaneAngle);
  if (!fMCReactionPlaneAngle == 0) cout<<" MC Reaction Plane Angle = "<<  fMCReactionPlaneAngle << endl;
  else cout<<" MC Reaction Plane Angle = unknown "<< endl;

  cout<<" iGoodTracks = "<<iGoodTracks<<endl;
  cout<<" # of selected tracks for int. flow  = "<<iSelParticlesInt<<endl;
  cout<<" # of selected tracks for diff. flow = "<<iSelParticlesDiff<<endl;  
  cout << "# " << ++fCount << " events processed" << endl;

  return pEvent;
}

//-----------------------------------------------------------------------   
AliFlowEventSimple*  AliFlowEventSimpleMaker::FillTracks(AliESDEvent* anInput, AliMCEvent* anInputMc, AliCFManager* intCFManager, AliCFManager* diffCFManager, Int_t anOption)
{
  //fills the event with tracks from the ESD and kinematics from the MC info via the track label

  
  if (!(anOption ==0 || anOption ==1)) {
    cout<<"WRONG OPTION IN AliFlowEventSimpleMaker::FillTracks(AliESDEvent* anInput, AliMCEvent* anInputMc, Int_t anOption)"<<endl;
    exit(1);
  }

  Int_t iNumberOfInputTracks = anInput->GetNumberOfTracks() ;

  Int_t iNumberOfInputTracksMC = anInputMc->GetNumberOfTracks() ;
  if (iNumberOfInputTracksMC==-1) {
    cout<<"Skipping Event -- No MC information available for this event"<<endl;
    return 0;
  }

  AliFlowEventSimple* pEvent = new AliFlowEventSimple(10);
    
  Int_t iN = iNumberOfInputTracks; //maximum number of tracks in AliFlowEventSimple
  Int_t iGoodTracks = 0;           //number of good tracks
  Int_t itrkN = 0;                 //track counter
  Int_t iSelParticlesDiff = 0;     //number of tracks selected for Diff
  Int_t iSelParticlesInt = 0;      //number of tracks selected for Int

 
  //normal loop
  while (iGoodTracks < iN && itrkN < iNumberOfInputTracks) {
    AliESDtrack* pParticle = anInput->GetTrack(itrkN);   //get input particle
    //get Label
    Int_t iLabel = pParticle->GetLabel();
    //match to mc particle
    AliMCParticle* pMcParticle = anInputMc->GetTrack(TMath::Abs(iLabel));
    
    //check
    if (TMath::Abs(pParticle->GetLabel())!=pMcParticle->Label()) cout<<"pParticle->GetLabel()!=pMcParticle->Label() "<<pParticle->GetLabel()<<"  "<<pMcParticle->Label()<<endl;
    
    //make new AliFlowTrackSimple
    AliFlowTrackSimple* pTrack = new AliFlowTrackSimple();
    if(anOption == 0) { //take the PID from the MC & the kinematics from the ESD
      pTrack->SetPt(pParticle->Pt() );
      pTrack->SetEta(pParticle->Eta() );
      pTrack->SetPhi(pParticle->Phi() );
    }
    else if (anOption == 1) { //take the PID and kinematics from the MC
      pTrack->SetPt(pMcParticle->Pt() );
      pTrack->SetEta(pMcParticle->Eta() );
      pTrack->SetPhi(pMcParticle->Phi() );
    }
    else { cout<<"Not a valid option"<<endl; }

    //check if pParticle passes the cuts
    if(anOption == 0) { 
      //cout<<"take the PID from the MC & the kinematics from the ESD"<<endl;
      if (intCFManager->CheckParticleCuts(AliCFManager::kPartGenCuts,pMcParticle,"mcGenCuts1") && 
	  intCFManager->CheckParticleCuts(AliCFManager::kPartRecCuts,pParticle)) {  
	pTrack->SetForRPSelection(kTRUE); }
      if (diffCFManager->CheckParticleCuts(AliCFManager::kPartGenCuts,pMcParticle,"mcGenCuts2") &&
	  diffCFManager->CheckParticleCuts(AliCFManager::kPartRecCuts,pParticle)) {  
	pTrack->SetForPOISelection(kTRUE);}
    }
    else if (anOption == 1) { 
      //cout<<"take the PID and kinematics from the MC"<<endl;
      if (intCFManager->CheckParticleCuts(AliCFManager::kPartGenCuts,pMcParticle)) {  
	pTrack->SetForRPSelection(kTRUE); }
      if (diffCFManager->CheckParticleCuts(AliCFManager::kPartGenCuts,pMcParticle)) {  
	pTrack->SetForPOISelection(kTRUE);}
    }
    else { cout<<"Not a valid option"<<endl; }
      
    //check if any bits are set
    TBits bFlowBits = pTrack->GetFlowBits();
    if (bFlowBits.CountBits() ==0) {
      delete pTrack; } //track will not be used anymore
    else {
      pEvent->TrackCollection()->Add(pTrack) ; 
      iGoodTracks++;  
    
      if (pTrack->InRPSelection())
	{ iSelParticlesInt++; }
      if (pTrack->InPOISelection())
	{ iSelParticlesDiff++; }
      	     
    }
    
    itrkN++; 
  }
  
  pEvent-> SetEventNSelTracksRP(iSelParticlesInt);  
  pEvent->SetNumberOfTracks(iGoodTracks);

  cout << " Number of MC input tracks = " << iNumberOfInputTracksMC << endl;

  pEvent->SetMCReactionPlaneAngle(fMCReactionPlaneAngle);
  if (!fMCReactionPlaneAngle == 0) cout<<" MC Reaction Plane Angle = "<<  fMCReactionPlaneAngle << endl;
  else cout<<" MC Reaction Plane Angle = unknown "<< endl;

  cout<<" iGoodTracks = "<<iGoodTracks<<endl;
  cout<<" # of selected tracks for int. flow  = "<<iSelParticlesInt<<endl;
  cout<<" # of selected tracks for diff. flow = "<<iSelParticlesDiff<<endl;  
  cout << "# " << ++fCount << " events processed" << endl;

  return pEvent;
}

