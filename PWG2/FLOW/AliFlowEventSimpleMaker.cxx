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

// AliFlowEventSimpleMaker:
// Class to fill the AliFlowEventSimple
// with AliFlowTrackSimple objects
// Has fill methods for TTree, AliMCEvent, AliESDEvent and AliAODEvent
// author: N. van der Kolk (kolk@nikhef.nl)



ClassImp(AliFlowEventSimpleMaker)
//----------------------------------------------------------------------- 
AliFlowEventSimpleMaker::AliFlowEventSimpleMaker():
  fEvent(0),
  fTrack(0),
  fParticle(0)
{

  //constructor

}


//-----------------------------------------------------------------------   
AliFlowEventSimpleMaker::~AliFlowEventSimpleMaker()
{
  //desstructor
}


//-----------------------------------------------------------------------   
AliFlowEventSimple* AliFlowEventSimpleMaker::FillTracks(TTree* fInput)
{
  //fills the event from a TTree of kinematic.root files
  Bool_t  fDoubleLoop = kFALSE; 

  Int_t fNumberOfInputTracks = fInput->GetEntries() ;
  //cerr<<"fNumberOfInputTracks = "<<fNumberOfInputTracks<<endl;
  fParticle = new TParticle();
  fInput->SetBranchAddress("Particles",&fParticle);  
  //  fEvent = new AliFlowEventSimple(fNumberOfInputTracks);
  fEvent = new AliFlowEventSimple(10);
  //cerr<<fEvent<<" fEvent "<<endl;
  
  Int_t fN = fNumberOfInputTracks;
  //  Int_t fN = 576; //multiplicity for chi=1.5
  //  Int_t fN = 256; //multiplicity for chi=1
  //  Int_t fN = 164; //multiplicity for chi=0.8
  Int_t fGoodTracks = 0;
  Int_t ftrkN = 0;
  Int_t fSelParticlesDiff = 0;
  Int_t fSelParticlesInt = 0;
  
  if (fDoubleLoop)
    {                   //double loop
      while (fGoodTracks < fN*2 && ftrkN < fNumberOfInputTracks) 
	{
	  fInput->GetEntry(ftrkN);   //get input particle
	  //cut on tracks
	  if(TMath::Abs(fParticle->Eta()) < 0.9) 
	    {
	      //	      Int_t fLoop = floor(2.*fParticle->Pt())+2;
	      //	      for(Int_t d=0;d<fLoop;d++) 
	      for(Int_t d=0;d<2;d++) 
		{
		  if(
		     TMath::Abs(fParticle->GetPdgCode()) == 211
		     //	      TMath::Abs(fParticle->GetPdgCode()) == 211 ||
		     //	      TMath::Abs(fParticle->GetPdgCode()) == 321 ||
		     //	      TMath::Abs(fParticle->GetPdgCode()) == 2212
		     )
		    {
		      fTrack = new AliFlowTrackSimple();
		      fTrack->SetPt(fParticle->Pt() );
		      fTrack->SetEta(fParticle->Eta() );
		      fTrack->SetPhi(fParticle->Phi() );
		      fTrack->SetForIntegratedFlow(kTRUE);
		      fTrack->SetForDifferentialFlow(kTRUE);

		      if (fTrack->UseForIntegratedFlow())
			{ fSelParticlesInt++; }
		      if (fTrack->UseForDifferentialFlow())
			{ fSelParticlesDiff++; }
		      fGoodTracks++;
		      fEvent->TrackCollection()->Add(fTrack);
		    }
			/*
		  else if(
			  TMath::Abs(fParticle->GetPdgCode()) == 2212
			  )
		    {
		      fTrack = new AliFlowTrackSimple();
		      fTrack->SetPt(fParticle->Pt() );
		      fTrack->SetEta(fParticle->Eta() );
		      fTrack->SetPhi(fParticle->Phi() );
		      fTrack->SetForIntegratedFlow(kFALSE);
		      fTrack->SetForDifferentialFlow(kTRUE);

		      if (fTrack->UseForIntegratedFlow())
			{ fSelParticlesInt++; }
		      if (fTrack->UseForDifferentialFlow())
			{ fSelParticlesDiff++; }
		      fGoodTracks++;
		      fEvent->TrackCollection()->Add(fTrack);     
		    }
			*/
		}
	    }
	  ftrkN++; 
	}
    }

  else {                                  //normal loop
    while (fGoodTracks < fN && ftrkN < fNumberOfInputTracks) {
      fInput->GetEntry(ftrkN);   //get input particle
      //cut on tracks
      if (TMath::Abs(fParticle->Eta()) < 0.2)
	{
	  if(
	     TMath::Abs(fParticle->GetPdgCode()) == 211
	     //	      TMath::Abs(fParticle->GetPdgCode()) == 211 ||
	     //	      TMath::Abs(fParticle->GetPdgCode()) == 321 ||
	     //	      TMath::Abs(fParticle->GetPdgCode()) == 2212
	     )
	    {
	      fTrack = new AliFlowTrackSimple() ;
	      fTrack->SetPt(fParticle->Pt() );
	      fTrack->SetEta(fParticle->Eta() );
	      fTrack->SetPhi(fParticle->Phi() );
	      fTrack->SetForIntegratedFlow(kTRUE);
	      fTrack->SetForDifferentialFlow(kTRUE);

	      if (fTrack->UseForIntegratedFlow())
		{ fSelParticlesInt++; }
	      if (fTrack->UseForDifferentialFlow())
		{ fSelParticlesDiff++; }
	      fGoodTracks++;
	      fEvent->TrackCollection()->Add(fTrack) ;  	     
	    }
	  /*	  else if(
		  TMath::Abs(fParticle->GetPdgCode()) == 211
		  )
	    {
	      fTrack = new AliFlowTrackSimple();
	      fTrack->SetPt(fParticle->Pt() );
	      fTrack->SetEta(fParticle->Eta() );
	      fTrack->SetPhi(fParticle->Phi() );
	      fTrack->SetForIntegratedFlow(kFALSE);
	      fTrack->SetForDifferentialFlow(kTRUE);

	      if (fTrack->UseForIntegratedFlow())
		{ fSelParticlesInt++; }
	      if (fTrack->UseForDifferentialFlow())
		{ fSelParticlesDiff++; }
	      fGoodTracks++;
	      fEvent->TrackCollection()->Add(fTrack);  	     
	    }
	  */
	}
      
      ftrkN++; 
    }
  }

  fEvent-> SetEventNSelTracksIntFlow(fSelParticlesInt);  
  fEvent->SetNumberOfTracks(fGoodTracks);
  cout<<" fGoodTracks = "<<fGoodTracks<<endl;
  cout << "  fSelectedTracksInt = " << fSelParticlesInt << endl;  
  return fEvent;
  
}

//-----------------------------------------------------------------------   
AliFlowEventSimple* AliFlowEventSimpleMaker::FillTracks(AliMCEvent* fInput)
{
  //Fills the event from the MC kinematic information
  
  Int_t fNumberOfInputTracks = fInput->GetNumberOfTracks() ;
  cerr<<"fInput->GetNumberOfTracks() = "<<fNumberOfInputTracks<<endl;
 
  fEvent = new AliFlowEventSimple(10);
    
  //Int_t fN = 256; //multiplicity for chi=1
  Int_t fN = fNumberOfInputTracks;
  Int_t fGoodTracks = 0;
  Int_t ftrkN = 0;
  Int_t fSelParticlesDiff = 0;
  Int_t fSelParticlesInt = 0;

   
  //normal loop
  while (fGoodTracks < fN && ftrkN < fNumberOfInputTracks) {
    AliMCParticle* fParticle = fInput->GetTrack(ftrkN);   //get input particle
    //cut on tracks
    if (TMath::Abs(fParticle->Eta()) < 0.2)
      {
	if(
	   TMath::Abs(fParticle->Particle()->GetPdgCode()) == 211
	   //	      TMath::Abs(fParticle->Particle()->GetPdgCode()) == 211 ||
	   //	      TMath::Abs(fParticle->Particle()->GetPdgCode()) == 321 ||
	   //	      TMath::Abs(fParticle->Particle()->GetPdgCode()) == 2212
	   )
	  {
	    fTrack = new AliFlowTrackSimple() ;
	    fTrack->SetPt(fParticle->Pt() );
	    fTrack->SetEta(fParticle->Eta() );
	    fTrack->SetPhi(fParticle->Phi() );
	    fTrack->SetForIntegratedFlow(kTRUE);
	    fTrack->SetForDifferentialFlow(kTRUE);

	    if (fTrack->UseForIntegratedFlow())
	      { fSelParticlesInt++; }
	    if (fTrack->UseForDifferentialFlow())
	      { fSelParticlesDiff++; }
	    fGoodTracks++;
	    fEvent->TrackCollection()->Add(fTrack) ;  	     
	  }
	  /*	  else if(
		  TMath::Abs(fParticle->Particle()->GetPdgCode()) == 211
		  )
	    {
	      fTrack = new AliFlowTrackSimple();
	      fTrack->SetPt(fParticle->Pt() );
	      fTrack->SetEta(fParticle->Eta() );
	      fTrack->SetPhi(fParticle->Phi() );
	      fTrack->SetForIntegratedFlow(kFALSE);
	      fTrack->SetForDifferentialFlow(kTRUE);

	      if (fTrack->UseForIntegratedFlow())
		{ fSelParticlesInt++; }
	      if (fTrack->UseForDifferentialFlow())
		{ fSelParticlesDiff++; }
	      fGoodTracks++;
	      fEvent->TrackCollection()->Add(fTrack);  	     
	    }
	  */
      }
      
    ftrkN++; 
  }
  
  fEvent-> SetEventNSelTracksIntFlow(fSelParticlesInt);  
  fEvent->SetNumberOfTracks(fGoodTracks);
  cout<<" fGoodTracks = "<<fGoodTracks<<endl;
  cout << "  fSelectedTracksInt = " << fSelParticlesInt << endl;  
  return fEvent;
  

}


//-----------------------------------------------------------------------   
AliFlowEventSimple* AliFlowEventSimpleMaker::FillTracks(AliESDEvent* fInput)
{
  //Fills the event from the ESD
  
  Int_t fNumberOfInputTracks = fInput->GetNumberOfTracks() ;
  cerr<<"fInput->GetNumberOfTracks() = "<<fNumberOfInputTracks<<endl;
  
  fEvent = new AliFlowEventSimple(10);
    
  //Int_t fN = 256; //multiplicity for chi=1
  Int_t fN = fNumberOfInputTracks;
  Int_t fGoodTracks = 0;
  Int_t ftrkN = 0;
  Int_t fSelParticlesDiff = 0;
  Int_t fSelParticlesInt = 0;


  //normal loop
  while (fGoodTracks < fN && ftrkN < fNumberOfInputTracks) {
    AliESDtrack* fParticle = fInput->GetTrack(ftrkN);   //get input particle
    //cut on tracks
    if (TMath::Abs(fParticle->Eta()) < 0.2)
      {
	fTrack = new AliFlowTrackSimple() ;
	fTrack->SetPt(fParticle->Pt() );
	fTrack->SetEta(fParticle->Eta() );
	fTrack->SetPhi(fParticle->Phi() );
	fTrack->SetForIntegratedFlow(kTRUE);
	fTrack->SetForDifferentialFlow(kTRUE);

	if (fTrack->UseForIntegratedFlow())
	  { fSelParticlesInt++; }
	if (fTrack->UseForDifferentialFlow())
	  { fSelParticlesDiff++; }
	fGoodTracks++;
	fEvent->TrackCollection()->Add(fTrack) ;  	     
      }
      
    ftrkN++; 
  }
  
  fEvent-> SetEventNSelTracksIntFlow(fSelParticlesInt);  
  fEvent->SetNumberOfTracks(fGoodTracks);
  cout<<" fGoodTracks = "<<fGoodTracks<<endl;
  cout << "  fSelectedTracksInt = " << fSelParticlesInt << endl;  
  return fEvent;


}



//-----------------------------------------------------------------------   
AliFlowEventSimple* AliFlowEventSimpleMaker::FillTracks(AliAODEvent* fInput)
{
  //Fills the event from the AOD
  
  Int_t fNumberOfInputTracks = fInput->GetNumberOfTracks() ;
  cerr<<"fInput->GetNumberOfTracks() = "<<fNumberOfInputTracks<<endl;
  
  fEvent = new AliFlowEventSimple(10);
    
  //Int_t fN = 256; //multiplicity for chi=1
  Int_t fN = fNumberOfInputTracks;
  Int_t fGoodTracks = 0;
  Int_t ftrkN = 0;
  Int_t fSelParticlesDiff = 0;
  Int_t fSelParticlesInt = 0;

  
  //normal loop
  while (fGoodTracks < fN && ftrkN < fNumberOfInputTracks) {
    AliAODTrack* fParticle = fInput->GetTrack(ftrkN);   //get input particle
    //cut on tracks
    if (TMath::Abs(fParticle->Eta()) < 0.2)
      {
	fTrack = new AliFlowTrackSimple() ;
	fTrack->SetPt(fParticle->Pt() );
	fTrack->SetEta(fParticle->Eta() );
	fTrack->SetPhi(fParticle->Phi() );
	fTrack->SetForIntegratedFlow(kTRUE);
	fTrack->SetForDifferentialFlow(kTRUE);

	if (fTrack->UseForIntegratedFlow())
	  { fSelParticlesInt++; }
	if (fTrack->UseForDifferentialFlow())
	  { fSelParticlesDiff++; }
	fGoodTracks++;
	fEvent->TrackCollection()->Add(fTrack) ;  	     
      }
      
    ftrkN++; 
  }
  
  fEvent-> SetEventNSelTracksIntFlow(fSelParticlesInt);  
  fEvent->SetNumberOfTracks(fGoodTracks);
  cout<<" fGoodTracks = "<<fGoodTracks<<endl;
  cout << "  fSelectedTracksInt = " << fSelParticlesInt << endl;  
  return fEvent;
  
}
//-----------------------------------------------------------------------   
AliFlowEventSimple*  AliFlowEventSimpleMaker::FillTracks(AliESDEvent* fInput, AliMCEvent* fInputMc, Int_t fOption)
{
  //fills the event with tracks from the ESD and kinematics from the MC info via the track label

  if (!(fOption ==0 || fOption ==1)) {
    cout<<"WRONG OPTION IN AliFlowEventSimpleMaker::FillTracks(AliESDEvent* fInput, AliMCEvent* fInputMc, Int_t fOption)"<<endl;
    exit(1);
  }

  Int_t fNumberOfInputTracks = fInput->GetNumberOfTracks() ;
  cerr<<"fInput->GetNumberOfTracks() = "<<fNumberOfInputTracks<<endl;
  
  fEvent = new AliFlowEventSimple(10);
    
  //Int_t fN = 256; //multiplicity for chi=1
  Int_t fN = fNumberOfInputTracks;
  Int_t fGoodTracks = 0;
  Int_t ftrkN = 0;
  Int_t fSelParticlesDiff = 0;
  Int_t fSelParticlesInt = 0;

  //normal loop
  while (fGoodTracks < fN && ftrkN < fNumberOfInputTracks) {
    AliESDtrack* fParticle = fInput->GetTrack(ftrkN);   //get input particle
    //get Label
    Int_t fLabel = fParticle->GetLabel();
    //match to mc particle
    AliMCParticle* fMcParticle = fInputMc->GetTrack(TMath::Abs(fLabel));
    
    //check
    if (TMath::Abs(fParticle->GetLabel())!=fMcParticle->Label()) cout<<"fParticle->GetLabel()!=fMcParticle->Label() "<<fParticle->GetLabel()<<"  "<<fMcParticle->Label()<<endl;
    
    //cut on tracks
    if (TMath::Abs(fParticle->Eta()) < 0.2)
      {
	if(
	   TMath::Abs(fMcParticle->Particle()->GetPdgCode()) == 211 //pions
	   //	      TMath::Abs(fMcParticle->Particle()->GetPdgCode()) == 211 ||
	   //	      TMath::Abs(fMcParticle->Particle()->GetPdgCode()) == 321 ||
	   //	      TMath::Abs(fMcParticle->Particle()->GetPdgCode()) == 2212
	   )
	  {
	    if(fOption == 0) { //take the PID from the MC & the kinematics from the ESD
	      fTrack = new AliFlowTrackSimple() ;
	      fTrack->SetPt(fParticle->Pt() );
	      fTrack->SetEta(fParticle->Eta() );
	      fTrack->SetPhi(fParticle->Phi() );
	      fTrack->SetForIntegratedFlow(kTRUE);
	      fTrack->SetForDifferentialFlow(kTRUE);
	    }
	    else if (fOption == 1) { //take the PID and kinematics from the MC
	      fTrack = new AliFlowTrackSimple() ;
	      fTrack->SetPt(fMcParticle->Pt() );
	      fTrack->SetEta(fMcParticle->Eta() );
	      fTrack->SetPhi(fMcParticle->Phi() );
	      fTrack->SetForIntegratedFlow(kTRUE);
	      fTrack->SetForDifferentialFlow(kTRUE);
	    }
	    else { cout<<"Not a valid option"<<endl; }
	    if (fTrack->UseForIntegratedFlow())
	      { fSelParticlesInt++; }
	    if (fTrack->UseForDifferentialFlow())
	      { fSelParticlesDiff++; }
	    fGoodTracks++;
	    fEvent->TrackCollection()->Add(fTrack) ;  	     
	  }
      }
    ftrkN++; 
  }
  
  fEvent-> SetEventNSelTracksIntFlow(fSelParticlesInt);  
  fEvent->SetNumberOfTracks(fGoodTracks);
  cout<<" fGoodTracks = "<<fGoodTracks<<endl;
  cout << "  fSelectedTracksInt = " << fSelParticlesInt << endl;  
  return fEvent;


}



/*
$Log$
*/ 

