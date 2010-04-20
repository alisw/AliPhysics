/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

/*$Id$

*/

/////////////////////////////////////////////////////////////
//
//   Filtering task to be run on STEERING level
//   1. Select randomly the event to be kept fully  0 cuttof - fKeepFraction - 
//      gRandmom->Rndm() used. 
//   2. Pt() cut used - fPtCut 
//
// /////////////////////////////////////////////////////////////

#include <TTree.h>
#include <TChain.h>

#include "AliLog.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliAnalysisTaskFilter.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskFilterSTEER.h"
#include "TRandom.h"


ClassImp(AliAnalysisTaskFilterSTEER)


//________________________________________________________________________
AliAnalysisTaskFilterSTEER::AliAnalysisTaskFilterSTEER():
AliAnalysisTaskFilter(),
fPtCut(0),
fKeepFraction(0),
fESDInput(0),
fESDfriendInput(0)
{
	// Constructor
	
	// Define input and output slots here
	// Input slot #0 works with a TChain
	DefineInput(0, TChain::Class());
	// Output slot #0 writes into a TTree
	//DefineOutput(0,TTree::Class());  
}

//________________________________________________________________________
AliAnalysisTaskFilterSTEER::AliAnalysisTaskFilterSTEER(const char* name, Double_t ptCut, Double_t fractionCut):
AliAnalysisTaskFilter(name),
fPtCut(ptCut),
fKeepFraction(fractionCut),
fESDInput(0),
fESDfriendInput(0)
{
	// Constructor
	
	// Define input and output slots here
	// Input slot #0 works with a TChain
	DefineInput(0, TChain::Class());
	// Output slot #0 writes into a TTree
	//DefineOutput(0,TTree::Class());  
}

//________________________________________________________________________
AliAnalysisTaskFilterSTEER::~AliAnalysisTaskFilterSTEER()
{

	// dtor

}  

//________________________________________________________________________
void AliAnalysisTaskFilterSTEER::Init()
{
	// Initialization
	
	return;
}

//________________________________________________________________________
void AliAnalysisTaskFilterSTEER::UserCreateOutputObjects()
{
	//
	// Create the output container
	//
	
	return;
}

//________________________________________________________________________
void AliAnalysisTaskFilterSTEER::UserExec(Option_t */*option*/)
{

	AliInfo("Filling Friends");

	fESDInput = dynamic_cast<AliESDEvent*>(InputEvent()); // get the input ESD
	fESDfriendInput = (AliESDfriend*)(fESDInput->FindListObject("AliESDfriend"));
	if(!fESDInput) {
		printf("AliAnalysisTaskFilterSTEER::Exec(): no ESD \n");
		return;
	} 
	if(!fESDfriendInput) {
		printf("AliAnalysisTaskFilterSTEER::Exec(): no ESDfriend \n");
		return;
	} 
	// attach ESDfriend
	
	AliESDfriend* esdFriendOutput = (AliESDfriend*)ESDfriend();  
	AliInfo(Form("Number of ESD tracks in input = %d ",fESDInput->GetNumberOfTracks()));
	AliInfo(Form("Number of tracks in input friends = %d ",fESDfriendInput->GetNumberOfTracks()));
	AliInfo(Form("Number of tracks in output friendsNew before filtering = %d ",esdFriendOutput->GetNumberOfTracks()));
	
	//
	// select all for fraction of events - event selected randomly
	// event number 0 always stored - to define the friends layout
	//
	AliESDfriendTrack* tNull = new AliESDfriendTrack();
	if (fESDInput->GetEventNumberInFile()==0 || gRandom->Rndm()<fKeepFraction){
	  //select all tracks for slected fraction
	  AliInfo(Form("Keeping Event number  %d",fESDInput->GetEventNumberInFile()));
	  for (Int_t i = 0; i< fESDInput->GetNumberOfTracks(); i++){
	    // keep friend
	    AliESDfriendTrack* tOld = (AliESDfriendTrack*)fESDfriendInput->GetTrack(i);
	    if (tOld)  AddFriendTrackAt(tOld,i);
	    else  AddFriendTrackAt(tNull,i);
	  }
	}
	//
	//
	//
	for (Int_t i = 0; i< fESDInput->GetNumberOfTracks(); i++){
	  // keep friend
	  AliESDfriendTrack* tOld = (AliESDfriendTrack*)fESDfriendInput->GetTrack(i);
	  Bool_t isOK=kFALSE;
	  if (tOld) 
	    if (tOld->GetTPCOut()) 
	      if (TMath::Abs(tOld->GetTPCOut()->Pt())>fPtCut) isOK=kTRUE;
	  if (isOK){
	    AliInfo(Form("Keeping %d-th track",i));
	    AddFriendTrackAt(tOld,i);
	  }else{
	    AddFriendTrackAt(tNull,i);	    
	  }
	}			 
	AliInfo(Form("Number of tracks in output friendsNew after filtering = %d ",esdFriendOutput->GetNumberOfTracks()));
	AliInfo(Form("Number of tracks in output friendsNew after filtering with GetEntries() = %d ",esdFriendOutput->GetEntriesInTracks()));
	return;
}

//________________________________________________________________________
void AliAnalysisTaskFilterSTEER::Terminate(Option_t */*option*/)
{
	// Terminate analysis
	//
	AliDebug(2,"AliAnalysisTaskFilterSTEER: Terminate() \n");
	
	return;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskFilterSTEER::UserSelectESDfriendForCurrentEvent()
{
  // 
  // Selecting or discarding current event
  //
  return kTRUE;	
}
