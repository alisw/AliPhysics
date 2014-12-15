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
//   Filtering task:  
//   Selection of only 1% of the events for which to keep  
//   the ESD friend
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
#include "AliAnalysisTaskFilterSteer.h"
#include "TRandom.h"


ClassImp(AliAnalysisTaskFilterSteer)


//________________________________________________________________________
AliAnalysisTaskFilterSteer::AliAnalysisTaskFilterSteer():
AliAnalysisTaskFilter(),
fFraction(0.01),
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
AliAnalysisTaskFilterSteer::AliAnalysisTaskFilterSteer(const char* name):
AliAnalysisTaskFilter(name),
fFraction(0.01),
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
AliAnalysisTaskFilterSteer::~AliAnalysisTaskFilterSteer()
{

	// dtor

}  

//________________________________________________________________________
void AliAnalysisTaskFilterSteer::Init()
{
	// Initialization
	
	return;
}

//________________________________________________________________________
void AliAnalysisTaskFilterSteer::UserCreateOutputObjects()
{
	//
	// Create the output container
	//
	
	return;
}

//________________________________________________________________________
void AliAnalysisTaskFilterSteer::UserExec(Option_t */*option*/)
{

	// 
	// Applying random selection of the events

	fESDInput = dynamic_cast<AliESDEvent*>(InputEvent()); // get the input ESD
	fESDfriendInput = InputFriend();  // get the input friend
	if(!fESDInput) {
		printf("AliAnalysisTaskFilterSteer::Exec(): no ESD \n");
		return;
	} 
	if(!fESDfriendInput) {
		printf("AliAnalysisTaskFilterSteer::Exec(): no ESDfriend \n");
		return;
	} 

	// attach ESDfriend
	
	AliESDfriend* esdFriendOutput = (AliESDfriend*)ESDfriend();  
	AliDebug(3,Form("Number of ESD tracks in input = %d ",fESDInput->GetNumberOfTracks()));
	AliDebug(3,Form("Number of tracks in input friends = %d ",fESDfriendInput->GetNumberOfTracks()));
	AliDebug(3,Form("Number of tracks in output friendsNew before filtering = %d ",esdFriendOutput->GetNumberOfTracks()));
	
	//
	//  keeping all the tracks for the randomly "fFraction" of the total number of events
	//

	for (Int_t i = 0; i< fESDInput->GetNumberOfTracks(); i++){
		AliESDfriendTrack* tOld = (AliESDfriendTrack*)fESDfriendInput->GetTrack(i);
		AddFriendTrackAt(tOld,i);
	}			 
	AliDebug(2,Form("Number of tracks in output friendsNew after filtering = %d ",esdFriendOutput->GetNumberOfTracks()));
	AliDebug(2,Form("Number of tracks in output friendsNew after filtering with GetEntries() = %d ",esdFriendOutput->GetEntriesInTracks()));

	return;
}

//________________________________________________________________________
void AliAnalysisTaskFilterSteer::Terminate(Option_t */*option*/)
{
	// Terminate analysis
	//
	AliDebug(2,"AliAnalysisTaskFilterSteer: Terminate() \n");
	
	return;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskFilterSteer::UserSelectESDfriendForCurrentEvent()
{
	// 
	// Selecting or discarding current event
	//
	Double_t number = gRandom->Rndm();
	if (number<fFraction){
		// keeping event
		AliDebug(2,Form("*****************Selecting event (number = %f)",number));
		return kTRUE;	
	}
	AliDebug(2,Form("*****************Skipping event (number = %f)",number));
	return kFALSE;
}
