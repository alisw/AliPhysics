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

/*$Id$*/

/////////////////////////////////////////////////////////////
//
//   Test task
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
#include "AliAnalysisTaskFilterFriendSecond.h"


ClassImp(AliAnalysisTaskFilterFriendSecond)


//________________________________________________________________________
AliAnalysisTaskFilterFriendSecond::AliAnalysisTaskFilterFriendSecond():
AliAnalysisTaskFilter(),
fESDInput(0),
fESDfriendInput(0)
{
	// Constructor
	
	// Define input and output slots here
	// Input slot #0 works with a TChain
	DefineInput(0, TChain::Class());
	// Output slot #0 writes into a TTree
	//	DefineOutput(0,TTree::Class());  
}

//________________________________________________________________________
AliAnalysisTaskFilterFriendSecond::AliAnalysisTaskFilterFriendSecond(const char* name):
AliAnalysisTaskFilter(name),
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
AliAnalysisTaskFilterFriendSecond::~AliAnalysisTaskFilterFriendSecond()
{

	// dtor

}  

//________________________________________________________________________
void AliAnalysisTaskFilterFriendSecond::Init()
{
	// Initialization
	
	return;
}

//________________________________________________________________________
void AliAnalysisTaskFilterFriendSecond::UserCreateOutputObjects()
{
	//
	// Create the output container
	//
	
	return;
}

//________________________________________________________________________
void AliAnalysisTaskFilterFriendSecond::UserExec(Option_t */*option*/)
{

	//
	// Filtering Friends
	//

	fESDInput = dynamic_cast<AliESDEvent*>(InputEvent()); // get the input ESD
	fESDfriendInput = InputFriend();  // get the input friend
	if(!fESDInput) {
		printf("AliAnalysisTaskFilterFriendSecond::Exec(): no ESD \n");
		return;
	} 
	if(!fESDfriendInput) {
		printf("AliAnalysisTaskFilterFriendSecond::Exec(): no ESDfriend \n");
		return;
	} 
	// attach ESDfriend
	
	AliESDfriend* esdFriendOutput = (AliESDfriend*)ESDfriend();  
	AliDebug(3,Form("Number of ESD tracks in input = %d ",fESDInput->GetNumberOfTracks()));
	AliDebug(3,Form("Number of tracks in input friends = %d ",fESDfriendInput->GetNumberOfTracks()));
	AliDebug(3,Form("Number of tracks in output friendsNew before filtering = %d ",esdFriendOutput->GetNumberOfTracks()));
	
	for (Int_t i = 0; i< fESDInput->GetNumberOfTracks(); i++){
		if (i%3 == 0){
			// keep friend
			AliDebug(2,Form("Keeping %d-th track",i));
			AliESDfriendTrack* tOld = (AliESDfriendTrack*)fESDfriendInput->GetTrack(i);
			// debugging printouts
			AliDebug(3,Form("1P of the %d-th track = %f",i,tOld->Get1P()));
			AddFriendTrackAt(tOld,i);
			//	tOld->Dump();
		}
		else {
			//discard friend 
			SkipFriendTrackAt(i);
		}
	} 
	AliDebug(2,Form("Number of tracks in output friendsNew after filtering with GetEntries() = %d ",esdFriendOutput->GetEntriesInTracks()));
	return;
}

//________________________________________________________________________
void AliAnalysisTaskFilterFriendSecond::Terminate(Option_t */*option*/)
{
	// Terminate analysis
	//
	AliDebug(2,"AliAnalysisTaskFilterFriendSecond: Terminate() \n");
	
	return;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskFilterFriendSecond::UserSelectESDfriendForCurrentEvent()
{
	// 
	// Selecting or discarding current event
	//

	/*		
	fESDInput = dynamic_cast<AliESDEvent*>(InputEvent()); // get the input ESD
	if ((fESDInput->GetNumberOfTracks())%2 == 1){
		AliDebug(2,"****************Selecting event");
		return kTRUE;
	}
	
	AliDebug(2,"****************Discarding event");	
	return kFALSE;
	*/
		AliDebug(2,"****************Selecting event");
	return kTRUE;
}
