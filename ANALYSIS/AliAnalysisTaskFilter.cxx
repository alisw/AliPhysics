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

/* $Id$ */

//////////////////////////////////////////////////////////////////////////
//
//  Base class for filtering friends
//
//////////////////////////////////////////////////////////////////////////
 
#include <TChain.h>
#include <TFile.h>
#include <TList.h>

#include "AliAnalysisTaskFilter.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataSlot.h"
#include "AliESDEvent.h"
#include "AliESD.h"
#include "AliVEvent.h"
#include "AliESDHandler.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliESDfriend.h"
#include "AliESDfriendTrack.h"


ClassImp(AliAnalysisTaskFilter)

////////////////////////////////////////////////////////////////////////

AliAnalysisTaskFilter::AliAnalysisTaskFilter():
	AliAnalysisTask(),
	fDebug(0),
	fEntry(0),
	fInputEvent(0x0),
	fInputHandler(0x0),
	fOutputESDfriend(0x0),
	fTreeEF(0x0),
	fInputESDfriend(0x0)
{
	//
	// Default constructor
	//
}

//______________________________________________________________________

AliAnalysisTaskFilter::AliAnalysisTaskFilter(const char* name):
	AliAnalysisTask(name, "AnalysisTaskFilter"),
	fDebug(0),
	fEntry(0),
	fInputEvent(0x0),
	fInputHandler(0x0),
	fOutputESDfriend(0x0),
	fTreeEF(0x0),
	fInputESDfriend(0x0)
{
	//
	// Default constructor
	//

	DefineInput (0, TChain::Class());
	DefineOutput(0,  TTree::Class());
}

//______________________________________________________________________

AliAnalysisTaskFilter::AliAnalysisTaskFilter(const AliAnalysisTaskFilter& obj):
	AliAnalysisTask(obj),
	fDebug(0),
	fEntry(0),
	fInputEvent(0x0),
	fInputHandler(0x0),
	fOutputESDfriend(0x0),
	fTreeEF(0x0),
	fInputESDfriend(0x0)
{
	//
	// Copy constructor
	//

	fDebug        = obj.fDebug;
	fEntry        = obj.fEntry;
	fInputEvent   = obj.fInputEvent;
	fInputHandler = obj.fInputHandler;
	fOutputESDfriend = obj.fOutputESDfriend;
	fTreeEF        = obj.fTreeEF;    
	fInputESDfriend = obj.fInputESDfriend;
}


//______________________________________________________________________

AliAnalysisTaskFilter& AliAnalysisTaskFilter::operator=(const AliAnalysisTaskFilter& other)
{
	//
	// Assignment
	//

	AliAnalysisTask::operator=(other);
	fDebug        = other.fDebug;
	fEntry        = other.fEntry;
	fInputEvent   = other.fInputEvent;
	fInputHandler = other.fInputHandler;
	fOutputESDfriend = other.fOutputESDfriend;
	fTreeEF        = other.fTreeEF;    
	fInputESDfriend = other.fInputESDfriend;
	return *this;
}


//______________________________________________________________________

void AliAnalysisTaskFilter::ConnectInputData(Option_t* /*option*/)
{
	//
	// Connect the input data
	//

	if (fDebug > 1) printf("AnalysisTaskFilter::ConnectInputData() \n");
	fInputHandler = (AliInputEventHandler*) 
		((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
	if (fInputHandler) {
		fInputEvent = fInputHandler->GetEvent();
		if (fInputEvent){
			fInputESDfriend = (AliESDfriend*)(fInputEvent->FindListObject("AliESDfriend"));
			if (!fInputESDfriend){
				AliError("No friend found");
			}
		}
		else {
			AliError("No Input Event found, the friend will remain empty");
		}
	} 
	else {
		AliError("No Input Event Handler connected") ; 
		return ; 
	}
}

//______________________________________________________________________

void AliAnalysisTaskFilter::CreateOutputObjects()
{
	//
	// Create the output container
	//

	if (fDebug > 1) printf("AnalysisTaskFilter::CreateOutPutData() \n");

	AliESDHandler* handler = (AliESDHandler*) ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
    
	if (handler) {
		fTreeEF = handler->GetTree();
	}
	else {
		AliWarning("No AOD Event Handler connected.") ; 
	}

	UserCreateOutputObjects();
}

//______________________________________________________________________

void AliAnalysisTaskFilter::Exec(Option_t* option)
{
	//
	// Exec analysis of one event
	//

	if (fDebug > 1) AliInfo("AliAnalysisTaskFilter::Exec() \n");

	if( fInputHandler ) {
		fEntry = fInputHandler->GetReadEntry();
	}
    
	   
	if ( !((Entry()-1)%100) && fDebug > 0) {
		AliInfo(Form("%s ----> Processing event # %lld", CurrentFileName(), Entry()));
	}
	AliESDHandler* handler = (AliESDHandler*)((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
	
	if (UserSelectESDfriendForCurrentEvent()){
		// Call the user analysis only if the event was selected   
		handler->SelectEventForFriends();
		fOutputESDfriend   = handler->GetESDfriend();
		UserExec(option);
		// copy the VZERO friend only if it is not already there
		if (fOutputESDfriend->GetVZEROfriend() == 0x0){
			AliDebug(2,"Copying VZERO friend object");
			AliESDVZEROfriend* vZEROfriend = fInputESDfriend->GetVZEROfriend();
			fOutputESDfriend->SetVZEROfriend(vZEROfriend);
		}
	}
	else {
		// Event not selected
		AliDebug(2,"The event was not selected");
	}

	AliAnalysisDataSlot *out0 = GetOutputSlot(0);
	if (out0 && out0->IsConnected()) PostData(0, fTreeEF);    
}

//______________________________________________________________________

const char* AliAnalysisTaskFilter::CurrentFileName()
{
	// Returns the current file name    
	if( fInputHandler ){
		return fInputHandler->GetTree()->GetCurrentFile()->GetName();
	}
	else return "";
}

//______________________________________________________________________

void AliAnalysisTaskFilter::AddFriendTrackAt(AliESDfriendTrack* t, Int_t index)
{
	//
	// Adds the friend track at the i-th position in the TClonesArray
	// of the ESD friend tracks
	//

	AliESDfriendTrack* currentTrack = (AliESDfriendTrack*)fOutputESDfriend->GetTrack(index);
	if(currentTrack){
		if (currentTrack->TestSkipBit()){
			AliDebug(2,Form("Friend at index %d already there but dummy - the skip bit will be set to FALSE", index));
			t->SetSkipBit(kFALSE);
		}
		else{
			AliDebug(2,Form("Friend at index %d already there and not dummy", index));
			return;
		}
	}
	else{
		AliDebug(2,Form("Track at %d not there yet ",index));
	}
	AliDebug(2,Form("Adding track at %d",index));
	fOutputESDfriend->AddTrackAt(t,index);
	return;
}

//______________________________________________________________________

void AliAnalysisTaskFilter::SkipFriendTrackAt(Int_t index)
{
	//
	// Skip the friend track at the i-th position in the TClonesArray
	// of the ESD friend tracks
	//

	AliESDfriendTrack* currentTrack = (AliESDfriendTrack*)fOutputESDfriend->GetTrack(index);
	if (currentTrack){
		AliDebug(2,Form("Track already there (no matter what validity) at %d, keeping it as it is", index));
	}
	else {
		AliDebug(2,Form("Adding NULL track at %d, and setting skip bit to TRUE",index));
		AliESDfriendTrack* tNull = new AliESDfriendTrack();
		tNull->SetSkipBit(kTRUE);
		fOutputESDfriend->AddTrackAt(tNull,index);
	}
	return;
}
