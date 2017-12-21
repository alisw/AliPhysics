/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. 	*
 *																			*
 * Authors: Friederike Bock													*
 * Version 1.0																*
 *																			*
 * Permission to use, copy, modify and distribute this software and its	 	*
 * documentation strictly for non-commercial purposes is hereby granted	 	*
 * without fee, provided that the above copyright notice appears in all	 	*
 * copies and that both the copyright notice and this permission notice	 	*
 * appear in the supporting documentation. The authors make no claims		*
 * about the suitability of this software for any purpose. It is			*
 * provided "as is" without express or implied warranty.					*
 **************************************************************************/

////////////////////////////////////////////////
//--------------------------------------------- 
// Class reconstructing primary electrons
//---------------------------------------------
////////////////////////////////////////////////


#include "AliPrimaryPionSelector.h"
#include "AliPrimaryPionCuts.h"
#include "AliAODv0.h"
#include "AliESDv0.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "TVector.h"
#include "AliESDtrack.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "TChain.h"
#include "AliMCEvent.h"

class iostream;

using namespace std;



ClassImp(AliPrimaryPionSelector)

//________________________________________________________________________
AliPrimaryPionSelector::AliPrimaryPionSelector(const char *name) : AliAnalysisTaskSE(name),
    fPionCuts(0),
    fPosPionsIndex(),
    fNegPionsIndex(),
    fEventIsSelected(kFALSE)
{
    // Default constructor
    DefineInput(0, TChain::Class());
}

//________________________________________________________________________
AliPrimaryPionSelector::~AliPrimaryPionSelector()
{
    // default deconstructor

   
}
//________________________________________________________________________
void AliPrimaryPionSelector::Init()
{
    // Initialize function to be called once before analysis

    if(fPionCuts == 0){
    //  fPionCuts=AliConversionCuts::GetStandardCuts2010PbPb();
        AliError("No Cut Selection initialized");
    }

}

//________________________________________________________________________
void AliPrimaryPionSelector::UserCreateOutputObjects()
{
    // Create User Output Objects
}

//________________________________________________________________________
void AliPrimaryPionSelector::UserExec(Option_t *){
    // User Exec
	fEventIsSelected=ProcessEvent(fInputEvent,fMCEvent);
}

//________________________________________________________________________
Bool_t AliPrimaryPionSelector::ProcessEvent(AliVEvent *inputEvent,AliMCEvent *mcEvent)
{
	//Reset the index

	fPosPionsIndex.clear();
	fNegPionsIndex.clear();


	fInputEvent=inputEvent;
	fMCEvent=mcEvent;

	if(!fInputEvent){
		AliError("No Input event");
		return kFALSE;
	}

	if(!fPionCuts){AliError("No ConversionCuts");return kFALSE;}


	if(fInputEvent->IsA()==AliESDEvent::Class()){
		ProcessESDs();
	}

	//if(fInputEvent->IsA()==AliAODEvent::Class()){
	//GetAODConversionGammas();
	//}


	return kTRUE;
}

///________________________________________________________________________
Bool_t AliPrimaryPionSelector::ProcessESDs(){
	// Process ESD V0s for conversion photon reconstruction
	AliESDEvent *fESDEvent=dynamic_cast<AliESDEvent*>(fInputEvent);
	if(fESDEvent){
		for(Int_t currentTrackIndex=0;currentTrackIndex<fESDEvent->GetNumberOfTracks();currentTrackIndex++){
			AliESDtrack *fCurrentTrack = (AliESDtrack*)(fESDEvent->GetTrack(currentTrackIndex));
			if(!fCurrentTrack){
				printf("Requested Track does not exist");
				continue;
			}
			if (  fPionCuts->PionIsSelected( fCurrentTrack ) ) {
				if( fCurrentTrack->GetSign() > 0.0 ){
					fPosPionsIndex.push_back(currentTrackIndex);
				} else {
					fNegPionsIndex.push_back(currentTrackIndex);
				}
			}
		}
	}
	return kTRUE;
}


//________________________________________________________________________
void AliPrimaryPionSelector::Terminate(Option_t *)
{

}
