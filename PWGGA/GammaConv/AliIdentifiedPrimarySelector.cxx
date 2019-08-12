/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. 	*
 *																			*
 * Authors: Friederike Bock
 * 2019. Transformed to Identified Charged Particles by A. Marin													*
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
// Class reconstructing identified primary charged particles
//---------------------------------------------
////////////////////////////////////////////////


#include "AliIdentifiedPrimarySelector.h"
#include "AliIdentifiedPrimaryCuts.h"
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



ClassImp(AliIdentifiedPrimarySelector)

//________________________________________________________________________
AliIdentifiedPrimarySelector::AliIdentifiedPrimarySelector(const char *name) : AliAnalysisTaskSE(name),
    fIdentifiedCuts(0),
    fPosIdentifiedIndex(),
    fNegIdentifiedIndex(),
    fEventIsSelected(kFALSE)
{
    // Default constructor
    DefineInput(0, TChain::Class());
}

//________________________________________________________________________
AliIdentifiedPrimarySelector::~AliIdentifiedPrimarySelector()
{
    // default deconstructor

   
}
//________________________________________________________________________
void AliIdentifiedPrimarySelector::Init()
{
    // Initialize function to be called once before analysis

    if(fIdentifiedCuts == 0){
    //  fKaonCuts=AliConversionCuts::GetStandardCuts2010PbPb();
        AliError("No Cut Selection initialized");
    }

}

//________________________________________________________________________
void AliIdentifiedPrimarySelector::UserCreateOutputObjects()
{
    // Create User Output Objects
}

//________________________________________________________________________
void AliIdentifiedPrimarySelector::UserExec(Option_t *){
    // User Exec
	fEventIsSelected=ProcessEvent(fInputEvent,fMCEvent);
}

//________________________________________________________________________
Bool_t AliIdentifiedPrimarySelector::ProcessEvent(AliVEvent *inputEvent,AliMCEvent *mcEvent)
{
	//Reset the index

	fPosIdentifiedIndex.clear();
	fNegIdentifiedIndex.clear();


	fInputEvent=inputEvent;
	fMCEvent=mcEvent;

	if(!fInputEvent){
		AliError("No Input event");
		return kFALSE;
	}

	if(!fIdentifiedCuts){AliError("No Identified charged particle cuts");return kFALSE;}


	if(fInputEvent->IsA()==AliESDEvent::Class()){
		ProcessESDs();
	}

    if(fInputEvent->IsA()==AliAODEvent::Class()){
      ProcessAODs();
    }


	return kTRUE;
}

///________________________________________________________________________
Bool_t AliIdentifiedPrimarySelector::ProcessESDs(){
	// Process ESD for Identified Primary charged particles
	AliESDEvent *fESDEvent=dynamic_cast<AliESDEvent*>(fInputEvent);
	if(fESDEvent){
		for(Int_t currentTrackIndex=0;currentTrackIndex<fESDEvent->GetNumberOfTracks();currentTrackIndex++){
            AliESDtrack *fCurrentTrack = dynamic_cast<AliESDtrack*> (fESDEvent->GetTrack(currentTrackIndex));
			if(!fCurrentTrack){
				printf("Requested Track does not exist");
				continue;
			}

			if (  fIdentifiedCuts->IdentifiedIsSelected( fCurrentTrack ) ) {
				if( fCurrentTrack->GetSign() > 0.0 ){
					fPosIdentifiedIndex.push_back(currentTrackIndex);
				} else {
					fNegIdentifiedIndex.push_back(currentTrackIndex);
				}
			}
		}
	}
	return kTRUE;
}

///________________________________________________________________________
Bool_t AliIdentifiedPrimarySelector::ProcessAODs(){
    // process AOD identified primary charged particles
    AliAODEvent *fAODEvent=dynamic_cast<AliAODEvent*>(fInputEvent);
    if(fAODEvent){
        for(Int_t currentTrackIndex=0;currentTrackIndex<fAODEvent->GetNumberOfTracks();currentTrackIndex++){
            AliAODTrack *fCurrentTrack = dynamic_cast<AliAODTrack*> (fAODEvent->GetTrack(currentTrackIndex));
            if(!fCurrentTrack){
                printf("Requested Track does not exist");
                continue;
            }

            Float_t sign = ( fCurrentTrack ? fCurrentTrack->Charge() : fCurrentTrack->GetSign() );
            if (  fIdentifiedCuts->IdentifiedIsSelectedAOD( fCurrentTrack ) ) {
                if( sign > 0.0 ){
                    fPosIdentifiedIndex.push_back(currentTrackIndex);
                } else {
                    fNegIdentifiedIndex.push_back(currentTrackIndex);
                }
            }
        }
    }
    return kTRUE;
}

//________________________________________________________________________
void AliIdentifiedPrimarySelector::Terminate(Option_t *)
{

}
