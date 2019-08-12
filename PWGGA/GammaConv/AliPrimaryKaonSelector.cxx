/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. 	*
 *																			*
 * Authors: Friederike Bock
 * 2019. Transformed to Kaons by A. Marin													*
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
// Class reconstructing primary kaons
//---------------------------------------------
////////////////////////////////////////////////


#include "AliPrimaryKaonSelector.h"
#include "AliPrimaryKaonCuts.h"
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



ClassImp(AliPrimaryKaonSelector)

//________________________________________________________________________
AliPrimaryKaonSelector::AliPrimaryKaonSelector(const char *name) : AliAnalysisTaskSE(name),
    fKaonCuts(0),
    fPosKaonsIndex(),
    fNegKaonsIndex(),
    fEventIsSelected(kFALSE)
{
    // Default constructor
    DefineInput(0, TChain::Class());
}

//________________________________________________________________________
AliPrimaryKaonSelector::~AliPrimaryKaonSelector()
{
    // default deconstructor

   
}
//________________________________________________________________________
void AliPrimaryKaonSelector::Init()
{
    // Initialize function to be called once before analysis

    if(fKaonCuts == 0){
    //  fKaonCuts=AliConversionCuts::GetStandardCuts2010PbPb();
        AliError("No Cut Selection initialized");
    }

}

//________________________________________________________________________
void AliPrimaryKaonSelector::UserCreateOutputObjects()
{
    // Create User Output Objects
}

//________________________________________________________________________
void AliPrimaryKaonSelector::UserExec(Option_t *){
    // User Exec
	fEventIsSelected=ProcessEvent(fInputEvent,fMCEvent);
}

//________________________________________________________________________
Bool_t AliPrimaryKaonSelector::ProcessEvent(AliVEvent *inputEvent,AliMCEvent *mcEvent)
{
	//Reset the index

	fPosKaonsIndex.clear();
	fNegKaonsIndex.clear();


	fInputEvent=inputEvent;
	fMCEvent=mcEvent;

	if(!fInputEvent){
		AliError("No Input event");
		return kFALSE;
	}

	if(!fKaonCuts){AliError("No ConversionCuts");return kFALSE;}


	if(fInputEvent->IsA()==AliESDEvent::Class()){
		ProcessESDs();
	}

    if(fInputEvent->IsA()==AliAODEvent::Class()){
      ProcessAODs();
    }


	return kTRUE;
}

///________________________________________________________________________
Bool_t AliPrimaryKaonSelector::ProcessESDs(){
	// Process ESD V0s for conversion photon reconstruction
	AliESDEvent *fESDEvent=dynamic_cast<AliESDEvent*>(fInputEvent);
	if(fESDEvent){
		for(Int_t currentTrackIndex=0;currentTrackIndex<fESDEvent->GetNumberOfTracks();currentTrackIndex++){
            AliESDtrack *fCurrentTrack = dynamic_cast<AliESDtrack*> (fESDEvent->GetTrack(currentTrackIndex));
			if(!fCurrentTrack){
				printf("Requested Track does not exist");
				continue;
			}

			if (  fKaonCuts->KaonIsSelected( fCurrentTrack ) ) {
				if( fCurrentTrack->GetSign() > 0.0 ){
					fPosKaonsIndex.push_back(currentTrackIndex);
				} else {
					fNegKaonsIndex.push_back(currentTrackIndex);
				}
			}
		}
	}
	return kTRUE;
}

///________________________________________________________________________
Bool_t AliPrimaryKaonSelector::ProcessAODs(){
    // process AOD primary kaons
    AliAODEvent *fAODEvent=dynamic_cast<AliAODEvent*>(fInputEvent);
    if(fAODEvent){
        for(Int_t currentTrackIndex=0;currentTrackIndex<fAODEvent->GetNumberOfTracks();currentTrackIndex++){
            AliAODTrack *fCurrentTrack = dynamic_cast<AliAODTrack*> (fAODEvent->GetTrack(currentTrackIndex));
            if(!fCurrentTrack){
                printf("Requested Track does not exist");
                continue;
            }

            Float_t sign = ( fCurrentTrack ? fCurrentTrack->Charge() : fCurrentTrack->GetSign() );
            if (  fKaonCuts->KaonIsSelectedAOD( fCurrentTrack ) ) {
                if( sign > 0.0 ){
                    fPosKaonsIndex.push_back(currentTrackIndex);
                } else {
                    fNegKaonsIndex.push_back(currentTrackIndex);
                }
            }
        }
    }
    return kTRUE;
}

//________________________________________________________________________
void AliPrimaryKaonSelector::Terminate(Option_t *)
{

}
