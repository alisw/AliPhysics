/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. 	*
 *																			*
 * Authors: Friederike Bock													*
 * Converted to deuteron by A. Marin
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
// Class reconstructing primary deuterons
//---------------------------------------------
////////////////////////////////////////////////


#include "AliPrimaryDeuteronSelector.h"
#include "AliPrimaryDeuteronCuts.h"
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



ClassImp(AliPrimaryDeuteronSelector)

//________________________________________________________________________
AliPrimaryDeuteronSelector::AliPrimaryDeuteronSelector(const char *name) : AliAnalysisTaskSE(name),
    fDeuteronCuts(0),
    fPosDeuteronsIndex(),
    fNegDeuteronsIndex(),
    fEventIsSelected(kFALSE)
{
    // Default constructor
    DefineInput(0, TChain::Class());
}

//________________________________________________________________________
AliPrimaryDeuteronSelector::~AliPrimaryDeuteronSelector()
{
    // default deconstructor

   
}
//________________________________________________________________________
void AliPrimaryDeuteronSelector::Init()
{
    // Initialize function to be called once before analysis

    if(fDeuteronCuts == 0){
    //  fDeuteronCuts=AliConversionCuts::GetStandardCuts2010PbPb();
        AliError("No Cut Selection initialized");
    }

}

//________________________________________________________________________
void AliPrimaryDeuteronSelector::UserCreateOutputObjects()
{
    // Create User Output Objects
}

//________________________________________________________________________
void AliPrimaryDeuteronSelector::UserExec(Option_t *){
    // User Exec
	fEventIsSelected=ProcessEvent(fInputEvent,fMCEvent);
}

//________________________________________________________________________
Bool_t AliPrimaryDeuteronSelector::ProcessEvent(AliVEvent *inputEvent,AliMCEvent *mcEvent)
{
	//Reset the index

	fPosDeuteronsIndex.clear();
	fNegDeuteronsIndex.clear();


	fInputEvent=inputEvent;
	fMCEvent=mcEvent;

	if(!fInputEvent){
		AliError("No Input event");
		return kFALSE;
	}

	if(!fDeuteronCuts){AliError("No ConversionCuts");return kFALSE;}


	if(fInputEvent->IsA()==AliESDEvent::Class()){
		ProcessESDs();
	}

    if(fInputEvent->IsA()==AliAODEvent::Class()){
      ProcessAODs();
    }


	return kTRUE;
}

///________________________________________________________________________
Bool_t AliPrimaryDeuteronSelector::ProcessESDs(){
	// Process ESD V0s for conversion photon reconstruction
	AliESDEvent *fESDEvent=dynamic_cast<AliESDEvent*>(fInputEvent);
	if(fESDEvent){
		for(Int_t currentTrackIndex=0;currentTrackIndex<fESDEvent->GetNumberOfTracks();currentTrackIndex++){
            AliESDtrack *fCurrentTrack = dynamic_cast<AliESDtrack*> (fESDEvent->GetTrack(currentTrackIndex));
			if(!fCurrentTrack){
				printf("Requested Track does not exist");
				continue;
			}

			if (  fDeuteronCuts->DeuteronIsSelected( fCurrentTrack ) ) {
				if( fCurrentTrack->GetSign() > 0.0 ){
					fPosDeuteronsIndex.push_back(currentTrackIndex);
				} else {
					fNegDeuteronsIndex.push_back(currentTrackIndex);
				}
			}
		}
	}
	return kTRUE;
}

///________________________________________________________________________
Bool_t AliPrimaryDeuteronSelector::ProcessAODs(){
    // process AOD primary pions
    AliAODEvent *fAODEvent=dynamic_cast<AliAODEvent*>(fInputEvent);
    if(fAODEvent){
        for(Int_t currentTrackIndex=0;currentTrackIndex<fAODEvent->GetNumberOfTracks();currentTrackIndex++){
            AliAODTrack *fCurrentTrack = dynamic_cast<AliAODTrack*> (fAODEvent->GetTrack(currentTrackIndex));
            if(!fCurrentTrack){
                printf("Requested Track does not exist");
                continue;
            }

            Float_t sign = ( fCurrentTrack ? fCurrentTrack->Charge() : fCurrentTrack->GetSign() );
            if (  fDeuteronCuts->DeuteronIsSelectedAOD( fCurrentTrack ) ) {
                if( sign > 0.0 ){
                    fPosDeuteronsIndex.push_back(currentTrackIndex);
                } else {
                    fNegDeuteronsIndex.push_back(currentTrackIndex);
                }
            }
        }
    }
    return kTRUE;
}

//________________________________________________________________________
void AliPrimaryDeuteronSelector::Terminate(Option_t *)
{

}
