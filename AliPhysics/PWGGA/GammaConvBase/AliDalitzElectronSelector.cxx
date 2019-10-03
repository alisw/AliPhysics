/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Pedro Gonzales                        *
 * Version 1.0                                                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

////////////////////////////////////////////////
//---------------------------------------------
// Class reconstructing primary electrons
//---------------------------------------------
////////////////////////////////////////////////


#include "AliDalitzElectronSelector.h"
#include "AliDalitzElectronCuts.h"
#include "AliAODv0.h"
#include "AliESDv0.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "TVector.h"
#include "AliESDtrack.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "TChain.h"
#include "AliDalitzAODESD.h"

class iostream;

using namespace std;



ClassImp(AliDalitzElectronSelector)

//________________________________________________________________________
AliDalitzElectronSelector::AliDalitzElectronSelector(const char *name) : AliAnalysisTaskSE(name),
  fElectronCuts(0),
  fPositronsIndex(),
  fElectronsIndex(),
  fEventIsSelected(kFALSE)
{
  // Default constructor
  DefineInput(0, TChain::Class());
}

//________________________________________________________________________
AliDalitzElectronSelector::~AliDalitzElectronSelector()
{
  // default deconstructor
}
//________________________________________________________________________
void AliDalitzElectronSelector::Init()
{
  // Initialize function to be called once before analysis

  if(fElectronCuts == 0){
  //  fElectronCuts=AliConversionCuts::GetStandardCuts2010PbPb();
      AliError("No Cut Selection initialized");
  }
}

//________________________________________________________________________
void AliDalitzElectronSelector::UserCreateOutputObjects()
{
  // Create User Output Objects
}

//________________________________________________________________________
void AliDalitzElectronSelector::UserExec(Option_t *){
  // User Exec
  fEventIsSelected=ProcessEvent(fInputEvent,fMCEvent);
}

//________________________________________________________________________
Bool_t AliDalitzElectronSelector::ProcessEvent(AliVEvent *inputEvent,AliMCEvent *mcEvent)
{
  //Reset the index

  fPositronsIndex.clear();
  fElectronsIndex.clear();


  fInputEvent=inputEvent;
  fMCEvent=mcEvent;

  if(!fInputEvent){
    AliError("No Input event");
    return kFALSE;
  }

  if(!fElectronCuts){AliError("No ConversionCuts");return kFALSE;}


  //if(fInputEvent->IsA()==AliESDEvent::Class()){
    ProcessESDs();
  //}

  //if(fInputEvent->IsA()==AliAODEvent::Class()){
  //GetAODConversionGammas();
  //}
  return kTRUE;
}

///________________________________________________________________________
Bool_t AliDalitzElectronSelector::ProcessESDs(){
  // Process ESD V0s for conversion photon reconstruction
    AliAODEvent *fAODEvent=0;
    AliESDEvent *fESDEvent=0;
    Int_t NumberOfTracks=0;

    if(fInputEvent->IsA()==AliESDEvent::Class()){
        fESDEvent=dynamic_cast<AliESDEvent*>(fInputEvent);
        NumberOfTracks=fESDEvent->GetNumberOfTracks();
    }
    else {fAODEvent=dynamic_cast<AliAODEvent*>(fInputEvent);
        NumberOfTracks=fAODEvent->GetNumberOfTracks();
    }
  //if(fESDEvent){
    for(Int_t currentTrackIndex=0;currentTrackIndex<NumberOfTracks;currentTrackIndex++){
      AliDalitzAODESD *fCurrentTrack =0;

      if(fESDEvent){
            fCurrentTrack= new AliDalitzAODESD((AliESDtrack*)(fESDEvent->GetTrack(currentTrackIndex)));
            fCurrentTrack->ComputeImpactParameter();
      }
      else {fCurrentTrack= new AliDalitzAODESD((AliAODTrack*)(fAODEvent->GetTrack(currentTrackIndex)));
        fCurrentTrack->ComputeImpactParameter(fAODEvent->GetPrimaryVertex(),fAODEvent->GetMagneticField());
      }

      if(!fCurrentTrack){
        printf("Requested Track does not exist");
        continue;
      }
      if (  fElectronCuts->ElectronIsSelected( fCurrentTrack ) ) {
        if( fCurrentTrack->GetSignG() > 0.0 ){
          fPositronsIndex.push_back(currentTrackIndex);
        } else {
          fElectronsIndex.push_back(currentTrackIndex);
        }
      }
        delete fCurrentTrack;
    }
  //}
  return kTRUE;
}

//________________________________________________________________________
void AliDalitzElectronSelector::Terminate(Option_t *)
{

}
