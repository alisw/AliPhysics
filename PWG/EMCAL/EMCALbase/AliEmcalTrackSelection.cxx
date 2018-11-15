/************************************************************************************
 * Copyright (C) 2017, Copyright Holders of the ALICE Collaboration                 *
 * All rights reserved.                                                             *
 *                                                                                  *
 * Redistribution and use in source and binary forms, with or without               *
 * modification, are permitted provided that the following conditions are met:      *
 *     * Redistributions of source code must retain the above copyright             *
 *       notice, this list of conditions and the following disclaimer.              *
 *     * Redistributions in binary form must reproduce the above copyright          *
 *       notice, this list of conditions and the following disclaimer in the        *
 *       documentation and/or other materials provided with the distribution.       *
 *     * Neither the name of the <organization> nor the                             *
 *       names of its contributors may be used to endorse or promote products       *
 *       derived from this software without specific prior written permission.      *
 *                                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND  *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED    *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE           *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY              *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES       *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND      *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT       *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                     *
 ************************************************************************************/
#include <TObjArray.h>
#include <TClonesArray.h>
#include "AliESDtrackCuts.h"
#include "AliEmcalESDtrackCutsWrapper.h"
#include "AliEmcalVCutsWrapper.h"
#include "AliEmcalTrackSelection.h"
#include "AliLog.h"
#include "AliVCuts.h"
#include "AliVTrack.h"
#include "AliVEvent.h"
#include <iostream>

/// \cond CLASSIMP
ClassImp(AliEmcalManagedObject)
ClassImp(AliEmcalTrackSelection)
/// \endcond

using namespace PWG::EMCAL;

AliEmcalTrackSelection::AliEmcalTrackSelection() :
	TObject(),
	fListOfTracks(NULL),
	fListOfCuts(NULL),
	fSelectionModeAny(kFALSE)
{
}

AliEmcalTrackSelection::AliEmcalTrackSelection(const AliEmcalTrackSelection& ref):
	TObject(ref),
	fListOfTracks(NULL),
	fListOfCuts(NULL),
	fSelectionModeAny(kFALSE)
{
	if(ref.fListOfTracks) fListOfTracks = new TObjArray(*(ref.fListOfTracks));
	if(ref.fListOfCuts){
	  fListOfCuts = new TObjArray;
	  fListOfCuts->SetOwner(true); // Ownership handled object-by-object by the smart pointer
	  for(auto cutIter : *(ref.fListOfCuts))
	    fListOfCuts->Add(new AliEmcalManagedObject(*(static_cast<AliEmcalManagedObject *>(cutIter))));
	}
}

AliEmcalTrackSelection& AliEmcalTrackSelection::operator=(const AliEmcalTrackSelection& ref) {
	TObject::operator=(ref);
	if(this != &ref){
		this->~AliEmcalTrackSelection();
		if(ref.fListOfTracks) fListOfTracks = new TObjArray(*(ref.fListOfTracks));
		if(ref.fListOfCuts){
		  fListOfCuts = new TObjArray;
		  fListOfCuts->SetOwner(true);  // Ownership handled object-by-object by the smart pointer
		  for(auto cutIter : *(ref.fListOfCuts))
		    fListOfCuts->Add(new AliEmcalManagedObject(*(static_cast<AliEmcalManagedObject *>(cutIter))));
		} else fListOfCuts = NULL;
	}
	return *this;
}

AliEmcalTrackSelection::~AliEmcalTrackSelection() {
	if(fListOfTracks) delete fListOfTracks;
	if(fListOfCuts) delete fListOfCuts;
}

void AliEmcalTrackSelection::AddTrackCuts(AliVCuts *cuts){
  AliInfoStream() << "Adding track cuts " << cuts->GetName() << " of type " << cuts->IsA()->GetName() << std::endl;
  if(!fListOfCuts){
    fListOfCuts = new TObjArray;
    fListOfCuts->SetOwner(true);
  }
  if(cuts) {
    // Special treatment for AliESDtrackCuts:
    // As the function IsSelected is not properly implemented for AliAODTracks
    // a wrapper needs to be used, which handles the expected behaviour for
    // both AliESDtracks and AliAODTracks
    AliVCuts *mycuts = cuts;
    if(AliESDtrackCuts *esdcuts = dynamic_cast<AliESDtrackCuts *>(cuts)) mycuts = new PWG::EMCAL::AliEmcalESDtrackCutsWrapper(esdcuts->GetName(), esdcuts);
    // Convert to AliEmcalCutBase
    fListOfCuts->Add(new AliEmcalManagedObject(new PWG::EMCAL::AliEmcalVCutsWrapper(mycuts), true));
  } 
}

void AliEmcalTrackSelection::AddTrackCuts(PWG::EMCAL::AliEmcalCutBase *cuts) {
  AliInfoStream() << "Adding track cuts " << cuts->GetName() << " of type " << cuts->IsA()->GetName() << std::endl;
  if(!fListOfCuts){
    fListOfCuts = new TObjArray;
    fListOfCuts->SetOwner(true);
  }
  if(cuts) {
    fListOfCuts->Add(new AliEmcalManagedObject(cuts));
  }
}

void AliEmcalTrackSelection::AddTrackCuts(TObjArray *cuts){
  if(!cuts) {
    AliErrorStream() << "Not setting cuts since cut array is null" << std::endl; 
    return;
  }
  for(auto c : *cuts){
    PWG::EMCAL::AliEmcalCutBase *emccuts = dynamic_cast<PWG::EMCAL::AliEmcalCutBase*>(c);
    if(emccuts){
      AddTrackCuts(emccuts);
    } else {
      AliVCuts *vcuts = dynamic_cast<AliVCuts *>(c);
      if(vcuts) {
        AddTrackCuts(vcuts);
      } else {
        AliErrorStream() << "Object not inheriting from AliVCuts - not added to track selection" << std::endl;
      }
    }
  }
}

Int_t AliEmcalTrackSelection::GetNumberOfCutObjects() const {
  if(!fListOfCuts) return 0;
  return fListOfCuts->GetEntries();
}

PWG::EMCAL::AliEmcalCutBase* AliEmcalTrackSelection::GetTrackCuts(Int_t icut) {
  if(!fListOfCuts) return NULL;
  if(icut < fListOfCuts->GetEntries()){
    AliEmcalManagedObject *ptr = static_cast<AliEmcalManagedObject *>(fListOfCuts->At(icut));
    return static_cast<PWG::EMCAL::AliEmcalCutBase *>(ptr->GetObject());
  }

  return NULL;
}

TObjArray* AliEmcalTrackSelection::GetAcceptedTracks(const TClonesArray* const tracks)
{
  if (!fListOfTracks) {
    fListOfTracks = new TObjArray;
    fListOfTracks->SetOwner(kTRUE);
  }
  else {
    fListOfTracks->Clear();
  }

  for(auto mytrack : *tracks) {
    fListOfTracks->AddLast(new PWG::EMCAL::AliEmcalTrackSelResultPtr(IsTrackAccepted(static_cast<AliVTrack *>(mytrack))));
  }
  return fListOfTracks;
}

TObjArray* AliEmcalTrackSelection::GetAcceptedTracks(const AliVEvent* const event)
{
  if (!fListOfTracks) {
    fListOfTracks = new TObjArray;
    fListOfTracks->SetOwner(kTRUE);
  }
  else {
    fListOfTracks->Clear();
  }

  for(int itrk = 0; itrk < event->GetNumberOfTracks(); itrk++){
    fListOfTracks->AddLast(new PWG::EMCAL::AliEmcalTrackSelResultPtr(IsTrackAccepted(static_cast<AliVTrack*>(event->GetTrack(itrk)))));
  }
  return fListOfTracks;
}

AliEmcalManagedObject::AliEmcalManagedObject():
    TObject(),
    fOwner(false),
    fManagedObject(nullptr)
{

}

AliEmcalManagedObject::AliEmcalManagedObject(TObject *managedObject, Bool_t owner):
    TObject(),
    fOwner(owner),
    fManagedObject(managedObject)
{
}

AliEmcalManagedObject::AliEmcalManagedObject(const AliEmcalManagedObject &ref):
    TObject(ref),
    fOwner(false),
    fManagedObject(ref.fManagedObject)
{
}

AliEmcalManagedObject &AliEmcalManagedObject::operator=(const AliEmcalManagedObject &ref){
  TObject::operator=(ref);

  if(this != &ref){
    Cleanup();
    fOwner = false;
    fManagedObject = ref.fManagedObject;
  }
  return *this;
}

void AliEmcalManagedObject::Cleanup(){
  if(fManagedObject && fOwner) delete fManagedObject;
  fManagedObject = nullptr;
}
