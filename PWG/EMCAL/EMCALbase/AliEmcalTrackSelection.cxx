/**************************************************************************
 * Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
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
#include "iostream"

/// \cond CLASSIMP
ClassImp(AliEmcalManagedObject)
ClassImp(AliEmcalTrackSelection)
/// \endcond

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
  AliInfoStream() << "Adding trackc cuts " << cuts->GetName() << " of type " << cuts->IsA()->GetName() << std::endl;
  if(!fListOfCuts){
    fListOfCuts = new TObjArray;
    fListOfCuts->SetOwner(true);
  }
  if(cuts) {
    fListOfCuts->Add(new AliEmcalManagedObject(cuts));
  }
}

void AliEmcalTrackSelection::AddTrackCuts(TObjArray *cuts){
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
