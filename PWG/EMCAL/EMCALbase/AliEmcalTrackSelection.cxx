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
#include <AliEmcalTrackSelection.h>
#include <TObjArray.h>
#include <TClonesArray.h>
#include <AliESDtrackCuts.h>
#include <AliEmcalESDtrackCutsWrapper.h>
#include <AliLog.h>
#include <AliVCuts.h>
#include <AliVTrack.h>
#include <AliVEvent.h>
#include <iostream>

/// \cond CLASSIMP
ClassImp(AliEmcalManagedObject)
ClassImp(AliEmcalTrackSelection)
/// \endcond

AliEmcalTrackSelection::AliEmcalTrackSelection() :
	TObject(),
	fListOfTracks(NULL),
  fListOfTrackBitmaps(NULL),
  fTrackBitmap(64),
	fListOfCuts(NULL),
	fSelectionModeAny(kFALSE)
{
}

AliEmcalTrackSelection::AliEmcalTrackSelection(const AliEmcalTrackSelection& ref):
	TObject(ref),
	fListOfTracks(NULL),
	fListOfTrackBitmaps(NULL),
	fTrackBitmap(64),
	fListOfCuts(NULL),
	fSelectionModeAny(kFALSE)
{
	if(ref.fListOfTracks) fListOfTracks = new TObjArray(*(ref.fListOfTracks));
	if(ref.fListOfTrackBitmaps) fListOfTrackBitmaps = new TClonesArray(*(ref.fListOfTrackBitmaps));
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
		if(ref.fListOfTrackBitmaps) fListOfTrackBitmaps = new TClonesArray(*(ref.fListOfTrackBitmaps));
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
	if(fListOfTrackBitmaps) delete fListOfTrackBitmaps;
	if(fListOfCuts) delete fListOfCuts;
}

void AliEmcalTrackSelection::AddTrackCuts(AliVCuts *cuts){
  AliInfoStream() << "Adding trackc cuts " << cuts->GetName() << " of type " << cuts->IsA()->GetName() << std::endl;
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
    fListOfCuts->Add(new AliEmcalManagedObject(mycuts, true));
  } 
}

void AliEmcalTrackSelection::AddTrackCuts(TObjArray *cuts){
  for(auto c : *cuts){
    AliVCuts *cuts = dynamic_cast<AliVCuts*>(c);
    if(cuts){
      AddTrackCuts(cuts);
    } else {
      AliErrorStream() << "Object not inheriting from AliVCuts - not added to track selection" << std::endl;
    }
  }
}

Int_t AliEmcalTrackSelection::GetNumberOfCutObjects() const {
  if(!fListOfCuts) return 0;
  return fListOfCuts->GetEntries();
}

AliVCuts* AliEmcalTrackSelection::GetTrackCuts(Int_t icut) {
  if(!fListOfCuts) return NULL;
  if(icut < fListOfCuts->GetEntries()){
    AliEmcalManagedObject *ptr = static_cast<AliEmcalManagedObject *>(fListOfCuts->At(icut));
    return static_cast<AliVCuts *>(ptr->GetObject());
  }

  return NULL;
}

TObjArray* AliEmcalTrackSelection::GetAcceptedTracks(const TClonesArray* const tracks)
{
  if (!fListOfTracks) {
    fListOfTracks = new TObjArray;
  }
  else {
    fListOfTracks->Clear();
  }

  if (!fListOfTrackBitmaps) {
    fListOfTrackBitmaps = new TClonesArray("TBits", 1000);
    fListOfTrackBitmaps->SetOwner(kTRUE);
  }
  else {
    fListOfTrackBitmaps->Delete();
  }

  TIter next(tracks);
  AliVTrack* track = 0;
  Int_t i = 0;
  while((track = static_cast<AliVTrack*>(next()))) {
    if (IsTrackAccepted(track)) {
      fListOfTracks->AddLast(track);
    }
    else {
      fListOfTracks->AddLast(0);
    }
    new ((*fListOfTrackBitmaps)[i]) TBits(fTrackBitmap);
    i++;
  }
  return fListOfTracks;
}

TObjArray* AliEmcalTrackSelection::GetAcceptedTracks(const AliVEvent* const event)
{
  if (!fListOfTracks) {
    fListOfTracks = new TObjArray;
  }
  else {
    fListOfTracks->Clear();
  }

  if (!fListOfTrackBitmaps) {
    fListOfTrackBitmaps = new TClonesArray("TBits", 1000);
    fListOfTrackBitmaps->SetOwner(kTRUE);
  }
  else {
    fListOfTrackBitmaps->Delete();
  }

  for(int itrk = 0; itrk < event->GetNumberOfTracks(); itrk++){
    AliVTrack *trk = static_cast<AliVTrack*>(event->GetTrack(itrk));
    if (IsTrackAccepted(trk)) {
      fListOfTracks->AddLast(trk);
    }
    else {
      fListOfTracks->AddLast(trk);
    }
    new ((*fListOfTrackBitmaps)[itrk]) TBits(fTrackBitmap);
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
