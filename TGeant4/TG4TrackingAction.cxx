// $Id$
// Category: event
//
// Author: I.Hrivnacova
//
// Class TG4TrackingAction
// -----------------------
// See the class description in the header file.

#include "TG4TrackingAction.h"
#include "TG4StepManager.h"
#include "TG4VSensitiveDetector.h"
#include "TG4SDServices.h"
#include "TG4Globals.h"

//_____________________________________________________________________________
TG4TrackingAction::TG4TrackingAction() {
//
}

//_____________________________________________________________________________
TG4TrackingAction::TG4TrackingAction(const TG4TrackingAction& right) {
//
  TG4Globals::Exception("TG4TrackingAction is protected from copying.");
}

//_____________________________________________________________________________
TG4TrackingAction::~TG4TrackingAction() {
//
}

// operators

//_____________________________________________________________________________
TG4TrackingAction& 
TG4TrackingAction::operator=(const TG4TrackingAction &right)
{
  // check assignement to self
  if (this == &right) return *this;
  
  TG4Globals::Exception("TG4TrackingAction is protected from assigning.");

  return *this;
}

// public methods

//_____________________________________________________________________________
void TG4TrackingAction::PreUserTrackingAction(const G4Track* track)
{
// Called by G4 kernel before starting tracking.
// ---

  // set step manager status
  TG4StepManager* stepManager = TG4StepManager::Instance();
  stepManager->SetStep((G4Track*)track, kVertex);
  
  // call pre-tracking action of derived class
  PreTrackingAction(track);

  // let sensitive detector process vertex step
  // (this ensures compatibility with G3 that
  // makes first step of zero length)
   
  G4VPhysicalVolume* pv = stepManager->GetCurrentPhysicalVolume();
  
  if (!pv) {
    G4String text = "TG4TrackingAction::PreUserTrackingAction: \n";
    text = text + "   Cannot locate track vertex."; 
    TG4Globals::Exception(text);
  }  
  
#ifdef TGEANT4_DEBUG
  TG4VSensitiveDetector* tsd
    = TG4SDServices::Instance()
         ->GetSensitiveDetector(
	      pv->GetLogicalVolume()->GetSensitiveDetector());

  if (tsd) tsd->UserProcessHits((G4Track*)track, 0);
#else
  TG4VSensitiveDetector* tsd
    = (TG4VSensitiveDetector*) pv->GetLogicalVolume()->GetSensitiveDetector();

  if (tsd) tsd->UserProcessHits((G4Track*)track, 0);
#endif  
}

//_____________________________________________________________________________
void TG4TrackingAction::PostUserTrackingAction(const G4Track* track)
{
// Called by G4 kernel after finishing tracking.
// ---

  // call post-tracking action of derived class
  PostTrackingAction(track);
}

