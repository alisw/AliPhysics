// $Id$
// Category: event
//
// See the class description in the header file.

#include "TG4SteppingAction.h"
#include "TG4StepManager.h"
#include "TG4VSensitiveDetector.h"
#include "TG4Globals.h"

TG4SteppingAction::TG4SteppingAction() {
//
}

TG4SteppingAction::TG4SteppingAction(const TG4SteppingAction& right) {
//
  TG4Globals::Exception("TG4SteppingAction is protected from copying.");
}

TG4SteppingAction::~TG4SteppingAction() {
//
}

// operators

TG4SteppingAction& 
TG4SteppingAction::operator=(const TG4SteppingAction &right)
{
  // check assignement to self
  if (this == &right) return *this;
  
  TG4Globals::Exception("TG4SteppingAction is protected from assigning.");

  return *this;
}

// public methods

void TG4SteppingAction::UserSteppingAction(const G4Step* step)
{
// Called by G4 kernel at the end of each step.
// ---

  // call stepping action of derived class
  SteppingAction(step);

  // let sensitive detector process boundary step
  // if crossing geometry border
  // (this ensures compatibility with G3 that
  // makes boundary step of zero length)

  if (step->GetPostStepPoint()->GetStepStatus() == fGeomBoundary &&
      step->GetTrack()->GetTrackStatus() == fAlive &&
      step->GetTrack()->GetNextVolume() != 0) {

    G4VSensitiveDetector* sd
      = step->GetPostStepPoint()->GetPhysicalVolume()
          ->GetLogicalVolume()->GetSensitiveDetector();

    if (sd) {
      TG4VSensitiveDetector* tsd = dynamic_cast<TG4VSensitiveDetector*>(sd);
      if (tsd) 
        tsd->ProcessHitsOnBoundary((G4Step*)step);
      else {
        G4String text = "TG4SteppingAction:::UserSteppingAction: \n";
        text = text + "   Unknown sensitive detector type"; 
        TG4Globals::Exception(text);
      }   	
    } 
  }  
}

