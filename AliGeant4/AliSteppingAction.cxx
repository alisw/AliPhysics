// $Id$
// Category: event
//
// See the class description in the header file.

#include "AliSteppingAction.h"
#include "AliSteppingActionMessenger.h"
#include "AliGlobals.h"

#include <G4Track.hh>
#include <G4SteppingManager.hh>

const G4double AliSteppingAction::fgkTolerance = 1e-6*mm;

AliSteppingAction::AliSteppingAction()
  : fKeptStepPoint(G4ThreeVector())
{
//
  fMessenger = new AliSteppingActionMessenger(this);
}

AliSteppingAction::AliSteppingAction(const AliSteppingAction& right) {
//
  AliGlobals::Exception("AliSteppingAction is protected from copying.");
}

AliSteppingAction::~AliSteppingAction() {
//
  delete fMessenger;
}

// operators

AliSteppingAction& 
AliSteppingAction::operator=(const AliSteppingAction &right)
{
  // check assignement to self
  if (this == &right) return *this;
  
  AliGlobals::Exception("AliSteppingAction is protected from assigning.");

  return *this;
}

// public methods

void AliSteppingAction::SteppingAction(const G4Step* step)
{
// After processing the given number of steps (kCheckNofSteps)
// the particle position is compared with the previus one
// - in case the distance is less than fgkTolerance value
// the verbose mode is switched on, particle is let 
// to process a small number of steps (kMaxNofLoopSteps)
// and then stopped and killed.
// ---

  G4Track* track = step->GetTrack();  

  // reset parameters at beginning of tracking
  G4int stepNumber = track->GetCurrentStepNumber();
  if (stepNumber == 1) {
    fKeptStepPoint = G4ThreeVector();
    return;
  }  
    
  if (stepNumber % kCheckNofSteps == 0) {  
    // detect looping track
    G4ThreeVector newStepPoint = step->GetPreStepPoint()->GetPosition();
    G4double trajectory = (newStepPoint-fKeptStepPoint).mag();
    G4bool kill = false;
    if (trajectory < fgkTolerance) {

      // print looping info
      if (fLoopVerboseLevel > 0) {
        G4cout << "*** Particle is looping. ***" << G4endl;
	if (fStandardVerboseLevel == 0) PrintTrackInfo(track);
      }	
      kill = true;
    }
    
    if (kill) {

      // set loop verbose level 
      fpSteppingManager->SetVerboseLevel(fLoopVerboseLevel);
      
      fLoopStepCounter++;
    }  
    fKeptStepPoint = newStepPoint;
  }  
}
