// $Id$
// Category: event
//
// Author: I. Hrivnacova
//
// Class AliSteppingAction
// -----------------------
// See the class description in the header file.

#include "AliSteppingAction.h"
#include "AliRun.h"

#include <G4Track.hh>
#include <G4SteppingManager.hh>

#include <math.h>

//_____________________________________________________________________________
AliSteppingAction::AliSteppingAction()
  : TG4SteppingAction(),
     fMessenger(this) {
//
}

//_____________________________________________________________________________
AliSteppingAction::~AliSteppingAction() {
//
}

// public methods

//_____________________________________________________________________________
void AliSteppingAction::SteppingAction(const G4Step* step)
{
// Stops particle if it gets outside of user defined tracking region.
// ---

  G4ThreeVector position 
    = step->GetPostStepPoint()->GetPosition();

  if (position.mag()    > gAlice->TrackingRmax() ||
      abs(position.z()) > gAlice->TrackingZmax()) {
 
    // print looping info
    if (fLoopVerboseLevel > 0) {
      G4cout << "*** Particle has reached user defined tracking region. ***" 
             << G4endl;
      if (fStandardVerboseLevel == 0) PrintTrackInfo(step->GetTrack());
    }  

    // stop the track
    step->GetTrack()->SetTrackStatus(fStopAndKill);      
  }
}
