// $Id$
// Category: event
//
// See the class description in the header file.

#include "AliSteppingAction.h"
#include "AliSteppingActionMessenger.h"
#include "AliGlobals.h"

#include <G4Track.hh>
#include <G4SteppingManager.hh>

const G4double AliSteppingAction::fgTolerance = 1e-6*mm;

AliSteppingAction::AliSteppingAction()
  : fKeptStepPoint(G4ThreeVector()),
    fLoopVerboseLevel(1),
    fStandardVerboseLevel(0),
    fLoopStepCounter(0)
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

// private methods

void AliSteppingAction::PrintTrackInfo(const G4Track* track) const
{
// Prints the track info
// - taken from private G4TrackingManager::Verbose()
// and the standard header for verbose tracking
// - taken from G4SteppingVerbose::TrackingStarted().
// ---

  // print track info
  G4cout << endl;
  G4cout << "*******************************************************"
         << "**************************************************"
         << endl;
  G4cout << "* G4Track Information: " 
         << "  Particle = " << track->GetDefinition()->GetParticleName() 
         << "," 
	 << "   Track ID = " << track->GetTrackID() 
         << "," 
	 << "   Parent ID = " << track->GetParentID() 
         << endl;
  G4cout << "*******************************************************"
         << "**************************************************"
         << endl;
  G4cout << endl;
      
  // print header    
#ifdef G4_USE_G4BESTUNIT_FOR_VERBOSE
    G4cout << setw( 5) << "Step#"  << " "
           << setw( 8) << "X"      << "     "
	   << setw( 8) << "Y"      << "     "  
	   << setw( 8) << "Z"      << "     "
	   << setw( 9) << "KineE"  << "     "
	   << setw( 8) << "dE"     << "     "  
	   << setw(12) << "StepLeng"   << " "  
	   << setw(12) << "TrackLeng"  << " "
	   << setw(12) << "NextVolume" << " "
	   << setw( 8) << "ProcName"   << endl;	     
#else
    G4cout << setw( 5) << "Step#"      << " "
	   << setw( 8) << "X(mm)"      << " "
	   << setw( 8) << "Y(mm)"      << " "  
	   << setw( 8) << "Z(mm)"      << " "
	   << setw( 9) << "KinE(MeV)"  << " "
	   << setw( 8) << "dE(MeV)"    << " "  
	   << setw( 8) << "StepLeng"   << " "  
	   << setw( 9) << "TrackLeng"  << " "
	   << setw(11) << "NextVolume" << " "
	   << setw( 8) << "ProcName"   << endl;	     
#endif
}

// public methods

void AliSteppingAction::UserSteppingAction(const G4Step* step)
{
// After processing the given number of steps (kCheckNofSteps)
// the particle position is compared with the previus one
// - in case the distance is less than fgTolerance value
// the verbose mode is switched on, particle is let 
// to process a small number of steps (kMaxNofLoopSteps)
// and then stopped and killed.
// ---

  G4Track* track = step->GetTrack();  

  // reset parameters at beginning of tracking
  G4int stepNumber = track->GetCurrentStepNumber();
  if (stepNumber == 0) {
    fKeptStepPoint = G4ThreeVector();
    fLoopStepCounter = 0;
    fStandardVerboseLevel = fpSteppingManager->GetverboseLevel();
    return;
  }  
    
  if (fLoopStepCounter) {
    // count steps after detecting looping track
    fLoopStepCounter++;
    if (fLoopStepCounter == kMaxNofLoopSteps) {

      // stop the looping track
      track->SetTrackStatus(fStopAndKill);      

      // reset back parameters
      fpSteppingManager->SetVerboseLevel(fStandardVerboseLevel);
      fKeptStepPoint = G4ThreeVector();
      fLoopStepCounter = 0;
    }  
  }  

  if (stepNumber % kCheckNofSteps == 0) {  
    // detect looping track
    G4ThreeVector newStepPoint = step->GetPreStepPoint()->GetPosition();
    G4double trajectory = (newStepPoint-fKeptStepPoint).mag();
    if (trajectory < fgTolerance) {

      // print looping info
      if (fLoopVerboseLevel > 0) {
        G4cout << "*** Particle is looping. ***" << endl;
	if (fStandardVerboseLevel == 0) PrintTrackInfo(track);
      }	
      // set loop verbose level 
      fpSteppingManager->SetVerboseLevel(fLoopVerboseLevel);
      
      fLoopStepCounter++;
    }  
    fKeptStepPoint = newStepPoint;
  }  
}
