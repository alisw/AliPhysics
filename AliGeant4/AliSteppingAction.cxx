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
  G4cout << G4endl;
  G4cout << "*******************************************************"
         << "**************************************************"
         << G4endl;
  G4cout << "* G4Track Information: " 
         << "  Particle = " << track->GetDefinition()->GetParticleName() 
         << "," 
	 << "   Track ID = " << track->GetTrackID() 
         << "," 
	 << "   Parent ID = " << track->GetParentID() 
         << G4endl;
  G4cout << "*******************************************************"
         << "**************************************************"
         << G4endl;
  G4cout << G4endl;
      
  // print header    
#ifdef G4_USE_G4BESTUNIT_FOR_VERBOSE
    G4cout << G4std::setw( 5) << "Step#"  << " "
           << G4std::setw( 8) << "X"      << "     "
	   << G4std::setw( 8) << "Y"      << "     "  
	   << G4std::setw( 8) << "Z"      << "     "
	   << G4std::setw( 9) << "KineE"  << "     "
	   << G4std::setw( 8) << "dE"     << "     "  
	   << G4std::setw(12) << "StepLeng"   << " "  
	   << G4std::setw(12) << "TrackLeng"  << " "
	   << G4std::setw(12) << "NextVolume" << " "
	   << G4std::setw( 8) << "ProcName"   << G4endl;	     
#else
    G4cout << G4std::setw( 5) << "Step#"      << " "
	   << G4std::setw( 8) << "X(mm)"      << " "
	   << G4std::setw( 8) << "Y(mm)"      << " "  
	   << G4std::setw( 8) << "Z(mm)"      << " "
	   << G4std::setw( 9) << "KinE(MeV)"  << " "
	   << G4std::setw( 8) << "dE(MeV)"    << " "  
	   << G4std::setw( 8) << "StepLeng"   << " "  
	   << G4std::setw( 9) << "TrackLeng"  << " "
	   << G4std::setw(11) << "NextVolume" << " "
	   << G4std::setw( 8) << "ProcName"   << G4endl;	     
#endif
}

// public methods

void AliSteppingAction::UserSteppingAction(const G4Step* step)
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
    G4bool kill = false;
    if (trajectory < fgkTolerance) {

      // print looping info
      if (fLoopVerboseLevel > 0) {
        G4cout << "*** Particle is looping. ***" << G4endl;
	if (fStandardVerboseLevel == 0) PrintTrackInfo(track);
      }	
      kill = true;
    }
    
    if (stepNumber> kMaxNofSteps) { 
      
      // print looping info
      if (fLoopVerboseLevel > 0) {
        G4cout << "*** Particle reached max step number ("
	       << kMaxNofSteps << "). ***" << G4endl;
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
