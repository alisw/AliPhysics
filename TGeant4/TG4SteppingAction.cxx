// $Id$
// Category: event
//
// See the class description in the header file.

#include "TG4SteppingAction.h"
#include "TG4VSensitiveDetector.h"
#include "TG4Globals.h"

#include <G4Track.hh>
#include <G4SteppingManager.hh>

TG4SteppingAction::TG4SteppingAction() 
  : fMaxNofSteps(kMaxNofSteps),
    fStandardVerboseLevel(0),
    fLoopVerboseLevel(1),
    fLoopStepCounter(0)
 {
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

// protected methods

void TG4SteppingAction::PrintTrackInfo(const G4Track* track) const
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

void TG4SteppingAction::UserSteppingAction(const G4Step* step)
{
// Called by G4 kernel at the end of each step.
// ---

  G4Track* track = step->GetTrack();  

  // reset parameters at beginning of tracking
  G4int stepNumber = track->GetCurrentStepNumber();
  if (stepNumber == 1) { 
    fStandardVerboseLevel = fpSteppingManager->GetverboseLevel();
    fLoopStepCounter = 0;
  }  
  else if (fLoopStepCounter) {
      // count steps after detecting looping track
      fLoopStepCounter++;
      if (fLoopStepCounter == kMaxNofLoopSteps) {

        // stop the looping track
        track->SetTrackStatus(fStopAndKill);      

        // reset back parameters
        fpSteppingManager->SetVerboseLevel(fStandardVerboseLevel);
        fLoopStepCounter = 0;
      } 
    }  
    else if (stepNumber> fMaxNofSteps) { 
      
       // print looping info
       if (fLoopVerboseLevel > 0) {
          G4cout << "*** Particle reached max step number ("
	         << fMaxNofSteps << "). ***" << G4endl;
	  if (fStandardVerboseLevel == 0) PrintTrackInfo(track);
       }	

       // set loop verbose level 
       fpSteppingManager->SetVerboseLevel(fLoopVerboseLevel);
      
       // start looping counter
       fLoopStepCounter++;
    }
      
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

