// $Id$
// Category: event
//
// See the class description in the header file.

#include "AliEventAction.h"
#include "AliEventActionMessenger.h"
#include "AliRun.h"
#include "AliTrackingAction.h"
#include "AliSensitiveDetector.h"
#include "AliGlobals.h"

#include "TG4GeometryManager.h"

#include <G4Event.hh>
#include <G4TrajectoryContainer.hh>
#include <G4Trajectory.hh>
#include <G4VVisManager.hh>
#include <G4UImanager.hh>
#include <G4LogicalVolumeStore.hh>
#include <G4VSensitiveDetector.hh>

AliEventAction::AliEventAction()
  : fVerboseLevel(1), 
    fDrawFlag("CHARGED")
{
//
  fMessenger = new AliEventActionMessenger(this);
}

AliEventAction::AliEventAction(const AliEventAction& right) {
//
  AliGlobals::Exception("AliEventAction is protected from copying.");
}

AliEventAction::~AliEventAction() {
//
  delete fMessenger;
}

// operators

AliEventAction& AliEventAction::operator=(const AliEventAction &right)
{
  // check assignement to self
  if (this == &right) return *this;
  
  AliGlobals::Exception("AliEventAction is protected from assigning.");

  return *this;
}

// private methods

void AliEventAction::DisplayEvent(const G4Event* event) const
{
// Draws trajectories.
// ---


  // trajectories processing
  G4TrajectoryContainer* trajectoryContainer 
    = event->GetTrajectoryContainer();

  G4int nofTrajectories = 0;
  if (trajectoryContainer)
  { nofTrajectories = trajectoryContainer->entries(); }
  
  if (fVerboseLevel>0) {
    G4cout << "    " << nofTrajectories; 
    G4cout << " trajectories stored." << G4endl;
  }  

  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager && nofTrajectories>0)
  {
    G4UImanager::GetUIpointer()->ApplyCommand("/vis~/draw/current");

    for (G4int i=0; i<nofTrajectories; i++)
    { 
      G4VTrajectory* vtrajectory = (*(event->GetTrajectoryContainer()))[i];
      G4Trajectory* trajectory = dynamic_cast<G4Trajectory*>(vtrajectory);
      if (!trajectory) {
        AliGlobals::Exception(
	  "AliEventAction::DisplayEvent: Unknown trajectory type.");
      }
      if ( (fDrawFlag == "ALL") ||
          ((fDrawFlag == "CHARGED") && (trajectory->GetCharge() != 0.))){
	 trajectory->DrawTrajectory(50); 
	    // the argument number defines the size of the step points
	    // use 2000 to make step points well visible
      }	
    }      
    G4UImanager::GetUIpointer()->ApplyCommand("/vis~/show/view");
  }  
}

// public methods

void AliEventAction::BeginOfEventAction(const G4Event* event)
{
// Called by G4 kernel at the beginning of event.
// ---

  G4int eventID = event->GetEventID();

  // reset the counters (primary tracks, saved tracks)
  AliTrackingAction* trackingAction
    = AliTrackingAction::Instance();
  trackingAction->PrepareNewEvent();   

  if (fVerboseLevel>0)
    G4cout << ">>> Event " << event->GetEventID() << G4endl;
}

void AliEventAction::EndOfEventAction(const G4Event* event)
{
// Called by G4 kernel at the end of event.
// ---

  // save the last primary track store of
  // the current event
  AliTrackingAction* trackingAction
    = AliTrackingAction::Instance();
  trackingAction->SaveAndDestroyTrack();   

  if (fVerboseLevel>0) {
    G4int nofPrimaryTracks = trackingAction->GetNofPrimaryTracks();
    G4cout  << "    " << nofPrimaryTracks << 
               " primary tracks processed." << G4endl;
  }	       

  // display event
  DisplayEvent(event);

  // aliroot
  // store event header data
  gAlice->GetHeader()->SetEvent(event->GetEventID());
  gAlice->GetHeader()->SetNvertex(event->GetNumberOfPrimaryVertex());
  gAlice->GetHeader()->SetNprimary(trackingAction->GetNofPrimaryTracks());
  gAlice->GetHeader()->SetNtrack(trackingAction->GetNofSavedTracks());

  gAlice->FinishEvent();    
}



