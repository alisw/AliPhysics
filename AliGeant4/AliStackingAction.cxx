// $Id$
// Category: event
//
// See the class description in the header file.

#include "AliStackingAction.h"
#include "AliStackingActionMessenger.h"
#include "AliTrackingAction.h"
#include "AliGlobals.h"

#include <G4Track.hh>
#include <G4TrackStack.hh>
#include <G4StackedTrack.hh>
#include <G4StackManager.hh>
#include <G4NeutrinoE.hh>
#include <G4NeutrinoMu.hh>
#include <G4NeutrinoTau.hh>
#include <G4AntiNeutrinoE.hh>
#include <G4AntiNeutrinoMu.hh>
#include <G4AntiNeutrinoTau.hh>

AliStackingAction::AliStackingAction()
  : fStage(0), 
    fVerboseLevel(0),
    fSavePrimaries(true),
    fTrackingAction(0) 
{
// 
  fPrimaryStack = new G4TrackStack();
  fMessenger = new AliStackingActionMessenger(this);
}

AliStackingAction::AliStackingAction(const AliStackingAction& right) {
//
  AliGlobals::Exception("AliStackingAction is protected from copying.");
}

AliStackingAction::~AliStackingAction() {
// 
  delete fPrimaryStack;
  delete fMessenger; 
}

// operators

AliStackingAction& 
AliStackingAction::operator=(const AliStackingAction &right)
{
  // check assignement to self
  if (this == &right) return *this;
  
  AliGlobals::Exception("AliStackingAction is protected from assigning.");

  return *this;
}

// public methods

G4ClassificationOfNewTrack 
AliStackingAction::ClassifyNewTrack(const G4Track* track)
{
// Classifies the new track.
// ---

  G4ClassificationOfNewTrack classification;
  if (fStage == 0) { 
    // move all primaries to PrimaryStack
    G4Track* nonconstTrack = (G4Track*)track;
    G4StackedTrack* newTrack = new G4StackedTrack(nonconstTrack);
    fPrimaryStack->PushToStack(newTrack);  
    classification = fPostpone;

    // save primary particle info
    // (secondary particles are stored 
    //  by AlTrackingAction::PreUserTrackingAction() method)
    if (fSavePrimaries)
      fTrackingAction->SaveTrack(track);
  }  
  else {
     // exclude neutrinos
    G4ParticleDefinition* particle = track->GetDefinition();
    if( particle == G4NeutrinoE::NeutrinoEDefinition() ||
        particle == G4NeutrinoMu::NeutrinoMuDefinition() ||
        particle == G4NeutrinoTau::NeutrinoTauDefinition() ||
        particle == G4AntiNeutrinoE::AntiNeutrinoEDefinition() ||
        particle == G4AntiNeutrinoMu::AntiNeutrinoMuDefinition() ||
        particle == G4AntiNeutrinoTau::AntiNeutrinoTauDefinition()) {

        return fKill;	 
     }	
     
     G4int parentID = track->GetParentID();
     if (parentID ==0) { 
        return fUrgent; 
     }
     else { 
        return fWaiting; 
     }
  }
  return classification;
}

void AliStackingAction::NewStage()
{
// Called by G4 kernel at the new stage of stacking.
// ---

  fStage++;
  if (fVerboseLevel>0) 
  {
    G4cout << "AliStackingAction::NewStage " << fStage 
           << " has been started." << G4endl;
  }

  G4int nofUrgent = stackManager->GetNUrgentTrack();
  if (nofUrgent == 0)
  {
    G4int nofPrimary = fPrimaryStack->GetNTrack();
    if (nofPrimary>0)
    { 
       G4StackedTrack* stackedTrack
         = fPrimaryStack->PopFromStack();
       G4Track* primaryTrack
         = stackedTrack->GetTrack();
       delete stackedTrack;
       stackManager->PushOneTrack(primaryTrack);
     }
  }
}
    
void AliStackingAction::ClearPrimaryStack()
{
// Clears the primary stack.
// ---

  stackManager->ClearPostponeStack();
}

void AliStackingAction::PrepareNewEvent()
{
// Called by G4 kernel at the beginning of event.
// ---

  fStage = 0;
  ClearPrimaryStack();
  fTrackingAction = AliTrackingAction::Instance();
  fSavePrimaries = fTrackingAction->GetSavePrimaries();
}


