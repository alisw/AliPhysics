// $Id$ //
// Category: physics
//
// Author: I. Hrivnacova
//
// Class TG4VSpecialControls
// -------------------------
// See the class description in the header file.

#include "TG4SpecialControls.h"
#include "TG4GeometryServices.h"
#include "TG4Limits.h"

#include <G4StepStatus.hh>
#include <G4ProcessManager.hh>
#include <G4ProcessVector.hh>

//_____________________________________________________________________________
TG4SpecialControls::TG4SpecialControls(const G4String& processName)
  : G4VProcess(processName),
    TG4Verbose("specialControls"),
    fSwitchControls(kUnswitch),
    fLastTrackID(0) {
//    
   verboseLevel = VerboseLevel();
   if (VerboseLevel() >0 ) {
     G4cout << GetProcessName() << " is created "<< G4endl;
   }
}

//_____________________________________________________________________________
TG4SpecialControls::TG4SpecialControls(const TG4SpecialControls& right)     
  : TG4Verbose("specialControls") {
// 
  TG4Globals::Exception(
    "TG4SpecialControls is protected from copying.");
}

//_____________________________________________________________________________
TG4SpecialControls::~TG4SpecialControls() {
//
}

// operators

//_____________________________________________________________________________
TG4SpecialControls& TG4SpecialControls::operator=(
                                          const TG4SpecialControls& right)
{
  // check assignement to self
  if (this == &right) return *this;

  TG4Globals::Exception(
    "TG4SpecialControls is protected from assigning.");
    
  return *this;  
} 

// private methods   

//_____________________________________________________________________________
void TG4SpecialControls::Reset()
{
// Resets the buffers to the initial state.
// ---
    			
  fSwitchControls = kUnswitch;

  // clear buffers
  fSwitchedProcesses.clear();
  fSwitchedControls.clear();
}

// public methods   
          
//_____________________________________________________________________________
G4double TG4SpecialControls::PostStepGetPhysicalInteractionLength(
                           const G4Track& track, G4double previousStepSize,
			   G4ForceCondition* condition)
{
// Returns the Step-size (actual length) which is allowed 
// by this process.
// ---

  *condition = NotForced;

  if (track.GetTrackID() != fLastTrackID) {
    // new track
    Reset();
    fLastTrackID = track.GetTrackID();
  }  

  G4double proposedStep = DBL_MAX;
  //G4double minStep = (1.0e-9)*m;
  G4double minStep = 0.;
    // must be greater than DBL_MIN - so that particle can get out of
    // the boundary 
    // proposedStep = 0.; causes navigator to fall into panic 
  
  G4StepStatus status     
    = track.GetStep()->GetPreStepPoint()->GetStepStatus();

  // get limits
#ifdef TGEANT4_DEBUG
  TG4Limits* limits 
     = TG4GeometryServices::Instance()
         ->GetLimits(track.GetVolume()->GetLogicalVolume()->GetUserLimits()); 

  if (!limits) {
    G4String text = "TG4VSpecialControls::PostStepGetPhysicalInteractionLength:\n";
    text = text + "    " + track.GetVolume()->GetLogicalVolume()->GetName();
    text = text + " has not limits.";
    TG4Globals::Exception(text);
  }  
#else  
  TG4Limits* limits 
    = (TG4Limits*) track.GetVolume()->GetLogicalVolume()->GetUserLimits();
#endif    

  if (fSwitchControls != kUnswitch) {
    if (status == fGeomBoundary) {
      if  (limits->IsControl()) {
        // particle is exiting a logical volume with special controls
        // and entering another logical volume with special controls 
	proposedStep = minStep;
        fSwitchControls = kReswitch;
        if (VerboseLevel() > 1) { 
	  G4cout << "kReswitch" << G4endl;
	}  
      }
      else {
        // particle is exiting a logical volume with special controls
        // and entering a logical volume without special controls 
	proposedStep = minStep;
        fSwitchControls = kUnswitch;
        if (VerboseLevel() > 1) { 
	  G4cout << "kUnswitch" << G4endl;
	}  
      }
    }
  }
  else if (limits->IsControl()) {
       // particle is entering a logical volume with special controls
       // that have not yet been set
       proposedStep = minStep;
       fSwitchControls = kSwitch;
       if (VerboseLevel() > 1) { 
         G4cout << "kSwitch" << G4endl;
       }	 
  }  
  return proposedStep;
}

//_____________________________________________________________________________
G4VParticleChange* TG4SpecialControls::PostStepDoIt(
                      const G4Track& track, const G4Step& step)
{
// Changes processes activation of the current track
// according to the current user limits.
// ---

  G4ProcessManager* processManager
    = track.GetDefinition()->GetProcessManager();
  G4ProcessVector* processVector = processManager->GetProcessList();

  if ((fSwitchControls==kUnswitch) || (fSwitchControls==kReswitch)) {
  
    // set processes activation back
    for (G4int i=0; i<fSwitchedProcesses.length(); i++) {
      if (VerboseLevel() > 1) {
        G4cout << "Reset process activation back in " 
  	       << track.GetVolume()->GetName() 
               << G4endl;
      }
      processManager
        ->SetProcessActivation(fSwitchedProcesses[i],fSwitchedControls[i]);
    }
    fSwitchedProcesses.clear();
    fSwitchedControls.clear();
  }

  if ((fSwitchControls==kSwitch) ||  (fSwitchControls==kReswitch)) {

    // set TG4Limits processes controls
    TG4Limits* limits 
    = (TG4Limits*) track.GetVolume()->GetLogicalVolume()->GetUserLimits();

    for (G4int i=0; i<processVector->length(); i++) {

      TG4G3ControlValue control = limits->GetControl((*processVector)[i]);
      G4bool activation = processManager->GetProcessActivation(i);

      if (control != kUnset && ! TG4Globals::Compare(activation, control)) {

        // store the current processes controls
        if (VerboseLevel() > 1) {
	  G4cout << "Something goes to fSwitchedProcesses" << G4endl;
	}  
        fSwitchedProcesses.insert((*processVector)[i]);
        fSwitchedControls.push_back(activation);

        // set new process activation
        if (control == kInActivate) {
          if (VerboseLevel() > 1) {
            G4cout << "Set process inactivation for " 
                   << (*processVector)[i]->GetProcessName() << " in " 
    	           << track.GetVolume()->GetName() 
	           << G4endl;
          }
          processManager->SetProcessActivation(i,false);
        }  
        else {
	  // ((control == kActivate) || (control == kActivate2)) 
          if (VerboseLevel() > 1) {
            G4cout << "Set process activation for " 
                   << (*processVector)[i]->GetProcessName() << " in " 
	           << track.GetVolume()->GetName() 
	           << G4endl;
          }
          processManager->SetProcessActivation(i,true);
        }
      }	 
    }
  }
    
  // processManager->DumpInfo();      	
  aParticleChange.Initialize(track);
  aParticleChange.SetStatusChange(fAlive);
  return &aParticleChange;
}

