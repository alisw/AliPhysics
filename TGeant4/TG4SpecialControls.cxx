// $Id$ //
// Category: physics
//
// See the class description in the header file.

#include "TG4SpecialControls.h"
#include "TG4Limits.h"

#include <G4StepStatus.hh>
#include <G4ProcessManager.hh>
#include <G4ProcessVector.hh>

TG4SpecialControls::TG4SpecialControls(const G4String& aName)
  : G4VProcess(aName),
    fSwitchControls(kUnswitch)
{
   // verboseLevel = 1;
   if (verboseLevel>0) {
     G4cout << GetProcessName() << " is created "<< G4endl;
   }
}

TG4SpecialControls::TG4SpecialControls(const TG4SpecialControls& right) {
// 
  TG4Globals::Exception(
    "TG4SpecialControls is protected from copying.");
}

TG4SpecialControls::~TG4SpecialControls() {
//
}

// operators

TG4SpecialControls& TG4SpecialControls::operator=(
                                          const TG4SpecialControls& right)
{
  // check assignement to self
  if (this == &right) return *this;

  TG4Globals::Exception(
    "TG4SpecialControls is protected from assigning.");
    
  return *this;  
} 

// public methods   
          
G4double TG4SpecialControls::PostStepGetPhysicalInteractionLength(
                           const G4Track& track, G4double previousStepSize,
			   G4ForceCondition* condition)
{
// Returns the Step-size (actual length) which is allowed 
// by this process.
// ---

  *condition = NotForced;

  G4double proposedStep = DBL_MAX;
  G4double minStep = (1.0e-9)*m;
    // must be greater than DBL_MIN - so that particle can get out of
    // the boundary 
    // proposedStep = 0.; causes navigator to fall into panic 
  
  G4StepStatus status     
    = track.GetStep()->GetPreStepPoint()->GetStepStatus();
  TG4Limits* limits 
    = (TG4Limits*) track.GetVolume()->GetLogicalVolume()->GetUserLimits();

  if (fSwitchControls != kUnswitch) {
    if (status == fGeomBoundary) {
      if  ((limits) && (limits->IsControl())) {
        // particle is exiting a logical volume with special controls
        // and entering another logical volume with special controls 
	proposedStep = minStep;
        fSwitchControls = kReswitch;
        if (verboseLevel>0) G4cout << "kReswitch" << G4endl;
      }
      else {
        // particle is exiting a logical volume with special controls
        // and entering a logical volume without special controls 
	proposedStep = minStep;
        fSwitchControls = kUnswitch;
        if (verboseLevel>0) G4cout << "kUnswitch" << G4endl;
      }
    }
  }
  else if ((limits) && (limits->IsControl())) {
       // particle is entering a logical volume with special controls
       // that have not yet been set
       proposedStep = minStep;
       fSwitchControls = kSwitch;
       if (verboseLevel>0) G4cout << "kSwitch" << G4endl;
  }  
  return proposedStep;
}

G4VParticleChange* TG4SpecialControls::PostStepDoIt(
                      const G4Track& track, const G4Step& step)
{
// Changes processes activation of the current track
// according to the current user limits.
// ---

  TG4Limits* limits 
    = (TG4Limits*) track.GetVolume()->GetLogicalVolume()->GetUserLimits();

  G4ProcessManager* processManager
    = track.GetDefinition()->GetProcessManager();
  G4ProcessVector* processVector = processManager->GetProcessList();
  // processManager->DumpInfo();

  if ((fSwitchControls==kUnswitch) || (fSwitchControls==kReswitch)) {
    // set processes activation back
    for (G4int i=0; i<fSwitchedProcesses.entries(); i++) {
      if (verboseLevel>0) {
        G4cout << "Reset process activation back in" 
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
    for (G4int i=0; i<processManager->GetProcessListLength(); i++) {
      G4int control = limits->GetControl((*processVector)[i]);
      if (control != kUnset) {
        // store the current processes controls;
        fSwitchedProcesses.insert((*processVector)[i]);
        //fSwitchedControls.insert(processManager->GetProcessActivation(i));
        fSwitchedControls.push_back(processManager->GetProcessActivation(i));
        if (control == kInActivate) {
          if (verboseLevel>0) {
            G4cout << "Set process inactivation for " 
                   << (*processVector)[i]->GetProcessName() << " in " 
    	           << track.GetVolume()->GetName() 
	           << G4endl;
          }
          processManager->SetProcessActivation(i,false);
        }  
        else {
	  // ((control == kActivate) || (control == kActivate2)) 
          if (verboseLevel>0) {
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

