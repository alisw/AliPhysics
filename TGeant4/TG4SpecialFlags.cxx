// $Id$ //
// Category: physics
//
// See the class description in the header file.

#include "TG4SpecialFlags.h"
#include "TG4Limits.h"

#include <G4StepStatus.hh>
#include <G4ProcessManager.hh>
#include <G4ProcessVector.hh>

TG4SpecialFlags::TG4SpecialFlags(const G4String& aName)
  : G4VProcess(aName),
    fSwitchFlags(kUnswitch)
{
   // verboseLevel = 1;
   if (verboseLevel>0) {
     G4cout << GetProcessName() << " is created "<< G4endl;
   }
}

TG4SpecialFlags::TG4SpecialFlags(const TG4SpecialFlags& right) {
// 
  TG4Globals::Exception(
    "TG4SpecialFlags is protected from copying.");
}

TG4SpecialFlags::~TG4SpecialFlags() {
//
}

// operators

TG4SpecialFlags& TG4SpecialFlags::operator=(const TG4SpecialFlags& right)
{
  // check assignement to self
  if (this == &right) return *this;

  TG4Globals::Exception(
    "TG4SpecialFlags is protected from assigning.");
    
  return *this;  
} 

// public methods   
          
G4double TG4SpecialFlags::PostStepGetPhysicalInteractionLength(
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

  if (fSwitchFlags != kUnswitch) {
    if (status == fGeomBoundary) {
      if  ((limits) && (limits->IsFlag())) {
        // particle is exiting a logical volume with special flags
        // and entering another logical volume with special flags 
	proposedStep = minStep;
        fSwitchFlags = kReswitch;
        if (verboseLevel>0) G4cout << "kReswitch" << G4endl;
      }
      else {
        // particle is exiting a logical volume with special flags
        // and entering a logical volume without special flags 
	proposedStep = minStep;
        fSwitchFlags = kUnswitch;
        if (verboseLevel>0) G4cout << "kUnswitch" << G4endl;
      }
    }
  }
  else if ((limits) && (limits->IsFlag())) {
       // particle is entering a logical volume with special flags
       // that have not yet been set
       proposedStep = minStep;
       fSwitchFlags = kSwitch;
       if (verboseLevel>0) G4cout << "kSwitch" << G4endl;
  }  
  return proposedStep;
}

G4VParticleChange* TG4SpecialFlags::PostStepDoIt(
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

  if ((fSwitchFlags==kUnswitch) || (fSwitchFlags==kReswitch)) {
    // set processes activation back
    for (G4int i=0; i<fSwitchedProcesses.entries(); i++) {
      if (verboseLevel>0) {
        G4cout << "Reset process activation back in" 
  	       << track.GetVolume()->GetName() 
               << G4endl;
      }
      processManager
        ->SetProcessActivation(fSwitchedProcesses[i],fSwitchedFlags[i]);
    }
    fSwitchedProcesses.clear();
    fSwitchedFlags.clear();
  }

  if ((fSwitchFlags==kSwitch) ||  (fSwitchFlags==kReswitch)) {
    // set TG4Limits processes flags
    for (G4int i=0; i<processManager->GetProcessListLength(); i++) {
      G4int flag = limits->GetFlag((*processVector)[i]);
      if (flag != kUnset) {
        // store the current processes flags;
        fSwitchedProcesses.insert((*processVector)[i]);
        //fSwitchedFlags.insert(processManager->GetProcessActivation(i));
        fSwitchedFlags.push_back(processManager->GetProcessActivation(i));
        if (flag == kInActivate) {
          if (verboseLevel>0) {
            G4cout << "Set process inactivation for " 
                   << (*processVector)[i]->GetProcessName() << " in " 
    	           << track.GetVolume()->GetName() 
	           << G4endl;
          }
          processManager->SetProcessActivation(i,false);
        }  
        else {
	  // ((flag == kActivate) || (flag == kActivate2)) 
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

