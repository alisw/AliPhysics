// $Id$
// Category: run
//
// See the class description in the header file.

#include "TG4Messenger.h"
#include "TG4GeometryManager.h"
#include "TG4StepManager.h"
#include "TG4PhysicsManager.h"

#include <G4UIcmdWithoutParameter.hh>
#include <G4UIcmdWithABool.hh>

TG4Messenger::TG4Messenger(TG4GeometryManager* geometryManager, 
                           TG4PhysicsManager* physicsManager, 
			   TG4StepManager* stepManager)
  : fGeometryManager(geometryManager),
    fPhysicsManager(physicsManager),
    fStepManager(stepManager)
{ 
//
  fSetEMCmd
     = new G4UIcmdWithABool("/g4mc/setEM", this);
  fSetEMCmd->SetGuidance("Set electromagnetic processes.");
  fSetEMCmd->SetParameterName("EMControl", false);
  fSetEMCmd->AvailableForStates(PreInit);

  fSetOpticalCmd
     = new G4UIcmdWithABool("/g4mc/setOptical", this);
  fSetOpticalCmd->SetGuidance("Set Cerenkov and optical processes.");
  fSetOpticalCmd->SetParameterName("OpticalControl", false);
  fSetOpticalCmd->AvailableForStates(PreInit);

  fSetHadronCmd
     = new G4UIcmdWithABool("/g4mc/setHadron", this);
  fSetHadronCmd->SetGuidance("Set hadron processes.");
  fSetHadronCmd->SetParameterName("HadronControl", false);
  fSetHadronCmd->AvailableForStates(PreInit);

  fSetSpecialCutsCmd
     = new G4UIcmdWithABool("/g4mc/setSpecialCuts", this);
  fSetSpecialCutsCmd->SetGuidance("Set special cuts process.");
  fSetSpecialCutsCmd
    ->SetGuidance("!! Support for this option is under development.");
  fSetSpecialCutsCmd->SetParameterName("SpecialCutsControl", false);
  fSetSpecialCutsCmd->AvailableForStates(PreInit);

  fSetSpecialControlsCmd
     = new G4UIcmdWithABool("/g4mc/setSpecialControls", this);
  fSetSpecialControlsCmd->SetGuidance("Set special controls process.");
  fSetSpecialControlsCmd
    ->SetGuidance("!! Support for this option is under development.");
  fSetSpecialControlsCmd->SetParameterName("SpecialFlagsControl", false);
  fSetSpecialControlsCmd->AvailableForStates(PreInit);

  fProcessActivationCmd
     = new G4UIcmdWithoutParameter("/g4mc/setProcessActivation", this);
  fProcessActivationCmd->SetGuidance("Activate/inactivate physics processes.");
  fProcessActivationCmd->AvailableForStates(Idle);
}

TG4Messenger::TG4Messenger(){
//
} 

TG4Messenger::TG4Messenger(const TG4Messenger& right) {
// 
  TG4Globals::Exception("TG4Messenger is protected from copying.");
}

TG4Messenger::~TG4Messenger() {
//
  delete fSetEMCmd;
  delete fSetOpticalCmd;
  delete fSetHadronCmd;
  delete fSetSpecialCutsCmd;
  delete fSetSpecialControlsCmd;
  delete fProcessActivationCmd;
}

// operators

TG4Messenger& TG4Messenger::operator=(const TG4Messenger& right)
{
  // check assignement to self
  if (this == &right) return *this;

  TG4Globals::Exception("TG4Messenger is protected from assigning.");
    
  return *this;  
}    
          
// public methods

void TG4Messenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
// Applies command to the associated object.
// ---

  if (command == fSetEMCmd) {
    fPhysicsManager
      ->SetEMPhysics(fSetOpticalCmd->GetNewBoolValue(newValue)); 
  }    
  else if (command == fSetOpticalCmd) {
    fPhysicsManager
      ->SetOpticalPhysics(fSetOpticalCmd->GetNewBoolValue(newValue)); 
  }    
  else if (command == fSetHadronCmd) {
    fPhysicsManager
      ->SetHadronPhysics(fSetHadronCmd->GetNewBoolValue(newValue)); 
  }    
  else if (command == fSetSpecialCutsCmd) {
    fPhysicsManager
      ->SetSpecialCutsPhysics(fSetSpecialCutsCmd->GetNewBoolValue(newValue)); 
  }    
  else if (command == fSetSpecialControlsCmd) {
    fPhysicsManager
      ->SetSpecialControlsPhysics(
          fSetSpecialControlsCmd->GetNewBoolValue(newValue)); 
  }    
  if (command == fProcessActivationCmd) {
    fPhysicsManager->SetProcessActivation();
  }  
}
