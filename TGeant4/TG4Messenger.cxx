// $Id$
// Category: run
//
// See the class description in the header file.

#include "TG4Messenger.h"
#include "TG4GeometryManager.h"
#include "TG4StepManager.h"
#include "TG4PhysicsManager.h"

#include <G4UIcmdWithoutParameter.hh>

TG4Messenger::TG4Messenger(TG4GeometryManager* geometryManager, 
       TG4PhysicsManager* physicsManager, TG4StepManager* stepManager)
  : fGeometryManager(geometryManager),
    fPhysicsManager(physicsManager),
    fStepManager(stepManager)
{ 
//
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

  if (command == fProcessActivationCmd) {
    fPhysicsManager->SetProcessActivation();
  }  
}
