// $Id$
// Category: physics
//
// Author: I. Hrivnacova
//
// Class TG4PhysicsConstructorSpecialControls
// ------------------------------------------
// See the class description in the header file.

#include "TG4PhysicsConstructorSpecialControls.h"
#include "TG4SpecialControls.h"
#include "TG4G3PhysicsManager.h"
#include "TG4G3ControlVector.h"

#include <G4ParticleDefinition.hh>
#include <G4ProcessManager.hh>
#include <G4VProcess.hh>

//_____________________________________________________________________________
TG4PhysicsConstructorSpecialControls::TG4PhysicsConstructorSpecialControls(
                                     const G4String& name)
  : TG4VPhysicsConstructor(name) {
//
}

//_____________________________________________________________________________
TG4PhysicsConstructorSpecialControls::TG4PhysicsConstructorSpecialControls(
				     G4int verboseLevel, 
                                     const G4String& name)
  : TG4VPhysicsConstructor(name, verboseLevel) {
//
}

//_____________________________________________________________________________
TG4PhysicsConstructorSpecialControls::~TG4PhysicsConstructorSpecialControls() {
//
}

// protected methods

//_____________________________________________________________________________
void TG4PhysicsConstructorSpecialControls::ConstructParticle()
{
// The particles are constructed in the 
// TG4ModularPhysicsList.
// ---
}

//_____________________________________________________________________________
void TG4PhysicsConstructorSpecialControls::ConstructProcess()
{
// Adds TG4SpecialControls "process" that activates
// the control process controls defined in TG4Limits.
// ---

  TG4G3PhysicsManager* g3PhysicsManager 
    = TG4G3PhysicsManager::Instance();

  if (g3PhysicsManager->IsSpecialControls()) {
    TG4boolVector* isControlVector 
      = g3PhysicsManager->GetIsControlVector(); 

    theParticleIterator->reset();
    while ((*theParticleIterator)())
    {
      G4ParticleDefinition* particle = theParticleIterator->value();
      TG4G3ParticleWSP particleWSP 
        = g3PhysicsManager->GetG3ParticleWSP(particle);
      //G4String name;
      //GetG3ParticleWSPName(particleWSP, name);

      // special process is set in case
      // the special control is set by TG4Limits    
      if ((particleWSP !=kNofParticlesWSP) && 
          ((*isControlVector)[particleWSP])) {
        // check if process already exists
	G4String processName = "specialControl";
	G4VProcess* process = g3PhysicsManager->FindProcess(processName);
	if (!process) {
          process = new TG4SpecialControls(processName);
	}  
        //particle->GetProcessManager()->AddProcess(process, 0, -1, 1);
        particle->GetProcessManager()->AddDiscreteProcess(process);
      }
    }

    if (VerboseLevel() > 0) {
      G4cout << "### Special Controls constructed. " << G4endl;
      G4cout << "    Special controls process is defined for: " << G4endl
             << "    ";
      for (G4int i=0; i<kNofParticlesWSP; i++) {
        if ((*isControlVector)[i]) 
	  G4cout << g3PhysicsManager->GetG3ParticleWSPName(i) << " ";
      }  
      G4cout << G4endl;
    }  
  }
}

