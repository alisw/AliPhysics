// $Id$
// Category: physics
//
// See the class description in the header file.

#include "TG4PhysicsConstructorSpecialControls.h"
#include "TG4SpecialControls.h"
#include "TG4G3PhysicsManager.h"
#include "TG4G3ControlVector.h"

#include <G4ParticleDefinition.hh>
#include <G4ProcessManager.hh>
#include <G4VProcess.hh>

TG4PhysicsConstructorSpecialControls::TG4PhysicsConstructorSpecialControls(
                                     const G4String& name)
  : G4VPhysicsConstructor(name)
{
//
  SetVerboseLevel(1);
}

TG4PhysicsConstructorSpecialControls::~TG4PhysicsConstructorSpecialControls() {
//
}

// protected methods

void TG4PhysicsConstructorSpecialControls::ConstructParticle()
{
// The particles are constructed in the 
// TG4ModularPhysicsList.
// ---
}

void TG4PhysicsConstructorSpecialControls::ConstructProcess()
{
// Adds TG4SpecialControls "process" that activates
// the control process controls defined in TG4Limits.
// ---

  TG4G3PhysicsManager* g3PhysicsManager 
    = TG4G3PhysicsManager::Instance();

  if (g3PhysicsManager->IsSpecialControls())
  {
    G4cout << "IsSpecialControls started" << G4endl;
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

    if (verboseLevel>0) {
      G4cout << "TG4PhysicsConstructorSpecialControls::ConstructProcess: " << G4endl;
      G4cout << "   Special controls process is defined for: " << G4endl
             << "   ";
      for (G4int i=0; i<kNofParticlesWSP; i++) {
        G4String name;
        g3PhysicsManager->GetG3ParticleWSPName(i, name);
        if ((*isControlVector)[i]) G4cout << name << " ";
      }  
      G4cout << G4endl;
    }  
  }
}

