// $Id$
// Category: physics
//
// See the class description in the header file.

#include "TG4PhysicsConstructorSpecialCuts.h"
#include "TG4G3PhysicsManager.h"
#include "TG4SpecialCuts.h"

#include <G4ParticleDefinition.hh>
#include <G4ProcessManager.hh>
#include <G4VProcess.hh>


TG4PhysicsConstructorSpecialCuts::TG4PhysicsConstructorSpecialCuts(
                                     const G4String& name)
  : G4VPhysicsConstructor(name)
{
//
  SetVerboseLevel(1);
}

TG4PhysicsConstructorSpecialCuts::~TG4PhysicsConstructorSpecialCuts() {
//
}

// protected methods

void TG4PhysicsConstructorSpecialCuts::ConstructParticle()
{
// The particles are constructed in the 
// TG4ModularPhysicsList.
// ---
}

void TG4PhysicsConstructorSpecialCuts::ConstructProcess()
{
// Adds TG4SpecialCuts "process" that activates
// the kinetic energy cuts defined in 
// the vector of cuts (PhysicsManager::fCutVector) or in TG4Limits.
// ---

  TG4G3PhysicsManager* g3PhysicsManager 
    = TG4G3PhysicsManager::Instance();

  if (g3PhysicsManager->IsSpecialCuts())
  {
    TG4G3CutVector* cutVector
      = g3PhysicsManager->GetCutVector(); 
    TG4boolVector* isCutVector 
      = g3PhysicsManager->GetIsCutVector(); 

    theParticleIterator->reset();
    while ((*theParticleIterator)())
    {
      G4ParticleDefinition* particle = theParticleIterator->value();
      TG4G3ParticleWSP particleWSP 
        = g3PhysicsManager->GetG3ParticleWSP(particle);
      G4String name;
      g3PhysicsManager->GetG3ParticleWSPName(particleWSP, name);
      
      // uncomment this to see all particles "WSP"
      //G4cout << "Iterating particle: " 
      //       << particle->GetParticleName() << " " << particleWSP << " "
      //       << name << G4endl;

      // special process is created in case
      // cutVector (vector of kinetic energy cuts) is set
      // or the special cut is set by TG4Limits
      if ((particleWSP !=kNofParticlesWSP) && 
          ((*isCutVector)[particleWSP])) {
        // check if process already exists
	G4String processName = "specialCutFor" + name;
	G4VProcess* process = g3PhysicsManager->FindProcess(processName);
	if (!process) {
          process = new TG4SpecialCuts(particleWSP, cutVector, processName);
	}  
        //particle->GetProcessManager()->AddProcess(process, 0, -1, 1);
        particle->GetProcessManager()->AddDiscreteProcess(process);
      }
    }

    if (verboseLevel>0) {
      G4cout << "TG4PhysicsList::ConstructSpecialCuts: " << G4endl;
      if (cutVector)
        G4cout << "   Global kinetic energy cuts are set." << G4endl;
      G4cout << "   Special cuts process is defined for: " << G4endl 
             << "   ";
      for (G4int i=0; i<kAny; i++) {
        G4String name;
        g3PhysicsManager->GetG3ParticleWSPName(i, name);
        if ((*isCutVector)[i]) G4cout << name << " ";
      }  
      G4cout << G4endl;
    }  
  }
}


