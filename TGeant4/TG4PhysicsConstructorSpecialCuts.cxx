// $Id$
// Category: physics
//
// Author: I. Hrivnacova
//
// Class TG4PhysicsConstructorSpecialCuts
// --------------------------------------
// See the class description in the header file.

#include "TG4PhysicsConstructorSpecialCuts.h"
#include "TG4G3PhysicsManager.h"
#include "TG4SpecialCutsForGamma.h"
#include "TG4SpecialCutsForElectron.h"
#include "TG4SpecialCutsForEplus.h"
#include "TG4SpecialCutsForChargedHadron.h"
#include "TG4SpecialCutsForNeutralHadron.h"
#include "TG4SpecialCutsForMuon.h"
#include "TG4SpecialCutsForOther.h"

#include <G4ParticleDefinition.hh>
#include <G4ProcessManager.hh>
#include <G4VProcess.hh>


//_____________________________________________________________________________
TG4PhysicsConstructorSpecialCuts::TG4PhysicsConstructorSpecialCuts(
                                     const G4String& name)
  : G4VPhysicsConstructor(name)
{
//
  SetVerboseLevel(1);
}

//_____________________________________________________________________________
TG4PhysicsConstructorSpecialCuts::~TG4PhysicsConstructorSpecialCuts() {
//
}

// protected methods

//_____________________________________________________________________________
void TG4PhysicsConstructorSpecialCuts::ConstructParticle()
{
// The particles are constructed in the 
// TG4ModularPhysicsList.
// ---
}

//_____________________________________________________________________________
void TG4PhysicsConstructorSpecialCuts::ConstructProcess()
{
// Adds TG4SpecialCuts "process" that activates
// the kinetic energy cuts defined in 
// the vector of cuts (PhysicsManager::fCutVector) or in TG4Limits.
// ---

  TG4G3PhysicsManager* g3PhysicsManager = TG4G3PhysicsManager::Instance();

  if (g3PhysicsManager->IsSpecialCuts())
  {
    TG4boolVector* isCutVector = g3PhysicsManager->GetIsCutVector(); 

    theParticleIterator->reset();
    while ((*theParticleIterator)())
    {
      G4ParticleDefinition* particle = theParticleIterator->value();
      TG4G3ParticleWSP particleWSP 
        = g3PhysicsManager->GetG3ParticleWSP(particle);
      G4String name =
        g3PhysicsManager->GetG3ParticleWSPName(particleWSP);
      
      // uncomment this to see all particles "WSP"
      //G4cout << "Iterating particle: " 
      //       << particle->GetParticleName() << " " << particleWSP << " "
      //       << name << G4endl;

      if ((particleWSP !=kNofParticlesWSP) && 
          ((*isCutVector)[particleWSP])) {
        // special process is created in case
        // cutVector (vector of kinetic energy cuts) is set
        // or the special cut is set by TG4Limits
  
        // check if process already exists
        G4String processName = "specialCutFor" + name;
	G4VProcess* process = g3PhysicsManager->FindProcess(processName);
	if (!process) {
          switch (particleWSP) {
            case kGamma:
              process = new TG4SpecialCutsForGamma(processName);
	      break;
            case kElectron:
              process = new TG4SpecialCutsForElectron(processName);
	      break;
            case kEplus:  
              process = new TG4SpecialCutsForEplus(processName);
	      break;
            case kChargedHadron:  
              process = new TG4SpecialCutsForChargedHadron(processName);
	      break;
            case kNeutralHadron:  
              process = new TG4SpecialCutsForNeutralHadron(processName);
	      break;
            case kMuon:  
              process = new TG4SpecialCutsForMuon(processName);
	      break;
            case kAny:
              process = new TG4SpecialCutsForOther(processName);
	      break;
	  }  
        }
        // add process to particle
        particle->GetProcessManager()->AddDiscreteProcess(process);
      }
    }

    if (verboseLevel>0) {
      G4cout << "###  Special Cuts constructed. " << G4endl;
      G4cout << "     Special cuts process is defined for: " << G4endl 
             << "     ";
      for (G4int i=0; i<kAny; i++) {
        if ((*isCutVector)[i]) 
	  G4cout << g3PhysicsManager->GetG3ParticleWSPName(i) << " ";
      }  
      G4cout << G4endl;
    }  
  }
}


