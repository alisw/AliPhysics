// $Id$
// Category: physics
//
// Author: I. Hrivnacova
//
// Class TG4PhysicsConstructorEM
// ----------0000---------------
// See the class description in the header file.
// According to corresponding part of:
// ExN04PhysicsList.cc,v 1.7 1999/12/15 14:49:26 gunter
// GEANT4 tag Name: geant4-01-01

#include "TG4PhysicsConstructorEM.h"
#include "TG4ProcessControlMap.h"
#include "TG4ProcessMCMap.h"
#include "TG4G3Control.h"

#include <G4ParticleDefinition.hh>
#include <G4ProcessManager.hh>
#include <G4ComptonScattering.hh>
#include <G4GammaConversion.hh>
#include <G4PhotoElectricEffect.hh>
#include <G4MultipleScattering.hh>
#include <G4eIonisation.hh>
#include <G4eBremsstrahlung.hh>
#include <G4eplusAnnihilation.hh>
#include <G4MuIonisation.hh>
#include <G4MuBremsstrahlung.hh>
#include <G4MuPairProduction.hh>
#include <G4hIonisation.hh>

//_____________________________________________________________________________
TG4PhysicsConstructorEM::TG4PhysicsConstructorEM(const G4String& name)
  : G4VPhysicsConstructor(name)
{
//
  SetVerboseLevel(1);
}

//_____________________________________________________________________________
TG4PhysicsConstructorEM::~TG4PhysicsConstructorEM() {
//
}

// protected methods

//_____________________________________________________________________________
void TG4PhysicsConstructorEM::ConstructParticle()
{
// The particles are constructed in the 
// TG4ModularPhysicsList.
// ---
}

//_____________________________________________________________________________
void TG4PhysicsConstructorEM::ConstructProcess()
{
// Constructs electromagnetic processes.
// ---

  TG4ProcessControlMap* controlMap = TG4ProcessControlMap::Instance();
  TG4ProcessMCMap* mcMap = TG4ProcessMCMap::Instance();

  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
     
    if (particleName == "gamma") {
      // construct processes for gamma
      G4VProcess* theGammaConversion     = new G4GammaConversion();
      G4VProcess* theComptonScattering   = new G4ComptonScattering();
      G4VProcess* thePhotoElectricEffect = new G4PhotoElectricEffect();
      
      // add processes
      pmanager->AddDiscreteProcess(theGammaConversion);
      pmanager->AddDiscreteProcess(theComptonScattering);      
      pmanager->AddDiscreteProcess(thePhotoElectricEffect);
      
      // map to G3 controls
      controlMap->Add(theGammaConversion, kPAIR); 
      controlMap->Add(theComptonScattering, kCOMP); 
      controlMap->Add(thePhotoElectricEffect, kPHOT); 

      // map to AliMCProcess codes
      mcMap->Add(theGammaConversion, kPPair); 
      mcMap->Add(theComptonScattering, kPCompton); 
      mcMap->Add(thePhotoElectricEffect, kPPhotoelectric); 
    } 
    else if (particleName == "e-") {
      // construct processes for electron
      G4VProcess* theeminusMultipleScattering = new G4MultipleScattering();
      G4VProcess* theeminusIonisation = new G4eIonisation();
      G4VProcess* theeminusBremsstrahlung = new G4eBremsstrahlung();
      
      // add processes
      pmanager->AddProcess(theeminusMultipleScattering);
      pmanager->AddProcess(theeminusIonisation);
      pmanager->AddProcess(theeminusBremsstrahlung);      

      // set ordering for AlongStepDoIt
      pmanager->SetProcessOrdering(theeminusMultipleScattering, idxAlongStep,  1);
      pmanager->SetProcessOrdering(theeminusIonisation, idxAlongStep,  2);

      // set ordering for PostStepDoIt
      pmanager->SetProcessOrdering(theeminusMultipleScattering, idxPostStep, 1);
      pmanager->SetProcessOrdering(theeminusIonisation, idxPostStep, 2);
      pmanager->SetProcessOrdering(theeminusBremsstrahlung, idxPostStep, 3);

      // map to G3 controls
      controlMap->Add(theeminusMultipleScattering, kMULS); 
      controlMap->Add(theeminusIonisation, kLOSS); 
      controlMap->Add(theeminusBremsstrahlung, kBREM); 

      // map to AliMCProcess codes
      mcMap->Add(theeminusMultipleScattering, kPMultipleScattering); 
      mcMap->Add(theeminusIonisation, kPEnergyLoss); 
      mcMap->Add(theeminusBremsstrahlung, kPBrem); 
    } 
    else if (particleName == "e+") {
      // Construct processes for positron
      G4VProcess* theeplusMultipleScattering = new G4MultipleScattering();
      G4VProcess* theeplusIonisation = new G4eIonisation();
      G4VProcess* theeplusBremsstrahlung = new G4eBremsstrahlung();
      G4VProcess* theeplusAnnihilation = new G4eplusAnnihilation();

      // add processes
      pmanager->AddProcess(theeplusMultipleScattering);
      pmanager->AddProcess(theeplusIonisation);
      pmanager->AddProcess(theeplusBremsstrahlung);
      pmanager->AddProcess(theeplusAnnihilation);

      // set ordering for AtRestDoIt
      pmanager->SetProcessOrderingToFirst(theeplusAnnihilation, idxAtRest);

      // set ordering for AlongStepDoIt
      pmanager->SetProcessOrdering(theeplusMultipleScattering, idxAlongStep,  1);
      pmanager->SetProcessOrdering(theeplusIonisation, idxAlongStep,  2);

      // set ordering for PostStepDoIt
      pmanager->SetProcessOrdering(theeplusMultipleScattering, idxPostStep, 1);
      pmanager->SetProcessOrdering(theeplusIonisation, idxPostStep, 2);
      pmanager->SetProcessOrdering(theeplusBremsstrahlung, idxPostStep, 3);
      pmanager->SetProcessOrdering(theeplusAnnihilation, idxPostStep, 4);
  
      // map to G3 controls
      controlMap->Add(theeplusMultipleScattering, kMULS); 
      controlMap->Add(theeplusIonisation, kLOSS); 
      controlMap->Add(theeplusBremsstrahlung, kBREM); 
      controlMap->Add(theeplusAnnihilation, kANNI); 
      
      // map to AliMCProcess codes
      mcMap->Add(theeplusMultipleScattering, kPMultipleScattering); 
      mcMap->Add(theeplusIonisation, kPEnergyLoss); 
      mcMap->Add(theeplusBremsstrahlung, kPBrem); 
      mcMap->Add(theeplusAnnihilation, kPAnnihilation); 
    } 
    else if (particleName == "mu+" || 
             particleName == "mu-") {
      // construct processes for muons
      G4VProcess* aMultipleScattering = new G4MultipleScattering();
      G4VProcess* aBremsstrahlung = new G4MuBremsstrahlung();
      G4VProcess* aPairProduction = new G4MuPairProduction();
      G4VProcess* anIonisation = new G4MuIonisation();

      // add processes
      pmanager->AddProcess(anIonisation);
      pmanager->AddProcess(aMultipleScattering);
      pmanager->AddProcess(aBremsstrahlung);
      pmanager->AddProcess(aPairProduction);

      // set ordering for AlongStepDoIt
      pmanager->SetProcessOrdering(aMultipleScattering, idxAlongStep,  1);
      pmanager->SetProcessOrdering(anIonisation, idxAlongStep,  2);

      // set ordering for PostStepDoIt
      pmanager->SetProcessOrdering(aMultipleScattering, idxPostStep, 1);
      pmanager->SetProcessOrdering(anIonisation, idxPostStep, 2);
      pmanager->SetProcessOrdering(aBremsstrahlung, idxPostStep, 3);
      pmanager->SetProcessOrdering(aPairProduction, idxPostStep, 4);
     
      // map to G3 controls
      controlMap->Add(aMultipleScattering, kMULS); 
      controlMap->Add(aBremsstrahlung, kBREM); 
      controlMap->Add(anIonisation, kLOSS); 
      controlMap->Add(aPairProduction, kPAIR); 

      // map to AliMCProcess codes
      mcMap->Add(aMultipleScattering, kPMultipleScattering); 
      mcMap->Add(aBremsstrahlung, kPBrem); 
      mcMap->Add(anIonisation, kPEnergyLoss); 
      mcMap->Add(aPairProduction, kPPair); 
    } 
    else if (particleName == "GenericIon") {
      // construct processes for ions
      G4VProcess* aionIonization = new G4hIonisation;
      G4VProcess* aMultipleScattering = new G4MultipleScattering();
      pmanager->AddProcess(aionIonization);
      pmanager->AddProcess(aMultipleScattering);

      // set ordering for AlongStepDoIt
      pmanager->SetProcessOrdering(aMultipleScattering, idxAlongStep,  1);
      pmanager->SetProcessOrdering(aionIonization, idxAlongStep,  2);

      // set ordering for PostStepDoIt
      pmanager->SetProcessOrdering(aMultipleScattering, idxPostStep, 1);
      pmanager->SetProcessOrdering(aionIonization, idxPostStep, 2);

      // map to G3 controls
      controlMap->Add(aionIonization, kLOSS); 
      controlMap->Add(aMultipleScattering, kMULS); 

      // map to AliMCProcess codes
      mcMap->Add(aionIonization, kPEnergyLoss); 
      mcMap->Add(aMultipleScattering, kPMultipleScattering); 
    } 
    else if ((!particle->IsShortLived()) &&
	    (particle->GetPDGCharge() != 0.0) && 
	    (particle->GetParticleName() != "chargedgeantino")) {
      // construct processes for all others charged particles 
      // except geantino
      G4VProcess* aMultipleScattering = new G4MultipleScattering();
      G4VProcess* anIonisation = new G4hIonisation();

      // add processes
      pmanager->AddProcess(anIonisation);
      pmanager->AddProcess(aMultipleScattering);

      // set ordering for AlongStepDoIt
      pmanager->SetProcessOrdering(aMultipleScattering, idxAlongStep,  1);
      pmanager->SetProcessOrdering(anIonisation, idxAlongStep,  2);

      // set ordering for PostStepDoIt
      pmanager->SetProcessOrdering(aMultipleScattering, idxPostStep, 1);
      pmanager->SetProcessOrdering(anIonisation, idxPostStep, 2);

      // map to G3 controls
      controlMap->Add(anIonisation, kLOSS); 
      controlMap->Add(aMultipleScattering, kMULS); 

      // map to AliMCProcess codes
      mcMap->Add(anIonisation, kPEnergyLoss); 
      mcMap->Add(aMultipleScattering, kPMultipleScattering); 
    }
  }
  if (verboseLevel>0)
    G4cout << "### Electromagnetic physics constructed." << G4endl;
}
