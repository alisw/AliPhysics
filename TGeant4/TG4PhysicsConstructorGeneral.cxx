// $Id$
// Category: physics
//
// Author: I. Hrivnacova
//
// Class TG4PhysicsConstructorGeneral
// ------------------------------
// See the class description in the header file.
// According to ExN04IonPhysics.cc,v 1.1.2.1 2001/06/28 19:07:37 gunter Exp 
// GEANT4 tag Name: geant4-03-02

#include "TG4PhysicsConstructorGeneral.h"
#include "TG4ProcessControlMap.h"
#include "TG4ProcessMCMap.h"
#include "TG4G3Control.h"
#include "TG4ExtDecayer.h"
#include "AliDecayer.h"
#include "AliMC.h"

#include <G4ParticleDefinition.hh>
#include <G4ProcessManager.hh>
#include <G4ChargedGeantino.hh>
#include <G4Geantino.hh>

//_____________________________________________________________________________
TG4PhysicsConstructorGeneral::TG4PhysicsConstructorGeneral(const G4String& name)
  : TG4VPhysicsConstructor(name) {
//
}

//_____________________________________________________________________________
TG4PhysicsConstructorGeneral::TG4PhysicsConstructorGeneral(
						   G4int verboseLevel,
                                                   const G4String& name)
  : TG4VPhysicsConstructor(name, verboseLevel) {
//
}

//_____________________________________________________________________________
TG4PhysicsConstructorGeneral::~TG4PhysicsConstructorGeneral() {
//
}

// protected methods

//_____________________________________________________________________________
void TG4PhysicsConstructorGeneral::ConstructParticle()
{
// Instantiates particles.
// ---

  // pseudo-particles
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();  
}

//_____________________________________________________________________________
void TG4PhysicsConstructorGeneral::ConstructProcess()
{
// Constructs electromagnetic processes for e+.
// ---

  // Set external decayer
  AliDecayer* aliDecayer = gMC->Decayer(); 
  if (aliDecayer) {
    TG4ExtDecayer* tg4Decayer = new TG4ExtDecayer(aliDecayer);
       // the tg4Decayer is deleted in G4Decay destructor
    tg4Decayer->VerboseLevel(VerboseLevel());   
    fDecayProcess.SetExtDecayer(tg4Decayer);
    
    if (VerboseLevel() > 0) { 
      G4cout << "### External decayer is set" << G4endl;
    }  
  } 

  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if (fDecayProcess.IsApplicable(*particle)) { 
      pmanager ->AddProcess(&fDecayProcess);
      // set ordering for PostStepDoIt and AtRestDoIt
      pmanager ->SetProcessOrdering(&fDecayProcess, idxPostStep);
      pmanager ->SetProcessOrdering(&fDecayProcess, idxAtRest);
    }
  }
  
  // map to G3 controls
  TG4ProcessControlMap* processMap = TG4ProcessControlMap::Instance();
  processMap->Add(&fDecayProcess, kDCAY); 

  if (VerboseLevel() > 0) {
    G4cout << "### General physics constructed." << G4endl;
  }  
}
