// $Id$
// Category: physics
//
// Author: I. Hrivnacova
//
// Class TG4PhysicsConstructorIon
// ------------------------------
// See the class description in the header file.
// According to ExN04IonPhysics.cc,v 1.1.2.1 2001/06/28 19:07:37 gunter Exp 
// GEANT4 tag Name: geant4-03-02

#include "TG4PhysicsConstructorIon.h"
#include "TG4ProcessControlMap.h"
#include "TG4ProcessMCMap.h"
#include "TG4Globals.h"

#include <G4ParticleDefinition.hh>
#include <G4ProcessManager.hh>
#include "G4IonConstructor.hh"

//_____________________________________________________________________________
TG4PhysicsConstructorIon::TG4PhysicsConstructorIon(const G4String& name)
  : TG4VPhysicsConstructor(name),
    fSetEM(true),
    fSetHadron(true) {
//
}

//_____________________________________________________________________________
TG4PhysicsConstructorIon::TG4PhysicsConstructorIon(G4int verboseLevel,
                                                   G4bool setEM,
						   G4bool setHadron,
						   const G4String& name)
  : TG4VPhysicsConstructor(name, verboseLevel),
    fSetEM(setEM),
    fSetHadron(setHadron) {
//
}

//_____________________________________________________________________________
TG4PhysicsConstructorIon::TG4PhysicsConstructorIon(
                                     const TG4PhysicsConstructorIon& right)
{
//
  TG4Globals::Exception("TG4PhysicsConstructorIon is protected from copying.");
}

//_____________________________________________________________________________
TG4PhysicsConstructorIon::~TG4PhysicsConstructorIon() {
//
}

// operators

//_____________________________________________________________________________
TG4PhysicsConstructorIon& 
TG4PhysicsConstructorIon::operator=(const TG4PhysicsConstructorIon &right)
{
  // check assignement to self
  if (this == &right) return *this;
  
  TG4Globals::Exception(
    "TG4PhysicsConstructorIon is protected from assigning.");

  return *this;
}


// private methods

//_____________________________________________________________________________
void TG4PhysicsConstructorIon::ConstructHadProcessForGenericIon()
{
// Constructs electromagnetic processes for generic ion.
// ---

  // add process
  G4ProcessManager* pManager = G4GenericIon::GenericIon()->GetProcessManager();
  pManager->AddDiscreteProcess(&fElasticProcess);

  // map to G3 controls
  TG4ProcessControlMap* controlMap = TG4ProcessControlMap::Instance();
  controlMap->Add(&fElasticProcess, kHADR); 

  // map to AliMCProcess codes
  TG4ProcessMCMap* mcMap = TG4ProcessMCMap::Instance();
  mcMap->Add(&fElasticProcess, kPHElastic); 
}

//_____________________________________________________________________________
void TG4PhysicsConstructorIon::ConstructHadProcessForDeuteron()
{
// Constructs electromagnetic processes for deuteron.
// ---

  // add process
  G4ProcessManager* pManager = G4Deuteron::Deuteron()->GetProcessManager();
  pManager->AddDiscreteProcess(&fElasticProcess);

  fDeuteronModel = new G4LEDeuteronInelastic();
  fDeuteronProcess.RegisterMe(fDeuteronModel);
  pManager->AddDiscreteProcess(&fDeuteronProcess);

  // map to G3 controls
  TG4ProcessControlMap* controlMap = TG4ProcessControlMap::Instance();
  controlMap->Add(&fDeuteronProcess, kHADR); 

  // map to AliMCProcess codes
  TG4ProcessMCMap* mcMap = TG4ProcessMCMap::Instance();
  mcMap->Add(&fDeuteronProcess, kPHInhelastic); 
}

//_____________________________________________________________________________
void TG4PhysicsConstructorIon::ConstructHadProcessForTriton()
{
// Constructs electromagnetic processes for triton.
// ---

  // add process
  G4ProcessManager* pManager = G4Triton::Triton()->GetProcessManager();
  pManager->AddDiscreteProcess(&fElasticProcess);

  fTritonModel = new G4LETritonInelastic();
  fTritonProcess.RegisterMe(fTritonModel);
  pManager->AddDiscreteProcess(&fTritonProcess);

  // map to G3 controls
  TG4ProcessControlMap* controlMap = TG4ProcessControlMap::Instance();
  controlMap->Add(&fTritonProcess, kHADR); 

  // map to AliMCProcess codes
  TG4ProcessMCMap* mcMap = TG4ProcessMCMap::Instance();
  mcMap->Add(&fTritonProcess, kPHInhelastic); 
}

//_____________________________________________________________________________
void TG4PhysicsConstructorIon::ConstructHadProcessForAlpha()
{
// Constructs electromagnetic processes for alpha.
// ---

  // add process
  G4ProcessManager* pManager = G4Alpha::Alpha()->GetProcessManager();
  pManager->AddDiscreteProcess(&fElasticProcess);

  fAlphaModel = new G4LEAlphaInelastic();
  fAlphaProcess.RegisterMe(fAlphaModel);
  pManager->AddDiscreteProcess(&fAlphaProcess);

  // map to G3 controls
  TG4ProcessControlMap* controlMap = TG4ProcessControlMap::Instance();
  controlMap->Add(&fAlphaProcess, kHADR); 

  // map to AliMCProcess codes
  TG4ProcessMCMap* mcMap = TG4ProcessMCMap::Instance();
  mcMap->Add(&fAlphaProcess, kPHInhelastic); 
}

//_____________________________________________________________________________
void TG4PhysicsConstructorIon::ConstructHadProcessForHe3()
{
// Constructs electromagnetic processes for He3.
// ---

  // add process
  G4ProcessManager* pManager = G4He3::He3()->GetProcessManager();
  pManager->AddDiscreteProcess(&fElasticProcess);
}

//_____________________________________________________________________________
void TG4PhysicsConstructorIon::ConstructEMProcessForGenericIon()
{
// Constructs electromagnetic processes for generic ion.
// ---

  // add process
  G4ProcessManager* pManager = G4GenericIon::GenericIon()->GetProcessManager();
  pManager->AddProcess(&fIonIonisation, ordInActive, 2, 2);
  pManager->AddProcess(&fIonMultipleScattering);

  // set ordering
  pManager->SetProcessOrdering(&fIonMultipleScattering, idxAlongStep,  1);
  pManager->SetProcessOrdering(&fIonMultipleScattering, idxPostStep,  1);

  // map to G3 controls
  TG4ProcessControlMap* controlMap = TG4ProcessControlMap::Instance();
  controlMap->Add(&fIonIonisation, kLOSS); 
  controlMap->Add(&fIonMultipleScattering, kMULS); 

  // map to AliMCProcess codes
  TG4ProcessMCMap* mcMap = TG4ProcessMCMap::Instance();
  mcMap->Add(&fIonIonisation, kPEnergyLoss); 
  mcMap->Add(&fIonMultipleScattering, kPMultipleScattering); 
}

//_____________________________________________________________________________
void TG4PhysicsConstructorIon::ConstructEMProcessForDeuteron()
{
// Constructs electromagnetic processes for deuteron.
// ---

  // add process
  G4ProcessManager* pManager = G4Deuteron::Deuteron()->GetProcessManager();

  pManager->AddProcess(&fDeuteronIonisation, ordInActive, 2, 2);
  pManager->AddProcess(&fDeuteronMultipleScattering);

  // set ordering
  pManager->SetProcessOrdering(&fDeuteronMultipleScattering, idxAlongStep,  1);
  pManager->SetProcessOrdering(&fDeuteronMultipleScattering, idxPostStep,  1);

  // map to G3 controls
  TG4ProcessControlMap* controlMap = TG4ProcessControlMap::Instance();
  controlMap->Add(&fDeuteronIonisation, kLOSS); 
  controlMap->Add(&fDeuteronMultipleScattering, kMULS); 

  // map to AliMCProcess codes
  TG4ProcessMCMap* mcMap = TG4ProcessMCMap::Instance();
  mcMap->Add(&fDeuteronIonisation, kPEnergyLoss); 
  mcMap->Add(&fDeuteronMultipleScattering, kPMultipleScattering); 
}

//_____________________________________________________________________________
void TG4PhysicsConstructorIon::ConstructEMProcessForTriton()
{
// Constructs electromagnetic processes for triton.
// ---

  // add process
  G4ProcessManager* pManager = G4Triton::Triton()->GetProcessManager();

  pManager->AddProcess(&fTritonIonisation, ordInActive, 2, 2);
  pManager->AddProcess(&fTritonMultipleScattering);

  // set ordering
  pManager->SetProcessOrdering(&fTritonMultipleScattering, idxAlongStep,  1);
  pManager->SetProcessOrdering(&fTritonMultipleScattering, idxPostStep,  1);

  // map to G3 controls
  TG4ProcessControlMap* controlMap = TG4ProcessControlMap::Instance();
  controlMap->Add(&fTritonIonisation, kLOSS); 
  controlMap->Add(&fTritonMultipleScattering, kMULS); 

  // map to AliMCProcess codes
  TG4ProcessMCMap* mcMap = TG4ProcessMCMap::Instance();
  mcMap->Add(&fTritonIonisation, kPEnergyLoss); 
  mcMap->Add(&fTritonMultipleScattering, kPMultipleScattering); 
}

//_____________________________________________________________________________
void TG4PhysicsConstructorIon::ConstructEMProcessForAlpha()
{
// Constructs electromagnetic processes for alpha.
// ---

  // add process
  G4ProcessManager* pManager = G4Alpha::Alpha()->GetProcessManager();

  pManager->AddProcess(&fAlphaIonisation, ordInActive, 2, 2);
  pManager->AddProcess(&fAlphaMultipleScattering);

  // set ordering
  pManager->SetProcessOrdering(&fAlphaMultipleScattering, idxAlongStep,  1);
  pManager->SetProcessOrdering(&fAlphaMultipleScattering, idxPostStep,  1);

  // map to G3 controls
  TG4ProcessControlMap* controlMap = TG4ProcessControlMap::Instance();
  controlMap->Add(&fAlphaIonisation, kLOSS); 
  controlMap->Add(&fAlphaMultipleScattering, kMULS); 

  // map to AliMCProcess codes
  TG4ProcessMCMap* mcMap = TG4ProcessMCMap::Instance();
  mcMap->Add(&fAlphaIonisation, kPEnergyLoss); 
  mcMap->Add(&fAlphaMultipleScattering, kPMultipleScattering); 
}

//_____________________________________________________________________________
void TG4PhysicsConstructorIon::ConstructEMProcessForHe3()
{
// Constructs electromagnetic processes for He3.
// ---

  // add process
  G4ProcessManager* pManager = G4He3::He3()->GetProcessManager();
  pManager->AddProcess(&fHe3Ionisation, ordInActive, 2, 2);
  pManager->AddProcess(&fHe3MultipleScattering);

  // set ordering
  pManager->SetProcessOrdering(&fHe3MultipleScattering, idxAlongStep,  1);
  pManager->SetProcessOrdering(&fHe3MultipleScattering, idxPostStep,  1);

  // map to G3 controls
  TG4ProcessControlMap* controlMap = TG4ProcessControlMap::Instance();
  controlMap->Add(&fHe3Ionisation, kLOSS); 
  controlMap->Add(&fHe3MultipleScattering, kMULS); 

  // map to AliMCProcess codes
  TG4ProcessMCMap* mcMap = TG4ProcessMCMap::Instance();
  mcMap->Add(&fHe3Ionisation, kPEnergyLoss); 
  mcMap->Add(&fHe3MultipleScattering, kPMultipleScattering); 
}


// protected methods

//_____________________________________________________________________________
void TG4PhysicsConstructorIon::ConstructParticle()
{
// Instantiates particles.
// ---

  //  Construct light ions
  G4IonConstructor pConstructor;
  pConstructor.ConstructParticle();  
}

//_____________________________________________________________________________
void TG4PhysicsConstructorIon::ConstructProcess()
{
// Constructs electromagnetic processes for e+.
// ---

  if (fSetHadron) {
    // Elastic Process
    fElasticModel = new G4LElastic();
    fElasticProcess.RegisterMe(fElasticModel);

    // Hadron processes
    ConstructHadProcessForGenericIon();
    ConstructHadProcessForDeuteron();
    ConstructHadProcessForTriton();
    ConstructHadProcessForAlpha();
    ConstructHadProcessForHe3();

    if (VerboseLevel() > 1) {
      G4cout << "### Ion EM physics constructed." << G4endl;
    }  
  }  

  if (fSetEM) {
    // EM processes
    ConstructEMProcessForGenericIon();
    ConstructEMProcessForDeuteron();
    ConstructEMProcessForTriton();
    ConstructEMProcessForAlpha();
    ConstructEMProcessForHe3();

    if (VerboseLevel() > 1) {
      G4cout << "### Ion hadron physics constructed." << G4endl;
    }  
  }  

  if (VerboseLevel() > 0) {
    G4cout << "### Ion physics constructed." << G4endl;
  }  
}
