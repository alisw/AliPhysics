// $Id$
// Category: physics
//
// Author: I. Hrivnacova
//
// Class TG4PhysicsConstructorHadron
// ---------------------------------
// See the class description in the header file.
// According to ExN04HadronPhysics.cc,v 1.1.2.1 2001/06/28 19:07:37 gunter Exp
// GEANT4 tag Name: geant4-03-02

#include "TG4PhysicsConstructorHadron.h"
#include "TG4ProcessControlMap.h"
#include "TG4ProcessMCMap.h"
#include "TG4Globals.h"

#include <G4ParticleDefinition.hh>
#include <G4ProcessManager.hh>
#include <G4MesonConstructor.hh>
#include <G4BaryonConstructor.hh>
#include <G4ShortLivedConstructor.hh>

//_____________________________________________________________________________
TG4PhysicsConstructorHadron::TG4PhysicsConstructorHadron(const G4String& name)
  : G4VPhysicsConstructor(name)
{
//
  SetVerboseLevel(1);
}

//_____________________________________________________________________________
TG4PhysicsConstructorHadron::TG4PhysicsConstructorHadron(
                                     const TG4PhysicsConstructorHadron& right)
{
//
  TG4Globals::Exception("TG4PhysicsConstructorHadron is protected from copying.");
}

//_____________________________________________________________________________
TG4PhysicsConstructorHadron::~TG4PhysicsConstructorHadron() {
//
  for (G4int i=0; i<fOtherProcesses.size(); i++) delete fOtherProcesses[i];
}

// operators

//_____________________________________________________________________________
TG4PhysicsConstructorHadron& 
TG4PhysicsConstructorHadron::operator=(const TG4PhysicsConstructorHadron &right)
{
  // check assignement to self
  if (this == &right) return *this;
  
  TG4Globals::Exception(
    "TG4PhysicsConstructorHadron is protected from assigning.");

  return *this;
}


// private methods

//_____________________________________________________________________________
void TG4PhysicsConstructorHadron::ConstructProcessForPionPlus()
{
// Construct processes for pi+.
// ---

  // add process
  G4ProcessManager* pManager = G4PionPlus::PionPlus()->GetProcessManager();
  pManager->AddDiscreteProcess(&fElasticProcess);

  fLEPionPlusModel = new G4LEPionPlusInelastic();
  fHEPionPlusModel = new G4HEPionPlusInelastic();
  fPionPlusInelastic.RegisterMe(fLEPionPlusModel);
  fPionPlusInelastic.RegisterMe(fHEPionPlusModel);
  pManager->AddDiscreteProcess(&fPionPlusInelastic);

  pManager->AddProcess(&fPionPlusIonisation, ordInActive, 2, 2);
  pManager->AddProcess(&fPionPlusMult);

  // set ordering
  pManager->SetProcessOrdering(&fPionPlusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(&fPionPlusMult, idxPostStep, 1);

  // map to G3 controls
  TG4ProcessControlMap* controlMap = TG4ProcessControlMap::Instance();
  controlMap->Add(&fElasticProcess, kHADR); 
  controlMap->Add(&fPionPlusInelastic, kHADR); 
  controlMap->Add(&fPionPlusIonisation, kLOSS); 
  controlMap->Add(&fPionPlusMult, kMULS); 

  // map to G3 ALIMCProcess codes
  TG4ProcessMCMap* mcMap = TG4ProcessMCMap::Instance();
  mcMap->Add(&fElasticProcess, kPHElastic); 
  mcMap->Add(&fPionPlusInelastic, kPHInhelastic); 
  mcMap->Add(&fPionPlusIonisation, kPEnergyLoss); 
  mcMap->Add(&fPionPlusMult, kPMultipleScattering); 
}

//_____________________________________________________________________________
void TG4PhysicsConstructorHadron::ConstructProcessForPionMinus()
{
// Construct processes for pi-.
// ---

  // add process & set ordering
  G4ProcessManager* pManager = G4PionMinus::PionMinus()->GetProcessManager();
  pManager->AddDiscreteProcess(&fElasticProcess);

  fLEPionMinusModel = new G4LEPionMinusInelastic();
  fHEPionMinusModel = new G4HEPionMinusInelastic();
  fPionMinusInelastic.RegisterMe(fLEPionMinusModel);
  fPionMinusInelastic.RegisterMe(fHEPionMinusModel);
  pManager->AddDiscreteProcess(&fPionMinusInelastic);

  pManager->AddProcess(&fPionMinusIonisation, ordInActive, 2, 2);

  pManager->AddProcess(&fPionMinusMult);
  pManager->SetProcessOrdering(&fPionMinusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(&fPionMinusMult, idxPostStep, 1);

  pManager->AddRestProcess(&fPionMinusAbsorption, ordDefault);

  // map to G3 controls
  TG4ProcessControlMap* controlMap = TG4ProcessControlMap::Instance();
  controlMap->Add(&fPionMinusInelastic, kHADR); 
  controlMap->Add(&fPionMinusIonisation, kLOSS); 
  controlMap->Add(&fPionMinusMult, kMULS); 
  controlMap->Add(&fPionMinusAbsorption, kHADR); 

  // map to G3 ALIMCProcess codes
  TG4ProcessMCMap* mcMap = TG4ProcessMCMap::Instance();
  mcMap->Add(&fPionMinusInelastic, kPHInhelastic); 
  mcMap->Add(&fPionMinusIonisation, kPEnergyLoss); 
  mcMap->Add(&fPionMinusMult, kPMultipleScattering); 
  mcMap->Add(&fPionMinusAbsorption, kPNuclearAbsorption); 
}

//_____________________________________________________________________________
void TG4PhysicsConstructorHadron::ConstructProcessForKaonPlus()
{
// Construct processes for K+.
// ---

  // add process
  G4ProcessManager* pManager = G4KaonPlus::KaonPlus()->GetProcessManager();
  pManager->AddDiscreteProcess(&fElasticProcess);

  fLEKaonPlusModel = new G4LEKaonPlusInelastic();
  fHEKaonPlusModel = new G4HEKaonPlusInelastic();
  fKaonPlusInelastic.RegisterMe(fLEKaonPlusModel);
  fKaonPlusInelastic.RegisterMe(fHEKaonPlusModel);
  pManager->AddDiscreteProcess(&fKaonPlusInelastic);

  pManager->AddProcess(&fKaonPlusIonisation, ordInActive, 2, 2);
  pManager->AddProcess(&fKaonPlusMult);

  // set ordering
  pManager->SetProcessOrdering(&fKaonPlusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(&fKaonPlusMult, idxPostStep, 1);

  // map to G3 controls
  TG4ProcessControlMap* controlMap = TG4ProcessControlMap::Instance();
  controlMap->Add(&fKaonPlusInelastic, kHADR); 
  controlMap->Add(&fKaonPlusIonisation, kLOSS); 
  controlMap->Add(&fKaonPlusMult, kMULS); 

  // map to G3 ALIMCProcess codes
  TG4ProcessMCMap* mcMap = TG4ProcessMCMap::Instance();
  mcMap->Add(&fKaonPlusInelastic, kPHInhelastic); 
  mcMap->Add(&fKaonPlusIonisation, kPEnergyLoss); 
  mcMap->Add(&fKaonPlusMult, kPMultipleScattering); 
}

//_____________________________________________________________________________
void TG4PhysicsConstructorHadron::ConstructProcessForKaonMinus()
{
// Construct processes for K-.
// ---

  // add process & set ordering
  G4ProcessManager* pManager = G4KaonMinus::KaonMinus()->GetProcessManager();
  pManager->AddDiscreteProcess(&fElasticProcess);

  fLEKaonMinusModel = new G4LEKaonMinusInelastic();
  fHEKaonMinusModel = new G4HEKaonMinusInelastic();
  fKaonMinusInelastic.RegisterMe(fLEKaonMinusModel);
  fKaonMinusInelastic.RegisterMe(fHEKaonMinusModel);
  pManager->AddDiscreteProcess(&fKaonMinusInelastic);

  pManager->AddProcess(&fKaonMinusIonisation, ordInActive, 2, 2);

  pManager->AddProcess(&fKaonMinusMult);
  pManager->SetProcessOrdering(&fKaonMinusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(&fKaonMinusMult, idxPostStep, 1);

  pManager->AddRestProcess(&fKaonMinusAbsorption, ordDefault);

  // map to G3 controls
  TG4ProcessControlMap* controlMap = TG4ProcessControlMap::Instance();
  controlMap->Add(&fKaonMinusInelastic, kHADR); 
  controlMap->Add(&fKaonMinusIonisation, kLOSS); 
  controlMap->Add(&fKaonMinusMult, kMULS); 
  controlMap->Add(&fKaonMinusAbsorption, kHADR); 

  // map to G3 ALIMCProcess codes
  TG4ProcessMCMap* mcMap = TG4ProcessMCMap::Instance();
  mcMap->Add(&fKaonMinusInelastic, kPHInhelastic); 
  mcMap->Add(&fKaonMinusIonisation, kPEnergyLoss); 
  mcMap->Add(&fKaonMinusMult, kPMultipleScattering); 
  mcMap->Add(&fKaonMinusAbsorption, kPNuclearAbsorption); 
}

//_____________________________________________________________________________
void TG4PhysicsConstructorHadron::ConstructProcessForKaonZeroLong()
{
// Construct processes for K0L.
// ---

  // add process
  G4ProcessManager* pManager 
    = G4KaonZeroLong::KaonZeroLong()->GetProcessManager();
  pManager->AddDiscreteProcess(&fElasticProcess);

  fLEKaonZeroLModel = new G4LEKaonZeroLInelastic();
  fHEKaonZeroLModel = new G4HEKaonZeroInelastic();
  fKaonZeroLInelastic.RegisterMe(fLEKaonZeroLModel);
  fKaonZeroLInelastic.RegisterMe(fHEKaonZeroLModel);
  pManager->AddDiscreteProcess(&fKaonZeroLInelastic);

  // map to G3 controls
  TG4ProcessControlMap* controlMap = TG4ProcessControlMap::Instance();
  controlMap->Add(&fKaonZeroLInelastic, kHADR); 

  // map to G3 ALIMCProcess codes
  TG4ProcessMCMap* mcMap = TG4ProcessMCMap::Instance();
  mcMap->Add(&fKaonZeroLInelastic, kPHInhelastic); 
}

//_____________________________________________________________________________
void TG4PhysicsConstructorHadron::ConstructProcessForKaonZeroShort()
{
// Construct processes for K0S.
// ---

  // add process
  G4ProcessManager* pManager 
    = G4KaonZeroShort::KaonZeroShort()->GetProcessManager();
  pManager->AddDiscreteProcess(&fElasticProcess);

  fLEKaonZeroSModel = new G4LEKaonZeroSInelastic();
  fHEKaonZeroSModel = new G4HEKaonZeroInelastic();
  fKaonZeroSInelastic.RegisterMe(fLEKaonZeroSModel);
  fKaonZeroSInelastic.RegisterMe(fHEKaonZeroSModel);
  pManager->AddDiscreteProcess(&fKaonZeroSInelastic);

  // map to G3 controls
  TG4ProcessControlMap* controlMap = TG4ProcessControlMap::Instance();
  controlMap->Add(&fKaonZeroSInelastic, kHADR); 

  // map to G3 ALIMCProcess codes
  TG4ProcessMCMap* mcMap = TG4ProcessMCMap::Instance();
  mcMap->Add(&fKaonZeroSInelastic, kPHInhelastic); 
}

//_____________________________________________________________________________
void TG4PhysicsConstructorHadron::ConstructProcessForProton()
{
// Construct processes for proton.
// ---

  // add process
  G4ProcessManager* pManager = G4Proton::Proton()->GetProcessManager();
  pManager->AddDiscreteProcess(&fElasticProcess);

  fLEProtonModel = new G4LEProtonInelastic();
  fHEProtonModel = new G4HEProtonInelastic();
  fProtonInelastic.RegisterMe(fLEProtonModel);
  fProtonInelastic.RegisterMe(fHEProtonModel);
  pManager->AddDiscreteProcess(&fProtonInelastic);

  pManager->AddProcess(&fProtonIonisation, ordInActive, 2, 2);
  pManager->AddProcess(&fProtonMult);

  // set ordering
  pManager->SetProcessOrdering(&fProtonMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(&fProtonMult, idxPostStep, 1);

  // map to G3 controls
  TG4ProcessControlMap* controlMap = TG4ProcessControlMap::Instance();
  controlMap->Add(&fProtonInelastic, kHADR); 
  controlMap->Add(&fProtonIonisation, kLOSS); 
  controlMap->Add(&fProtonMult, kMULS); 

  // map to G3 ALIMCProcess codes
  TG4ProcessMCMap* mcMap = TG4ProcessMCMap::Instance();
  mcMap->Add(&fProtonInelastic, kPHInhelastic); 
  mcMap->Add(&fProtonIonisation, kPEnergyLoss); 
  mcMap->Add(&fProtonMult, kPMultipleScattering); 
}

//_____________________________________________________________________________
void TG4PhysicsConstructorHadron::ConstructProcessForAntiProton()
{
// Construct processes for anti-proton.
// ---

  // add process & set ordering
  G4ProcessManager* pManager = G4AntiProton::AntiProton()->GetProcessManager();
  pManager->AddDiscreteProcess(&fElasticProcess);

  fLEAntiProtonModel = new G4LEAntiProtonInelastic();
  fHEAntiProtonModel = new G4HEAntiProtonInelastic();
  fAntiProtonInelastic.RegisterMe(fLEAntiProtonModel);
  fAntiProtonInelastic.RegisterMe(fHEAntiProtonModel);
  pManager->AddDiscreteProcess(&fAntiProtonInelastic);

  pManager->AddProcess(&fAntiProtonIonisation, ordInActive, 2, 2);

  pManager->AddProcess(&fAntiProtonMult);
  pManager->SetProcessOrdering(&fAntiProtonMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(&fAntiProtonMult, idxPostStep, 1);

  pManager->AddRestProcess(&fAntiProtonAnnihilation);

  // map to G3 controls
  TG4ProcessControlMap* controlMap = TG4ProcessControlMap::Instance();
  controlMap->Add(&fAntiProtonInelastic, kHADR); 
  controlMap->Add(&fAntiProtonIonisation, kLOSS); 
  controlMap->Add(&fAntiProtonMult, kMULS); 
  controlMap->Add(&fAntiProtonAnnihilation, kHADR); 

  // map to G3 ALIMCProcess codes
  TG4ProcessMCMap* mcMap = TG4ProcessMCMap::Instance();
  mcMap->Add(&fAntiProtonInelastic, kPHInhelastic); 
  mcMap->Add(&fAntiProtonIonisation, kPEnergyLoss); 
  mcMap->Add(&fAntiProtonMult, kPMultipleScattering); 
  mcMap->Add(&fAntiProtonAnnihilation, kPPbarAnnihilation); 
}

//_____________________________________________________________________________
void TG4PhysicsConstructorHadron::ConstructProcessForNeutron()
{
// Construct processes for neutron.
// ---

  // add process
  G4ProcessManager* pManager = G4Neutron::Neutron()->GetProcessManager();
  pManager->AddDiscreteProcess(&fElasticProcess);

  fLENeutronModel = new G4LENeutronInelastic();
  fHENeutronModel = new G4HENeutronInelastic();
  fNeutronInelastic.RegisterMe(fLENeutronModel);
  fNeutronInelastic.RegisterMe(fHENeutronModel);
  pManager->AddDiscreteProcess(&fNeutronInelastic);

  //fNeutronFissionModel = new G4LFission();
  //fNeutronFission.RegisterMe(fNeutronFissionModel);
  //pManager->AddDiscreteProcess(&NeutronFission);

  //fNeutronCaptureModel = new G4LCapture();
  //fNeutronCapture.RegisterMe(fNeutronCaptureModel);
  //pManager->AddDiscreteProcess(&fNeutronCapture);

  // map to G3 controls
  TG4ProcessControlMap* controlMap = TG4ProcessControlMap::Instance();
  controlMap->Add(&fNeutronInelastic, kHADR); 

  // map to G3 ALIMCProcess codes
  TG4ProcessMCMap* mcMap = TG4ProcessMCMap::Instance();
  mcMap->Add(&fNeutronInelastic, kPHInhelastic); 
}

//_____________________________________________________________________________
void TG4PhysicsConstructorHadron::ConstructProcessForAntiNeutron()
{
// Construct processes for anti-neutron.
// ---

  // add process
  G4ProcessManager* pManager = G4AntiNeutron::AntiNeutron()->GetProcessManager();
  pManager->AddDiscreteProcess(&fElasticProcess);

  fLEAntiNeutronModel = new G4LEAntiNeutronInelastic();
  fHEAntiNeutronModel = new G4HEAntiNeutronInelastic();
  fAntiNeutronInelastic.RegisterMe(fLEAntiNeutronModel);
  fAntiNeutronInelastic.RegisterMe(fHEAntiNeutronModel);
  pManager->AddDiscreteProcess(&fAntiNeutronInelastic);

  pManager->AddRestProcess(&fAntiNeutronAnnihilation);

  // map to G3 controls
  TG4ProcessControlMap* controlMap = TG4ProcessControlMap::Instance();
  controlMap->Add(&fAntiNeutronInelastic, kHADR); 
  controlMap->Add(&fAntiNeutronAnnihilation, kHADR); 

  // map to G3 ALIMCProcess codes
  TG4ProcessMCMap* mcMap = TG4ProcessMCMap::Instance();
  mcMap->Add(&fAntiNeutronInelastic, kPHInhelastic); 
  mcMap->Add(&fAntiNeutronAnnihilation, kPNoProcess);   
}

//_____________________________________________________________________________
void TG4PhysicsConstructorHadron::ConstructProcessForLambda()
{
// Construct processes for Lambda.
// ---

  // add process
  G4ProcessManager* pManager = G4Lambda::Lambda()->GetProcessManager();
  pManager->AddDiscreteProcess(&fElasticProcess);

  fLELambdaModel = new G4LELambdaInelastic();
  fHELambdaModel = new G4HELambdaInelastic();
  fLambdaInelastic.RegisterMe(fLELambdaModel);
  fLambdaInelastic.RegisterMe(fHELambdaModel);
  pManager->AddDiscreteProcess(&fLambdaInelastic);

  // map to G3 controls
  TG4ProcessControlMap* controlMap = TG4ProcessControlMap::Instance();
  controlMap->Add(&fLambdaInelastic, kHADR); 

  // map to G3 ALIMCProcess codes
  TG4ProcessMCMap* mcMap = TG4ProcessMCMap::Instance();
  mcMap->Add(&fLambdaInelastic, kPHInhelastic); 
}

//_____________________________________________________________________________
void TG4PhysicsConstructorHadron::ConstructProcessForAntiLambda()
{
// Construct processes for anti-Lambda.
// ---

  // add process
  G4ProcessManager* pManager = G4AntiLambda::AntiLambda()->GetProcessManager();
  pManager->AddDiscreteProcess(&fElasticProcess);

  fLEAntiLambdaModel = new G4LEAntiLambdaInelastic();
  fHEAntiLambdaModel = new G4HEAntiLambdaInelastic();
  fAntiLambdaInelastic.RegisterMe(fLEAntiLambdaModel);
  fAntiLambdaInelastic.RegisterMe(fHEAntiLambdaModel);
  pManager->AddDiscreteProcess(&fAntiLambdaInelastic);

  // map to G3 controls
  TG4ProcessControlMap* controlMap = TG4ProcessControlMap::Instance();
  controlMap->Add(&fAntiLambdaInelastic, kHADR); 

  // map to G3 ALIMCProcess codes
  TG4ProcessMCMap* mcMap = TG4ProcessMCMap::Instance();
  mcMap->Add(&fAntiLambdaInelastic, kPHInhelastic); 
}

//_____________________________________________________________________________
void TG4PhysicsConstructorHadron::ConstructProcessForSigmaMinus()
{
// Construct processes for Sigma-.
// ---

  // add process & set ordering
  G4ProcessManager* pManager = G4SigmaMinus::SigmaMinus()->GetProcessManager();
  pManager->AddDiscreteProcess(&fElasticProcess);

  fLESigmaMinusModel = new G4LESigmaMinusInelastic();
  fHESigmaMinusModel = new G4HESigmaMinusInelastic();
  fSigmaMinusInelastic.RegisterMe(fLESigmaMinusModel);
  fSigmaMinusInelastic.RegisterMe(fHESigmaMinusModel);
  pManager->AddDiscreteProcess(&fSigmaMinusInelastic);

  pManager->AddProcess(&fSigmaMinusIonisation, ordInActive, 2, 2);

  pManager->AddProcess(&fSigmaMinusMult);
  pManager->SetProcessOrdering(&fSigmaMinusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(&fSigmaMinusMult, idxPostStep, 1);

  // map to G3 controls
  TG4ProcessControlMap* controlMap = TG4ProcessControlMap::Instance();
  controlMap->Add(&fSigmaMinusInelastic, kHADR); 
  controlMap->Add(&fSigmaMinusIonisation, kLOSS); 
  controlMap->Add(&fSigmaMinusMult, kMULS); 

  // map to G3 ALIMCProcess codes
  TG4ProcessMCMap* mcMap = TG4ProcessMCMap::Instance();
  mcMap->Add(&fSigmaMinusInelastic, kPHInhelastic); 
  mcMap->Add(&fSigmaMinusIonisation, kPEnergyLoss); 
  mcMap->Add(&fSigmaMinusMult, kPMultipleScattering); 
}

//_____________________________________________________________________________
void TG4PhysicsConstructorHadron::ConstructProcessForAntiSigmaMinus()
{
// Construct processes for anti-Sigma-.
// ---

  // add process & set ordering
  G4ProcessManager* pManager 
    = G4AntiSigmaMinus::AntiSigmaMinus()->GetProcessManager();
  pManager->AddDiscreteProcess(&fElasticProcess);

  fLEAntiSigmaMinusModel = new G4LEAntiSigmaMinusInelastic();
  fHEAntiSigmaMinusModel = new G4HEAntiSigmaMinusInelastic();
  fAntiSigmaMinusInelastic.RegisterMe(fLEAntiSigmaMinusModel);
  fAntiSigmaMinusInelastic.RegisterMe(fHEAntiSigmaMinusModel);
  pManager->AddDiscreteProcess(&fAntiSigmaMinusInelastic);

  pManager->AddProcess(&fAntiSigmaMinusIonisation, ordInActive, 2, 2);

  pManager->AddProcess(&fAntiSigmaMinusMult);
  pManager->SetProcessOrdering(&fAntiSigmaMinusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(&fAntiSigmaMinusMult, idxPostStep, 1);

  // map to G3 controls
  TG4ProcessControlMap* controlMap = TG4ProcessControlMap::Instance();
  controlMap->Add(&fAntiSigmaMinusInelastic, kHADR); 
  controlMap->Add(&fAntiSigmaMinusIonisation, kLOSS); 
  controlMap->Add(&fAntiSigmaMinusMult, kMULS); 

  // map to G3 ALIMCProcess codes
  TG4ProcessMCMap* mcMap = TG4ProcessMCMap::Instance();
  mcMap->Add(&fAntiSigmaMinusInelastic, kPHInhelastic); 
  mcMap->Add(&fAntiSigmaMinusIonisation, kPEnergyLoss); 
  mcMap->Add(&fAntiSigmaMinusMult, kPMultipleScattering); 
}

//_____________________________________________________________________________
void TG4PhysicsConstructorHadron::ConstructProcessForSigmaPlus()
{
// Construct processes for Sigma+.
// ---

  // add process & set ordering
  G4ProcessManager* pManager = G4SigmaPlus::SigmaPlus()->GetProcessManager();
  pManager->AddDiscreteProcess(&fElasticProcess);

  fLESigmaPlusModel = new G4LESigmaPlusInelastic();
  fHESigmaPlusModel = new G4HESigmaPlusInelastic();
  fSigmaPlusInelastic.RegisterMe(fLESigmaPlusModel);
  fSigmaPlusInelastic.RegisterMe(fHESigmaPlusModel);
  pManager->AddDiscreteProcess(&fSigmaPlusInelastic);

  pManager->AddProcess(&fSigmaPlusIonisation, ordInActive, 2, 2);

  pManager->AddProcess(&fSigmaPlusMult);
  pManager->SetProcessOrdering(&fSigmaPlusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(&fSigmaPlusMult, idxPostStep, 1);

  // map to G3 controls
  TG4ProcessControlMap* controlMap = TG4ProcessControlMap::Instance();
  controlMap->Add(&fSigmaPlusInelastic, kHADR); 
  controlMap->Add(&fSigmaPlusIonisation, kLOSS); 
  controlMap->Add(&fSigmaPlusMult, kMULS); 

  // map to G3 ALIMCProcess codes
  TG4ProcessMCMap* mcMap = TG4ProcessMCMap::Instance();
  mcMap->Add(&fSigmaPlusInelastic, kPHInhelastic); 
  mcMap->Add(&fSigmaPlusIonisation, kPEnergyLoss); 
  mcMap->Add(&fSigmaPlusMult, kPMultipleScattering); 
}

//_____________________________________________________________________________
void TG4PhysicsConstructorHadron::ConstructProcessForAntiSigmaPlus()
{
// Construct processes for anti-Sigma+.
// ---

  // add process & set ordering
  G4ProcessManager* pManager 
    = G4AntiSigmaPlus::AntiSigmaPlus()->GetProcessManager();
  pManager->AddDiscreteProcess(&fElasticProcess);

  fLEAntiSigmaPlusModel = new G4LEAntiSigmaPlusInelastic();
  fHEAntiSigmaPlusModel = new G4HEAntiSigmaPlusInelastic();
  fAntiSigmaPlusInelastic.RegisterMe(fLEAntiSigmaPlusModel);
  fAntiSigmaPlusInelastic.RegisterMe(fHEAntiSigmaPlusModel);
  pManager->AddDiscreteProcess(&fAntiSigmaPlusInelastic);

  pManager->AddProcess(&fAntiSigmaPlusIonisation, ordInActive, 2, 2);

  pManager->AddProcess(&fAntiSigmaPlusMult);
  pManager->SetProcessOrdering(&fAntiSigmaPlusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(&fAntiSigmaPlusMult, idxPostStep, 1);

  // map to G3 controls
  TG4ProcessControlMap* controlMap = TG4ProcessControlMap::Instance();
  controlMap->Add(&fAntiSigmaPlusInelastic, kHADR); 
  controlMap->Add(&fAntiSigmaPlusIonisation, kLOSS); 
  controlMap->Add(&fAntiSigmaPlusMult, kMULS); 

  // map to G3 ALIMCProcess codes
  TG4ProcessMCMap* mcMap = TG4ProcessMCMap::Instance();
  mcMap->Add(&fAntiSigmaPlusInelastic, kPHInhelastic); 
  mcMap->Add(&fAntiSigmaPlusIonisation, kPEnergyLoss); 
  mcMap->Add(&fAntiSigmaPlusMult, kPMultipleScattering); 
}

//_____________________________________________________________________________
void TG4PhysicsConstructorHadron::ConstructProcessForXiMinus()
{
// Construct processes for Xi-.
// ---

  // add process & set ordering
  G4ProcessManager* pManager = G4XiMinus::XiMinus()->GetProcessManager();
  pManager->AddDiscreteProcess(&fElasticProcess);

  fLEXiMinusModel = new G4LEXiMinusInelastic();
  fHEXiMinusModel = new G4HEXiMinusInelastic();
  fXiMinusInelastic.RegisterMe(fLEXiMinusModel);
  fXiMinusInelastic.RegisterMe(fHEXiMinusModel);
  pManager->AddDiscreteProcess(&fXiMinusInelastic);

  pManager->AddProcess(&fXiMinusIonisation, ordInActive, 2, 2);

  pManager->AddProcess(&fXiMinusMult);
  pManager->SetProcessOrdering(&fXiMinusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(&fXiMinusMult, idxPostStep, 1);

  // map to G3 controls
  TG4ProcessControlMap* controlMap = TG4ProcessControlMap::Instance();
  controlMap->Add(&fXiMinusInelastic, kHADR); 
  controlMap->Add(&fXiMinusIonisation, kLOSS); 
  controlMap->Add(&fXiMinusMult, kMULS); 

  // map to G3 ALIMCProcess codes
  TG4ProcessMCMap* mcMap = TG4ProcessMCMap::Instance();
  mcMap->Add(&fXiMinusInelastic, kPHInhelastic); 
  mcMap->Add(&fXiMinusIonisation, kPEnergyLoss); 
  mcMap->Add(&fXiMinusMult, kPMultipleScattering); 
}

//_____________________________________________________________________________
void TG4PhysicsConstructorHadron::ConstructProcessForAntiXiMinus()
{
// Construct processes for anti-Xi-.
// ---

  // add process & set ordering
  G4ProcessManager* pManager 
    = G4AntiXiMinus::AntiXiMinus()->GetProcessManager();
  pManager->AddDiscreteProcess(&fElasticProcess);

  fLEAntiXiMinusModel = new G4LEAntiXiMinusInelastic();
  fHEAntiXiMinusModel = new G4HEAntiXiMinusInelastic();
  fAntiXiMinusInelastic.RegisterMe(fLEAntiXiMinusModel);
  fAntiXiMinusInelastic.RegisterMe(fHEAntiXiMinusModel);
  pManager->AddDiscreteProcess(&fAntiXiMinusInelastic);

  pManager->AddProcess(&fAntiXiMinusIonisation, ordInActive, 2, 2);

  pManager->AddProcess(&fAntiXiMinusMult);
  pManager->SetProcessOrdering(&fAntiXiMinusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(&fAntiXiMinusMult, idxPostStep, 1);

  // map to G3 controls
  TG4ProcessControlMap* controlMap = TG4ProcessControlMap::Instance();
  controlMap->Add(&fAntiXiMinusInelastic, kHADR); 
  controlMap->Add(&fAntiXiMinusIonisation, kLOSS); 
  controlMap->Add(&fAntiXiMinusMult, kMULS); 

  // map to G3 ALIMCProcess codes
  TG4ProcessMCMap* mcMap = TG4ProcessMCMap::Instance();
  mcMap->Add(&fAntiXiMinusInelastic, kPHInhelastic); 
  mcMap->Add(&fAntiXiMinusIonisation, kPEnergyLoss); 
  mcMap->Add(&fAntiXiMinusMult, kPMultipleScattering); 
}

//_____________________________________________________________________________
void TG4PhysicsConstructorHadron::ConstructProcessForXiZero()
{
// Construct processes for Xi0.
// ---

  // add process
  G4ProcessManager* pManager = G4XiZero::XiZero()->GetProcessManager();
  pManager->AddDiscreteProcess(&fElasticProcess);

  fLEXiZeroModel = new G4LEXiZeroInelastic();
  fHEXiZeroModel = new G4HEXiZeroInelastic();
  fXiZeroInelastic.RegisterMe(fLEXiZeroModel);
  fXiZeroInelastic.RegisterMe(fHEXiZeroModel);
  pManager->AddDiscreteProcess(&fXiZeroInelastic);

  // map to G3 controls
  TG4ProcessControlMap* controlMap = TG4ProcessControlMap::Instance();
  controlMap->Add(&fXiZeroInelastic, kHADR); 

  // map to G3 ALIMCProcess codes
  TG4ProcessMCMap* mcMap = TG4ProcessMCMap::Instance();
  mcMap->Add(&fXiZeroInelastic, kPHInhelastic); 
}

//_____________________________________________________________________________
void TG4PhysicsConstructorHadron::ConstructProcessForAntiXiZero()
{
// Construct processes for anti-Xi0.
// ---

  // add process
  G4ProcessManager* pManager = G4AntiXiZero::AntiXiZero()->GetProcessManager();
  pManager->AddDiscreteProcess(&fElasticProcess);

  fLEAntiXiZeroModel = new G4LEAntiXiZeroInelastic();
  fHEAntiXiZeroModel = new G4HEAntiXiZeroInelastic();
  fAntiXiZeroInelastic.RegisterMe(fLEAntiXiZeroModel);
  fAntiXiZeroInelastic.RegisterMe(fHEAntiXiZeroModel);
  pManager->AddDiscreteProcess(&fAntiXiZeroInelastic);

  // map to G3 controls
  TG4ProcessControlMap* controlMap = TG4ProcessControlMap::Instance();
  controlMap->Add(&fAntiXiZeroInelastic, kHADR); 

  // map to G3 ALIMCProcess codes
  TG4ProcessMCMap* mcMap = TG4ProcessMCMap::Instance();
  mcMap->Add(&fAntiXiZeroInelastic, kPHInhelastic); 
}

//_____________________________________________________________________________
void TG4PhysicsConstructorHadron::ConstructProcessForOmegaMinus()
{
// Construct processes for Omega-.
// ---

  // add process & set ordering
  G4ProcessManager* pManager = G4OmegaMinus::OmegaMinus()->GetProcessManager();
  pManager->AddDiscreteProcess(&fElasticProcess);

  fLEOmegaMinusModel = new G4LEOmegaMinusInelastic();
  fHEOmegaMinusModel = new G4HEOmegaMinusInelastic();
  fOmegaMinusInelastic.RegisterMe(fLEOmegaMinusModel);
  fOmegaMinusInelastic.RegisterMe(fHEOmegaMinusModel);
  pManager->AddDiscreteProcess(&fOmegaMinusInelastic);

  pManager->AddProcess(&fOmegaMinusIonisation, ordInActive, 2, 2);

  pManager->AddProcess(&fOmegaMinusMult);
  pManager->SetProcessOrdering(&fOmegaMinusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(&fOmegaMinusMult, idxPostStep, 1);

  // map to G3 controls
  TG4ProcessControlMap* controlMap = TG4ProcessControlMap::Instance();
  controlMap->Add(&fOmegaMinusInelastic, kHADR); 
  controlMap->Add(&fOmegaMinusIonisation, kLOSS); 
  controlMap->Add(&fOmegaMinusMult, kMULS); 

  // map to G3 ALIMCProcess codes
  TG4ProcessMCMap* mcMap = TG4ProcessMCMap::Instance();
  mcMap->Add(&fOmegaMinusInelastic, kPHInhelastic); 
  mcMap->Add(&fOmegaMinusIonisation, kPEnergyLoss); 
  mcMap->Add(&fOmegaMinusMult, kPMultipleScattering); 
}

//_____________________________________________________________________________
void TG4PhysicsConstructorHadron::ConstructProcessForAntiOmegaMinus()
{
// Construct processes for pi+.
// ---

  // add process & set ordering
  G4ProcessManager* pManager 
    = G4AntiOmegaMinus::AntiOmegaMinus()->GetProcessManager();
  pManager->AddDiscreteProcess(&fElasticProcess);

  fLEAntiOmegaMinusModel = new G4LEAntiOmegaMinusInelastic();
  fHEAntiOmegaMinusModel = new G4HEAntiOmegaMinusInelastic();
  fAntiOmegaMinusInelastic.RegisterMe(fLEAntiOmegaMinusModel);
  fAntiOmegaMinusInelastic.RegisterMe(fHEAntiOmegaMinusModel);
  pManager->AddDiscreteProcess(&fAntiOmegaMinusInelastic);

  pManager->AddProcess(&fAntiOmegaMinusIonisation, ordInActive, 2, 2);

  pManager->AddProcess(&fAntiOmegaMinusMult);
  pManager->SetProcessOrdering(&fAntiOmegaMinusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(&fAntiOmegaMinusMult, idxPostStep, 1);

  // map to G3 controls
  TG4ProcessControlMap* controlMap = TG4ProcessControlMap::Instance();
  controlMap->Add(&fAntiOmegaMinusInelastic, kHADR); 
  controlMap->Add(&fAntiOmegaMinusIonisation, kLOSS); 
  controlMap->Add(&fAntiOmegaMinusMult, kMULS); 

  // map to G3 ALIMCProcess codes
  TG4ProcessMCMap* mcMap = TG4ProcessMCMap::Instance();
  mcMap->Add(&fAntiOmegaMinusInelastic, kPHInhelastic); 
  mcMap->Add(&fAntiOmegaMinusIonisation, kPEnergyLoss); 
  mcMap->Add(&fAntiOmegaMinusMult, kPMultipleScattering); 
}

//_____________________________________________________________________________
void TG4PhysicsConstructorHadron::ConstructProcessForOther()
{
// Construct processes for other hadrons.
// ---

  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pManager = particle->GetProcessManager();
    G4int nofAlongStepProcesses 
      = pManager->GetAlongStepProcessVector()->length();

    if ( !particle->IsShortLived()       &&
         particle->GetPDGCharge() != 0.0 &&
         nofAlongStepProcesses == 1      && 
         particle->GetParticleName() != "chargedgeantino") {
	
      // create processes
      G4VProcess* aMultipleScattering = new G4MultipleScattering();
      G4VProcess* anIonisation = new G4hIonisation();

      // add processes
      pManager->AddProcess(anIonisation, ordInActive, 2, 2);
      pManager->AddProcess(aMultipleScattering);

      // set ordering
      pManager->SetProcessOrdering(aMultipleScattering, idxAlongStep, 1);
      pManager->SetProcessOrdering(aMultipleScattering, idxPostStep, 1);

      // map to G3 controls
      TG4ProcessControlMap* controlMap = TG4ProcessControlMap::Instance();
      controlMap->Add(anIonisation, kLOSS); 
      controlMap->Add(aMultipleScattering, kMULS); 

      // map to AliMCProcess codes
      TG4ProcessMCMap* mcMap = TG4ProcessMCMap::Instance();
      mcMap->Add(anIonisation, kPEnergyLoss); 
      mcMap->Add(aMultipleScattering, kPMultipleScattering); 
      
      // keep for deleting 
      fOtherProcesses.push_back(aMultipleScattering);
      fOtherProcesses.push_back(anIonisation);
    }
  }    
}

// protected methods

//_____________________________________________________________________________
void TG4PhysicsConstructorHadron::ConstructParticle()
{
// Construct all hadrons.
// ---

  //  Construct all mesons
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  //  Construct all barions
  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  //  Construct  resonances and quarks
  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  
}

//_____________________________________________________________________________
void TG4PhysicsConstructorHadron::ConstructProcess()
{
// Construct all hadronic processes.
// ---

  // Elastic Process
  fElasticModel = new G4LElastic();
  fElasticProcess.RegisterMe(fElasticModel);
  
  ConstructProcessForPionPlus();
  ConstructProcessForPionMinus();
  ConstructProcessForKaonPlus();
  ConstructProcessForKaonMinus();
  ConstructProcessForKaonZeroLong();
  ConstructProcessForKaonZeroShort();
  ConstructProcessForProton();
  ConstructProcessForAntiProton();
  ConstructProcessForNeutron();
  ConstructProcessForAntiNeutron();
  ConstructProcessForLambda();
  ConstructProcessForAntiLambda();
  ConstructProcessForSigmaMinus();
  ConstructProcessForAntiSigmaMinus();
  ConstructProcessForSigmaPlus();
  ConstructProcessForAntiSigmaPlus();
  ConstructProcessForXiMinus();
  ConstructProcessForAntiXiMinus();
  ConstructProcessForXiZero();
  ConstructProcessForAntiXiZero();
  ConstructProcessForOmegaMinus();
  ConstructProcessForAntiOmegaMinus();
  ConstructProcessForOther();

  if (verboseLevel>0)
     G4cout << "### Hadron physics constructed." << G4endl;
}

