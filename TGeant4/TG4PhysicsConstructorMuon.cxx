// $Id$
// Category: physics
//
// Author: I. Hrivnacova
//
// Class TG4PhysicsConstructorMuon
// -----------------------------
// See the class description in the header file.
// According to ExN04MuonPhysics.cc,v 1.2.2.1 2001/06/28 19:07:37 gunter Exp 
// GEANT4 tag Name: geant4-03-02

#include "TG4PhysicsConstructorMuon.h"
#include "TG4ProcessControlMap.h"
#include "TG4ProcessMCMap.h"
#include "TG4G3Control.h"

#include <G4ParticleDefinition.hh>
#include <G4ProcessManager.hh>
#include <G4MuonPlus.hh>
#include <G4MuonMinus.hh>
#include <G4TauMinus.hh>
#include <G4TauPlus.hh>
#include <G4NeutrinoTau.hh>
#include <G4AntiNeutrinoTau.hh>
#include <G4NeutrinoMu.hh>
#include <G4AntiNeutrinoMu.hh>

//TBR
#include <G4MuIonisation.hh>
#include <G4MuBremsstrahlung.hh>
#include <G4MuPairProduction.hh>
#include <G4hIonisation.hh>

//_____________________________________________________________________________
TG4PhysicsConstructorMuon::TG4PhysicsConstructorMuon(const G4String& name)
  : G4VPhysicsConstructor(name)
{
//
  SetVerboseLevel(1);
}

//_____________________________________________________________________________
TG4PhysicsConstructorMuon::~TG4PhysicsConstructorMuon() {
//
}

// private methods

//_____________________________________________________________________________
void TG4PhysicsConstructorMuon::ConstructProcessForMuonPlus()
{
// Constructs electromagnetic processes for mu+.
// ---
  
  // add processes
  G4ProcessManager* pManager = G4MuonPlus::MuonPlus()->GetProcessManager();
  pManager->AddProcess(&fMuPlusIonisation, ordInActive,2, 2);
  pManager->AddDiscreteProcess(&fMuPlusBremsstrahlung);
  pManager->AddDiscreteProcess(&fMuPlusPairProduction);
  pManager->AddProcess(&fMuPlusMultipleScattering);

  // set ordering
  pManager->SetProcessOrdering(&fMuPlusMultipleScattering, idxAlongStep,  1);
  pManager->SetProcessOrdering(&fMuPlusMultipleScattering, idxPostStep,  1);

  // map to G3 controls
  TG4ProcessControlMap* controlMap = TG4ProcessControlMap::Instance();
  controlMap->Add(&fMuPlusIonisation, kLOSS); 
  controlMap->Add(&fMuPlusBremsstrahlung, kBREM); 
  controlMap->Add(&fMuPlusPairProduction, kPAIR); 
  controlMap->Add(&fMuPlusMultipleScattering, kMULS); 


  // map to AliMCProcess codes
  TG4ProcessMCMap* mcMap = TG4ProcessMCMap::Instance();
  mcMap->Add(&fMuPlusIonisation, kPEnergyLoss); 
  mcMap->Add(&fMuPlusBremsstrahlung, kPBrem); 
  mcMap->Add(&fMuPlusPairProduction, kPPair); 
  mcMap->Add(&fMuPlusMultipleScattering, kPMultipleScattering); 
}  

//_____________________________________________________________________________
void TG4PhysicsConstructorMuon::ConstructProcessForMuonMinus()
{
// Constructs electromagnetic processes for mu-.
// ---
  
  // add processes & set ordering
  G4ProcessManager* pManager = G4MuonMinus::MuonMinus()->GetProcessManager();
  pManager->AddProcess(&fMuMinusIonisation, ordInActive,2, 2);
  pManager->AddDiscreteProcess(&fMuMinusBremsstrahlung);
  pManager->AddDiscreteProcess(&fMuMinusPairProduction);
  pManager->AddProcess(&fMuMinusMultipleScattering);

  pManager->SetProcessOrdering(&fMuMinusMultipleScattering, idxAlongStep,  1);
  pManager->SetProcessOrdering(&fMuMinusMultipleScattering, idxPostStep,  1);

  pManager->AddRestProcess(&fMuMinusCaptureAtRest);

  // map to G3 controls
  TG4ProcessControlMap* controlMap = TG4ProcessControlMap::Instance();
  controlMap->Add(&fMuMinusIonisation, kLOSS); 
  controlMap->Add(&fMuMinusBremsstrahlung, kBREM); 
  controlMap->Add(&fMuMinusPairProduction, kPAIR); 
  controlMap->Add(&fMuMinusMultipleScattering, kMULS); 
  controlMap->Add(&fMuMinusCaptureAtRest, kMUNU); 

  // map to AliMCProcess codes
  TG4ProcessMCMap* mcMap = TG4ProcessMCMap::Instance();
  mcMap->Add(&fMuMinusIonisation, kPEnergyLoss); 
  mcMap->Add(&fMuMinusBremsstrahlung, kPBrem); 
  mcMap->Add(&fMuMinusPairProduction, kPPair); 
  mcMap->Add(&fMuMinusMultipleScattering, kPMultipleScattering); 
  mcMap->Add(&fMuMinusCaptureAtRest, kPMuonNuclear); 
}  

//_____________________________________________________________________________
void TG4PhysicsConstructorMuon::ConstructProcessForTauPlus()
{
// Constructs electromagnetic processes for tau+.
// ---
  
  // add processes
  G4ProcessManager* pManager = G4TauPlus::TauPlus()->GetProcessManager();
  pManager->AddProcess(&fTauPlusIonisation, ordInActive,2, 2);
  pManager->AddProcess(&fTauPlusMultipleScattering);

  // set ordering
  pManager->SetProcessOrdering(&fTauPlusMultipleScattering, idxAlongStep,  1);
  pManager->SetProcessOrdering(&fTauPlusMultipleScattering, idxPostStep,  1);

  // map to G3 controls
  TG4ProcessControlMap* controlMap = TG4ProcessControlMap::Instance();
  controlMap->Add(&fTauPlusIonisation, kLOSS); 
  controlMap->Add(&fTauPlusMultipleScattering, kMULS); 


  // map to AliMCProcess codes
  TG4ProcessMCMap* mcMap = TG4ProcessMCMap::Instance();
  mcMap->Add(&fTauPlusIonisation, kPEnergyLoss); 
  mcMap->Add(&fTauPlusMultipleScattering, kPMultipleScattering); 
}  

//_____________________________________________________________________________
void TG4PhysicsConstructorMuon::ConstructProcessForTauMinus()
{
// Constructs electromagnetic processes for tau-.
// ---
  
  // add processes
  G4ProcessManager* pManager = G4TauMinus::TauMinus()->GetProcessManager();
  pManager->AddProcess(&fTauMinusIonisation, ordInActive,2, 2);
  pManager->AddProcess(&fTauMinusMultipleScattering);

  // set ordering
  pManager->SetProcessOrdering(&fTauMinusMultipleScattering, idxAlongStep,  1);
  pManager->SetProcessOrdering(&fTauMinusMultipleScattering, idxPostStep,  1);

  // map to G3 controls
  TG4ProcessControlMap* controlMap = TG4ProcessControlMap::Instance();
  controlMap->Add(&fTauMinusIonisation, kLOSS); 
  controlMap->Add(&fTauMinusMultipleScattering, kMULS); 


  // map to AliMCProcess codes
  TG4ProcessMCMap* mcMap = TG4ProcessMCMap::Instance();
  mcMap->Add(&fTauMinusIonisation, kPEnergyLoss); 
  mcMap->Add(&fTauMinusMultipleScattering, kPMultipleScattering); 
}  


// protected methods

//_____________________________________________________________________________
void TG4PhysicsConstructorMuon::ConstructParticle()
{
// Instantiates particles.
// ---

  // Mu
  G4MuonPlus::MuonPlusDefinition();
  G4MuonMinus::MuonMinusDefinition();
  G4NeutrinoMu::NeutrinoMuDefinition();
  G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();

  // Tau
  G4TauMinus::TauMinusDefinition();
  G4TauPlus::TauPlusDefinition();
  G4NeutrinoTau::NeutrinoTauDefinition();
  G4AntiNeutrinoTau::AntiNeutrinoTauDefinition();
}

//_____________________________________________________________________________
void TG4PhysicsConstructorMuon::ConstructProcess()
{
// Constructs electromagnetic processes for muons.
// ---

  ConstructProcessForMuonPlus();
  ConstructProcessForMuonMinus();
  ConstructProcessForTauPlus();
  ConstructProcessForTauMinus();

  if (verboseLevel>0)
    G4cout << "### Muon physics constructed." << G4endl;
}
