// $Id$
// Category: run
//
// Author: I. Hrivnacova
//
// Class AliPrimaryGeneratorAction
// -------------------------------
// See the class description in the header file.

#include "AliPrimaryGeneratorAction.h"
#include "AliPrimaryGeneratorMessenger.h"
#include "AliParticleGun.h"
#include "AliGunParticle.h"
#include "AliGlobals.h"

#include "TG4G3Units.h"
#include "TG4PrimaryGeneratorAction.h"
#include "TG4TrackingAction.h"

#include <G4Event.hh>
#include <G4ParticleTable.hh>
#include <G4ParticleDefinition.hh>

#include <Randomize.hh>

#include <TParticle.h>
#include <TVirtualMCApplication.h>

//_____________________________________________________________________________
AliPrimaryGeneratorAction::AliPrimaryGeneratorAction()
  : AliVerbose("primaryGeneratorAction"),
    fGenerator(kStack),
    fNofGunParticles(1),
    fParticleGun(),
    fMessenger(this) {
//
}

//_____________________________________________________________________________
AliPrimaryGeneratorAction::~AliPrimaryGeneratorAction() {
//
}

// private methods

//_____________________________________________________________________________
void AliPrimaryGeneratorAction::ConstructGenerator()
{
// Constructs selected generator.
// ---

  switch (fGenerator) { 
    case kGun:
      // gun is constructed interactively      
      return;
    case kGeantino: 
      ConstructGeantinoGenerator(false);
      return;
    case kChargedGeantino:  
      ConstructGeantinoGenerator(true);
      return;
    case kStack:
      return;
  }
}   
      
//_____________________________________________________________________________
void AliPrimaryGeneratorAction::ConstructGeantinoGenerator(G4bool isCharged)
{
// Geantino with random momentum direction
// (the default generator).
// ---

  // reset gun
  fParticleGun.Reset();     

  G4ParticleTable* particleTable 
    = G4ParticleTable::GetParticleTable();

  for (G4int i=0; i< fNofGunParticles; i++)
  {
    G4ParticleDefinition* particleDef = 0; 
    if (!isCharged)
      particleDef = particleTable->FindParticle("geantino");
    else  
      particleDef = particleTable->FindParticle("chargedgeantino");

    if (!particleDef) {
      G4String text = "AliPrimaryGeneratorAction::GenerateGeantino:\n";
      text = text + "   G4ParticleTable::FindParticle() failed.";
      AliGlobals::Exception(text);
    }  
      
    G4double rn[3];
    RandFlat::shootArray(3,rn);
    G4double px=rn[0];
    G4double py=rn[1];
    G4double pz=rn[2];
    G4ThreeVector momentumDir(px, py, pz);

    G4double energy = 1.*GeV;
    G4ThreeVector position(0.,0.,0.);
    G4double time = 0.;
    G4ThreeVector polarization(0.,0.,0.);
    
    AliGunParticle * gunParticle
      = new AliGunParticle(particleDef, momentumDir, energy, position, time, 
              polarization);

    fParticleGun.AddParticle(gunParticle);     
  } 
  if (VerboseLevel() > 1) { 
    G4cout << "Geantino generator has been built." << G4endl; 
  }
} 
            
// public methods

//_____________________________________________________________________________
void AliPrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
// Generates primary particles by the selected generator.
// ---

  if (fGenerator == kStack)  {

    // Use MC stack
    TG4PrimaryGeneratorAction action;
    action.GeneratePrimaries(event);

    // Do not save primary particles
    // (they would be stored twice)
    TG4TrackingAction* trackingAction
      =  TG4TrackingAction::Instance();
    if (trackingAction) trackingAction->SetSavePrimaries(false);
  }  
  else {
  
    // Begin of Event
    TVirtualMCApplication::Instance()->BeginEvent();

    // Construct particle gun
    ConstructGenerator();
  
    // Generate primary vertices
    fParticleGun.GeneratePrimaryVertex(event);

    // Primary particles have to be saved in stack
    TG4TrackingAction* trackingAction
      =  TG4TrackingAction::Instance();
    if (trackingAction) trackingAction->SetSavePrimaries(true);
  }  
}

//_____________________________________________________________________________
void AliPrimaryGeneratorAction::SetGenerator(AliPrimaryGenerator generator)
{ 
// Sets generator.
// ---

  fGenerator = generator; 
}

//_____________________________________________________________________________
void AliPrimaryGeneratorAction::SetNofGunParticles(G4int nofParticles)
{ 
// Sets number of primary particles.
// This method is applied only to "gun" type generators
// (and not to AliGenerator from AliRoot).
// ---

  fNofGunParticles = nofParticles;
}
