// $Id$
// Category: run
//
// See the class description in the header file.

#include "AliPrimaryGeneratorAction.h"
#include "AliPrimaryGeneratorMessenger.h"
#include "AliTrackingAction.h"
#include "AliParticleGun.h"
#include "AliGunParticle.h"
#include "AliGlobals.h"
#include "AliRun.h"
#include "AliGenerator.h"

#include "TG3Units.h"

#include <G4Event.hh>
#include <G4ParticleTable.hh>
#include <G4ParticleDefinition.hh>

#include <Randomize.hh>

#include <TParticle.h>

AliPrimaryGeneratorAction::AliPrimaryGeneratorAction()
  : fGenerator(kAliGenerator),
    fNofGunParticles(1),
    fVerboseLevel(0)
{
//
  fParticleGun = new AliParticleGun();
  fMessenger = new AliPrimaryGeneratorMessenger(this);
}

AliPrimaryGeneratorAction::AliPrimaryGeneratorAction(
                                    const AliPrimaryGeneratorAction& right) {
//				    
  AliGlobals::Exception(
    "AliPrimaryGeneratorAction is protected from copying.");
}

AliPrimaryGeneratorAction::~AliPrimaryGeneratorAction() {
//
  delete fMessenger;
  delete fParticleGun;
}

// operators

AliPrimaryGeneratorAction& 
AliPrimaryGeneratorAction::operator=(const AliPrimaryGeneratorAction &right)
{
  // check assignement to self
  if (this == &right) return *this;
  
  AliGlobals::Exception(
    "AliPrimaryGeneratorAction is protected from assigning.");

  return *this;
}

// private methods

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
    case kAliGenerator:
      ConstructAliGenerator();
      return;
  }
}   
      
void AliPrimaryGeneratorAction::ConstructGeantinoGenerator(G4bool isCharged)
{
// Geantino with random momentum direction
// (the default generator).
// ---

  // reset gun
  fParticleGun->Reset();     

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

    fParticleGun->AddParticle(gunParticle);     
  } 
  if (fVerboseLevel>1) { 
    G4cout << "Geantino generator has been built." << G4endl; 
  }
} 
            
void AliPrimaryGeneratorAction::ConstructAliGenerator()
{
// Generator from AliRoot
// AliGenerator::Generate() fills the AliRun::fParticles array.
// ---

  // check if AliGenerator is set
  AliGenerator* generator = gAlice->Generator();
  if (!generator) {
    G4String text = "AliPrimaryGeneratorAction::ConstructGenerator:\n";
    text = text + "   No AliGenerator is defined in gAlice.";
    AliGlobals::Exception(text);
  }  
  // fill AliRun::fParticleMap array 
  generator->Generate();
}

void AliPrimaryGeneratorAction::GenerateAliGeneratorPrimaries(G4Event* event)
{
// Creates a new G4PrimaryVertex objects for each TParticle
// in fParticles array.
// ---

  G4PrimaryVertex* previousVertex = 0;
  G4ThreeVector previousPosition = G4ThreeVector(); 
  G4double previousTime = 0.; 

  G4int nofParticles = gAlice->GetNtrack();
  // add verbose
  //G4cout << " nofParticles: " <<  nofParticles << G4endl;
  for( G4int i=0; i<nofParticles; i++ ) {    
  
    // get the particle from AliRun stack
    TParticle* particle = gAlice->Particle(i);

    // get particle definition from G4ParticleTable
    G4int pdgEncoding = particle->GetPdgCode();
    G4ParticleTable* particleTable 
      = G4ParticleTable::GetParticleTable();                
    G4ParticleDefinition* particleDefinition = 0;    
    if (pdgEncoding != 0) 
      particleDefinition = particleTable->FindParticle(pdgEncoding);
    else {
      G4String name = particle->GetName();
      if (name == "Rootino")	
        particleDefinition = particleTable->FindParticle("geantino");
    }	
  
    if (particleDefinition==0) {
      G4cout << "pdgEncoding: " << pdgEncoding << G4endl;
      G4String text = 
        "AliPrimaryGeneratorAction::GenerateAliGeneratorPrimaries:\n";
      text = text + "   G4ParticleTable::FindParticle() failed.";
      AliGlobals::Exception(text);
    }	

    // get/create vertex
    G4ThreeVector position 
      = G4ThreeVector(particle->Vx()*TG3Units::Length(),
                      particle->Vy()*TG3Units::Length(),
		      particle->Vz()*TG3Units::Length());
    G4double time = particle->T()*TG3Units::Time(); 
    G4PrimaryVertex* vertex;
    if ( i==0 || position != previousPosition || time != previousTime ) {   
      // create a new vertex 
      // in case position and time of gun particle are different from 
      // previous values
      // (vertex objects are destroyed in G4EventManager::ProcessOneEvent()
      // when event is deleted)  
      vertex = new G4PrimaryVertex(position, time);
      event->AddPrimaryVertex(vertex);

      previousVertex = vertex;
      previousPosition = position;
      previousTime = time;
    }
    else 
      vertex = previousVertex;

    // create a primary particle and add it to the vertex
    // (primaryParticle objects are destroyed in G4EventManager::ProcessOneEvent()
    // when event and then vertex is deleted)
    G4double px = particle->Px()*TG3Units::Energy();
    G4double py = particle->Py()*TG3Units::Energy();
    G4double pz = particle->Pz()*TG3Units::Energy();
    G4PrimaryParticle* primaryParticle 
      = new G4PrimaryParticle(particleDefinition, px, py, pz);
      
    // set polarization
    TVector3 polarization;
    particle->GetPolarisation(polarization);
    primaryParticle
      ->SetPolarization(polarization.X(), polarization.Y(), polarization.Z());    
    
    // add primary particle to the vertex
    vertex->SetPrimary(primaryParticle);
  }
}

// public methods

void AliPrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
// Generates primary particles by the selected generator.
// ---

  // reset AliRun
  gAlice->BeginEvent();

  // fill particle gun/particle array
  ConstructGenerator();
  
  // generate primary vertices
  if (fGenerator == kAliGenerator)  {
    // use AliGenerator if set
    GenerateAliGeneratorPrimaries(event);

    // do not save primary particles
    // (they would be stored twice)
    AliTrackingAction* trackingAction
      =  AliTrackingAction::Instance();
    trackingAction->SetSavePrimaries(false);
  }  
  else {
    // use particle gun otherwise
    fParticleGun->GeneratePrimaryVertex(event);

    // primary particles have to be saved
    AliTrackingAction* trackingAction
      =  AliTrackingAction::Instance();
    trackingAction->SetSavePrimaries(true);
  }  
}

void AliPrimaryGeneratorAction::SetGenerator(AliPrimaryGenerator generator)
{ 
// Sets generator.
// ---

  fGenerator = generator; 
}

void AliPrimaryGeneratorAction::SetNofGunParticles(G4int nofParticles)
{ 
// Sets number of primary particles.
// This method is applied only to "gun" type generators
// (and not to AliGenerator from AliRoot).
// ---

  fNofGunParticles = nofParticles;
}






