// $Id$
// Category: event
//
// See the class description in the header file.

#include "AliParticleGun.h"
#include "AliParticleGunMessenger.h"
#include "AliGlobals.h"

#include <G4PrimaryParticle.hh>
#include <G4Event.hh>
#include <G4ios.hh>


AliParticleGun::AliParticleGun() {
//
  fMessenger = new AliParticleGunMessenger(this);
}

AliParticleGun::AliParticleGun(const AliParticleGun& right)
{
  // particles vector
  fGunParticlesVector.clearAndDestroy();
  for (G4int i=0; i<right.fGunParticlesVector.entries(); i++) {
    AliGunParticle* rhsParticle = right.fGunParticlesVector[i];
    AliGunParticle* newParticle = new AliGunParticle(*rhsParticle);
    fGunParticlesVector.insert(newParticle);
  }  

  fMessenger = new AliParticleGunMessenger(this);
}

AliParticleGun::~AliParticleGun() {
//
  fGunParticlesVector.clearAndDestroy();
  delete fMessenger;
}

// operators

AliParticleGun& AliParticleGun::operator=(const AliParticleGun& right)
{
  // check assignement to self
  if (this == &right) return *this;

  // particles vector
  fGunParticlesVector.clearAndDestroy();
  for (G4int i=0; i<right.fGunParticlesVector.entries(); i++) {
    AliGunParticle* rhsParticle = right.fGunParticlesVector[i];
    AliGunParticle* newParticle = new AliGunParticle(*rhsParticle);
    fGunParticlesVector.insert(newParticle);
  }
  
  return *this;  
}
  
// public methods

void AliParticleGun::AddParticle(AliGunParticle* particle)
{ 
// Adds particle.
// ---

  fGunParticlesVector.insert(particle); 
}

void AliParticleGun::RemoveParticle(G4int iParticle)
{ 
// Removes particle.
// ---

  AliGunParticle* particle = fGunParticlesVector[iParticle];
  fGunParticlesVector.remove(particle); 
  delete particle;  
}

void AliParticleGun::GeneratePrimaryVertex(G4Event* event)
{
// Generates primary vertices.
// ---

  G4PrimaryVertex* previousVertex = 0;
  G4ThreeVector previousPosition = G4ThreeVector(); 
  G4double previousTime = 0.; 

  G4int nofGunParticles = fGunParticlesVector.entries();
  for( G4int i=0; i<nofGunParticles; i++ )
  {    
    AliGunParticle* particle = fGunParticlesVector[i];

    G4ParticleDefinition* particleDefinition
      = particle->GetParticleDefinition();
    if (particleDefinition==0) {
      AliGlobals::Exception(
        "AliParticleGun::GeneratePrimaryVertex: Unknown particle definition.");
    }	

    G4ThreeVector position = particle->GetPosition(); 
    G4double time = particle->GetTime(); 
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
    G4double px = particle->GetMomentum().x();
    G4double py = particle->GetMomentum().y();
    G4double pz = particle->GetMomentum().z();
    G4PrimaryParticle* primaryParticle 
      = new G4PrimaryParticle(particleDefinition, px, py, pz);
      
    G4double mass =  particleDefinition->GetPDGMass();
    primaryParticle->SetMass(mass);
    primaryParticle->SetPolarization(particle->GetPolarization().x(),
                                     particle->GetPolarization().y(),
				     particle->GetPolarization().z());    
    
    vertex->SetPrimary(primaryParticle);
  }
  
  // delete gun particles
  fGunParticlesVector.clearAndDestroy();

  // add verbose
  G4cout << "AliParticleGun::GeneratePrimaryVertex:" << G4endl;
  G4cout << "   " 
         << event->GetNumberOfPrimaryVertex() << " of primary vertices,"
         << "   " << nofGunParticles << " of primary particles " << G4endl;  
}

void AliParticleGun::Reset()
{ 
// Resets the particle gun.
// ---

  fGunParticlesVector.clearAndDestroy(); 
}

void AliParticleGun::List()
{
// Lists the particle gun.
// ---

  G4int nofGunParticles = fGunParticlesVector.entries();

  G4cout << "Particle Gun: " << G4endl;

  if (nofGunParticles==0)
  { G4cout << "   No particles are defined." << G4endl; }
  else
  {
    for (G4int i=0; i<nofGunParticles; i++)
    {
      G4cout << i << " th particle properties: " << G4endl;
      G4cout << "============================" << G4endl;
      fGunParticlesVector[i]->Print();
    }
  }
}

