// $Id$
// Category: event
//
// Author: I. Hrivnacova
//
// Class AliParticleGun
// --------------------
// See the class description in the header file.

#include "AliParticleGun.h"
#include "AliParticleGunMessenger.h"
#include "AliGunParticle.h"
#include "AliGlobals.h"

#include <G4ParticleDefinition.hh>
#include <G4PrimaryParticle.hh>
#include <G4PrimaryVertex.hh>
#include <G4Event.hh>

//_____________________________________________________________________________
AliParticleGun::AliParticleGun() 
  : G4VPrimaryGenerator(),
    AliVerbose("particleGun"),
    fMessenger(this) {
//
}

//_____________________________________________________________________________
AliParticleGun::AliParticleGun(const AliParticleGun& right)
  : G4VPrimaryGenerator(right),
    AliVerbose("particleGun"),
    fMessenger(this)
{
  // copy stuff
  *this = right;
}

//_____________________________________________________________________________
AliParticleGun::~AliParticleGun() {
//
  // delete gun particles
  Reset();
}

// operators

//_____________________________________________________________________________
AliParticleGun& AliParticleGun::operator=(const AliParticleGun& right)
{
  // check assignement to self
  if (this == &right) return *this;

  // base class assignment
  this->G4VPrimaryGenerator::operator=(right);

  // delete gun particles
  Reset();

  // copy  gun particles
  GunParticleConstIterator it;
  for (it = (right.fGunParticleVector).begin();
       it != (right.fGunParticleVector).end(); it++)
    fGunParticleVector.push_back(new AliGunParticle(**it));

  return *this;  
}
  
// public methods

//_____________________________________________________________________________
void AliParticleGun::AddParticle(AliGunParticle* particle)
{ 
// Adds particle.
// ---

  fGunParticleVector.push_back(particle); 
}

//_____________________________________________________________________________
void AliParticleGun::RemoveParticle(G4int iParticle)
{ 
// Removes particle.
// ---

  GunParticleIterator it = fGunParticleVector.begin();
  it += iParticle;
  //advance(it,iParticle);
   
  delete *it;
  fGunParticleVector.erase(it); 
}

//_____________________________________________________________________________
void AliParticleGun::GeneratePrimaryVertex(G4Event* event)
{
// Generates primary vertices.
// ---

  G4PrimaryVertex* previousVertex = 0;
  G4ThreeVector previousPosition = G4ThreeVector(); 
  G4double previousTime = 0.; 

  GunParticleIterator it;
  for (it = fGunParticleVector.begin(); it != fGunParticleVector.end(); it++) {    

    AliGunParticle* particle = *it;

    G4ParticleDefinition* particleDefinition
      = particle->GetParticleDefinition();
    if (particleDefinition==0) {
      AliGlobals::Exception(
        "AliParticleGun::GeneratePrimaryVertex: Unknown particle definition.");
    }	

    G4ThreeVector position = particle->GetPosition(); 
    G4double time = particle->GetTime(); 
    G4PrimaryVertex* vertex;
    if (it == fGunParticleVector.begin() || 
        position != previousPosition || time != previousTime ) {   
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
  
  if (VerboseLevel() > 0) {
    G4cout << "AliParticleGun::GeneratePrimaryVertex:" << G4endl;
    G4cout << "   " 
           << event->GetNumberOfPrimaryVertex()  << " of primary vertices,"
           << "   " << fGunParticleVector.size() << " of primary particles " 
  	   << G4endl;  
  }	   

  // delete gun particles
  Reset();
}

//_____________________________________________________________________________
void AliParticleGun::Reset()
{ 
// Resets the particle gun.
// ---

  GunParticleIterator it;
  for (it = fGunParticleVector.begin(); it != fGunParticleVector.end(); it++)
    delete *it;

  fGunParticleVector.clear(); 
}

//_____________________________________________________________________________
void AliParticleGun::List()
{
// Lists the particle gun.
// ---

  G4cout << "Particle Gun: " << G4endl;

  if (fGunParticleVector.size() == 0) 
    G4cout << "   No particles are defined." << G4endl; 
  else {
    G4int i = 0;
    GunParticleIterator it;
    for (it = fGunParticleVector.begin(); 
         it != fGunParticleVector.end(); it++) {    
	 
      G4cout << i++ << " th particle properties: " << G4endl;
      G4cout << "============================" << G4endl;
      (*it)->Print();
    }
  }
}

