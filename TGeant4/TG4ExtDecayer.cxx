// $Id$
// Category: physics
//
// Author: I. Hrivnacova
//
// Class TG4ExtDecayer
// -------------------
// See the class description in the header file.

#include "AliDecayer.h"

#include "TG4ExtDecayer.h"
#include "TG4ParticlesManager.h"
#include "TG4G3Units.h"

#include <G4DynamicParticle.hh>
#include <G4DecayProducts.hh>
#include <G4DecayTable.hh>
#include <G4ParticleTable.hh>
#include <G4DynamicParticle.hh>
#include <G4Track.hh>

#include <TParticle.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>

#include <math.h>

//_____________________________________________________________________________
TG4ExtDecayer::TG4ExtDecayer(AliDecayer* externalDecayer)
  : G4VExtDecayer("TG4ExtDecayer"),
    TG4Verbose("extDecayer"),
    fExternalDecayer(externalDecayer) {
//
  fDecayProductsArray = new  TClonesArray("TParticle", 1000);
  fParticlesManager = TG4ParticlesManager::Instance();
}

//_____________________________________________________________________________
TG4ExtDecayer::TG4ExtDecayer(const TG4ExtDecayer& right)
  : TG4Verbose("extDecayer") {
// 
  TG4Globals::Exception(
    "TG4ExtDecayer is protected from copying.");
}

//_____________________________________________________________________________
TG4ExtDecayer::~TG4ExtDecayer() {
//
  delete fDecayProductsArray;
}


// operators

//_____________________________________________________________________________
TG4ExtDecayer& TG4ExtDecayer::operator=(const TG4ExtDecayer& right)
{
  // check assignement to self
  if (this == &right) return *this;

  TG4Globals::Exception(
    "TG4ExtDecayer is protected from assigning.");
    
  return *this;  
} 

// public methods

//_____________________________________________________________________________
G4DecayProducts* TG4ExtDecayer::ImportDecayProducts(const G4Track& track)
{
  // check if external decayer is defined
  if (!fExternalDecayer) {
#ifdef G4VERBOSE
     G4cerr << "TG4ExtDecayer::ImportDecayProducts: " << G4endl
	    << " No fExternalDecayer is defined." << G4endl;
#endif
    return 0;
  }  
  
  // get particle momentum
  G4ThreeVector momentum = track.GetMomentum(); 
  G4double mag = momentum.mag();  
  TLorentzVector p;    
  p[0] = momentum.x()/mag;
  p[1] = momentum.x()/mag;
  p[2] = momentum.x()/mag;
  p[3] = mag/TG4G3Units::Energy();
  
  // get particle PDG
  // ask TG4ParticlesManager to get PDG encoding 
  // (in order to get PDG from extended TDatabasePDG
  // in case the standard PDG code is not defined)
  G4ParticleDefinition* particleDef = track.GetDefinition();
  G4int pdgEncoding = fParticlesManager->GetPDGEncodingFast(particleDef);

  // let AliDecayer decay the particle
  // and import the decay products
  fExternalDecayer->Decay(pdgEncoding, &p);
  G4int nofParticles
    = fExternalDecayer->ImportParticles(fDecayProductsArray);
  
  if (VerboseLevel()>1) {
    G4cout << "nofParticles: " <<  nofParticles << G4endl;
  }  

  // convert decay products TParticle type 
  // to G4DecayProducts  
  G4DecayProducts* decayProducts
    = new G4DecayProducts(*(track.GetDynamicParticle()));

  G4int counter = 0;
  for (G4int i=0; i<nofParticles; i++) {

    // get particle from TClonesArray
    TParticle* particle
      = fParticlesManager->GetParticle(fDecayProductsArray, i);
      
    G4int status = particle->GetStatusCode();
    G4int pdg = particle->GetPdgCode();
    if ( status>0 && status<11 && 
         abs(pdg)!=12 && abs(pdg)!=14 && abs(pdg)!=16 ) {
      // pass to tracking final particles only;
      // skip neutrinos

      if (VerboseLevel()>1) {
        G4cout << "  " << i << "th particle PDG: " << pdg << "   ";
      }
            
      // create G4DynamicParticle 
      G4DynamicParticle* dynamicParticle 
        = fParticlesManager->CreateDynamicParticle(particle);

      if (dynamicParticle) {

        if (VerboseLevel()>1) {
          G4cout << "  G4 particle name: " 
                 << dynamicParticle->GetDefinition()->GetParticleName()
	         << G4endl;
	}	 

        // add dynamicParticle to decayProducts
        decayProducts->PushProducts(dynamicParticle);
	
	counter++;
      }
    }       
  }			     
  if (VerboseLevel()>1) {
    G4cout << "nofParticles for tracking: " <<  counter << G4endl;
  }  
     
  return decayProducts;
}
    
