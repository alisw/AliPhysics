// $Id$
// Category: physics
//
// See the class description in the header file.

#include "TG4ParticlesManager.h"
#include "TG4G3Units.h"

#include <G4ParticleDefinition.hh>
#include <G4DynamicParticle.hh>
#include <G4ParticleTable.hh>

#include <TDatabasePDG.h>
#include <TParticle.h>
#include <TClonesArray.h>

TG4ParticlesManager* TG4ParticlesManager::fgInstance = 0;

TG4ParticlesManager::TG4ParticlesManager()
{ 
//
  if (fgInstance) {
    TG4Globals::Exception(
      "TG4ParticlesManager: attempt to create two instances of singleton.");
  }
      
  fgInstance = this;  
}

TG4ParticlesManager::TG4ParticlesManager(const TG4ParticlesManager& right) {
// 
  TG4Globals::Exception(
    "Attempt to copy TG4ParticlesManager singleton.");
}

TG4ParticlesManager::~TG4ParticlesManager() {
//
}

// operators

TG4ParticlesManager& 
TG4ParticlesManager::operator=(const TG4ParticlesManager& right)
{
  // check assignement to self
  if (this == &right) return *this;

  TG4Globals::Exception(
    "Attempt to assign TG4ParticlesManager singleton.");
    
  return *this;  
}    
          
// private methods

G4int TG4ParticlesManager::GetPDGEncoding(G4ParticleDefinition* particle)
{
// Returns the PDG code of particle;
// if standard PDG code is not defined the TDatabasePDG
// is used.
// ---

  // get PDG encoding from G4 particle definition
  G4int pdgEncoding = particle->GetPDGEncoding();

  if (pdgEncoding == 0) {
    // get PDG encoding from TDatabasePDG
  
    // get particle name from the name map
    G4String g4name = particle->GetParticleName();
    G4String tname = fParticleNameMap.GetSecond(g4name);
    if (tname == "Undefined") {
      G4String text = "TG4ParticlesManager::GetPDGEncoding: \n";
      text = text + "    Particle " + g4name;
      text = text + " was not found in the name map.";
      TG4Globals::Exception(text);
    }  
  
    // get particle from TDatabasePDG
    TDatabasePDG* pdgDB = TDatabasePDG::Instance();
    TParticlePDG* particle = pdgDB->GetParticle(tname);
    if (!particle) {
      G4String text = "TG4ParticlesManager::GetPDGEncoding: \n";
      text = text + "    Particle " + tname;
      text = text + " was not found in TDatabasePDG.";
      TG4Globals::Exception(text);
    }  
    
    // get PDG encoding
    pdgEncoding = particle->PdgCode();
  }
    
  return pdgEncoding;  
}  
     
G4int TG4ParticlesManager::GetPDGEncoding(G4String particleName)
{
// Returns the PDG code of particle sepcified by name.
// ---

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  G4ParticleDefinition* particle = 0;  
  particle = particleTable->FindParticle(particleName);
  if (!particle) {    
    G4String text = "TG4ParticlesManager::GetPDGEncoding:\n";
    text = text + "   G4ParticleTable::FindParticle() " + particleName;
    text = text + " failed.";
    TG4Globals::Exception(text);
  }	

  return GetPDGEncoding(particle);
}  
  
void  TG4ParticlesManager::MapParticles()
{
  // map G4 particle names to TDatabasePDG names
  // (the map is built only for particles that have not
  //  defined standard PDG encoding)
  
  fParticleNameMap.Add("deuteron","Deuteron");
  fParticleNameMap.Add("triton",  "Triton");
  fParticleNameMap.Add("alpha",   "Alpha");
  fParticleNameMap.Add("He3",     "HE3");
  fParticleNameMap.Add("opticalphoton","Cherenkov");
  // fParticleNameMap.Add("???","FeedbackPhoton");
  fParticleNameMap.Add("geantino", "Rootino");
  
  // map G4 particle names to TDatabasePDG encodings
  fParticlePDGMap.Add("deuteron", GetPDGEncoding("deuteron"));
  fParticlePDGMap.Add("triton", GetPDGEncoding("triton"));
  fParticlePDGMap.Add("alpha", GetPDGEncoding("alpha"));
  fParticlePDGMap.Add("He3", GetPDGEncoding("He3") );
  fParticlePDGMap.Add("opticalphoton", GetPDGEncoding("opticalphoton"));
  // fParticlePDGMap.Add("???","FeedbackPhoton");
  fParticleNameMap.Add("geantino", GetPDGEncoding("geantino"));

  // add verbose
  G4cout << "Particle maps have been filled." << G4endl;
  //fParticleNameMap.PrintAll();
  //fParticlePDGMap.PrintAll();
}    

// public methods

G4int TG4ParticlesManager::GetPDGEncodingFast(G4ParticleDefinition* particle)
{
// Returns the PDG code of particle;
// if standard PDG code is not defined the preregistred
// fParticlePDGMap is used.
// ---

  // get PDG encoding from G4 particle definition
  G4int pdgEncoding = particle->GetPDGEncoding();

  if (pdgEncoding == 0) {
    // use FParticlePDGMap if standard PDG code is not defined
    G4String name = particle->GetParticleName();
    pdgEncoding = fParticlePDGMap.GetSecond(name);
  }
    
  return pdgEncoding;  
}  
     
     
TParticle* TG4ParticlesManager::GetParticle(const TClonesArray* particles, 
                                          G4int index) const
{
// Retrives particle with given index from TClonesArray 
// and checks type.
// ---

  TObject* particleTObject = particles->UncheckedAt(index);      
  TParticle* particle
    = dynamic_cast<TParticle*>(particleTObject);

  // check particle type
  if (!particle) 
    TG4Globals::Exception(
      "TG4ParticlesManager::GetParticle: Unknown particle type");
  
  return particle;  
}     


G4ParticleDefinition* TG4ParticlesManager::GetParticleDefinition(
                               const TParticle* particle) const
{
// Returns G4 particle definition for given TParticle
// TO DO: replace with using particles name map    
// ---

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
      "TG4ParticlesManager::GetParticleDefinition:\n";
    text = text + "   G4ParticleTable::FindParticle() failed.";
    TG4Globals::Warning(text);
  }	
  
  return particleDefinition;
}

G4DynamicParticle* TG4ParticlesManager::CreateDynamicParticle(
                                   const TParticle* particle) const
{ 
// Creates G4DynamicParticle.
// ---

  // get particle properties
  G4ParticleDefinition* particleDefinition 
    = GetParticleDefinition(particle);    
  if (!particleDefinition) return 0;  
        
  G4ThreeVector momentum = GetParticleMomentum(particle);

  // create G4DynamicParticle
  G4DynamicParticle* dynamicParticle 
    = new G4DynamicParticle(particleDefinition, momentum);
  
  return dynamicParticle;
}

G4ThreeVector TG4ParticlesManager::GetParticlePosition(
                                   const TParticle* particle) const 
{
// Returns particle vertex position.
// ---

  G4ThreeVector position 
     = G4ThreeVector(particle->Vx()*TG4G3Units::Length(),
                     particle->Vy()*TG4G3Units::Length(),
		     particle->Vz()*TG4G3Units::Length());
  return position;
}  		     
			
			
G4ThreeVector TG4ParticlesManager::GetParticleMomentum(
                                   const TParticle* particle) const
{
// Returns particle momentum.
// ---
  G4ThreeVector momentum 
     = G4ThreeVector(particle->Px()*TG4G3Units::Energy(),
                     particle->Py()*TG4G3Units::Energy(),
		     particle->Pz()*TG4G3Units::Energy());
  return momentum;
}

