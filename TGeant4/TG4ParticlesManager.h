// $Id$
// Category: physics
//
// This class provides mapping between TDatabasePDG
// and Geant4 particles. 

#ifndef TG4_PARTICLES_MANAGER_H
#define TG4_PARTICLES_MANAGER_H

#include "TG4Globals.h"
#include "TG4NameMap.h"
#include "TG4IntMap.h"

#include <globals.hh>
#include <G4ThreeVector.hh>

#include <Rtypes.h>

class G4DynamicParticle;
class G4ParticleDefinition;

class TParticle;
class TClonesArray;

class TG4ParticlesManager
{
  public:
    TG4ParticlesManager();
    // --> protected
    // TG4ParticlesManager(const TG4ParticlesManager& right);
    virtual ~TG4ParticlesManager();

    // static access method
    static TG4ParticlesManager* Instance();
        
    // methods
    void MapParticles();

    // get methods
         // for G4 particle types   
    G4int GetPDGEncodingFast(G4ParticleDefinition* particle);

         // for Root particle types;
    TParticle* GetParticle(const TClonesArray* particles, G4int index) const;
    G4ParticleDefinition* GetParticleDefinition(
                           const TParticle* particle) const;
    G4DynamicParticle* CreateDynamicParticle(
                           const TParticle* particle) const;
    G4ThreeVector GetParticlePosition(
                           const TParticle* particle) const;
    G4ThreeVector GetParticleMomentum(
                           const TParticle* particle) const;        
    
  protected:
    TG4ParticlesManager(const TG4ParticlesManager& right);

    // operators
    TG4ParticlesManager& operator=(const TG4ParticlesManager& right);

  private:
    // methods
    G4int GetPDGEncoding(G4ParticleDefinition* particle);
    G4int GetPDGEncoding(G4String particleName);

    // static data members
    static TG4ParticlesManager*  fgInstance; //this instance
    
    // data members
    TG4NameMap  fParticleNameMap;  //the mapping between G4 particle names
                                   //and TDatabasePDG names 
    TG4IntMap   fParticlePDGMap;   //the mapping between G4 particle names
                                   //and TDatabasePDG codes
};

// inline methods

inline TG4ParticlesManager* TG4ParticlesManager::Instance() 
{ return fgInstance; }

#endif //TG4_PARTICLES_MANAGER_H

