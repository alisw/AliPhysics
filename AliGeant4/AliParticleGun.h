// $Id$
// Category: event
//
// Particle gun that can be interactively composed by a user.

#ifndef ALI_PARTICLE_GUN_H
#define ALI_PARTICLE_GUN_H

#include <G4VPrimaryGenerator.hh>
#include <globals.hh>

#include <g4rw/tpordvec.h>

class AliGunParticle;
class AliParticleGunMessenger;

class G4Event;

class AliParticleGun : public G4VPrimaryGenerator
{
  typedef G4RWTPtrOrderedVector<AliGunParticle>  AliGunParticleVector;

  public:
     AliParticleGun();
     AliParticleGun(const AliParticleGun &right);
     virtual ~AliParticleGun();

     // operators
     AliParticleGun& operator=(const AliParticleGun& right);

     // methods
     void AddParticle(AliGunParticle* particle);
     void RemoveParticle(G4int iParticle);
     virtual void GeneratePrimaryVertex(G4Event* evt);
     void Reset();
     void List();

     // get methods
     G4int GetNofGunParticles() const;
  
  private:
     // data members
     AliGunParticleVector      fGunParticlesVector; //vector of AliGunParticle
     AliParticleGunMessenger*  fMessenger;          //messenger
};

// inline methods

inline G4int AliParticleGun::GetNofGunParticles() const
{ return fGunParticlesVector.entries(); }

#endif //ALI_PARTICLE_GUN_H







