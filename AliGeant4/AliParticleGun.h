// $Id$
// Category: event
//
// Author: I. Hrivnacova
//
// Class AliParticleGun
// --------------------
// Particle gun that can be interactively composed by a user.

#ifndef ALI_PARTICLE_GUN_H
#define ALI_PARTICLE_GUN_H

#include <G4VPrimaryGenerator.hh>
#include <globals.hh>
#include <g4std/vector>

#include "AliParticleGunMessenger.h"

class AliGunParticle;

class G4Event;

class AliParticleGun : public G4VPrimaryGenerator
{
  typedef G4std::vector<AliGunParticle*>    GunParticleVector;
  typedef GunParticleVector::iterator       GunParticleIterator;
  typedef GunParticleVector::const_iterator GunParticleConstIterator;

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
     GunParticleVector        fGunParticleVector; //vector of AliGunParticle
     AliParticleGunMessenger  fMessenger;         //messenger
};

// inline methods

inline G4int AliParticleGun::GetNofGunParticles() const
{ return fGunParticleVector.size(); }

#endif //ALI_PARTICLE_GUN_H







