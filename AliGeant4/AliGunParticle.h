// $Id$
// Category: event
//
// Data type class that stores properties of a gun particle.
// Used in AliParticleGun.

#ifndef ALI_GUN_PARTICLE_H
#define ALI_GUN_PARTICLE_H

#include <G4ParticleMomentum.hh>
#include <G4ThreeVector.hh>
#include <globals.hh>

class G4ParticleDefinition;

class AliGunParticle
{
  public:
    AliGunParticle();
    AliGunParticle(G4ParticleDefinition* particleDef, G4ParticleMomentum momentum, 
       G4ThreeVector position, G4double time, G4ThreeVector polarization );
    AliGunParticle( G4ParticleDefinition* particleDef, G4ParticleMomentum momentumDir, 
      G4double kinEnergy, G4ThreeVector position, G4double time, 
      G4ThreeVector polarization );
    AliGunParticle(const AliGunParticle& right);
    ~AliGunParticle();     

    // operators
    AliGunParticle& operator=(const AliGunParticle& right);
    G4int operator==(const AliGunParticle& right) const;
    G4int operator!=(const AliGunParticle& right) const;

    // methods
    void Print() const;  

    // set methods
    void SetParticleDefinition(G4ParticleDefinition* particleDef);
    void SetMomentum(G4ParticleMomentum  momentum);
    void SetPosition(G4ThreeVector position);
    void SetTime(G4double time);
    void SetPolarization(G4ThreeVector polarization);
    void SetMomentumDirection(G4ParticleMomentum  momentumDir);
    void SetKineticEnergy(G4double kinEnergy);
      
    // get methods  
    G4ParticleDefinition* GetParticleDefinition() const;
    G4ParticleMomentum GetMomentum() const;
    G4ThreeVector GetPosition() const;
    G4double GetTime() const;
    G4ThreeVector GetPolarization() const;
    G4ParticleMomentum GetMomentumDirection() const;
    G4double GetKineticEnergy() const;
                         
  private:
    // data members
    G4ParticleDefinition*  fParticleDefinition; //G4ParticleDefinition
    G4ParticleMomentum     fParticleMomentum;   //G4ParticleMomentum
    G4ThreeVector          fPosition;           //position
    G4double               fTime;               //time
    G4ThreeVector          fPolarization;       //polarization
};

// inline methods

inline void AliGunParticle::SetParticleDefinition(G4ParticleDefinition* particleDef)
{ fParticleDefinition = particleDef; }

inline void AliGunParticle::SetMomentum(G4ParticleMomentum  momentum)
{ fParticleMomentum = momentum; }

inline void AliGunParticle::SetPosition(G4ThreeVector position)
{ fPosition = position; }

inline void AliGunParticle::SetTime(G4double time)
{ fTime = time; }

inline void AliGunParticle::SetPolarization(G4ThreeVector polarization)
{ fPolarization = polarization; }

inline void AliGunParticle::SetMomentumDirection(G4ParticleMomentum  momentumDir)
{ fParticleMomentum = fParticleMomentum.mag()*momentumDir.unit(); }

inline G4ParticleDefinition* AliGunParticle::GetParticleDefinition() const
{ return fParticleDefinition; }

inline G4ParticleMomentum AliGunParticle::GetMomentum() const
{ return fParticleMomentum; }

inline G4ThreeVector AliGunParticle::GetPosition() const
{ return fPosition; }

inline G4double AliGunParticle::GetTime() const
{ return fTime; }

inline G4ThreeVector AliGunParticle::GetPolarization() const
{ return fPolarization; }

inline G4ParticleMomentum AliGunParticle::GetMomentumDirection() const
{ return fParticleMomentum.unit(); }

#endif //ALI_GUN_PARTICLE_H   
   

