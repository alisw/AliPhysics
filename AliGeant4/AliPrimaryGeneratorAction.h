// $Id$
// Category: run
//
// Author: I. Hrivnacova
//
// Class AliPrimaryGeneratorAction
// -------------------------------
// Class that defines primary generator action. 
// Available primary generators (AliPrimaryGenerator):
//  kGun,               // gun (can be set interactively) 
//  kGeantino,          // geantino with random direction
//  kChargedGeantino,   // chargedgeantino with random direction
//  kStack              // AliGenerator from the MC stack

#ifndef ALI_PRIMARY_GENERATOR_ACTION_H
#define ALI_PRIMARY_GENERATOR_ACTION_H

#include "AliVerbose.h"
#include "AliPrimaryGenerator.h"
#include "AliPrimaryGeneratorMessenger.h"
#include "AliParticleGun.h"

#include <G4VUserPrimaryGeneratorAction.hh>
#include <globals.hh>

class AliParticleGun;
class G4ParticleGun;
class G4Event;

class AliPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction,
                                  public AliVerbose
{
  public:
    AliPrimaryGeneratorAction();
    virtual ~AliPrimaryGeneratorAction();

    // methods
    virtual void GeneratePrimaries(G4Event* event);
    
    // set methods
    void SetGenerator(AliPrimaryGenerator generator);
    void SetNofGunParticles(G4int nofParticles);

    // get methods
    AliPrimaryGenerator GetGenerator() const;
    G4int GetNofGunParticles() const;
    
  private:
    // methods
    void ConstructGenerator();
    void ConstructGeantinoGenerator(G4bool isCharged);
    void ConstructAliGenerator();
    void GenerateAliGeneratorPrimaries(G4Event* event);
    
    // data members
    AliPrimaryGenerator  fGenerator;       //selected AliPrimaryGenerator
    G4int                fNofGunParticles; //number of gun particles
    AliParticleGun       fParticleGun;     //AliParticleGun
    AliPrimaryGeneratorMessenger  fMessenger; //messenger
};

// inline methods

inline AliPrimaryGenerator AliPrimaryGeneratorAction::GetGenerator() const
{ return fGenerator; }

inline G4int AliPrimaryGeneratorAction::GetNofGunParticles() const
{ return fNofGunParticles; }

#endif //ALI_PRIMARY_GENERATOR_ACTION_H


