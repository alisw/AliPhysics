// $Id$
// Category: run
//
// Class that defines primary generator action. 
// Available primary generators (AliPrimaryGenerator):
//  kGun,               // gun (can be set interactively) 
//  kGeantino,          // geantino with random direction
//  kChargedGeantino,   // chargedgeantino with random direction
//  kAliGenerator       // AliGenerator from AliRoot

#ifndef ALI_PRIMARY_GENERATOR_ACTION_H
#define ALI_PRIMARY_GENERATOR_ACTION_H

#include "AliPrimaryGenerator.h"

#include <G4VUserPrimaryGeneratorAction.hh>
#include <globals.hh>

class AliParticleGun;
class AliPrimaryGeneratorMessenger;
class G4ParticleGun;
class G4Event;
class TClonesArray;

class AliPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    AliPrimaryGeneratorAction();
    // --> protected
    // AliPrimaryGeneratorAction(const AliPrimaryGeneratorAction& right);
    virtual ~AliPrimaryGeneratorAction();

    // methods
    virtual void GeneratePrimaries(G4Event* event);
    
    // set methods
    void SetGenerator(AliPrimaryGenerator generator);
    void SetNofGunParticles(G4int nofParticles);
    void SetVerboseLevel(G4int level);

    // get methods
    AliPrimaryGenerator GetGenerator() const;
    G4int GetNofGunParticles() const;
    G4int GetVerboseLevel() const;
    
  protected:
    AliPrimaryGeneratorAction(const AliPrimaryGeneratorAction& right);

    // operators
    AliPrimaryGeneratorAction& operator=(
                              const AliPrimaryGeneratorAction& right);

  private:
    // methods
    void ConstructGenerator();
    void ConstructGeantinoGenerator(G4bool isCharged);
    void ConstructAliGenerator();
    void GenerateAliGeneratorPrimaries(G4Event* event);
    
    // data members
    AliPrimaryGenerator  fGenerator;       //selected AliPrimaryGenerator
    G4int                fNofGunParticles; //number of gun particles
    G4int                fVerboseLevel;    //verbose level
    AliParticleGun*      fParticleGun;     //AliParticleGun
    TClonesArray*        fParticleArray;   //AliRun::fParticles
    AliPrimaryGeneratorMessenger*  fMessenger; //messenger
};

// inline methods

inline void AliPrimaryGeneratorAction::SetVerboseLevel(G4int level)
{ fVerboseLevel = level; }

inline AliPrimaryGenerator AliPrimaryGeneratorAction::GetGenerator() const
{ return fGenerator; }

inline G4int AliPrimaryGeneratorAction::GetNofGunParticles() const
{ return fNofGunParticles; }

inline G4int AliPrimaryGeneratorAction::GetVerboseLevel() const
{ return fVerboseLevel; }

#endif //ALI_PRIMARY_GENERATOR_ACTION_H


