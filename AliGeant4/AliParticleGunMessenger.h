// $Id$
// Category: event
//
// Messenger class that defines commands for AliParticleGun.

#ifndef ALI_PARTICLE_GUN_MESSENGER_H
#define ALI_PARTICLE_GUN_MESSENGER_H

#include <G4UImessenger.hh>
#include <globals.hh>

class AliParticleGun;
class AliGunParticle;
class G4ParticleTable;
class G4UIcommand;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger;
class G4UIcmdWith3Vector;
class G4UIcmdWith3VectorAndUnit;

class AliParticleGunMessenger: public G4UImessenger
{
  public:
    AliParticleGunMessenger(AliParticleGun* gun);
    // --> protected
    // AliParticleGunMessenger();
    // AliParticleGunMessenger(const AliParticleGunMessenger& right);
    virtual ~AliParticleGunMessenger();

    // methods
    virtual void SetNewValue(G4UIcommand* command, G4String newValues);
    virtual G4String GetCurrentValue(G4UIcommand* command);

  protected:
    AliParticleGunMessenger();
    AliParticleGunMessenger(const AliParticleGunMessenger& right);

    // operators
    AliParticleGunMessenger& operator=(
                            const AliParticleGunMessenger& right);

  private:
    // data members
    AliParticleGun*   fGun;           //associated class
    AliGunParticle*   fParticle;      //current AliGunParticle
    G4ParticleTable*  fParticleTable; //G4ParticleTable
 
    // commands data members
    G4UIdirectory*              fGunDirectory;      //command directory
    G4UIcmdWithoutParameter*    fListAvailableCmd;  //command: listAvailable
    G4UIcmdWithoutParameter*    fListCurrentCmd;    //command: listCurrent
    G4UIcmdWithAString*         fParticleCmd;       //command: particle
    G4UIcmdWith3VectorAndUnit*  fMomentumCmd;       //command: momentum
    G4UIcmdWith3VectorAndUnit*  fPositionCmd;       //command: position
    G4UIcmdWithADoubleAndUnit*  fTimeCmd;           //command: time
    G4UIcmdWith3Vector*         fPolarizationCmd;   //command: polarization
    G4UIcmdWith3Vector*         fDirectionCmd;      //command: direction 
    G4UIcmdWithADoubleAndUnit*  fKinEnergyCmd;      //command: kinEnergy
    G4UIcmdWithoutParameter*    fListCmd;           //command: list
    G4UIcmdWithoutParameter*    fAddParticleCmd;    //command: addParticle
    G4UIcmdWithAnInteger*       fRemoveParticleCmd; //command: removeParticle
    G4UIcmdWithoutParameter*    fResetCmd;          //command: reset
};

#endif //ALI_PARTICLE_GUN_MESSENGER_H


