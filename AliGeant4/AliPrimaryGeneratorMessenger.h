// $Id$
// Category: run
//
// Messenger class that defines commands for AliPrimaryGeneratorAction.

#ifndef ALI_PRIMARY_GENERATOR_MESSENGER_H
#define ALI_PRIMARY_GENERATOR_MESSENGER_H

#include <G4UImessenger.hh>
#include <globals.hh>

class AliPrimaryGeneratorAction;

class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;

class AliPrimaryGeneratorMessenger: public G4UImessenger
{
  public:
    AliPrimaryGeneratorMessenger(AliPrimaryGeneratorAction* primaryGenAction);
    // --> protected
    // AliPrimaryGeneratorMessenger();
    // AliPrimaryGeneratorMessenger(const AliPrimaryGeneratorMessenger& right);
    virtual ~AliPrimaryGeneratorMessenger();
    
    // methods
    void SetNewValue(G4UIcommand* command, G4String string);
    
  protected:
    AliPrimaryGeneratorMessenger();
    AliPrimaryGeneratorMessenger(const AliPrimaryGeneratorMessenger& right);

    // operators
    AliPrimaryGeneratorMessenger& operator=(
                                 const AliPrimaryGeneratorMessenger& right);

  private:
    // data members
    AliPrimaryGeneratorAction*  fPrimaryGenAction;   //associated class
    G4UIdirectory*              fPrimariesDirectory; //command directory
    G4UIcmdWithAString*         fGeneratorCmd;       //command: set
    G4UIcmdWithAnInteger*       fNofParticlesCmd;    //command: nofParticles
    G4UIcmdWithAnInteger*       fVerboseCmd;         //command: verbose
};

#endif //ALI_PRIMARY_GENERATOR_MESSENGER_H

