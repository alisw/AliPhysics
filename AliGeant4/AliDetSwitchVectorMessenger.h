// $Id$
// Category: geometry
//
// Author: I. Hrivnacova
//
// Class AliDetSwitchVectorMessenger
// ------------------------------------
// Messenger class that defines commands for AliDetSwitchVector.

#ifndef ALI_DET_SWITCH_VECTOR_MESSENGER_H
#define ALI_DET_SWITCH_VECTOR_MESSENGER_H

#include <G4UImessenger.hh>
#include <globals.hh>

class AliDetSwitchVector;

class G4UIcommand;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithoutParameter;

class AliDetSwitchVectorMessenger: public G4UImessenger
{
  public:
    AliDetSwitchVectorMessenger(AliDetSwitchVector* detSwitchVector);
    // --> protected
    // AliDetSwitchVectorMessenger();
    // AliDetSwitchVectorMessenger(const AliDetSwitchVectorMessenger& right);
    virtual ~AliDetSwitchVectorMessenger();

    // methods
    virtual void SetNewValue(G4UIcommand* command, G4String newValues);
    void Update();
    
  protected:
    AliDetSwitchVectorMessenger();
    AliDetSwitchVectorMessenger(const AliDetSwitchVectorMessenger& right);

    // operators
    AliDetSwitchVectorMessenger& 
    operator=(const AliDetSwitchVectorMessenger &right);
             
  private:
    // methods
    void SetGuidance();
    void SetCandidates();

    AliDetSwitchVector*  fDetSwitchVector; //associated class
    
    // commands data members
    G4UIcmdWithAString*         fSwitchOnCmd;     //command: switchOn 
    G4UIcmdWithAString*         fSwitchOffCmd;    //command: switchOn
    G4UIcmdWithoutParameter*    fListCmd;         //command: list
    G4UIcmdWithoutParameter*    fListAvailableCmd;//command: listAvailable
};

#endif //ALI_DET_SWITCH_VECTOR_MESSENGER_H

