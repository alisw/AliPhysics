// $Id$
// Category: run
//
// Messenger class that defines commands for AliRun.

#ifndef ALI_RUN_MESSENGER_H
#define ALI_RUN_MESSENGER_H 

#include <G4UImessenger.hh>
#include <globals.hh>

class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAnInteger;

class AliRunMessenger: public G4UImessenger
{
  public:
    AliRunMessenger();
    // --> protected
    // AliRunMessenger(const AliRunMessenger& right);
    virtual ~AliRunMessenger();
   
    // methods 
    virtual void SetNewValue(G4UIcommand* command, G4String string);
    
  protected:
    AliRunMessenger(const AliRunMessenger& right);

    // operators
    AliRunMessenger& operator=(const AliRunMessenger& right);

  private:
    // data members
    G4UIdirectory*            fRunDirectory;  //command directory
    G4UIcmdWithoutParameter*  fInitializeCmd; //command: initialize
    G4UIcmdWithAnInteger*     fBeamOnCmd;     //command: beamOn
    G4UIcmdWithoutParameter*  fLegoCmd;       //command: lego
};

#endif //ALI_RUN_MESSENGER_H
