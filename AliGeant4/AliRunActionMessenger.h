// $Id$
// Category: run
//
// Messenger class that defines commands for AliRunAction.

#ifndef ALI_RUN_ACTION_MESSENGER_H
#define ALI_RUN_ACTION_MESSENGER_H 

#include <G4UImessenger.hh>
#include <globals.hh>

class AliRunAction;

class G4UIdirectory;
class G4UIcmdWithAnInteger;

class AliRunActionMessenger: public G4UImessenger
{
  public:
    AliRunActionMessenger(AliRunAction* eventAction);
    // --> protected
    // AliRunActionMessenger();
    // AliRunActionMessenger(const AliRunActionMessenger& right);
    virtual ~AliRunActionMessenger();
   
    // methods 
    virtual void SetNewValue(G4UIcommand* command, G4String string);
    
  protected:
    AliRunActionMessenger();
    AliRunActionMessenger(const AliRunActionMessenger& right);

    // operators
    AliRunActionMessenger& operator=(
                            const AliRunActionMessenger& right);

  private:
    // data members
    AliRunAction*          fRunAction;          //associated class
    G4UIdirectory*         fRunActionDirectory; //command directory
    G4UIcmdWithAnInteger*  fVerboseCmd;         //command: verbose 
};

#endif //ALI_RUN_ACTION_MESSENGER_H
