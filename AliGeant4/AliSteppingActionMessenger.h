// $Id$
// Category: event
//
// Messenger class that defines commands for AliSteppingAction.

#ifndef ALI_STEPPING_ACTION_MESSENGER_H
#define ALI_STEPPING_ACTION_MESSENGER_H 

#include <G4UImessenger.hh>
#include <globals.hh>

class AliSteppingAction;

class G4UIdirectory;
class G4UIcmdWithAnInteger;

class AliSteppingActionMessenger: public G4UImessenger
{
  public:
    AliSteppingActionMessenger(AliSteppingAction* trackingAction);
    // --> protected
    // AliSteppingActionMessenger();
    // AliSteppingActionMessenger(const AliSteppingActionMessenger& right);
    virtual ~AliSteppingActionMessenger();
   
    // methods 
    virtual void SetNewValue(G4UIcommand* command, G4String string);
    
  protected:
    AliSteppingActionMessenger();
    AliSteppingActionMessenger(const AliSteppingActionMessenger& right);

    // operators
    AliSteppingActionMessenger& operator=(
                            const AliSteppingActionMessenger& right);

  private:
    // data members
    AliSteppingAction*     fSteppingAction; //associated class  
    G4UIcmdWithAnInteger*  fLoopVerboseCmd; //command: loopVerbose
};

#endif //ALI_STEPPING_ACTION_MESSENGER_H
