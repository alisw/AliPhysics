// $Id$
// Category: event
//
// Messenger class that defines commands for AliStackingAction.

#ifndef ALI_STACKING_ACTION_MESSENGER_H
#define ALI_STACKING_ACTION_MESSENGER_H 

#include <G4UImessenger.hh>
#include <globals.hh>

class AliStackingAction;

class G4UIdirectory;
class G4UIcmdWithAnInteger;
class G4UIcmdWithoutParameter;

class AliStackingActionMessenger: public G4UImessenger
{
  public:
    AliStackingActionMessenger(AliStackingAction* trackingAction);
    // --> protected
    // AliStackingActionMessenger();
    // AliStackingActionMessenger(const AliStackingActionMessenger& right);
    virtual ~AliStackingActionMessenger();
   
    // methods 
    virtual void SetNewValue(G4UIcommand* command, G4String string);
    
  protected:
    AliStackingActionMessenger();
    AliStackingActionMessenger(const AliStackingActionMessenger& right);

    // operators
    AliStackingActionMessenger& operator=(
                            const AliStackingActionMessenger& right);

  private:
    // data members
    AliStackingAction*        fStackingAction;    //associated class 
    G4UIdirectory*            fStackingDirectory; //command directory
    G4UIcmdWithoutParameter*  fClearStackCmd;     //command: clearStack
    G4UIcmdWithAnInteger*     fVerboseCmd;        //command: verbose
};

#endif //ALI_STACKING_ACTION_MESSENGER_H
