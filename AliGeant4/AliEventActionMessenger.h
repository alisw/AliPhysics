// $Id$
// Category: event
//
// Messenger class that defines commands for AliEventAction.

#ifndef ALI_EVENT_ACTION_MESSENGER_H
#define ALI_EVENT_ACTION_MESSENGER_H 

#include <G4UImessenger.hh>
#include <globals.hh>

class AliEventAction;

class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;

class AliEventActionMessenger: public G4UImessenger
{
  public:
    AliEventActionMessenger(AliEventAction* eventAction);
    // --> protected
    // AliEventActionMessenger();
    // AliEventActionMessenger(const AliEventActionMessenger& right);
    virtual ~AliEventActionMessenger();
   
    // methods 
    virtual void SetNewValue(G4UIcommand* command, G4String string);
    
  protected:
    AliEventActionMessenger();
    AliEventActionMessenger(const AliEventActionMessenger& right);

    // operators
    AliEventActionMessenger& operator=(
                            const AliEventActionMessenger& right);

  private:
    // data members
    AliEventAction*        fEventAction;    //associated class
    G4UIdirectory*         fEventDirectory; //command directory
    G4UIcmdWithAString*    fDrawTracksCmd;  //command: drawTracks
    G4UIcmdWithAnInteger*  fVerboseCmd;     //command: verbose
};

#endif //ALI_EVENT_ACTION_MESSENGER_H
