// $Id$
// Category: event
//
// Author: I. Hrivnacova
//
// Class AliTrackingActionMessenger
// --------------------------------
// Messenger class that defines commands for AliTrackingAction.
 
#ifndef ALI_TRACKING_ACTION_MESSENGER_H
#define ALI_TRACKING_ACTION_MESSENGER_H 

#include <G4UImessenger.hh>
#include <globals.hh>

class AliTrackingAction;

class G4UIdirectory;
class G4UIcmdWithAnInteger;

class AliTrackingActionMessenger: public G4UImessenger
{
  public:
    AliTrackingActionMessenger(AliTrackingAction* trackingAction);
    // --> protected
    // AliTrackingActionMessenger();
    // AliTrackingActionMessenger(const AliTrackingActionMessenger& right);
    virtual ~AliTrackingActionMessenger();
   
    // methods 
    virtual void SetNewValue(G4UIcommand* command, G4String string);
    
  protected:
    AliTrackingActionMessenger();
    AliTrackingActionMessenger(const AliTrackingActionMessenger& right);

    // operators
    AliTrackingActionMessenger& operator=(
                            const AliTrackingActionMessenger& right);

  private:
    // data members
    AliTrackingAction*     fTrackingAction;    //associated class 
    G4UIdirectory*         fTrackingDirectory; //command directory
    G4UIcmdWithAnInteger*  fVerboseCmd;        //command: verbose
    G4UIcmdWithAnInteger*  fNewVerboseCmd;     //command: newVerbose
    G4UIcmdWithAnInteger*  fNewVerboseTrackCmd;//command: newVerboseTrack
};

#endif //ALI_TRACKING_ACTION_MESSENGER_H
