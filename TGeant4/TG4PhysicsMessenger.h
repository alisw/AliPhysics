// $Id$
// Category: physics
//
// Author: I. Hrivnacova
//
// Class TG4PhysicsMessenger
// -------------------------
// Messenger class that defines commands for the physics manager.

#ifndef TG4_PHYSICS_MESSENGER_H
#define TG4_PHYSICS_MESSENGER_H 

#include <G4UImessenger.hh>
#include <globals.hh>

class TG4PhysicsManager;

class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithABool;
class G4UIcmdWithAString;

class TG4PhysicsMessenger: public G4UImessenger
{
  public:
    TG4PhysicsMessenger(TG4PhysicsManager* physicsManager); 
    // --> protected   
    // TG4PhysicsMessenger();
    // TG4PhysicsMessenger(const TG4PhysicsMessenger& right);
    virtual ~TG4PhysicsMessenger();
   
    // methods 
    virtual void SetNewValue(G4UIcommand* command, G4String string);
    
  protected:
    TG4PhysicsMessenger();  
    TG4PhysicsMessenger(const TG4PhysicsMessenger& right);

    // operators
    TG4PhysicsMessenger& operator=(const TG4PhysicsMessenger& right);

  private:
    // data members
    TG4PhysicsManager* fPhysicsManager;       //associated class
    G4UIdirectory*     fDirectory;            //command directory
    
    G4UIcmdWithABool*  fSetEMCmd;             //setEM command   
    G4UIcmdWithABool*  fSetMuonCmd;           //setMuon command   
    G4UIcmdWithABool*  fSetHadronCmd;         //setHadron command   
    G4UIcmdWithABool*  fSetOpticalCmd;        //setOptical command   
    G4UIcmdWithABool*  fSetSpecialCutsCmd;    //setSpecialCuts command   
    G4UIcmdWithABool*  fSetSpecialControlsCmd;//setSpecialControls command   
    G4UIcmdWithoutParameter*  fProcessActivationCmd; //.
                                              //setProcessActivation command    
    G4UIcmdWithoutParameter*  fPrintProcessMCMapCmd; //.
                                              //printProcessMCMap command
    G4UIcmdWithoutParameter*  fPrintProcessControlMapCmd; //.
                                              //printProcessControlsMap command
    G4UIcmdWithAString*       fPrintVolumeLimitsCmd; //.
                                              //printVolumeLimits command
    G4UIcmdWithoutParameter*  fPrintGeneralCutsCmd; //.
                                              //printGeneralCuts command
    G4UIcmdWithoutParameter*  fPrintGeneralControlsCmd; //.
                                              //printGeneralControls command
};

#endif //TG4_PHYSICS_MESSENGER_H
