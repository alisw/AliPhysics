// $Id$
// Category: physics
//
// Messenger class that defines commands for TG4PhysicsList

#ifndef TG4_PHYSICS_LIST_MESSENGER_H
#define TG4_PHYSICS_LIST_MESSENGER_H 

#include <G4UImessenger.hh>
#include <globals.hh>

class TG4PhysicsList;

class G4UIcmdWithABool;

class TG4PhysicsListMessenger: public G4UImessenger
{
  public:
    TG4PhysicsListMessenger(TG4PhysicsList* physicsList); 
    // --> protected                       
    // TG4PhysicsListMessenger(const TG4PhysicsListMessenger& right);
    virtual ~TG4PhysicsListMessenger();
   
    // methods 
    virtual void SetNewValue(G4UIcommand* command, G4String string);
    
  protected:
    TG4PhysicsListMessenger(const TG4PhysicsListMessenger& right);

    // operators
    TG4PhysicsListMessenger& operator=(
                            const TG4PhysicsListMessenger& right);

  private:
    // data members
    TG4PhysicsList*    fPhysicsList;        //associated class 
    G4UIcmdWithABool*  fSetOpticalCmd;      //setCerenkov command   
    G4UIcmdWithABool*  fSetHadronCmd;       //setHadron command   
    G4UIcmdWithABool*  fSetSpecialCutsCmd;  //setSpecialCuts command   
    G4UIcmdWithABool*  fSetSpecialFlagsCmd; //setSpecialFlags command   
};

#endif //TG4_PHYSICS_LIST_MESSENGER_H
