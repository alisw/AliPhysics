// $Id$
// Category: run
//
// Messenger class that defines commands for 
// geometry, physics and step managers

#ifndef TG4_MESSENGER_H
#define TG4_MESSENGER_H 

#include <G4UImessenger.hh>
#include <globals.hh>

class TG4GeometryManager;
class TG4PhysicsManager;
class TG4StepManager;

class G4UIcmdWithoutParameter;
class G4UIcmdWithABool;

class TG4Messenger: public G4UImessenger
{
  public:
    TG4Messenger(TG4GeometryManager* geometryManager, 
                 TG4PhysicsManager* physicsManager, 
		 TG4StepManager* stepManager);
    // --> protected   
    // TG4Messenger();
    // TG4Messenger(const TG4Messenger& right);
    virtual ~TG4Messenger();
   
    // methods 
    virtual void SetNewValue(G4UIcommand* command, G4String string);
    
  protected:
    TG4Messenger();  
    TG4Messenger(const TG4Messenger& right);

    // operators
    TG4Messenger& operator=(const TG4Messenger& right);

  private:
    // data members
    TG4GeometryManager*       fGeometryManager; //geometry manager
    TG4PhysicsManager*        fPhysicsManager;  //physics manager
    TG4StepManager*           fStepManager;     //step manager
    
    G4UIcmdWithABool*  fSetEMCmd;             //setEM command   
    G4UIcmdWithABool*  fSetOpticalCmd;        //setOptical command   
    G4UIcmdWithABool*  fSetHadronCmd;         //setHadron command   
    G4UIcmdWithABool*  fSetSpecialCutsCmd;    //setSpecialCuts command   
    G4UIcmdWithABool*  fSetSpecialControlsCmd;//setSpecialControls command   
    G4UIcmdWithoutParameter*  fProcessActivationCmd; //.
                                              //setProcessActivation command    
};

#endif //TG4_MESSENGER_H
