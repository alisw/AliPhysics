// $Id$
// Category: geometry
//
// Messenger class that defines command directory for each
// AliModuleConstruction instance.

#ifndef ALI_MODULE_CONSTRUCTION_MESSENGER_H
#define ALI_MODULE_CONSTRUCTION_MESSENGER_H

#include <G4UImessenger.hh>
#include <globals.hh>

class AliModuleConstruction;

class G4UIdirectory;
class G4UIcmdWithABool;
class G4UIcmdWithAString;
class G4UIcmdWithoutParameter;
class G4LogicalVolume;

class AliModuleConstructionMessenger: public G4UImessenger
{
  public:
    AliModuleConstructionMessenger(
       AliModuleConstruction* moduleConstruction, G4String moduleName);
    // --> protected   
    // AliModuleConstructionMessenger();
    // AliModuleConstructionMessenger(
    //       const AliModuleConstructionMessenger& right);
    virtual ~AliModuleConstructionMessenger();

    // methods
    virtual void SetNewValue(G4UIcommand* command, G4String newValues);

  protected:
    AliModuleConstructionMessenger();
    AliModuleConstructionMessenger(
             const AliModuleConstructionMessenger& right);

    // operators
    AliModuleConstructionMessenger& operator=(
             const AliModuleConstructionMessenger &right);
             
  private:
    // data members
    AliModuleConstruction*    fModuleConstruction;   //associated class
    G4UIdirectory*            fModuleDirectory;      //command directory
    G4UIcmdWithAString*       fSetFrameCmd;          //command: setFrame
    G4UIcmdWithoutParameter*  fListCmd;              //command: list
    G4UIcmdWithoutParameter*  fListLongCmd;          //command: listLong
    G4UIcmdWithAString*       fListDaughtersCmd;     //command: listDaughters
    G4UIcmdWithAString*       fListLongDaughtersCmd; //command: listLongDaughters  
    
#ifdef ALICE_VISUALIZE
    // commands data members
    G4LogicalVolume*          fCurrentVolume;          //current logical volume
    G4UIcmdWithAString*       fSetCurrentLVCmd;        //command: setVolume
    G4UIcmdWithABool*         fSetDetVisibilityCmd;    //command: setDetVisibility
    G4UIcmdWithABool*         fSetLVTreeVisibilityCmd; //command: setLVTreeVisibility   
    G4UIcmdWithABool*         fSetVolVisibilityCmd;    //command: setVolVisibility
    G4UIcmdWithAString*       fSetDetColourCmd;        //command: setDetColour
    G4UIcmdWithAString*       fSetLVTreeColourCmd;     //command: setLVTreeColour  
    G4UIcmdWithAString*       fSetVolColourCmd;        //command: setVolColour
#endif //ALICE_VISUALIZE
};

#endif //ALI_MODULE_CONSTRUCTION_MESSENGER_H

