// $Id$
// Category: geometry
//
// Author: I. Hrivnacova
//
// Class AliLVTreeMessenger
// ------------------------------------
// Messenger class that defines commands for AliLVTree.

#ifndef ALI_LV_TREE_MESSENGER_H
#define ALI_LV_TREE_MESSENGER_H

#include <G4UImessenger.hh>
#include <globals.hh>

class AliLVTree;

class G4UIdirectory;
class G4UIcmdWithABool;
class G4UIcmdWithAString;
class G4UIcmdWithoutParameter;
class G4LogicalVolume;

class AliLVTreeMessenger: public G4UImessenger
{
  public:
    AliLVTreeMessenger(AliLVTree* lvTree);
    // --> protected   
    // AliLVTreeMessenger();
    // AliLVTreeMessenger(const AliLVTreeMessenger& right);
    virtual ~AliLVTreeMessenger();

    // methods
    virtual void SetNewValue(G4UIcommand* command, G4String newValues);

  protected:
    AliLVTreeMessenger();
    AliLVTreeMessenger(const AliLVTreeMessenger& right);

    // operators
    AliLVTreeMessenger& operator=(const AliLVTreeMessenger &right);
             
  private:
    // data members
    AliLVTree*                fLVTree;               //associated class
    G4LogicalVolume*          fCurrentVolume;        //current logical volume
    G4UIdirectory*            fDirectory;            //command directory
    G4UIcmdWithAString*       fSetCurrentLVCmd;      //command: setVolume
    G4UIcmdWithoutParameter*  fListCmd;              //command: list
    G4UIcmdWithoutParameter*  fListLongCmd;          //command: listLong
    G4UIcmdWithAString*       fListDaughtersCmd;     //command: listDaughters
    G4UIcmdWithAString*       fListLongDaughtersCmd; //command: listLongDaughters  

#ifdef G4VIS_USE
    G4UIcmdWithABool*         fSetLVTreeVisibilityCmd; //command: setLVTreeVisibility   
    G4UIcmdWithABool*         fSetVolVisibilityCmd;    //command: setVolVisibility
    G4UIcmdWithAString*       fSetLVTreeColourCmd;     //command: setLVTreeColour  
    G4UIcmdWithAString*       fSetVolColourCmd;        //command: setVolColour
#endif //G4VIS_USE
};

#endif //ALI_LV_TREE_MESSENGER_H

