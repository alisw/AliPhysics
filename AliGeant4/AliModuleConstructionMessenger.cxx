// $Id$
// Category: geometry
//
// See the class description in the header file.

#include "AliModuleConstructionMessenger.h"
#include "AliModuleConstruction.h"
#include "AliGlobals.h"
#ifdef ALICE_VISUALIZE
#include "AliColourStore.h"
#endif

#include <G4UIdirectory.hh>
#include <G4UIcmdWithABool.hh>
#include <G4UIcmdWithAString.hh>
#include <G4UIcmdWithoutParameter.hh>

AliModuleConstructionMessenger::AliModuleConstructionMessenger(
   AliModuleConstruction* moduleConstruction, G4String moduleName)
 : fModuleConstruction(moduleConstruction)
{
//
  G4String dirName = "/aliDet/"; 
  dirName = dirName + moduleName + "/"; 
  fModuleDirectory = new G4UIdirectory(dirName);
  G4String guidance = "AlSubDetConstruction ";
  guidance = guidance + moduleName + " control commands.";
  fModuleDirectory->SetGuidance(guidance);

  G4String commandPath = dirName + "setFrame";
  fSetFrameCmd= new G4UIcmdWithAString(commandPath, this);
  fSetFrameCmd ->SetGuidance("Set detector frame");
  fSetFrameCmd->SetParameterName("frameName", false);
  fSetFrameCmd->AvailableForStates(PreInit, Idle);  
 
  commandPath = dirName + "list";
  fListCmd = new G4UIcmdWithoutParameter(commandPath, this);
  guidance = "List all logical volumes of ";
  guidance = guidance + moduleName + " detector";
  fListCmd->SetGuidance(guidance);
  fListCmd->AvailableForStates(PreInit,Idle);  
 
  commandPath = dirName + "listLong";
  fListLongCmd = new G4UIcmdWithoutParameter(commandPath, this);
  fListLongCmd
    ->SetGuidance("List all logical volumes and number of its physical volumes");
  guidance = "of " + moduleName + " detector";
  fListLongCmd->SetGuidance(guidance);
  fListLongCmd->AvailableForStates(PreInit,Idle);  
 
  commandPath = dirName + "listDaughters";
  fListDaughtersCmd = new G4UIcmdWithAString(commandPath, this);
  fListDaughtersCmd->SetGuidance("List daughters of the given logical volumes");
  fListDaughtersCmd->SetParameterName("lvName", false);
  fListDaughtersCmd->AvailableForStates(PreInit,Idle);  
 
  commandPath = dirName + "listLongDaughters";
  fListLongDaughtersCmd = new G4UIcmdWithAString(commandPath, this);
  fListLongDaughtersCmd
    ->SetGuidance("List daughters of the given logical volumes");
  fListLongDaughtersCmd->SetGuidance("and number of its physical volumes");
  fListLongDaughtersCmd->SetParameterName("lvName", false);
  fListLongDaughtersCmd->AvailableForStates(PreInit,Idle);  
 
#ifdef ALICE_VISUALIZE
  fCurrentVolume = 0;

  commandPath = dirName + "setVolume";
  fSetCurrentLVCmd = new G4UIcmdWithAString(commandPath, this);
  fSetCurrentLVCmd->SetGuidance("Set the current logical volume.");
  fSetCurrentLVCmd->SetParameterName("curVolume", false);
  fSetCurrentLVCmd->AvailableForStates(PreInit,Idle);  
 
  commandPath = dirName + "setVisibility";
  fSetDetVisibilityCmd = new G4UIcmdWithABool(commandPath, this);
  guidance = "Make ";
  guidance = guidance + moduleName + " detector visible/invisible.";
  fSetDetVisibilityCmd->SetGuidance(guidance);
  fSetDetVisibilityCmd->SetParameterName("detVisibility", false);
  fSetDetVisibilityCmd->AvailableForStates(PreInit,Idle);  
 
  commandPath = dirName + "setLVTreeVisibility";
  fSetLVTreeVisibilityCmd = new G4UIcmdWithABool(commandPath, this);
  fSetLVTreeVisibilityCmd 
    ->SetGuidance("Make current volume tree visible/invisible.");
  fSetLVTreeVisibilityCmd->SetParameterName("volVisibility", false);
  fSetLVTreeVisibilityCmd->AvailableForStates(PreInit,Idle);  
 
  commandPath = dirName + "setVolVisibility";
  fSetVolVisibilityCmd = new G4UIcmdWithABool(commandPath, this);
  fSetVolVisibilityCmd 
    ->SetGuidance("Make current volume visible/invisible.");
  fSetVolVisibilityCmd->SetParameterName("volVisibility", false);
  fSetVolVisibilityCmd->AvailableForStates(PreInit,Idle);  
 
  commandPath = dirName + "setColour";
  fSetDetColourCmd = new G4UIcmdWithAString(commandPath, this);
  AliColourStore* pColours = AliColourStore::Instance();
  guidance = "Set colour for all ";
  guidance = guidance + moduleName + " detector volumes.";
  fSetDetColourCmd->SetGuidance(guidance);
  fSetDetColourCmd->SetGuidance("Available colours:");
  fSetDetColourCmd->SetGuidance(pColours->GetColoursListWithCommas());
  fSetDetColourCmd->SetParameterName("detColour", false);
  G4String candidatesList = pColours->GetColoursList();  
  fSetDetColourCmd->SetCandidates(candidatesList);
  fSetDetColourCmd->AvailableForStates(PreInit,Idle);  

  commandPath = dirName + "setLVTreeColour";
  fSetLVTreeColourCmd = new G4UIcmdWithAString(commandPath, this);
  fSetLVTreeColourCmd->SetGuidance("Set colour for the current volume tree.");
  fSetLVTreeColourCmd->SetGuidance("Available colours:");
  fSetLVTreeColourCmd->SetGuidance(pColours->GetColoursListWithCommas());
  fSetLVTreeColourCmd->SetParameterName("volColour", false);
  candidatesList = pColours->GetColoursList();  
  fSetLVTreeColourCmd->SetCandidates(candidatesList);
  fSetLVTreeColourCmd->AvailableForStates(PreInit,Idle);  

  commandPath = dirName + "setVolColour";
  fSetVolColourCmd = new G4UIcmdWithAString(commandPath, this);
  fSetVolColourCmd->SetGuidance("Set colour for the current volume.");
  fSetVolColourCmd->SetGuidance("Available colours:");
  fSetVolColourCmd->SetGuidance(pColours->GetColoursListWithCommas());
  fSetVolColourCmd->SetParameterName("volColour", false);
  candidatesList = pColours->GetColoursList();  
  fSetVolColourCmd->SetCandidates(candidatesList);
  fSetVolColourCmd->AvailableForStates(PreInit,Idle);  
#endif //ALICE_VISUALIZE
}

AliModuleConstructionMessenger::AliModuleConstructionMessenger() {
//
}

AliModuleConstructionMessenger::AliModuleConstructionMessenger(
                                const AliModuleConstructionMessenger& right)
{
//
  AliGlobals::Exception(
    "AliModuleConstructionMessenger is protected from copying.");
}

AliModuleConstructionMessenger::~AliModuleConstructionMessenger()
{
//
  delete fModuleDirectory;
  delete fSetFrameCmd;  
  delete fListCmd;
  delete fListLongCmd;
  delete fListDaughtersCmd;
  delete fListLongDaughtersCmd;
#ifdef ALICE_VISUALIZE
  delete fSetDetVisibilityCmd;
  delete fSetDetColourCmd;
  delete fSetCurrentLVCmd;
  delete fSetVolVisibilityCmd;    
  delete fSetVolColourCmd;    
  delete fSetLVTreeVisibilityCmd;    
  delete fSetLVTreeColourCmd;    
#endif //ALICE_VISUALIZE
}

// operators

AliModuleConstructionMessenger& 
AliModuleConstructionMessenger::operator=(
                                const AliModuleConstructionMessenger& right)
{
  // check assignement to self
  if (this == &right) return *this;

  AliGlobals::Exception(
     "AliModuleConstructionMessenger is protected from assigning.");
    
  return *this;  
}    
          
// public methods
  
void AliModuleConstructionMessenger::SetNewValue(G4UIcommand* command,
                                                 G4String newValues)
{
// Applies command to the associated object.
// ---

  if (command == fSetFrameCmd) {
    fModuleConstruction->SetDetFrame(newValues);
  }
  else if (command == fListCmd) {
    fModuleConstruction->ListAllLVTree();
  }
  else if (command == fListLongCmd) {
    fModuleConstruction->ListAllLVTreeLong();
  }
  else if (command == fListDaughtersCmd) {
    fModuleConstruction->ListLVTree(newValues);
  }
  else if (command == fListLongDaughtersCmd) {
    fModuleConstruction->ListLVTreeLong(newValues);
  }
#ifdef ALICE_VISUALIZE
  if (command == fSetCurrentLVCmd) {
    fCurrentVolume = fModuleConstruction->FindLogicalVolume(newValues);
  }    
  else if (command == fSetDetVisibilityCmd) {
    fModuleConstruction
      ->SetDetVisibility(fSetDetVisibilityCmd->GetNewBoolValue(newValues)); 
  } 
  else if (command == fSetLVTreeVisibilityCmd) {
    fModuleConstruction
      ->SetLVTreeVisibility(fCurrentVolume,
          fSetVolVisibilityCmd->GetNewBoolValue(newValues)); 
  } 
  else if (command == fSetVolVisibilityCmd) {
    fModuleConstruction
      ->SetVolumeVisibility(fCurrentVolume,
          fSetVolVisibilityCmd->GetNewBoolValue(newValues)); 
  } 
  else if (command == fSetDetColourCmd) {
    fModuleConstruction
      ->SetDetColour(newValues);
  }     
  else if (command == fSetLVTreeColourCmd) {
    fModuleConstruction
      ->SetLVTreeColour(fCurrentVolume, newValues);
  }     
  else if (command == fSetVolColourCmd) {
    fModuleConstruction
      ->SetVolumeColour(fCurrentVolume, newValues);
  }     
#endif //ALICE_VISUALIZE
}

