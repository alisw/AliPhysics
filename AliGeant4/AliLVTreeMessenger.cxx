// $Id$
// Category: geometry
//
// Author: I. Hrivnacova
//
// Class AliLVTreeMessenger
// ------------------------------------
// See the class description in the header file.

#include "AliLVTreeMessenger.h"
#include "AliLVTree.h"
#include "AliGlobals.h"
#ifdef G4VIS_USE
#include "AliColourStore.h"
#endif

#include "TG4GeometryServices.h"

#include <G4UIdirectory.hh>
#include <G4UIcmdWithABool.hh>
#include <G4UIcmdWithAString.hh>
#include <G4UIcmdWithoutParameter.hh>

//_____________________________________________________________________________
AliLVTreeMessenger::AliLVTreeMessenger(AliLVTree* lvTree)
 : fLVTree(lvTree),
   fCurrentVolume(0)
{
//
  G4String dirName = "/aliTree/"; 
  fDirectory = new G4UIdirectory(dirName);
  fDirectory->SetGuidance("LV tree control commands.");

  G4String commandPath = dirName + "setVolume";
  fSetCurrentLVCmd = new G4UIcmdWithAString(commandPath, this);
  fSetCurrentLVCmd->SetGuidance("Set the current logical volume.");
  fSetCurrentLVCmd->SetParameterName("curVolume", false);
  fSetCurrentLVCmd->AvailableForStates(PreInit,Idle);  
 
  commandPath = dirName + "list";
  fListCmd = new G4UIcmdWithoutParameter(commandPath, this);
  G4String guidance = "List LV tree of the current volume. ";
  fListCmd->SetGuidance(guidance);
  fListCmd->AvailableForStates(PreInit,Idle);  
 
  commandPath = dirName + "listLong";
  fListLongCmd = new G4UIcmdWithoutParameter(commandPath, this);
  guidance = "List LV tree of the current volume \n";
  guidance = guidance + "including number of its daughters.";
  fListLongCmd->SetGuidance(guidance);
  fListLongCmd->AvailableForStates(PreInit,Idle);  
 
#ifdef G4VIS_USE
  commandPath = dirName + "setLVTreeVisibility";
  fSetLVTreeVisibilityCmd = new G4UIcmdWithABool(commandPath, this);
  fSetLVTreeVisibilityCmd 
    ->SetGuidance("Make current volume tree visible/invisible.");
  fSetLVTreeVisibilityCmd->SetParameterName("lvtreeVisibility", false);
  fSetLVTreeVisibilityCmd->AvailableForStates(PreInit,Idle);  
 
  commandPath = dirName + "setVolVisibility";
  fSetVolVisibilityCmd = new G4UIcmdWithABool(commandPath, this);
  fSetVolVisibilityCmd 
    ->SetGuidance("Make current volume visible/invisible.");
  fSetVolVisibilityCmd->SetParameterName("volVisibility", false);
  fSetVolVisibilityCmd->AvailableForStates(PreInit,Idle);  
 
  commandPath = dirName + "setLVTreeColour";
  fSetLVTreeColourCmd = new G4UIcmdWithAString(commandPath, this);
  fSetLVTreeColourCmd->SetGuidance("Set colour for the current volume tree.");
  fSetLVTreeColourCmd->SetGuidance("Available colours:");
  guidance = AliColourStore::Instance()->GetColoursListWithCommas();
  fSetLVTreeColourCmd->SetGuidance(guidance);
  fSetLVTreeColourCmd->SetParameterName("lvtreeColour", false);
  G4String candidatesList = AliColourStore::Instance()->GetColoursList();  
  fSetLVTreeColourCmd->SetCandidates(candidatesList);
  fSetLVTreeColourCmd->AvailableForStates(PreInit,Idle);  

  commandPath = dirName + "setVolColour";
  fSetVolColourCmd = new G4UIcmdWithAString(commandPath, this);
  fSetVolColourCmd->SetGuidance("Set colour for the current volume.");
  fSetVolColourCmd->SetGuidance("Available colours:");
  guidance = AliColourStore::Instance()->GetColoursListWithCommas();
  fSetVolColourCmd->SetGuidance(guidance);
  fSetVolColourCmd->SetParameterName("volColour", false);
  candidatesList = AliColourStore::Instance()->GetColoursList();  
  fSetVolColourCmd->SetCandidates(candidatesList);
  fSetVolColourCmd->AvailableForStates(PreInit,Idle);  
#endif //G4VIS_USE
}

//_____________________________________________________________________________
AliLVTreeMessenger::AliLVTreeMessenger() {
//
}

//_____________________________________________________________________________
AliLVTreeMessenger::AliLVTreeMessenger(const AliLVTreeMessenger& right)
{
//
  AliGlobals::Exception(
    "AliLVTreeMessenger is protected from copying.");
}

//_____________________________________________________________________________
AliLVTreeMessenger::~AliLVTreeMessenger()
{
//
  delete fDirectory;
  delete fSetCurrentLVCmd;
  delete fListCmd;
  delete fListLongCmd;
#ifdef G4VIS_USE
  delete fSetLVTreeVisibilityCmd;    
  delete fSetVolVisibilityCmd;    
  delete fSetLVTreeColourCmd;    
  delete fSetVolColourCmd;    
#endif //G4VIS_USE
}

// operators

//_____________________________________________________________________________
AliLVTreeMessenger& 
AliLVTreeMessenger::operator=(const AliLVTreeMessenger& right)
{
  // check assignement to self
  if (this == &right) return *this;

  AliGlobals::Exception(
     "AliLVTreeMessenger is protected from assigning.");
    
  return *this;  
}    
          
// public methods
  
//_____________________________________________________________________________
void AliLVTreeMessenger::SetNewValue(G4UIcommand* command,
                                     G4String newValues)
{
// Applies command to the associated object.
// ---

  G4String dirName = "/aliTree/"; 
  fDirectory = new G4UIdirectory(dirName);
  G4String guidance = "LV tree control commands ";
  fDirectory->SetGuidance(guidance);

  if (command == fSetCurrentLVCmd) {
    fCurrentVolume 
      = TG4GeometryServices::Instance()->FindLogicalVolume(newValues);
  }    
  else if (command == fListCmd) {
    fLVTree->List(fCurrentVolume);
  }
  else if (command == fListLongCmd) {
    fLVTree->ListLong(fCurrentVolume);
  }
#ifdef G4VIS_USE
  if (command == fSetLVTreeVisibilityCmd) {
    fLVTree
      ->SetLVTreeVisibility(fCurrentVolume,
          fSetVolVisibilityCmd->GetNewBoolValue(newValues)); 
  } 
  else if (command == fSetVolVisibilityCmd) {
    fLVTree
      ->SetVolumeVisibility(fCurrentVolume,
          fSetVolVisibilityCmd->GetNewBoolValue(newValues)); 
  } 
  else if (command == fSetLVTreeColourCmd) {
    fLVTree
      ->SetLVTreeColour(fCurrentVolume, newValues);
  }     
  else if (command == fSetVolColourCmd) {
    fLVTree
      ->SetVolumeColour(fCurrentVolume, newValues);
  }     
#endif //G4VIS_USE
}

