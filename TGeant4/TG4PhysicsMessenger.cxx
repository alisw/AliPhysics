// $Id$
// Category: physics
//
// Author: I. Hrivnacova
//
// Class TG4PhysicsMessenger
// -------------------------
// See the class description in the header file.

#include "TG4PhysicsMessenger.h"
#include "TG4GeometryServices.h"
#include "TG4PhysicsManager.h"
#include "TG4G3PhysicsManager.h"
#include "TG4G3CutVector.h"
#include "TG4G3ControlVector.h"
#include "TG4ProcessControlMap.h"
#include "TG4ProcessMCMap.h"

#include <G4UIdirectory.hh>
#include <G4UIcmdWithoutParameter.hh>
#include <G4UIcmdWithABool.hh>
#include <G4UIcmdWithAString.hh>

//_____________________________________________________________________________
TG4PhysicsMessenger::TG4PhysicsMessenger(TG4PhysicsManager* physicsManager)
  : fPhysicsManager(physicsManager)
{ 
//
  fDirectory = new G4UIdirectory("/tg4Physics/");
  fDirectory->SetGuidance("TGeant4 physics control commands.");

  fSetEMCmd
     = new G4UIcmdWithABool("/tg4Physics/setEM", this);
  fSetEMCmd->SetGuidance("Set electromagnetic physics.");
  fSetEMCmd->SetParameterName("EMControl", false);
  fSetEMCmd->AvailableForStates(PreInit);

  fSetMuonCmd
     = new G4UIcmdWithABool("/tg4Physics/setMuon", this);
  fSetMuonCmd->SetGuidance("Set muon physics.");
  fSetMuonCmd->SetParameterName("EMControl", false);
  fSetMuonCmd->AvailableForStates(PreInit);

  fSetHadronCmd
     = new G4UIcmdWithABool("/tg4Physics/setHadron", this);
  fSetHadronCmd->SetGuidance("Set hadron physics.");
  fSetHadronCmd->SetParameterName("HadronControl", false);
  fSetHadronCmd->AvailableForStates(PreInit);

  fSetOpticalCmd
     = new G4UIcmdWithABool("/tg4Physics/setOptical", this);
  fSetOpticalCmd->SetGuidance("Set Cerenkov and optical physics.");
  fSetOpticalCmd->SetParameterName("OpticalControl", false);
  fSetOpticalCmd->AvailableForStates(PreInit);

  fSetSpecialCutsCmd
     = new G4UIcmdWithABool("/tg4Physics/setSpecialCuts", this);
  fSetSpecialCutsCmd->SetGuidance("Set special cuts process.");
  fSetSpecialCutsCmd
    ->SetGuidance("!! Support for this option is under development.");
  fSetSpecialCutsCmd->SetParameterName("SpecialCutsControl", false);
  fSetSpecialCutsCmd->AvailableForStates(PreInit);

  fSetSpecialControlsCmd
     = new G4UIcmdWithABool("/tg4Physics/setSpecialControls", this);
  fSetSpecialControlsCmd->SetGuidance("Set special controls process.");
  fSetSpecialControlsCmd
    ->SetGuidance("!! Support for this option is under development.");
  fSetSpecialControlsCmd->SetParameterName("SpecialFlagsControl", false);
  fSetSpecialControlsCmd->AvailableForStates(PreInit);

  fProcessActivationCmd
     = new G4UIcmdWithoutParameter("/tg4Physics/setProcessActivation", this);
  fProcessActivationCmd->SetGuidance("Activate/inactivate physics processes.");
  fProcessActivationCmd->AvailableForStates(Idle);

  fPrintProcessMCMapCmd
     = new G4UIcmdWithoutParameter("/tg4Physics/printProcessMCMap", this);
  fPrintProcessMCMapCmd
    ->SetGuidance("Prints mapping of G4 processes to G3 controls.");
  fPrintProcessMCMapCmd->AvailableForStates(Idle);

  fPrintProcessControlMapCmd
     = new G4UIcmdWithoutParameter("/tg4Physics/printProcessControlMap", this);
  fPrintProcessControlMapCmd
    ->SetGuidance("Prints mapping of G4 processes to G3 controls.");
  fPrintProcessControlMapCmd->AvailableForStates(Idle);

  fPrintVolumeLimitsCmd
     = new G4UIcmdWithAString("/tg4Physics/printVolumeLimits", this);
  fPrintVolumeLimitsCmd
    ->SetGuidance("Prints the limits set to the specified volume.");
  fPrintVolumeLimitsCmd->SetParameterName("PrintVolumeLimits", false);
  fPrintVolumeLimitsCmd->AvailableForStates(Idle);

  fPrintGeneralCutsCmd
     = new G4UIcmdWithoutParameter("/tg4Physics/printGeneralCuts", this);
  fPrintGeneralCutsCmd
    ->SetGuidance("Prints the general G3 cuts.");
  fPrintGeneralCutsCmd->AvailableForStates(Idle);

  fPrintGeneralControlsCmd
     = new G4UIcmdWithoutParameter("/tg4Physics/printGeneralControls", this);
  fPrintGeneralControlsCmd
    ->SetGuidance("Prints the general G3 process controls.");
  fPrintGeneralControlsCmd->AvailableForStates(Idle);
}

//_____________________________________________________________________________
TG4PhysicsMessenger::TG4PhysicsMessenger(){
//
} 

//_____________________________________________________________________________
TG4PhysicsMessenger::TG4PhysicsMessenger(const TG4PhysicsMessenger& right) {
// 
  TG4Globals::Exception("TG4PhysicsMessenger is protected from copying.");
}

//_____________________________________________________________________________
TG4PhysicsMessenger::~TG4PhysicsMessenger() {
//

  delete fDirectory;
  delete fSetEMCmd;
  delete fSetMuonCmd;
  delete fSetHadronCmd;
  delete fSetOpticalCmd;
  delete fSetSpecialCutsCmd;
  delete fSetSpecialControlsCmd;
  delete fProcessActivationCmd;
  delete fPrintProcessMCMapCmd;
  delete fPrintProcessControlMapCmd;
  delete fPrintVolumeLimitsCmd;
  delete fPrintGeneralCutsCmd;
  delete fPrintGeneralControlsCmd;
}

// operators

//_____________________________________________________________________________
TG4PhysicsMessenger& TG4PhysicsMessenger::operator=(const TG4PhysicsMessenger& right)
{
  // check assignement to self
  if (this == &right) return *this;

  TG4Globals::Exception("TG4PhysicsMessenger is protected from assigning.");
    
  return *this;  
}    
          
// public methods

//_____________________________________________________________________________
void TG4PhysicsMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
// Applies command to the associated object.
// ---

  if (command == fSetEMCmd) {
    fPhysicsManager
      ->SetEMPhysics(fSetEMCmd->GetNewBoolValue(newValue)); 
  }    
  if (command == fSetMuonCmd) {
    fPhysicsManager
      ->SetMuonPhysics(fSetMuonCmd->GetNewBoolValue(newValue)); 
  }    
  else if (command == fSetHadronCmd) {
    fPhysicsManager
      ->SetHadronPhysics(fSetHadronCmd->GetNewBoolValue(newValue)); 
  }    
  else if (command == fSetOpticalCmd) {
    fPhysicsManager
      ->SetOpticalPhysics(fSetOpticalCmd->GetNewBoolValue(newValue)); 
  }    
  else if (command == fSetSpecialCutsCmd) {
    fPhysicsManager
      ->SetSpecialCutsPhysics(fSetSpecialCutsCmd->GetNewBoolValue(newValue)); 
  }    
  else if (command == fSetSpecialControlsCmd) {
    fPhysicsManager
      ->SetSpecialControlsPhysics(
          fSetSpecialControlsCmd->GetNewBoolValue(newValue)); 
  }    
  else if (command == fProcessActivationCmd) {
    fPhysicsManager->SetProcessActivation();
  }  
  else if (command == fPrintProcessMCMapCmd) {
    TG4ProcessMCMap::Instance()->PrintAll();
  }  
  else if (command == fPrintProcessControlMapCmd) {
    TG4ProcessControlMap::Instance()->PrintAll();
  }  
  else if (command == fPrintVolumeLimitsCmd) {
    TG4GeometryServices::Instance()->PrintVolumeLimits(newValue);
  }  
  else if (command == fPrintGeneralCutsCmd) {
    TG4G3PhysicsManager::Instance()->GetCutVector()->Print();
  }  
  else if (command == fPrintGeneralControlsCmd) {
    TG4G3PhysicsManager::Instance()->GetControlVector()->Print();
  }  
}
