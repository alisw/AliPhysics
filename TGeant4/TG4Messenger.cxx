// $Id$
// Category: run
//
// Author: I. Hrivnacova
//
// Class TG4Messenger
// ------------------
// See the class description in the header file.

#include "TG4Messenger.h"
#include "TG4GeometryManager.h"
#include "TG4GeometryServices.h"
#include "TG4StepManager.h"
#include "TG4PhysicsManager.h"
#include "TG4G3PhysicsManager.h"
#include "TG4G3CutVector.h"
#include "TG4G3ControlVector.h"
#include "TG4ProcessControlMap.h"
#include "TG4ProcessMCMap.h"

#include <G4UIcmdWithoutParameter.hh>
#include <G4UIcmdWithABool.hh>
#include <G4UIcmdWithAString.hh>

//_____________________________________________________________________________
TG4Messenger::TG4Messenger(TG4GeometryManager* geometryManager, 
                           TG4PhysicsManager* physicsManager, 
			   TG4StepManager* stepManager)
  : fGeometryManager(geometryManager),
    fPhysicsManager(physicsManager),
    fStepManager(stepManager)
{ 
//
  fSetEMCmd
     = new G4UIcmdWithABool("/g4mc/setEM", this);
  fSetEMCmd->SetGuidance("Set electromagnetic physics.");
  fSetEMCmd->SetParameterName("EMControl", false);
  fSetEMCmd->AvailableForStates(PreInit);

  fSetMuonCmd
     = new G4UIcmdWithABool("/g4mc/setMuon", this);
  fSetMuonCmd->SetGuidance("Set muon physics.");
  fSetMuonCmd->SetParameterName("EMControl", false);
  fSetMuonCmd->AvailableForStates(PreInit);

  fSetHadronCmd
     = new G4UIcmdWithABool("/g4mc/setHadron", this);
  fSetHadronCmd->SetGuidance("Set hadron physics.");
  fSetHadronCmd->SetParameterName("HadronControl", false);
  fSetHadronCmd->AvailableForStates(PreInit);

  fSetOpticalCmd
     = new G4UIcmdWithABool("/g4mc/setOptical", this);
  fSetOpticalCmd->SetGuidance("Set Cerenkov and optical physics.");
  fSetOpticalCmd->SetParameterName("OpticalControl", false);
  fSetOpticalCmd->AvailableForStates(PreInit);

  fSetSpecialCutsCmd
     = new G4UIcmdWithABool("/g4mc/setSpecialCuts", this);
  fSetSpecialCutsCmd->SetGuidance("Set special cuts process.");
  fSetSpecialCutsCmd
    ->SetGuidance("!! Support for this option is under development.");
  fSetSpecialCutsCmd->SetParameterName("SpecialCutsControl", false);
  fSetSpecialCutsCmd->AvailableForStates(PreInit);

  fSetSpecialControlsCmd
     = new G4UIcmdWithABool("/g4mc/setSpecialControls", this);
  fSetSpecialControlsCmd->SetGuidance("Set special controls process.");
  fSetSpecialControlsCmd
    ->SetGuidance("!! Support for this option is under development.");
  fSetSpecialControlsCmd->SetParameterName("SpecialFlagsControl", false);
  fSetSpecialControlsCmd->AvailableForStates(PreInit);

  fProcessActivationCmd
     = new G4UIcmdWithoutParameter("/g4mc/setProcessActivation", this);
  fProcessActivationCmd->SetGuidance("Activate/inactivate physics processes.");
  fProcessActivationCmd->AvailableForStates(Idle);

  fPrintProcessMCMapCmd
     = new G4UIcmdWithoutParameter("/g4mc/printProcessMCMap", this);
  fPrintProcessMCMapCmd
    ->SetGuidance("Prints mapping of G4 processes to G3 controls.");
  fPrintProcessMCMapCmd->AvailableForStates(Idle);

  fPrintProcessControlMapCmd
     = new G4UIcmdWithoutParameter("/g4mc/printProcessControlMap", this);
  fPrintProcessControlMapCmd
    ->SetGuidance("Prints mapping of G4 processes to G3 controls.");
  fPrintProcessControlMapCmd->AvailableForStates(Idle);

  fPrintVolumeLimitsCmd
     = new G4UIcmdWithAString("/g4mc/printVolumeLimits", this);
  fPrintVolumeLimitsCmd
    ->SetGuidance("Prints the limits set to the specified volume.");
  fPrintVolumeLimitsCmd->SetParameterName("PrintVolumeLimits", false);
  fPrintVolumeLimitsCmd->AvailableForStates(Idle);

  fPrintGeneralCutsCmd
     = new G4UIcmdWithoutParameter("/g4mc/printGeneralCuts", this);
  fPrintGeneralCutsCmd
    ->SetGuidance("Prints the general G3 cuts.");
  fPrintGeneralCutsCmd->AvailableForStates(Idle);

  fPrintGeneralControlsCmd
     = new G4UIcmdWithoutParameter("/g4mc/printGeneralControls", this);
  fPrintGeneralControlsCmd
    ->SetGuidance("Prints the general G3 process controls.");
  fPrintGeneralControlsCmd->AvailableForStates(Idle);
}

//_____________________________________________________________________________
TG4Messenger::TG4Messenger(){
//
} 

//_____________________________________________________________________________
TG4Messenger::TG4Messenger(const TG4Messenger& right) {
// 
  TG4Globals::Exception("TG4Messenger is protected from copying.");
}

//_____________________________________________________________________________
TG4Messenger::~TG4Messenger() {
//
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
TG4Messenger& TG4Messenger::operator=(const TG4Messenger& right)
{
  // check assignement to self
  if (this == &right) return *this;

  TG4Globals::Exception("TG4Messenger is protected from assigning.");
    
  return *this;  
}    
          
// public methods

//_____________________________________________________________________________
void TG4Messenger::SetNewValue(G4UIcommand* command, G4String newValue)
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
