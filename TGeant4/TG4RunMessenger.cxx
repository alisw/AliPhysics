// $Id$
// Category: run
//
// See the class description in the header file.

#include "TG4RunMessenger.h"
#include "TG4RunManager.h"
#include "TG4Globals.h"
#include "TG4UICmdWithAComplexString.h"

#include <G4UIdirectory.hh>
#include <G4UIcmdWithoutParameter.hh>
#include <G4UIcmdWithAString.hh>

TG4RunMessenger::TG4RunMessenger(TG4RunManager* runManager)
  : fRunManager(runManager)
{ 
//
  fDirectory = new G4UIdirectory("/g4mc/");
  fDirectory->SetGuidance("TGeant4 control commands.");

  fRootCmd = new G4UIcmdWithoutParameter("/g4mc/root", this);
  fRootCmd->SetGuidance("Switch to Root interactive shell.");
  fRootCmd->AvailableForStates(PreInit, Init, Idle, GeomClosed, EventProc);

  fRootMacroCmd = new G4UIcmdWithAString("/g4mc/rootMacro", this);
  fRootMacroCmd->SetGuidance("Process Root macro with given name (from file name.C)");
  fRootMacroCmd->SetParameterName("macroName", true);
  fRootMacroCmd->AvailableForStates(PreInit, Init, Idle, GeomClosed, EventProc);

  fRootCommandCmd = new TG4UICmdWithAComplexString("/g4mc/rootCmd", this);
  fRootCommandCmd->SetGuidance("Process Root command");
  fRootCommandCmd->SetParameterName("command", false);
  fRootCommandCmd->SetDefaultValue(" ");
  fRootCommandCmd->AvailableForStates(PreInit, Init, Idle, GeomClosed, EventProc);

  fG3DefaultsCmd = new G4UIcmdWithoutParameter("/g4mc/g3Defaults", this);
  fG3DefaultsCmd->SetGuidance("Set G3 default parameters (cut values,");
  fG3DefaultsCmd->SetGuidance("tracking media max step values, ...)");
  fG3DefaultsCmd->AvailableForStates(PreInit);
}

TG4RunMessenger::TG4RunMessenger(){
//
} 

TG4RunMessenger::TG4RunMessenger(const TG4RunMessenger& right) {
// 
  TG4Globals::Exception("TG4RunMessenger is protected from copying.");
}

TG4RunMessenger::~TG4RunMessenger() {
//
  delete fDirectory;
  delete fRootCmd;
  delete fRootMacroCmd;
  delete fRootCommandCmd;
  delete fG3DefaultsCmd;
}

// operators

TG4RunMessenger& TG4RunMessenger::operator=(const TG4RunMessenger& right)
{
  // check assignement to self
  if (this == &right) return *this;

  TG4Globals::Exception("TG4RunMessenger is protected from assigning.");
    
  return *this;  
}    
          
// public methods

void TG4RunMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
// Applies command to the associated object.
// ---

  if (command == fRootCmd) {
    fRunManager->StartRootUI(); 
  }
  else if (command == fRootMacroCmd) {  
    fRunManager->ProcessRootMacro(newValue); 
  }
  else if (command == fRootCommandCmd) {
    fRunManager->ProcessRootCommand(newValue); 
  }
  else if (command == fG3DefaultsCmd) {
    fRunManager->UseG3Defaults(); 
  }
}
