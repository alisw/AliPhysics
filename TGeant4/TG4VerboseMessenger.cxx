// $Id$
// Category: global
//
// Author: I. Hrivnacova
//
// Class TG4VerboseMessenger
// ------------------
// See the class description in the header file.

#include "TG4VerboseMessenger.h"
#include "TG4VVerbose.h"
#include "TG4Globals.h"

#include <G4UIcmdWithAnInteger.hh>
#include <G4UIdirectory.hh>
#include <G4UImanager.hh>
#include <G4UIcommandTree.hh>

//_____________________________________________________________________________
TG4VerboseMessenger::TG4VerboseMessenger(const G4String& directoryName)
  : fkDirectoryName(directoryName)
{ 
//

  fDirectory = new G4UIdirectory(directoryName);
  fDirectory->SetGuidance("TGeant4 verbose control commands.");

  // sets the given level to all verbose instances
  fGlobalVerboseCmd 
    = new G4UIcmdWithAnInteger(G4String(directoryName + "all"), this);
  G4String guidance("Set a given verbose level to all verbose instances.");
  fGlobalVerboseCmd->SetGuidance(guidance);
  fGlobalVerboseCmd->SetParameterName("GlobalVerbose", false);
  fGlobalVerboseCmd->AvailableForStates(PreInit, Init, Idle);
}


//_____________________________________________________________________________
void TG4VerboseMessenger::AddCommand(TG4VVerbose* verbose, 
                                     const G4String& cmdName)
{
//
//--

  G4UIcmdWithAnInteger* cmd 
    = new G4UIcmdWithAnInteger(G4String(fkDirectoryName + cmdName), this);

  fVerboseVector.push_back(verbose);
  fCommandVector.push_back(cmd);

  G4String guidance("Set  verbose level.");
  guidance.insert(4,cmdName);
  cmd->SetGuidance(guidance);

  G4String parameterName("Verbose");
  parameterName.insert(0,cmdName);
  cmd->SetParameterName(parameterName, false);
  
  cmd->AvailableForStates(PreInit, Init, Idle);
}

//_____________________________________________________________________________
TG4VerboseMessenger::TG4VerboseMessenger(){
//
} 

//_____________________________________________________________________________
TG4VerboseMessenger::TG4VerboseMessenger(const TG4VerboseMessenger& right) {
// 
  TG4Globals::Exception("TG4VerboseMessenger is protected from copying.");
}

//_____________________________________________________________________________
TG4VerboseMessenger::~TG4VerboseMessenger() {
//

  delete fDirectory;
  delete fGlobalVerboseCmd;

  for (G4int i=0; i<fCommandVector.size(); i++)
    delete fCommandVector[i];
}

// operators

//_____________________________________________________________________________
TG4VerboseMessenger& TG4VerboseMessenger::operator=(const TG4VerboseMessenger& right)
{
  // check assignement to self
  if (this == &right) return *this;

  TG4Globals::Exception("TG4VerboseMessenger is protected from assigning.");
    
  return *this;  
}    
          
// private methods

//_____________________________________________________________________________
void TG4VerboseMessenger::SetNewValueToAll(const G4String value) const
{
// Sets the value to all registered verbose instances.
// ---
   
   G4UIcommandTree* cmdTree
     = G4UImanager::GetUIpointer()->GetTree()->GetTree(fkDirectoryName);

   for (G4int i=0; i<cmdTree->GetCommandEntry(); i++) {
     if (cmdTree->GetCommand(i+1)->GetCommandName() != "all") {    
        // skip the first command in the tree ("all")
        cmdTree->GetCommand(i+1)->DoIt(value);     
     }	
   }  
}

// public methods

//_____________________________________________________________________________
void TG4VerboseMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
// Applies command to the associated object.
// ---

  if (command == fGlobalVerboseCmd) {
    G4cout << "SetNewValueToAll  " << G4endl;
    SetNewValueToAll(newValue); 
  }    
  for (G4int i=0; i<fCommandVector.size(); i++)  
    if (command == fCommandVector[i]) {
      fVerboseVector[i]
        ->VerboseLevel(fCommandVector[i]->GetNewIntValue(newValue)); 
    }    
}
