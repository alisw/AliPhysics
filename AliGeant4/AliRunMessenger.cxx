// $Id$
// Category: run
//
// See the class description in the header file.

#include "AliRunMessenger.h"
#include "AliFiles.h"
#include "AliGlobals.h"
#include "AliRun.h"

#include <G4UIdirectory.hh>
#include <G4UIcmdWithAnInteger.hh>
#include <G4UIcmdWithoutParameter.hh>

AliRunMessenger::AliRunMessenger()
{
// 
  fRunDirectory = new G4UIdirectory("/aliRun/");
  fRunDirectory->SetGuidance("AliRun control commands.");

  fInitializeCmd = new G4UIcmdWithoutParameter("/aliRun/initialize", this);
  fInitializeCmd->SetGuidance("Initialize AliRun");
  fInitializeCmd->AvailableForStates(PreInit);

  fBeamOnCmd = new G4UIcmdWithAnInteger("/aliRun/beamOn", this);
  fBeamOnCmd->SetGuidance("Run the specified number of events");
  fBeamOnCmd->SetParameterName("NofEvents", true);
  fBeamOnCmd->SetDefaultValue(1);
  fBeamOnCmd->SetRange("NofEvents >= 0");
  fBeamOnCmd->AvailableForStates(Idle);

  fLegoCmd = new G4UIcmdWithoutParameter("/aliRun/lego", this);
  fLegoCmd->SetGuidance("Lego run");
  fLegoCmd->AvailableForStates(Idle);
}

AliRunMessenger::AliRunMessenger(const AliRunMessenger& right) {
//
  AliGlobals::Exception("AliRunMessenger is protected from copying.");
}

AliRunMessenger::~AliRunMessenger() {
//
  delete fRunDirectory;
  delete fInitializeCmd;
  delete fBeamOnCmd;
  delete fLegoCmd;
}

// operators

AliRunMessenger& AliRunMessenger::operator=(const AliRunMessenger &right)
{
  // check assignement to self
  if (this == &right) return *this;
  
  AliGlobals::Exception("AliRunMessenger is protected from assigning.");

  return *this;
}

// public methods

void AliRunMessenger::SetNewValue(G4UIcommand* command, 
       G4String newValue)
{ 
// Applies command to the associated object.
// ---

  // test gAlice
  if (!gAlice) {
    AliGlobals::Exception(
      "AliRunMessenger: gAlice has not been instantiated yet.");
  }      

  if(command == fInitializeCmd) { 
    gAlice->Init(AliFiles::Config()); 
  }   
  else if(command == fBeamOnCmd) { 
    gAlice->Run(fBeamOnCmd->GetNewIntValue(newValue)); 
  }   
  else if(command == fLegoCmd) { 
    gAlice->RunLego(); 
  }   
}
