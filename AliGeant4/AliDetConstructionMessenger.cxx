// $Id$
// Category: geometry
//
// See the class description in the header file.

#include "AliDetConstructionMessenger.h"
#include "AliDetConstruction.h"
#include "AliGlobals.h"

#include <G4UIcmdWithABool.hh>

AliDetConstructionMessenger::AliDetConstructionMessenger(
                                AliDetConstruction* detConstruction)
  : fDetConstruction(detConstruction)
{
//
  fSetAllSensitiveCmd
    = new G4UIcmdWithABool("/aliDet/setAllSensitive", this);
  fSetAllSensitiveCmd 
    ->SetGuidance("If true: set all logical volumes sensitive.");
  fSetAllSensitiveCmd 
    ->SetGuidance("         (Each logical is volume associated with a sensitive");
  fSetAllSensitiveCmd 
    ->SetGuidance("          detector.)");
  fSetAllSensitiveCmd 
    ->SetGuidance("If false: only volumes defined with a sensitive tracking");
  fSetAllSensitiveCmd 
    ->SetGuidance("          medium are associated with a sensitive detector.");
  fSetAllSensitiveCmd->SetParameterName("sensitivity", false);
  fSetAllSensitiveCmd->AvailableForStates(PreInit);  

  fSetReadGeometryCmd 
    = new G4UIcmdWithABool("/aliDet/readGeometry", this);
  fSetReadGeometryCmd->SetGuidance("Read geometry from g3calls.dat files");
  fSetReadGeometryCmd->SetParameterName("readGeometry", false);
  fSetReadGeometryCmd->AvailableForStates(PreInit);  
 
  fSetWriteGeometryCmd 
    = new G4UIcmdWithABool("/aliDet/writeGeometry", this);
  fSetWriteGeometryCmd->SetGuidance("Write geometry to g3calls.dat file");
  fSetWriteGeometryCmd->SetParameterName("writeGeometry", false);
  fSetWriteGeometryCmd->AvailableForStates(PreInit);   
}

AliDetConstructionMessenger::AliDetConstructionMessenger() {
//
}

AliDetConstructionMessenger::AliDetConstructionMessenger(
                                const AliDetConstructionMessenger& right)
{
//
  AliGlobals::Exception(
    "AliDetConstructionMessenger is protected from copying.");
}

AliDetConstructionMessenger::~AliDetConstructionMessenger() {
//
  delete fSetAllSensitiveCmd;
  delete fSetReadGeometryCmd;
  delete fSetWriteGeometryCmd;
}

// operators

AliDetConstructionMessenger& 
AliDetConstructionMessenger::operator=(const AliDetConstructionMessenger& right)
{
  // check assignement to self
  if (this == &right) return *this;

  AliGlobals::Exception(
     "AliDetConstructionMessenger is protected from assigning.");
    
  return *this;  
}    
          
// public methods
  
void AliDetConstructionMessenger::SetNewValue(G4UIcommand* command, G4String newValues)
{
// Applies command to the associated object.
// ---

  if (command == fSetAllSensitiveCmd) {
    fDetConstruction->SetAllLVSensitive(
                         fSetAllSensitiveCmd->GetNewBoolValue(newValues));
  }
  else if (command == fSetReadGeometryCmd) {
    fDetConstruction->SetReadGeometry(newValues);
  }  
  else if (command == fSetWriteGeometryCmd) {
    fDetConstruction->SetWriteGeometry(newValues);
  }    
}

