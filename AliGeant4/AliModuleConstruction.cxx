// $Id$
// Category: geometry
//
// Author: I. Hrivnacova
//
// Class AliModuleConstruction
// ---------------------------
// See the class description in the header file.

#include "AliModuleConstruction.h"
#include "AliGlobals.h"
#include "AliFiles.h"
#include "AliRun.h"
#include "AliModule.h"

#include <G4UImanager.hh>

#include <TROOT.h> 
#include <TCint.h> 


//_____________________________________________________________________________
AliModuleConstruction::AliModuleConstruction(const G4String& moduleName,
                                             G4int version, 
			                     AliModuleType moduleType)
  : fAliModule(0),
    fModuleName(moduleName), 
    fType(moduleType),
    fVersion(version),
    fProcessConfig(true),
    fReadGeometry(false),
    fWriteGeometry(false),
    fDataFilePath("") {
//
}

//_____________________________________________________________________________
AliModuleConstruction::AliModuleConstruction()
  : fAliModule(0),
    fModuleName(""), 
    fType(kDetector),
    fVersion(-1),
    fProcessConfig(true),
    fReadGeometry(false),
    fWriteGeometry(false),
    fDataFilePath("") {
//
}

//_____________________________________________________________________________
AliModuleConstruction::AliModuleConstruction(const AliModuleConstruction& right)
{
//
  // copy stuff
  *this = right;
}

//_____________________________________________________________________________
AliModuleConstruction::~AliModuleConstruction()
{
//
}

// operators

//_____________________________________________________________________________
AliModuleConstruction& 
AliModuleConstruction::operator=(const AliModuleConstruction& right)
{    
  // check assignement to self
  if (this == &right) return *this;
  
  fAliModule = right.fAliModule;
  fModuleName = right.fModuleName; 
  fVersion = right.fVersion;
  fType = right.fType;
  fProcessConfig = right.fProcessConfig;
  fReadGeometry = right.fReadGeometry;
  fWriteGeometry = right.fWriteGeometry;
  fDataFilePath = right.fDataFilePath;

  return *this;
}

//_____________________________________________________________________________
G4int 
AliModuleConstruction::operator==(const AliModuleConstruction& right) const
{
//    
  return 0;
}

//_____________________________________________________________________________
G4int 
AliModuleConstruction::operator!=(const AliModuleConstruction& right) const
{
//    
  G4int returnValue = 1;
  if (*this == right) returnValue = 0; 
  
  return returnValue;
}

//_____________________________________________________________________________
void AliModuleConstruction::Configure()
{ 
// Executes the detector setup Root macro
// (extracted from AliRoot Config.C) and
// G4 macro.
// ---

  AliFiles* files = AliFiles::Instance();

  // filepaths and macro names 
  G4bool isStructure = (fType == kStructure);
  G4String rootFilePath 
    = files->GetRootMacroPath(GetDetName(), isStructure);
  G4String g4FilePath
    = files->GetG4MacroPath(GetDetName(), isStructure);
  fDataFilePath 
    = files->GetG3CallsDatPath(GetDetName(), fVersion, isStructure); 
  
  // load and execute aliroot config macro
  if (fProcessConfig) {
    gROOT->LoadMacro(rootFilePath);
    G4String macroName = files->GetDefaultMacroName();
    //macroName = macroName + "_" + GetDetName();
    macroName = macroName + "(";
    AliGlobals::AppendNumberToString(macroName, fVersion);
    macroName = macroName + ")";
    gInterpreter->ProcessLine(macroName);
  } 
  
  // process g4 config macro
  G4String command = "/control/execute ";
  G4UImanager* pUI = G4UImanager::GetUIpointer();  
  pUI->ApplyCommand(command + g4FilePath);
  
  // get AliModule created in Config.C macro
  fAliModule = gAlice->GetModule(GetDetName());
  if (!fAliModule) {
    G4String text = "AliModuleConstruction::Configure:\n";
    text = text + "    AliModule " + GetDetName();
    text = text + " has not been found in gAlice.";
    AliGlobals::Exception(text);
  }  
}

