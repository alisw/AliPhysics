// $Id$
// Category: geometry
//
// Author: I. Hrivnacova
//
// Class AliSingleModuleConstruction
// ---------------------------------
// See the class description in the header file.

#include "AliSingleModuleConstruction.h"
#include "AliGlobals.h"
#include "AliFiles.h"
#include "AliRun.h"
#include "AliModule.h"

#include "TG4GeometryManager.h"

#include <G4UImanager.hh>
//#include <G4Element.hh>

#include <TROOT.h> 
#include <TCint.h> 

G4VPhysicalVolume* AliSingleModuleConstruction::fgWorld = 0;

//_____________________________________________________________________________
AliSingleModuleConstruction::AliSingleModuleConstruction(
                                G4String moduleName, G4int version,
				AliModuleType moduleType)
  : AliModuleConstruction(moduleName),
    fVersion(version),
    fType(moduleType),
    fProcessConfig(true)
{
//
}

//_____________________________________________________________________________
AliSingleModuleConstruction::AliSingleModuleConstruction(
                                const AliSingleModuleConstruction& right)
  : AliModuleConstruction(right)				
{
//
  // copy stuff
  *this = right;
}

//_____________________________________________________________________________
AliSingleModuleConstruction::AliSingleModuleConstruction() {
//
}

//_____________________________________________________________________________
AliSingleModuleConstruction::~AliSingleModuleConstruction() {
//
}

// operators

//_____________________________________________________________________________
AliSingleModuleConstruction& 
AliSingleModuleConstruction::operator=(const AliSingleModuleConstruction& right)
{    
  // check assignement to self
  if (this == &right) return *this;
  
  // base class assignement
  AliModuleConstruction::operator=(right);
  
  fVersion = right.fVersion;
  fType = right.fType;
  fProcessConfig = right.fProcessConfig;
 
  return *this;
}

// public methods 

//_____________________________________________________________________________
void AliSingleModuleConstruction::Configure(const AliFiles& files)
{ 
// Executes the detector setup Root macro
// (extracted from AliRoot Config.C) and
// G4 macro.
// ---

  // filepaths and macro names 
  G4bool isStructure = (fType == kStructure);
  G4String rootFilePath 
    = files.GetRootMacroPath(fModuleName, isStructure);
  G4String g4FilePath
    = files.GetG4MacroPath(fModuleName, isStructure);
  fDataFilePath 
    = files.GetG3CallsDatPath(fModuleName, fVersion, isStructure); 
  
  // load and execute aliroot config macro
  if (fProcessConfig) {
    gROOT->LoadMacro(rootFilePath);
    G4String macroName = files.GetDefaultMacroName();
    //macroName = macroName + "_" + fModuleName;
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
  fAliModule = gAlice->GetModule(fModuleName);
  if (!fAliModule) {
    G4String text = "AliSingleModuleConstruction::Configure:\n";
    text = text + "    AliModule " + fModuleName;
    text = text + " has not been found in gAlice.";
    AliGlobals::Exception(text);
  }  
}

//_____________________________________________________________________________
void AliSingleModuleConstruction::Construct()
{ 
// Constructs geometry.
// ---

  // print default element table
  // const G4ElementTable* table = G4Element::GetElementTable();
  // G4cout << "Default elemnt table: " << G4endl;
  // for (G4int i=0; i<table->entries(); i++) {
  //   G4cout << *(*table)[i] << G4endl;
  // }  

  // Configure();

  // get geometry manager
  TG4GeometryManager* pGeometryManager = TG4GeometryManager::Instance();

  // register module name in the name map
  pGeometryManager->SetMapSecond(fAliModule->GetName());	

  if (fReadGeometry) {
    // create G3 geometry from g3calls.dat
    pGeometryManager->SetWriteGeometry(false);
    pGeometryManager->ReadG3Geometry(fDataFilePath);
  }
  else {
    // set geometry output stream for this module
    pGeometryManager->SetWriteGeometry(fWriteGeometry);
    if (fWriteGeometry) 
      pGeometryManager->OpenOutFile(fDataFilePath);

    // create geometry from AliRoot

    // construct materials
    fAliModule->CreateMaterials();

    // construct G3 geometry
    fAliModule->CreateGeometry();
        
    if (fWriteGeometry) 
      pGeometryManager->CloseOutFile();
  }  
  
  // construct G4 geometry
  G4VPhysicalVolume* world = pGeometryManager->CreateG4Geometry();
  if (!fgWorld) fgWorld = world; 
  
  // set the detector frame (envelope)
  // (without warning output if enevelope is not defined)
  SetDetFrame(false);

  // construct geometry for display
  fAliModule->BuildGeometry();

  // reset TG4GeometryManager 
  pGeometryManager->ClearG3Tables();

#ifdef ALICE_VISUALIZE
  if (GetDetFrame()) {
    // set visualization attributes
    // if detector envelope is defined
    SetDetVisibility(true);
    SetDetColour("Yellow");
  }  
#endif
}
