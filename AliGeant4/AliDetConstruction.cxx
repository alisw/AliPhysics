// $Id$
// Category: geometry
//
// Author: I. Hrivnacova
//
// Class AliDetConstruction
// ------------------------
// See the class description in the header file.

#include "AliDetConstruction.h"
#include "AliSingleModuleConstruction.h"
#include "AliDetSwitch.h"
#include "AliGlobals.h"
#include "AliRun.h"
#include "AliModule.h"

//_____________________________________________________________________________
AliDetConstruction::AliDetConstruction()
  : fTopVolumeName("ALIC")
{
  // initialize det switch vector: 
  // moduleName nofVersions defaultVersion PPRVersion [type isStandalone]     
        // det switch objects are deleted in
	// tbe base class (AliModulesCompositions) destructor

  AliDetSwitch* detSwitch;
  detSwitch = new AliDetSwitch("ABSO",   1, 0, 0, kStructure);
  AddDetSwitch(detSwitch); 
  detSwitch = new AliDetSwitch("DIPO",   3, 2, 2, kStructure, false);
  AddDetSwitch(detSwitch); 
  detSwitch = new AliDetSwitch("FRAME",  3, 2, 2, kStructure, false);
  AddDetSwitch(detSwitch); 
  detSwitch = new AliDetSwitch("HALL",   1, 0, 0, kStructure);
  AddDetSwitch(detSwitch); 
  detSwitch = new AliDetSwitch("MAG",    1, 0, 0, kStructure);
  AddDetSwitch(detSwitch); 
  detSwitch = new AliDetSwitch("PIPE",   5, 0, 0, kStructure);
  AddDetSwitch(detSwitch); 
  detSwitch = new AliDetSwitch("SHIL",   1, 0, 0, kStructure);
  AddDetSwitch(detSwitch); 
  detSwitch = new AliDetSwitch("CASTOR", 2, 1, 1);
  AddDetSwitch(detSwitch); 
  detSwitch = new AliDetSwitch("FMD",    2, 1, 0);
  AddDetSwitch(detSwitch); 
  detSwitch = new AliDetSwitch("ITS",    7, 5, 5);
  AddDetSwitch(detSwitch); 
  detSwitch = new AliDetSwitch("MUON",   2, 1, 1);
  AddDetSwitch(detSwitch); 
  detSwitch = new AliDetSwitch("PHOS",   5, 1, 1);
  AddDetSwitch(detSwitch); 
  detSwitch = new AliDetSwitch("PMD",    3, 1, 1);
  AddDetSwitch(detSwitch); 
  detSwitch = new AliDetSwitch("RICH",   3, 1, 1);
  AddDetSwitch(detSwitch); 
  detSwitch = new AliDetSwitch("START",  2, 1, 1);
  AddDetSwitch(detSwitch); 
  detSwitch = new AliDetSwitch("TOF",    5, 2, 2, kDetector, false);
  AddDetSwitch(detSwitch); 
  detSwitch = new AliDetSwitch("TPC",    4, 2, 2);
  AddDetSwitch(detSwitch); 
  detSwitch = new AliDetSwitch("TRD",    2, 1, 1, kDetector, false);
  AddDetSwitch(detSwitch); 
  detSwitch = new AliDetSwitch("ZDC",    2, 1, 1, kDetector, false);
  AddDetSwitch(detSwitch);  
}

//_____________________________________________________________________________
AliDetConstruction::AliDetConstruction(const AliDetConstruction& right)
  : AliModulesComposition(right)
{
  // AliModuleComposition is protected from copying
}  

//_____________________________________________________________________________
AliDetConstruction::~AliDetConstruction() {
//
}

// operators

//_____________________________________________________________________________
AliDetConstruction& 
AliDetConstruction::operator=(const AliDetConstruction& right)
{
  // check assignement to self
  if (this == &right) return *this;

  // base class assignement
  // AliModuleComposition is protected from assigning
  AliModulesComposition::operator=(right);

  return *this;  
}    
          
// private methods

//_____________________________________________________________________________
void AliDetConstruction::BuildDetectors()
{
// Create module constructions for AliModules 
// that have been created and registered by gAlice 
// ---

  TObjArray* pDetectors = gAlice->Detectors();
  TIter next(pDetectors);

  // the first AliModule is expected to be the top volume
  AliModule *module = (AliModule*)next();
  if (G4String(module->GetName()) != "BODY") {
    G4String text = "AliDetConstruction::BuildDetectors():\n";
    text = text + "    Instead of BODY - the first module ";
    text = text + module->GetName() + " has been found.";
    AliGlobals::Exception(text);
  }  
  AddSingleModuleConstruction("BODY", 0, kStructure);

  G4bool first = true;
  while ((module = (AliModule*)next())) {
    // register moduleConstruction in fDetSwitchVector
    // in order to keep availability of /AlDet/list command
    G4String modName = module->GetName();
    G4int modVersion = module->IsVersion();
    if (first)
      // skip registering of the top volume 
      first = false;
    else 
      SwitchDetOn(modName, modVersion);
 
    // all modules will be processed alltogether
    AddMoreModuleConstruction(modName, modVersion);

    G4cout << "Created module construction for " 
           << modName << "v" << modVersion << "." << G4endl;   
  }
  
  // do not process Config.C 
  // (it was processed when creating modules by gAlice)
  SetProcessConfigToModules(false); 	
}

//_____________________________________________________________________________
void AliDetConstruction::CreateDetectors()
{
// Creates AliModules and their module constructions 
// according to the fDetSwitchVector
// ---

  // add top volume (AliBODY) construction first
  AddSingleModuleConstruction("BODY", 0, kStructure);

  // add modules constructions
  for (G4int id=0; id<fDetSwitchVector.size(); id++)
  {
    G4String detName = fDetSwitchVector[id]->GetDetName();
    G4int version = fDetSwitchVector[id]->GetSwitchedVersion();
    G4bool isStandalone = fDetSwitchVector[id]->IsStandalone();
    AliModuleType type = fDetSwitchVector[id]->GetType();
    
    if (version > -1)
      if (isStandalone)
        AddSingleModuleConstruction(detName, version, type);
      else
        AddMoreModuleConstruction(detName, version, type);
  }    
}

//_____________________________________________________________________________
void AliDetConstruction::CheckDetDependencies()
{
// Checks modules dependencies.
// Dependent modules FRAME, TOF, TRD 
// TOF always requires FRAMEv1
// TRD can be built with both (??)
// ZDC requires DIPO
// ---

  // get switched versions of dependent modules
  G4int verTOF = GetDetSwitch("TOF")->GetSwitchedVersion(); 
  G4int verTRD = GetDetSwitch("TRD")->GetSwitchedVersion(); 
  G4int verZDC = GetDetSwitch("ZDC")->GetSwitchedVersion(); 
  G4int verFRAME = GetDetSwitch("FRAME")->GetSwitchedVersion(); 
  
  // check dependencies  
  if (verTOF > -1) {
    // TOF requires FRAMEv1 - obsolete? 
    if (verFRAME != 2) {
      GetDetSwitch("FRAME")->SwitchOn(2);
      G4String text = "AliDetConstruction::CheckDetDependencies: \n";
      text = text + "    Switched TOF requires FRAME v1.\n"; 
      text = text + "    The det switch for FRAME has been changed."; 
      AliGlobals::Warning(text);
    }  
  }
  if (verTRD > -1) {
    // TRD requires FRAME
    verFRAME = GetDetSwitch("FRAME")->GetSwitchedVersion(); 
    if (verFRAME < 0) {
      GetDetSwitch("FRAME")->SwitchOnDefault();
      G4String text = "AliDetConstruction::CheckDetDependencies: \n";
      text = text + "    Switched TRD requires FRAME.\n"; 
      text = text + "    The det switch for FRAME has been changed."; 
      AliGlobals::Warning(text);
    }  
  }  
  if (verZDC > 0) {
    // ZDC requires PIPE, ABSO, DIPO, SHIL 
    G4int verPIPE = GetDetSwitch("PIPE")->GetSwitchedVersion(); 
    G4int verABSO = GetDetSwitch("ABSO")->GetSwitchedVersion(); 
    G4int verDIPO = GetDetSwitch("DIPO")->GetSwitchedVersion(); 
    G4int verSHIL = GetDetSwitch("SHIL")->GetSwitchedVersion(); 
    if ( verPIPE != 1 || verABSO !=0 || verDIPO == -1 || verSHIL == -1) {
      G4String text = "AliDetConstruction::CheckDetDependencies: \n";
      text = text + "    Switched ZDC requires PIPE, ABSO, DIPO and SHIL.\n"; 
      if (verPIPE == -1) {
        GetDetSwitch("PIPE")->SwitchOnDefault();
        text = text + "    The det switch for PIPE has been changed.\n"; 
      }  
      if (verABSO == -1) {
        GetDetSwitch("ABSO")->SwitchOnDefault();
        text = text + "    The det switch for ABSO has been changed.\n"; 
      }  
      if (verDIPO == -1) {
        GetDetSwitch("DIPO")->SwitchOnDefault();
        text = text + "    The det switch for DIPO has been changed.\n"; 
      }  
      if (verSHIL == -1) {
        GetDetSwitch("SHIL")->SwitchOnDefault();
        text = text + "    The det switch for SHIL has been changed."; 
      }  
      AliGlobals::Warning(text);
    }  
  }    
}  

// public methods

//_____________________________________________________________________________
G4VPhysicalVolume* AliDetConstruction::Construct()
{
// Constructs geometry.
// This method is called by G4RunManager in initialization.
// ---

  if (gAlice->Modules()->GetLast() < 0) {
    // create geometry (including AliModules) according to 
    // the fDetSwitchVector
    CheckDetDependencies();
    CreateDetectors();
  }   
  else {
    // create geometry for AliModules 
    // that have been created and registered by gAlice 
    BuildDetectors();
  }  
  // construct modules geometry
  ConstructModules();

  return AliSingleModuleConstruction::GetWorld();      
}

