// $Id$
// Category: geometry
//
// See the class description in the header file.

#include "AliDetConstruction.h"
#include "AliSingleModuleConstruction.h"
#include "AliGlobals.h"
#include "AliRun.h"
#include "AliModule.h"

AliDetConstruction::AliDetConstruction()
  : fTopVolumeName("ALIC")
{
  // initialize det switch vector: 
  // moduleName nofVersions defaultVersion [type isStandalone] 

  AliDetSwitch* detSwitch;
  detSwitch = new AliDetSwitch("ABSO",   1, 0, kStructure);
  AddDetSwitch(detSwitch); 
  detSwitch = new AliDetSwitch("DIPO",   3, 2, kStructure);
  AddDetSwitch(detSwitch); 
  detSwitch = new AliDetSwitch("FRAME",  2, 1, kStructure, false);
  AddDetSwitch(detSwitch); 
  detSwitch = new AliDetSwitch("HALL",   1, 0, kStructure);
  AddDetSwitch(detSwitch); 
  detSwitch = new AliDetSwitch("MAG",    1, 0, kStructure);
  AddDetSwitch(detSwitch); 
  detSwitch = new AliDetSwitch("PIPE",   4, 0, kStructure);
  AddDetSwitch(detSwitch); 
  detSwitch = new AliDetSwitch("SHIL",   1, 0, kStructure);
  AddDetSwitch(detSwitch); 
  detSwitch = new AliDetSwitch("CASTOR", 2, 1);
  AddDetSwitch(detSwitch); 
  detSwitch = new AliDetSwitch("FMD",    2, 1);
  AddDetSwitch(detSwitch); 
  detSwitch = new AliDetSwitch("ITS",    6, 5);
  AddDetSwitch(detSwitch); 
  detSwitch = new AliDetSwitch("MUON",   2, 0);
  AddDetSwitch(detSwitch); 
  detSwitch = new AliDetSwitch("PHOS",   5, 1);
  AddDetSwitch(detSwitch); 
  detSwitch = new AliDetSwitch("PMD",    3, 0);
  AddDetSwitch(detSwitch); 
  detSwitch = new AliDetSwitch("RICH",   3, 1);
  AddDetSwitch(detSwitch); 
  detSwitch = new AliDetSwitch("START",  2, 0);
  AddDetSwitch(detSwitch); 
  detSwitch = new AliDetSwitch("TOF",    5, 1, kDetector, false);
  AddDetSwitch(detSwitch); 
  detSwitch = new AliDetSwitch("TPC",    4, 1);
  AddDetSwitch(detSwitch); 
  detSwitch = new AliDetSwitch("TRD",    2, 0, kDetector, false);
  AddDetSwitch(detSwitch); 
  detSwitch = new AliDetSwitch("ZDC",    2, 1);
  AddDetSwitch(detSwitch);  
}

AliDetConstruction::AliDetConstruction(const AliDetConstruction& right)
  : AliModulesComposition(right)
{
  // AliModuleComposition is protected from copying
}  

AliDetConstruction::~AliDetConstruction() {
//
}

// operators

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

void AliDetConstruction::CreateDetectors()
{
// Creates AliModules and their module constructions 
// according to the fDetSwitchVector
// ---

  // add top volume (AliBODY) construction first
  AddSingleModuleConstruction("BODY", 0, kStructure);

  // add modules constructions
  const G4RWTPtrOrderedVector<AliDetSwitch>& krDetSwitchVector 
    = GetDetSwitchVector();
  for (G4int id=0; id<krDetSwitchVector.entries(); id++)
  {
    G4String detName = krDetSwitchVector[id]->GetDetName();
    G4int version = krDetSwitchVector[id]->GetSwitchedVersion();
    G4bool isStandalone = krDetSwitchVector[id]->IsStandalone();
    AliModuleType type = krDetSwitchVector[id]->GetType();
    
    if (version > -1)
      if (isStandalone)
        AddSingleModuleConstruction(detName, version, type);
      else
        AddMoreModuleConstruction(detName, version, type);
  }    
}

void AliDetConstruction::CheckDetDependencies()
{
// Checks modules dependencies.
// Dependent modules FRAME, TOF, TRD 
// TOF always requires FRAMEv1
// TRD can be built with both (??)
// ---

  const G4RWTPtrOrderedVector<AliDetSwitch>& krDetSwitchVector 
    = GetDetSwitchVector();

  // get switched versions of dependent modules
  G4int nofDets = krDetSwitchVector.entries();
  G4int verFRAME = -1; 
  G4int verTOF = -1; 
  G4int verTRD = -1; 
  AliDetSwitch* detSwitchFRAME = 0;
  for (G4int id=0; id<nofDets; id++) {  
    G4String detName = krDetSwitchVector[id]->GetDetName();
    if (detName == "FRAME") { 
      verFRAME = krDetSwitchVector[id]->GetSwitchedVersion();  
      detSwitchFRAME = krDetSwitchVector[id];
    }  
    if (detName == "TOF")  
      verTOF = krDetSwitchVector[id]->GetSwitchedVersion();  
    if (detName == "TRD")  
      verTRD = krDetSwitchVector[id]->GetSwitchedVersion();  
  }
  
  // check dependencies  
  if (verTRD > -1 && verTOF > -1) {
    // both TRD and TOF 
    if (verFRAME != 1) {
      detSwitchFRAME->SwitchOn(1);
      G4String text = "AliDetConstruction::CheckDetDependencies: \n";
      text = text + "    Switched TOF and TRD require FRAME v1.\n"; 
      text = text + "    The det switch for FRAME has been changed."; 
      AliGlobals::Warning(text);
    }  
  }
  else if (verTRD > -1 && verTOF == -1)   {
    // only TRD
    if (verFRAME < 0) {
      detSwitchFRAME->SwitchOn(1);
      G4String text = "AliDetConstruction::CheckDetDependencies: \n";
      text = text + "    Switched TRD require FRAME.\n"; 
      text = text + "    The det switch for FRAME has been changed."; 
      AliGlobals::Warning(text);
    }  
  }  
  else if (verTRD == -1 && verTOF > -1)   {
    // only TOF
    if (verFRAME != 1) {
      detSwitchFRAME->SwitchOn(1);
      G4String text = "AliDetConstruction::CheckDetDependencies: \n";
      text = text + "    Switched TOF requires FRAME v1.\n"; 
      text = text + "    The det switch for FRAME has been changed."; 
      AliGlobals::Warning(text);
    }  
  }
/*  
  if (verTRD > -1 && verTOF > -1) {
    // both TRD and TOF 
    if (verTOF == 2 || verTOF == 3 || verTOF == 5 || verTOF == 6) {
      G4String text = "AliDetConstruction::CheckDetDependencies: \n";
      text = text + "    Switched TOF and TRD require different FRAME versions."; 
      AliGlobals::Exception(text);
    }  
    if (verFRAME != 0) {
      detSwitchFRAME->SwitchOn(0);
      G4String text = "AliDetConstruction::CheckDetDependencies: \n";
      text = text + "    Switched TOF and TRD require FRAME v0.\n"; 
      text = text + "    The det switch for FRAME has been changed."; 
      AliGlobals::Warning(text);
    }  
  }
  else if (verTRD > -1 && verTOF == -1)   {
    // only TRD
    if (verFRAME != 0) {
      detSwitchFRAME->SwitchOn(0);
      G4String text = "AliDetConstruction::CheckDetDependencies: \n";
      text = text + "    Switched TRD requires FRAME v0.\n"; 
      text = text + "    The det switch for FRAME has been changed."; 
      AliGlobals::Warning(text);
    }          
  }  
  else if (verTRD == -1 && verTOF > -1)   {
    // only TOF
    if ((verTOF == 0 || verTOF == 1 || verTOF == 4) && (verFRAME !=0)) {
      detSwitchFRAME->SwitchOn(0);
      G4String text = "AliDetConstruction::CheckDetDependencies: \n";
      text = text + "    Switched TOF requires FRAME v0.\n"; 
      text = text + "    The det switch for FRAME has been changed."; 
      AliGlobals::Warning(text);
    }
    if ((verTOF == 2 || verTOF == 3 || verTOF == 5 || verTOF == 6) &&
        (verFRAME != 1)) {
      detSwitchFRAME->SwitchOn(1);
      G4String text = "AliDetConstruction::CheckDetDependencies: \n";
      text = text + "    Switched TOF requires FRAME v1.\n"; 
      text = text + "    The det switch for FRAME has been changed."; 
      AliGlobals::Warning(text);
    }
  }
*/    
}  

// public methods

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

