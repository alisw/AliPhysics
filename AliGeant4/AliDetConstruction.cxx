// $Id$
// Category: geometry
//
// Author: I. Hrivnacova
//
// Class AliDetConstruction
// ------------------------
// See the class description in the header file.

#include "AliDetConstruction.h"
#include "AliDetSwitch.h"
#include "AliGlobals.h"
#include "AliFiles.h"
#include "AliRun.h"
#include "AliModule.h"

#include "TG4XMLGeometryGenerator.h"
#include "TG4GeometryServices.h"
#include "TG4LVTree.h"

#include <G4VPhysicalVolume.hh>

//_____________________________________________________________________________
AliDetConstruction::AliDetConstruction()
  : AliModulesComposition()
{
  // initialize det switch vector: 
  // moduleName nofVersions defaultVersion [type]
        // det switch objects are deleted in fDetSwitchVector destructor

  fDetSwitchVector.Add(new AliDetSwitch("MAG",    1, 0, kStructure));
  fDetSwitchVector.Add(new AliDetSwitch("ABSO",   1, 0, kStructure));
  fDetSwitchVector.Add(new AliDetSwitch("DIPO",   3, 2, kStructure));
  fDetSwitchVector.Add(new AliDetSwitch("FRAME",  3, 2, kStructure));
  fDetSwitchVector.Add(new AliDetSwitch("HALL",   1, 0, kStructure));
  fDetSwitchVector.Add(new AliDetSwitch("PIPE",   5, 0, kStructure));
  fDetSwitchVector.Add(new AliDetSwitch("SHIL",   2, 2, kStructure));
  fDetSwitchVector.Add(new AliDetSwitch("CRT",    1, 0));
  fDetSwitchVector.Add(new AliDetSwitch("EMCAL",  2, 1));
  fDetSwitchVector.Add(new AliDetSwitch("FMD",    2, 1));
  fDetSwitchVector.Add(new AliDetSwitch("ITS",    7, 5));
  fDetSwitchVector.Add(new AliDetSwitch("MUON",   2, 1));
  fDetSwitchVector.Add(new AliDetSwitch("PHOS",   2, 1));
  fDetSwitchVector.Add(new AliDetSwitch("PMD",    3, 1));
  fDetSwitchVector.Add(new AliDetSwitch("RICH",   3, 1));
  fDetSwitchVector.Add(new AliDetSwitch("START",  2, 1));
  fDetSwitchVector.Add(new AliDetSwitch("TOF",    5, 2));
  fDetSwitchVector.Add(new AliDetSwitch("TPC",    4, 2));
  fDetSwitchVector.Add(new AliDetSwitch("TRD",    2, 1));
  fDetSwitchVector.Add(new AliDetSwitch("ZDC",    3, 2));

  // update messenger
  fDetSwitchVector.UpdateMessenger();

  // instantiate LVtree browser
  TG4LVTree::Instance();
}

//_____________________________________________________________________________
AliDetConstruction::AliDetConstruction(const AliDetConstruction& right)
  : AliModulesComposition(right)
{
  // AliModuleComposition is protected from copying
}  

//_____________________________________________________________________________
AliDetConstruction::~AliDetConstruction() 
{
  // delete LVtree browser
  delete TG4LVTree::Instance();
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
  AddModule("BODY", 0, kStructure);

  G4bool first = true;
  while ((module = (AliModule*)next())) {
  
    // register moduleConstruction in fDetSwitchVector
    // in order to keep availability of /aliDet/list command
    G4String modName = module->GetName();
    G4int modVersion = module->IsVersion();
    if (first)
      // skip registering of the top volume 
      first = false;
    else 
      fDetSwitchVector.SwitchDetOn(modName, modVersion);
 
    // all modules will be processed alltogether
    AddModule(modName, modVersion, fDetSwitchVector.GetDetSwitch(modName)->GetType());

    if (VerboseLevel() > 0) {
      G4cout << "Created module construction for " 
             << modName << "v" << modVersion << "." << G4endl;   
    }	     
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
  AddModule("BODY", 0, kStructure);

  // add modules constructions
  for (G4int i=0; i<fDetSwitchVector.GetSize(); i++)
  {
    AliDetSwitch* detSwitch = fDetSwitchVector.GetDetSwitch(i);
    G4String detName = detSwitch->GetDetName();
    G4int version = detSwitch->GetSwitchedVersion();
    AliModuleType type = detSwitch->GetType();
    
    if (version > -1)
      AddModule(detName, version, type);
  }    
}

//_____________________________________________________________________________
void AliDetConstruction::CheckDependence(const G4String& master, 
                                         const G4String& slave)
{
// Checks modules dependence.
// If master is switch on and slave off, the default version
// of slave is switched on and a  warning is issued.
// ---

  AliDetSwitch* masterSwitch = fDetSwitchVector.GetDetSwitch(master);
  AliDetSwitch* slaveSwitch = fDetSwitchVector.GetDetSwitch(slave);

  if ( masterSwitch->GetSwitchedVersion() > -1 && 
       slaveSwitch->GetSwitchedVersion() < 0 ) {
     
    slaveSwitch->SwitchOnDefault();
    
    // warning
    G4String text = "AliDetConstruction::CheckDetDependence: \n";
    text = text + "    Switched " + master + " requires " + slave + ".\n"; 
    text = text + "    The det switch for " + slave + " has been changed."; 
    AliGlobals::Warning(text);
  }
}  
  
//_____________________________________________________________________________
void AliDetConstruction::CheckDetDependencies()
{
// Checks modules dependencies.
// ---

  CheckDependence("TOF", "FRAME");
  CheckDependence("TRD", "FRAME");
  CheckDependence("ZDC", "PIPE");
  CheckDependence("ZDC", "ABSO");
  CheckDependence("ZDC", "DIPO");
  CheckDependence("ZDC", "SHIL");
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

  return TG4GeometryServices::Instance()->GetWorld();      
}

//_____________________________________________________________________________
void AliDetConstruction::GenerateXMLGeometry() const 
{
// Generates XML geometry file from the top volume.
// The file name is set according the last switched detector
// registered in the det switch vector.
// ---

  G4VPhysicalVolume* world = TG4GeometryServices::Instance()->GetWorld();

  // XML filename
  // according to last switched detector
  G4String detName;
  G4String detVersion = "";
  G4int version = -1;
  for (G4int i=fDetSwitchVector.GetSize()-1; i>=0; i--) {
    version = fDetSwitchVector.GetDetSwitch(i)->GetSwitchedVersion();
    if (version > -1) {
      detName = fDetSwitchVector.GetDetSwitch(i)->GetDetName();
      AliGlobals::AppendNumberToString(detVersion,version); 
      break;
    }  
  }  
  G4String filePath 
    = AliFiles::Instance()->GetXMLFilePath(detName, version);
  
  // set top volume name
  G4String topName = world->GetName() + "_comp";
  
  // generate XML
  
  TG4XMLGeometryGenerator xml;
  xml.OpenFile(filePath);

  // generate materials 
  // not implemented
  // xml.GenerateMaterials(version, "today", "Generated from G4",
  //                       "v4", world->GetLogicalVolume());

  // generate volumes tree
  xml.GenerateSection("v6", detName, detVersion, "today", "Generated from Geant4",
                      topName, world->GetLogicalVolume());
  xml.CloseFile();
  
  if (VerboseLevel() > 0) {
    G4cout << "File " << detName << "v" << version << ".xml has been generated." 
           << G4endl;
  }	   
}  

