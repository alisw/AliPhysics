// $Id$
// Category: geometry
//
// Author: I. Hrivnacova
//
// Class AliModulesComposition
// ---------------------------
// See the class description in the header file.

#include "AliModulesComposition.h"
#include "AliSingleModuleConstruction.h"
#include "AliMoreModulesConstruction.h"
#include "AliDetSwitch.h"
#include "AliMagneticField.h"
#include "AliGlobals.h"
#include "AliFiles.h"

#include "TG4XMLGeometryGenerator.h"
#include "TG4GeometryManager.h"

#include <G4Material.hh>
#include <G4VPhysicalVolume.hh>

//_____________________________________________________________________________
AliModulesComposition::AliModulesComposition()
  : fReadGeometry(false),
    fWriteGeometry(false),
    fMagneticField(0),
    fMessenger(this)
{
//
  fMoreModulesConstruction = new AliMoreModulesConstruction();
}

//_____________________________________________________________________________
AliModulesComposition::AliModulesComposition(const AliModulesComposition& right)
  : fMessenger(this)
{
//
  AliGlobals::Exception("AliModulesComposition is protected from copying.");  
}

//_____________________________________________________________________________
AliModulesComposition::~AliModulesComposition() {
//   
  delete fMoreModulesConstruction;
  delete fMagneticField;
  
  // destroy det switch vector
  DetSwitchIterator it;
  for (it = fDetSwitchVector.begin(); it != fDetSwitchVector.end(); it++)
    delete *it; 
  
  // destroy det construction vector
  SingleModuleIterator itm;
  for (itm = fModuleConstructionVector.begin(); 
       itm != fModuleConstructionVector.end(); it++)
    delete *itm;
}

// operators

//_____________________________________________________________________________
AliModulesComposition& 
AliModulesComposition::operator=(const AliModulesComposition& right)
{
  // check assignement to self
  if (this == &right) return *this;
  
  AliGlobals::Exception("AliModulesComposition is protected from assigning.");  

  return *this;  
}    
          
// protected methods

//_____________________________________________________________________________
void AliModulesComposition::AddDetSwitch(AliDetSwitch* detSwitch)
{
// Adds detSwitch to the detSwitch vector.
// ---

  fDetSwitchVector.push_back(detSwitch);
  fMessenger.SetCandidates();
}  
  
//_____________________________________________________________________________
void AliModulesComposition::AddSingleModuleConstruction(
                               const G4String& name, G4int version, 
			       AliModuleType moduleType)
{
// Adds SingleModuleConstruction.
// ---

  fModuleConstructionVector
    .push_back(new AliSingleModuleConstruction(name, version, moduleType));
}  
		  
//_____________________________________________________________________________
void AliModulesComposition::AddMoreModuleConstruction(
                               const G4String& name, G4int version, 
			       AliModuleType moduleType)
{
// Adds module to MoreModulesConstruction (construction of dependent
// modules.)
// ---

  fMoreModulesConstruction->AddModule(name, version, moduleType);
}  
		  		  
//_____________________________________________________________________________
void AliModulesComposition::ConstructModules()
{
// Construct geometry of all modules (both standalone and dependent.)
// ---

  // set common options
  SetReadGeometryToModules(fReadGeometry);
  SetWriteGeometryToModules(fWriteGeometry);
  
  // configure single modules
  SingleModuleIterator it;
  for (it  = fModuleConstructionVector.begin(); 
       it != fModuleConstructionVector.end(); it++) {
       
    (*it)->Configure(*AliFiles::Instance());
    cout << "Module " << (*it)->GetDetName() << " configured." << endl;
  }  
     
  // configure dependent modules
  if (fMoreModulesConstruction->GetNofModules() > 0)
    fMoreModulesConstruction->Configure(*AliFiles::Instance());

  // construct single modules
  for (it  = fModuleConstructionVector.begin(); 
       it != fModuleConstructionVector.end(); it++) {

    G4cout << "Module " << (*it)->GetDetName()
           << " will be constructed now." << G4endl;
    (*it)->Construct();
  }  
    
  // construct dependent modules
  if (fMoreModulesConstruction->GetNofModules() > 0) {
    G4cout << "Dependent modules will be constructed now." << G4endl;
    fMoreModulesConstruction->Construct();
  }  
}  

//_____________________________________________________________________________
AliDetSwitch* AliModulesComposition::GetDetSwitch(
                                        const G4String& moduleName) const
{
// Returns the detector switch with given detector name.
// ---

  DetSwitchConstIterator it;
  for (it = fDetSwitchVector.begin(); it != fDetSwitchVector.end(); it++)
    if ((*it)->GetDetName() == moduleName) return *it; 

  G4String text = "AliModulesComposition::GetDetSwitch:\n";
  text = text + "Wrong detector name for " + moduleName;   
  AliGlobals::Exception(text);
  return 0;  
} 

//_____________________________________________________________________________
void AliModulesComposition::SetReadGeometryToModules(G4bool readGeometry)
{
// Sets readGeometry control to all modules.
// ---

  // single module constructions
  SingleModuleIterator it;
  for (it  = fModuleConstructionVector.begin(); 
       it != fModuleConstructionVector.end(); it++)        
    (*it)->SetReadGeometry(readGeometry);

  // more modules construction
  for (G4int i=0; i<fMoreModulesConstruction->GetNofModules(); i++) 
    fMoreModulesConstruction
      ->GetModuleConstruction(i)->SetReadGeometry(readGeometry);
}    
  
//_____________________________________________________________________________
void AliModulesComposition::SetWriteGeometryToModules(G4bool writeGeometry)
{
// Sets writeGeometry control to all modules.
// ---

  // single module constructions
  SingleModuleIterator it;
  for (it  = fModuleConstructionVector.begin(); 
       it != fModuleConstructionVector.end(); it++)        
    (*it)->SetWriteGeometry(writeGeometry);

  // more modules construction
  for (G4int i=0; i<fMoreModulesConstruction->GetNofModules(); i++) 
    fMoreModulesConstruction
      ->GetModuleConstruction(i)->SetWriteGeometry(writeGeometry);
}    

//_____________________________________________________________________________
void AliModulesComposition::SetProcessConfigToModules(G4bool processConfig)
{
// Sets processConfig control to all modules.
// ---

  // single module constructions
  SingleModuleIterator it;
  for (it  = fModuleConstructionVector.begin(); 
       it != fModuleConstructionVector.end(); it++)        
    (*it)->SetProcessConfig(processConfig);
  
  // more modules construction
  for (G4int i=0; i<fMoreModulesConstruction->GetNofModules(); i++) 
    fMoreModulesConstruction
      ->GetModuleConstruction(i)->SetProcessConfig(processConfig);
}    

// public methods

//_____________________________________________________________________________
void AliModulesComposition::SwitchDetOn(const G4String& moduleNameVer)
{ 
// Switchs on module specified by name and version.
// ---

  DetSwitchIterator it;

  if (moduleNameVer == "ALL") {
    for (it = fDetSwitchVector.begin(); it != fDetSwitchVector.end(); it++)
      (*it)->SwitchOnDefault(); 
  }
  else if (moduleNameVer == "PPR") {
    for (it = fDetSwitchVector.begin(); it != fDetSwitchVector.end(); it++)
      (*it)->SwitchOnPPR(); 
    AliFiles::Instance()->SetMacroName("ConfigPPR");
  }
  else if (moduleNameVer == "NONE") {
    for (it = fDetSwitchVector.begin(); it != fDetSwitchVector.end(); it++)
      (*it)->SwitchOff(); 
  }
  else {
    // get version number
    G4int len = moduleNameVer.length();
    G4String moduleName = moduleNameVer.substr(0, len-1);
    G4String version = moduleNameVer.substr(len-1, 1);
    G4int iVersion = AliGlobals::StringToInt(version);

    if (iVersion < 0) {
      // in case the version number is not provided
      // the default one is set
      SwitchDetOnDefault(moduleNameVer);
    }  
    else 
      SwitchDetOn(moduleName, iVersion);
  }
}

//_____________________________________________________________________________
void AliModulesComposition::SwitchDetOn(const G4String& moduleName, 
                                        G4int version)
{ 
// Switchs on module specified by name and version.
// ---

  GetDetSwitch(moduleName)->SwitchOn(version);
}

//_____________________________________________________________________________
void AliModulesComposition::SwitchDetOnDefault(const G4String& moduleName)
{ 
// Switchs on module specified by name with default version.
// ---

  GetDetSwitch(moduleName)->SwitchOnDefault();
}

//_____________________________________________________________________________
void AliModulesComposition::SwitchDetOnPPR(const G4String& moduleName)
{ 
// Switchs on module specified by name with PPR version.
// ---

  GetDetSwitch(moduleName)->SwitchOnPPR();
}

//_____________________________________________________________________________
void AliModulesComposition::SwitchDetOff(const G4String& moduleName)
{ 
// Switchs off module specified by name.
// ---

  if (moduleName == "ALL") {
    DetSwitchIterator it;
    for (it = fDetSwitchVector.begin(); it != fDetSwitchVector.end(); it++)
      (*it)->SwitchOff(); 
  }
  else 
    GetDetSwitch(moduleName)->SwitchOff();
}

//_____________________________________________________________________________
void AliModulesComposition::PrintSwitchedDets() const
{ 
// Lists switched detectors.
// ---

  G4String svList = GetSwitchedDetsList();
    
  G4cout << "Switched Alice detectors: " << G4endl;
  G4cout << "--------------------------" << G4endl;
  G4cout << svList << G4endl;
}

//_____________________________________________________________________________
void AliModulesComposition::PrintAvailableDets() const
{ 
// Lists available detectors.
// ---

  G4String avList = GetAvailableDetsList();
    
  G4cout << "Available Alice detectors: " << G4endl;
  G4cout << "---------------------------" << G4endl;
  G4cout << avList << G4endl;
}

//_____________________________________________________________________________
void AliModulesComposition::PrintMaterials() const
{
// Prints all materials.
// ---

  const G4MaterialTable* matTable = G4Material::GetMaterialTable();
  G4cout << *matTable;
}

//_____________________________________________________________________________
void AliModulesComposition::GenerateXMLGeometry() const 
{
// Generates XML geometry file from the top volume.
// The file name is set according the last switched detector
// registered in the det switch vector.
// ---

  G4VPhysicalVolume* world = AliSingleModuleConstruction::GetWorld();

  // XML filename
  // according to last switched detector
  G4String detName;
  G4String detVersion = "";
  G4int version = -1;
  for (G4int i=fDetSwitchVector.size()-1; i>=0; i--) {
    version = fDetSwitchVector[i]->GetSwitchedVersion();
    if (version > -1) {
      detName = fDetSwitchVector[i]->GetDetName();
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
  // not yet implemented
  // xml.GenerateMaterials(version, "today", "Generated from G4",
  //                     "v4", world->GetLogicalVolume());

  // generate volumes tree
  xml.GenerateSection(detName, detVersion, "today", "Generated from Geant4",
                      topName, world->GetLogicalVolume());
  xml.CloseFile();
  
  // set verbose
  G4cout << "File " << detName << "v" << version << ".xml has been generated." 
         << G4endl;
}  

//_____________________________________________________________________________
G4String AliModulesComposition::GetSwitchedDetsList() const
{ 
// Returns list of switched detectors.
// ---

  G4String svList = "";  
  G4int nofSwitchedDets = 0;
  DetSwitchConstIterator it;
  
  for (it = fDetSwitchVector.begin(); it != fDetSwitchVector.end(); it++) {
    G4int iVersion = (*it)->GetSwitchedVersion();
    if (iVersion > -1) {
      nofSwitchedDets++;
      G4String moduleNameVer = (*it)->GetDetName();
      AliGlobals::AppendNumberToString(moduleNameVer, iVersion);
      svList += moduleNameVer;
      svList += " "; 
    }
  }

  if (nofSwitchedDets == fDetSwitchVector.size()) svList = "ALL: " + svList;
  if (nofSwitchedDets == 0) svList = "NONE";   

  return svList;
}

//_____________________________________________________________________________
G4String AliModulesComposition::GetAvailableDetsList() const
{ 
// Returns list of available detectors.
// ---

  G4String svList = "";
  DetSwitchConstIterator it;
  
  for (it = fDetSwitchVector.begin(); it != fDetSwitchVector.end(); it++)
    for (G4int iv=0; iv<(*it)->GetNofVersions(); iv++) {
      G4String moduleNameVer = (*it)->GetDetName();
      AliGlobals::AppendNumberToString(moduleNameVer, iv);
      svList += moduleNameVer;
      svList += " ";
    } 

  return svList;
}

//_____________________________________________________________________________
G4String AliModulesComposition::GetAvailableDetsListWithCommas() const
{ 
// Returns list of available detectors with commas.
// ---

  G4String svList = "";
  G4int id =0;
  DetSwitchConstIterator it;

  for (it = fDetSwitchVector.begin(); it != fDetSwitchVector.end(); it++)
    for (G4int iv=0; iv<(*it)->GetNofVersions(); iv++) {
      G4String moduleNameVer = (*it)->GetDetName();
      AliGlobals::AppendNumberToString(moduleNameVer, iv);
      svList += moduleNameVer;
      if (iv < (*it)->GetNofVersions()-1)        svList += "/";
      else if (id++ < fDetSwitchVector.size()-1) svList += ", ";
    }

  return svList;
}

//_____________________________________________________________________________
G4String AliModulesComposition::GetDetNamesList() const
{ 
// Returns list of detector names.
// ---

  G4String svList = "";
  DetSwitchConstIterator it;
  
  for (it = fDetSwitchVector.begin(); it != fDetSwitchVector.end(); it++) {
    svList += (*it)->GetDetName();
    svList += " ";
  }

  return svList;
}

//_____________________________________________________________________________
G4String AliModulesComposition::GetDetNamesListWithCommas() const
{ 
// Returns list of detector names with commas.
// ---

  G4String svList = "";
  G4int id =0;
  DetSwitchConstIterator it;

  for (it = fDetSwitchVector.begin(); it != fDetSwitchVector.end(); it++) {
    svList += (*it)->GetDetName();
    if (id++ < fDetSwitchVector.size()-1) svList += ", ";
  }

  return svList;
}

//_____________________________________________________________________________
void AliModulesComposition::SetMagField(G4double fieldValue)
{
// Sets uniform magnetic field to specified value.
// ---

  // create fields if it does not exist
  if (!fMagneticField) fMagneticField = new AliMagneticField();
  
  // set value
  fMagneticField->SetFieldValue(fieldValue);
}

