// $Id$
// Category: geometry
//
// See the class description in the header file.

#include "AliModulesComposition.h"
#include "AliModulesCompositionMessenger.h"
#include "AliSingleModuleConstruction.h"
#include "AliMoreModulesConstruction.h"
#include "AliMagneticField.h"
#include "AliGlobals.h"

#include <G4UniformMagField.hh>
#include <G4FieldManager.hh>
#include <G4TransportationManager.hh>

AliModulesComposition::AliModulesComposition()
  : fAllLVSensitive(false),
    fReadGeometry(false),
    fWriteGeometry(false)    
{
//
  fMoreModulesConstruction = new AliMoreModulesConstruction();
  fMagneticField = new AliMagneticField();
  fMessenger = new AliModulesCompositionMessenger(this);
}

AliModulesComposition::AliModulesComposition(const AliModulesComposition& right)
{
//
  AliGlobals::Exception("AliModulesComposition is protected from copying.");  
}

AliModulesComposition::~AliModulesComposition() {
//   
  delete fMoreModulesConstruction;
  delete fMagneticField;
  delete fMessenger; 
  
  // destroy det switch vector
  fDetSwitchVector.clearAndDestroy();  
  
  // destroy det construction vector
  fModuleConstructionVector.clearAndDestroy();
}

// operators

AliModulesComposition& 
AliModulesComposition::operator=(const AliModulesComposition& right)
{
  // check assignement to self
  if (this == &right) return *this;
  
  AliGlobals::Exception("AliModulesComposition is protected from assigning.");  

  return *this;  
}    
          
// protected methods

void AliModulesComposition::AddDetSwitch(AliDetSwitch* detSwitch)
{
// Adds detSwitch to the detSwitch vector.
// ---

  fDetSwitchVector.insert(detSwitch);
  fMessenger->SetCandidates();
}  
  
void AliModulesComposition::AddSingleModuleConstruction(G4String moduleName, 
                               G4int version, AliModuleType moduleType)
{
// Adds SingleModuleConstruction.
// ---

  AliSingleModuleConstruction* moduleConstruction
    = new AliSingleModuleConstruction(moduleName, version, moduleType);
  fModuleConstructionVector.insert(moduleConstruction);
}  
		  
void AliModulesComposition::AddMoreModuleConstruction(G4String moduleName, 
                               G4int version, AliModuleType moduleType)
{
// Adds module to MoreModulesConstruction (construction of dependent
// modules.)
// ---

  fMoreModulesConstruction->AddModule(moduleName, version, moduleType);
}  
		  		  
void AliModulesComposition::ConstructModules()
{
// Construct geometry of all modules (both standalone and dependent.)
// ---

  // set common options
  SetReadGeometryToModules(fReadGeometry);
  SetWriteGeometryToModules(fWriteGeometry);
  SetAllLVSensitiveToModules(fAllLVSensitive);
     // common setAllLVSensitive is overridden by Config.in
     // macro
  
  // one module constructions
  G4int nofDets = fModuleConstructionVector.entries();
  for (G4int i=0; i<nofDets; i++) {
    G4cout << "Module " << fModuleConstructionVector[i]->GetDetName()
           << " will be constructed now." << endl;
    fModuleConstructionVector[i]->Construct();
  }  
    
  // more modules construction
  G4int nofModules = fMoreModulesConstruction->GetNofModules();
  if (nofModules>0) {
    G4cout << "Dependent modules will be constructed now." << endl;
    fMoreModulesConstruction->Construct();
  }  
}  

void AliModulesComposition::SetReadGeometryToModules(G4bool readGeometry)
{
// Sets readGeometry control to all modules.
// ---

  // single module constructions
  G4int nofDets = fModuleConstructionVector.entries();
  G4int i;
  for (i=0; i<nofDets; i++)
    fModuleConstructionVector[i]->SetReadGeometry(readGeometry);

  // more modules construction
  nofDets = fMoreModulesConstruction->GetNofModules();
  for (i=0; i<nofDets; i++) { 
    AliSingleModuleConstruction* moduleConstruction
      = fMoreModulesConstruction->GetModuleConstruction(i);
    moduleConstruction->SetReadGeometry(readGeometry);
  }  
}    
  
void AliModulesComposition::SetWriteGeometryToModules(G4bool writeGeometry)
{
// Sets writeGeometry control to all modules.
// ---

  // single module constructions
  G4int nofDets = fModuleConstructionVector.entries();
  G4int i;
  for (i=0; i<nofDets; i++)
    fModuleConstructionVector[i]->SetWriteGeometry(writeGeometry);

  // more modules construction
  nofDets = fMoreModulesConstruction->GetNofModules();
  for (i=0; i<nofDets; i++) { 
    AliSingleModuleConstruction* moduleConstruction
      = fMoreModulesConstruction->GetModuleConstruction(i);
    moduleConstruction->SetWriteGeometry(writeGeometry);
  }  
}    

void AliModulesComposition::SetAllLVSensitiveToModules(G4bool allSensitive)
{
// Sets setAllSensitive control to all modules.
// ---

  // single module constructions
  G4int nofDets = fModuleConstructionVector.entries();
  G4int i;
  for (i=0; i<nofDets; i++) {
    AliSingleModuleConstruction* moduleConstruction
      = dynamic_cast<AliSingleModuleConstruction*>(fModuleConstructionVector[i]);
    if (!moduleConstruction) {
      G4String text = "AliModulesComposition::SetProcessConfig: \n";
      text = text + "    Unknown module construction type.";
      AliGlobals::Exception(text);
    }      
    moduleConstruction->SetAllLVSensitive(allSensitive);
  }  

  // more modules construction
  nofDets = fMoreModulesConstruction->GetNofModules();
  for (i=0; i<nofDets; i++) { 
    AliSingleModuleConstruction* moduleConstruction
      = fMoreModulesConstruction->GetModuleConstruction(i);
    moduleConstruction->SetAllLVSensitive(allSensitive);
  }  
}    

void AliModulesComposition::SetProcessConfigToModules(G4bool processConfig)
{
// Sets processConfig control to all modules.
// ---

  // single module constructions
  G4int nofDets = fModuleConstructionVector.entries();
  G4int i;
  for (i=0; i<nofDets; i++) {
    AliSingleModuleConstruction* moduleConstruction
      = dynamic_cast<AliSingleModuleConstruction*>(fModuleConstructionVector[i]);
    if (!moduleConstruction) {
      G4String text = "AliModulesComposition::SetProcessConfig: \n";
      text = text + "    Unknown module construction type.";
      AliGlobals::Exception(text);
    }      
    moduleConstruction->SetProcessConfig(processConfig);
  }  
  
  // more modules construction
  nofDets = fMoreModulesConstruction->GetNofModules();
  for (i=0; i<nofDets; i++) { 
    AliSingleModuleConstruction* moduleConstruction
      = fMoreModulesConstruction->GetModuleConstruction(i);
    moduleConstruction->SetProcessConfig(processConfig);
  }    
}    

// public methods

void AliModulesComposition::SwitchDetOn(G4String moduleNameVer)
{ 
// Switchs on module specified by name and version.
// ---

  G4int nofDets = fDetSwitchVector.entries(); 
  if (moduleNameVer == "ALL") {
    for (G4int id=0; id<nofDets; id++) {
      G4int defaultVersion = fDetSwitchVector[id]->GetDefaultVersion();
      fDetSwitchVector[id]->SwitchOn(defaultVersion); 
    }  
  }
  else if (moduleNameVer == "NONE") {
    for (G4int id=0; id<nofDets; id++)
      fDetSwitchVector[id]->SwitchOff(); 
  }
  else {
    // get version number
    G4int len = moduleNameVer.length();
    G4String moduleName = moduleNameVer(0, len-1);
    G4String version = moduleNameVer(len-1, 1);
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

void AliModulesComposition::SwitchDetOn(G4String moduleName, G4int version)
{ 
// Switchs on module specified by name and version.
// ---

  G4int nofDets = fDetSwitchVector.entries(); 
  for (G4int id=0; id<nofDets; id++) {
    G4String iDetName = fDetSwitchVector[id]->GetDetName();
    if (iDetName == moduleName) {
      fDetSwitchVector[id]->SwitchOn(version);
      return;
    }
  }
  AliGlobals::Exception(
    "AliModulesComposition: Wrong detector name for " + moduleName + ".");
}

void AliModulesComposition::SwitchDetOnDefault(G4String moduleName)
{ 
// Switchs on module specified by name with default version.
// ---

  G4int nofDets = fDetSwitchVector.entries(); 
  for (G4int id=0; id<nofDets; id++) {
    G4String iDetName = fDetSwitchVector[id]->GetDetName();
    if (iDetName == moduleName) {
      fDetSwitchVector[id]->SwitchOnDefault();
      return;
    }
  }
  AliGlobals::Exception(
    "AliModulesComposition: Wrong detector name for " + moduleName + ".");
}

void AliModulesComposition::SwitchDetOff(G4String moduleName)
{ 
// Switchs off module specified by name.
// ---

  G4int nofDets = fDetSwitchVector.entries(); 
  if (moduleName == "ALL") {
    for (G4int id=0; id<nofDets; id++)
      fDetSwitchVector[id]->SwitchOff(); 
  }
  else {
    for (G4int id=0; id<nofDets; id++) {
      G4String iDetName = fDetSwitchVector[id]->GetDetName();
      if (iDetName == moduleName) { 
        fDetSwitchVector[id]->SwitchOff();
        return;
      }
    }
  }
  AliGlobals::Exception(
    "AliModulesComposition: Wrong detector name for " + moduleName + ".");
}

void AliModulesComposition::PrintSwitchedDets() const
{ 
// Lists switched detectors.
// ---

  G4String svList = GetSwitchedDetsList();
    
  G4cout << "Switched Alice detectors: " << endl;
  G4cout << "--------------------------" << endl;
  G4cout << svList << endl;
}

void AliModulesComposition::PrintAvailableDets() const
{ 
// Lists available detectors.
// ---

  G4String avList = GetAvailableDetsList();
    
  G4cout << "Available Alice detectors: " << endl;
  G4cout << "---------------------------" << endl;
  G4cout << avList << endl;
}

G4String AliModulesComposition::GetSwitchedDetsList() const
{ 
// Returns list of switched detectors.
// ---

  G4String svList = "";
  
  G4int nofDets = fDetSwitchVector.entries(); 
  G4int nofSwitchedDets = 0;
  for (G4int id=0; id<nofDets; id++) {
    G4int iVersion = fDetSwitchVector[id]->GetSwitchedVersion();
    if (iVersion > -1) {
      nofSwitchedDets++;
      G4String moduleNameVer = fDetSwitchVector[id]->GetDetName();
      AliGlobals::AppendNumberToString(moduleNameVer, iVersion);
      svList += moduleNameVer;
      svList += " "; 
    }
  }

  if (nofSwitchedDets ==  nofDets) svList = "ALL: " + svList;
  if (nofSwitchedDets ==  0)       svList = "NONE";   

  return svList;
}

const G4RWTPtrOrderedVector<AliDetSwitch>& 
AliModulesComposition::GetDetSwitchVector() const
{
// Returns detSwitch vector.
// ---

  //const AliDetSwitchVector& vector = fDetSwitchVector;
  const G4RWTPtrOrderedVector<AliDetSwitch>& vector = fDetSwitchVector;
  return vector;
}  

G4String AliModulesComposition::GetAvailableDetsList() const
{ 
// Returns list of available detectors.
// ---

  G4String svList = "";
  
  G4int nofDets = fDetSwitchVector.entries(); 
  for (G4int id=0; id<nofDets; id++) {
    G4int nofVersions = fDetSwitchVector[id]->GetNofVersions();
    for (G4int iv=0; iv<nofVersions; iv++) {
      G4String moduleNameVer = fDetSwitchVector[id]->GetDetName();
      AliGlobals::AppendNumberToString(moduleNameVer, iv);
      svList += moduleNameVer;
      svList += " ";
    }
  }

  return svList;
}

G4String AliModulesComposition::GetAvailableDetsListWithCommas() const
{ 
// Returns list of available detectors with commas.
// ---

  G4String svList = "";
  
  G4int nofDets = fDetSwitchVector.entries(); 
  for (G4int id=0; id<nofDets; id++) {
    G4int nofVersions = fDetSwitchVector[id]->GetNofVersions();
    for (G4int iv=0; iv<nofVersions; iv++) {
      G4String moduleNameVer = fDetSwitchVector[id]->GetDetName();
      AliGlobals::AppendNumberToString(moduleNameVer, iv);
      svList += moduleNameVer;
      if (iv<nofVersions-1)    svList += "/";
      else if (id < nofDets-1) svList += ", ";
    }
  }

  return svList;
}

G4String AliModulesComposition::GetDetNamesList() const
{ 
// Returns list of detector names.
// ---

  G4String svList = "";
  
  G4int nofDets = fDetSwitchVector.entries(); 
  for (G4int id=0; id<nofDets; id++) { 
    svList += fDetSwitchVector[id]->GetDetName();
    svList += " ";
  }

  return svList;
}

G4String AliModulesComposition::GetDetNamesListWithCommas() const
{ 
// Returns list of detector names with commas.
// ---

  G4String svList = "";
  
  G4int nofDets = fDetSwitchVector.entries(); 
  for (G4int id=0; id<nofDets; id++) { 
    svList += fDetSwitchVector[id]->GetDetName();
    if (id < nofDets-1) svList += ", ";
  }

  return svList;
}

void AliModulesComposition::SetMagField(G4double fieldValue)
{
// Sets uniform magnetic field to specified value.
// ---

  fMagneticField->SetFieldValue(fieldValue);
}

