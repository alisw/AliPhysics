// $Id$
// Category: geometry
//
// Author: I. Hrivnacova
//
// Class AliDetSwitchVector
// ---------------------------
// See the class description in the header file.

#include "AliDetSwitchVector.h"
#include "AliDetSwitch.h"
#include "AliGlobals.h"
#include "AliFiles.h"

//_____________________________________________________________________________
AliDetSwitchVector::AliDetSwitchVector()
  : fMessenger(this)
{
//
}

//_____________________________________________________________________________
AliDetSwitchVector::AliDetSwitchVector(const AliDetSwitchVector& right)
  : fMessenger(this)
{
//
  AliGlobals::Exception("AliDetSwitchVector is protected from copying.");  
}

//_____________________________________________________________________________
AliDetSwitchVector::~AliDetSwitchVector() {
//   
  // destroy det switch vector
  DetSwitchIterator it;
  for (it = fDetSwitchVector.begin(); it != fDetSwitchVector.end(); it++)
    delete *it; 
}

// operators

//_____________________________________________________________________________
AliDetSwitchVector& 
AliDetSwitchVector::operator=(const AliDetSwitchVector& right)
{
  // check assignement to self
  if (this == &right) return *this;
  
  AliGlobals::Exception("AliDetSwitchVector is protected from assigning.");  

  return *this;  
}    
          
// protected methods

//_____________________________________________________________________________
AliDetSwitch* 
AliDetSwitchVector::GetDetSwitch(const G4String& moduleName) const
{
// Returns the detector switch with given detector name.
// ---

  DetSwitchConstIterator it;
  for (it = fDetSwitchVector.begin(); it != fDetSwitchVector.end(); it++)
    if ((*it)->GetDetName() == moduleName) return *it; 

  G4String text = "AliDetSwitchVector::GetDetSwitch:\n";
  text = text + "Wrong detector name for " + moduleName;   
  AliGlobals::Exception(text);
  return 0;  
} 

// public methods

//_____________________________________________________________________________
void AliDetSwitchVector::Add(AliDetSwitch* detSwitch)
{
// Adds detSwitch to the detSwitch vector.
// ---

  fDetSwitchVector.push_back(detSwitch);
  fMessenger.Update();
}  
  
//_____________________________________________________________________________
void AliDetSwitchVector::SwitchDetOn(const G4String& moduleNameVer)
{ 
// Switchs on module specified by name and version.
// ---

  DetSwitchIterator it;

  if (moduleNameVer == "ALL") {
    for (it = fDetSwitchVector.begin(); it != fDetSwitchVector.end(); it++)
      (*it)->SwitchOnDefault(); 
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
void AliDetSwitchVector::SwitchDetOn(const G4String& moduleName, 
                                     G4int version)
{ 
// Switchs on module specified by name and version.
// ---

  GetDetSwitch(moduleName)->SwitchOn(version);
}

//_____________________________________________________________________________
void AliDetSwitchVector::SwitchDetOnDefault(const G4String& moduleName)
{ 
// Switchs on module specified by name with default version.
// ---

  GetDetSwitch(moduleName)->SwitchOnDefault();
}

//_____________________________________________________________________________
void AliDetSwitchVector::SwitchDetOff(const G4String& moduleName)
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
void AliDetSwitchVector::PrintSwitchedDets() const
{ 
// Lists switched detectors.
// ---

  G4String svList = GetSwitchedDetsList();
    
  G4cout << "Switched Alice detectors: " << G4endl;
  G4cout << "--------------------------" << G4endl;
  G4cout << svList << G4endl;
}

//_____________________________________________________________________________
void AliDetSwitchVector::PrintAvailableDets() const
{ 
// Lists available detectors.
// ---

  G4String avList = GetAvailableDetsList();
    
  G4cout << "Available Alice detectors: " << G4endl;
  G4cout << "---------------------------" << G4endl;
  G4cout << avList << G4endl;
}

//_____________________________________________________________________________
G4String AliDetSwitchVector::GetSwitchedDetsList() const
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

  if (nofSwitchedDets == G4int(fDetSwitchVector.size())) svList = "ALL: " + svList;
  if (nofSwitchedDets == 0) svList = "NONE";   

  return svList;
}

//_____________________________________________________________________________
G4String AliDetSwitchVector::GetAvailableDetsList() const
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
G4String AliDetSwitchVector::GetAvailableDetsListWithCommas() const
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
      if (iv < (*it)->GetNofVersions()-1)        
        svList += "/";
      else if (id++ < G4int(fDetSwitchVector.size())-1) 
        svList += ", ";
    }

  return svList;
}

//_____________________________________________________________________________
G4String AliDetSwitchVector::GetDetNamesList() const
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
G4String AliDetSwitchVector::GetDetNamesListWithCommas() const
{ 
// Returns list of detector names with commas.
// ---

  G4String svList = "";
  G4int id =0;
  DetSwitchConstIterator it;

  for (it = fDetSwitchVector.begin(); it != fDetSwitchVector.end(); it++) {
    svList += (*it)->GetDetName();
    if (id++ < G4int(fDetSwitchVector.size())-1)  
      svList += ", ";
  }

  return svList;
}

