// $Id$
// Category: global
//
// Author: I. Hrivnacova
//
// Class TG4ProcessControlMap
// --------------------------
// See the class description in the header file.

#include "TG4ProcessControlMap.h"
#include "TG4G3ControlVector.h"
#include "TG4Globals.h"

#include <G4VProcess.hh>
#include "g4std/iomanip"
#include "globals.hh"

TG4ProcessControlMap* TG4ProcessControlMap::fgInstance = 0;

//_____________________________________________________________________________
TG4ProcessControlMap::TG4ProcessControlMap() {
//
  if (fgInstance) {
    TG4Globals::Exception(
      "TG4ProcessControlMap: attempt to create two instances of singleton.");
  }
      
  fgInstance = this;  
}

//_____________________________________________________________________________
TG4ProcessControlMap::TG4ProcessControlMap(const TG4ProcessControlMap& right) {
//
  TG4Globals::Exception(    
    "Attempt to copy TG4ProcessControlMap singleton.");
}  

//_____________________________________________________________________________
TG4ProcessControlMap::~TG4ProcessControlMap() {
//
}

// operators

//_____________________________________________________________________________
TG4ProcessControlMap& TG4ProcessControlMap::operator=(const TG4ProcessControlMap& right)
{
  // check assignement to self
  if (this == &right) return *this;

  TG4Globals::Exception(
    "Attempt to assign TG4ProcessControlMap singleton.");
    
  return *this;  
}    
          
// private methods

//_____________________________________________________________________________
G4bool TG4ProcessControlMap::IsDefined(const G4String& processName)
{
// Returns true if the first is already in the map.
// ---

  if (fMap.find(processName) == fMap.end()) 
    return false;
  else                 
    return true;
}

// public methods

//_____________________________________________________________________________
G4bool TG4ProcessControlMap::Add(G4VProcess* process, TG4G3Control control)
{  
// Adds the pair to the map.
// ---

  if (!process) return false;

  return Add(process->GetProcessName(), control); 
}

//_____________________________________________________________________________
G4bool TG4ProcessControlMap::Add(G4String processName, TG4G3Control control)
{  
// Adds the pair to the map.
// ---

  if (!IsDefined(processName)) {
    // insert into map 
    // only in case it is not yet here
    fMap[processName] = control;
    return true;
  }
  return false;  
}

//_____________________________________________________________________________
void TG4ProcessControlMap::PrintAll() const
{
// Dumps all map.
// ---

  if (fMap.size()) {
    G4cout << "Dump of TG4ProcessControlMap - " << fMap.size() << " entries:" << G4endl;
    G4int counter = 0;
    for (MapConstIterator i=fMap.begin(); i != fMap.end(); i++) {
      G4String processName = (*i).first;
      TG4G3Control control = (*i).second;
      G4cout << "Map element " << G4std::setw(3) << counter++ << "   " 
             << processName << "   " 
	     << TG4G3ControlVector::GetControlName(control)
	     << G4endl;
    }
  }
}

//_____________________________________________________________________________
void TG4ProcessControlMap::Clear() 
{
// Clears the map.
// ---

  fMap.clear();
}  

//_____________________________________________________________________________
TG4G3Control 
TG4ProcessControlMap::GetControl(const G4VProcess* process)
{
// Returns the G3 process control for the process with a given name.
// ---

  if (!process) return kNoG3Controls;

  return GetControl(process->GetProcessName());
}

//_____________________________________________________________________________
TG4G3Control 
TG4ProcessControlMap::GetControl(const G4String& processName)
{
// Returns the G3 process control for the process with a given name.
// ---

  MapIterator i = fMap.find(processName);
  if (i == fMap.end()) 
    return kNoG3Controls;
  else                 
    return (*i).second;
}

//_____________________________________________________________________________
const G4String& 
TG4ProcessControlMap::GetControlName(const G4VProcess* process)
{
// Returns the G3 process control name for the process with a given name.
// ---

  if (!process) 
    return TG4G3ControlVector::GetControlName(kNoG3Controls);

  return GetControlName(process->GetProcessName());
}
	    
//_____________________________________________________________________________
const G4String& 
TG4ProcessControlMap::GetControlName(const G4String& processName)
{
// Returns the G3 process control name for the process with a given name.
// ---

  return TG4G3ControlVector::GetControlName(GetControl(processName));
}
	    
