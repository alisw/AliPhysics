// $Id$
// Category: global
//
// See the class description in the header file.

#include "TG4NameMap.h"
#include "TG4Globals.h"

#include "g4std/iomanip"
#include "globals.hh"

G4String TG4NameMap::fgUndefined = "Undefined";

typedef G4std::map<G4String, G4String, G4std::less<G4String> >
  :: iterator MapIterator;

TG4NameMap::TG4NameMap() 
  : fSecond(fgUndefined) {
//
}

TG4NameMap::TG4NameMap(const TG4NameMap& right) {
//
  TG4Globals::Exception("TG4NameMap is protected from copying.");
}  

TG4NameMap::~TG4NameMap() {
//
}

// operators

TG4NameMap& TG4NameMap::operator=(const TG4NameMap& right)
{
  // check assignement to self
  if (this == &right) return *this;

  TG4Globals::Exception("TG4NameMap is protected from assigning.");
    
  return *this;  
}    
          
// public methods

G4bool TG4NameMap::Add(const G4String& first, const G4String& second)
{  
// Adds names pair to the map.
// fSecond is not used in this add method.
// ---

  if (GetSecond(first) == fgUndefined) {
    // insert into map 
    // only in case it is not yet here
    fMap[first] = second;
    return true;
  }
  return false;
}


G4bool TG4NameMap::AddName(const G4String& name)
{  
// Adds name to the map.
// ---

  if (GetSecond(name) == fgUndefined) {
    // insert into map 
    // only in case it is not yet here
    fMap[name] = fSecond;
    return true;
  }
  return false;
}


const G4String& TG4NameMap::GetSecond(const G4String& name)
{
// Gets second name associated with given name.
// ---

  MapIterator i = fMap.find(name);
  if (i == fMap.end()) 
    return fgUndefined;
  else                 
    return (*i).second;
}


void TG4NameMap::PrintAll()
{
// Dumps all map.
// ---

  if (fMap.size()) {
    G4cout << "Dump of TG4NameMap - " << fMap.size() << " entries:" << G4endl;
    G4int counter = 0;
    for (MapIterator i=fMap.begin(); i != fMap.end(); i++) {
      const G4String& first  = (*i).first;
      const G4String& second = (*i).second;
      G4cout << "Map element " << G4std::setw(3) << counter++ << "   " 
             << first << "   " << second << G4endl;
    }
  }
}


void TG4NameMap::Clear() 
{
// Clears the map.
// ---

  if (fMap.size()>0){
    for (MapIterator i=fMap.begin(); i != fMap.end(); i++) {
      delete (*i).second;
    }
    fMap.clear();
  }
  fSecond = "Undefined";
}  
