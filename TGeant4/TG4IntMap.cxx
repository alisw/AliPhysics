// $Id$
// Category: global
//
// See the class description in the header file.

#include "TG4IntMap.h"
#include "TG4Globals.h"

#include "g4std/iomanip"
#include "globals.hh"

typedef G4std::map<G4String, G4int, G4std::less<G4String> >
  :: iterator IntMapIterator;

TG4IntMap::TG4IntMap(){
//
}

TG4IntMap::TG4IntMap(const TG4IntMap& right) {
//
  TG4Globals::Exception("TG4IntMap is protected from copying.");
}  

TG4IntMap::~TG4IntMap() {
//
}

// operators

TG4IntMap& TG4IntMap::operator=(const TG4IntMap& right)
{
  // check assignement to self
  if (this == &right) return *this;

  TG4Globals::Exception("TG4IntMap is protected from assigning.");
    
  return *this;  
}    
          
// private methods

G4bool TG4IntMap::IsDefined(const G4String& first)
{
// Returns true if the first is already in the map.
// ---

  IntMapIterator i = fMap.find(first);
  if (i == fMap.end()) 
    return false;
  else                 
    return true;
}

// public methods

G4bool TG4IntMap::Add(const G4String& first, G4int second)
{  
// Adds pair (name, int number) to the map.
// ---

  if (!IsDefined(first)) {
    // insert into map 
    // only in case it is not yet here
    fMap[first] = second;
    return true;
  }
  return false;
}

G4int TG4IntMap::GetSecond(const G4String& name)
{
// Gets second name associated with given name.
// ---

  IntMapIterator i = fMap.find(name);
  if (i == fMap.end()) {
    G4String text = "G4IntMap::GetSecond: ";
    text = text + name + " is not defined.";
    TG4Globals::Warning(text);
    return 0;
  }  
  else {                
    return (*i).second;
  }  
}

void TG4IntMap::PrintAll()
{
// Dumps all map.
// ---

  if (fMap.size()) {
    G4cout << "Dump of TG4IntMap - " << fMap.size() << " entries:" << G4endl;
    G4int counter = 0;
    for (IntMapIterator i=fMap.begin(); i != fMap.end(); i++) {
      const G4String& first  = (*i).first;
      G4int second = (*i).second;
      G4cout << "Map element " << G4std::setw(3) << counter++ << "   " 
             << first << "   " << second << G4endl;
    }
  }
}


void TG4IntMap::Clear() 
{
// Clears the map.
// ---

  if (fMap.size()>0) fMap.clear();
}  
