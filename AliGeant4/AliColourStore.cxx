// $Id$
// Category: visualization
//
// See the class description in the header file.

#include "AliColourStore.h"
#include "AliColour.h"
#include "AliGlobals.h"

#include <G4Element.hh>

// static data members

AliColourStore* AliColourStore::fgInstance = 0;

// lifecycle

AliColourStore::AliColourStore() {
//
  fColours.insert(AliColour("White",     1.0, 1.0, 1.0));    
  fColours.insert(AliColour("Black",     0.0, 0.0, 0.0));     
  fColours.insert(AliColour("Red",       1.0, 0.0, 0.0));   
  fColours.insert(AliColour("RoseDark",  1.0, 0.0, 0.5));  
  fColours.insert(AliColour("Green",     0.0, 1.0, 0.0));     
  fColours.insert(AliColour("Green2",    0.0, 1.0, 0.5));     
  fColours.insert(AliColour("GreenClair",0.5, 1.0, 0.0));
  fColours.insert(AliColour("Yellow",    1.0, 1.0, 0.0));     
  fColours.insert(AliColour("BlueDark",  0.0, 0.0, 1.0)); 
  fColours.insert(AliColour("BlueClair", 0.0, 1.0, 1.0)); 
  fColours.insert(AliColour("BlueClair2",0.0, 0.5, 1.0));
  fColours.insert(AliColour("Magenta",   1.0, 0.0, 1.0));    
  fColours.insert(AliColour("Magenta2",  0.5, 0.0, 1.0));   
  fColours.insert(AliColour("BrownClair",1.0, 0.5, 0.0));
  fColours.insert(AliColour("Gray",      0.3, 0.3, 0.3));    
  fColours.insert(AliColour("GrayClair", 0.6, 0.6, 0.6));
}

AliColourStore::AliColourStore(const AliColourStore& right) {
// 
  AliGlobals::Exception(
    "Attempt to copy AliColourStore singleton.");
}

AliColourStore::~AliColourStore() {
//
}

// operators

AliColourStore& AliColourStore::operator=(const AliColourStore& right)
{
  // check assignement to self
  if (this == &right) return *this;

  AliGlobals::Exception(
    "Attempt to assign AliColourStore singleton.");
    
  return *this;  
}    

// static methods
  
AliColourStore* AliColourStore::Instance() 
{
// Returns the singleton instance.
// Creates the instance if it does not exist.
// ---

  if (fgInstance == 0 ) {
    fgInstance = new AliColourStore();
  }
  
  return fgInstance;
}

// public methods

G4Colour AliColourStore::GetColour(G4String name) const
{
// Retrieves the colour by name.
// ---

  G4int nofCol = fColours.entries();
  for (G4int i=0; i<nofCol; i++)
  {
    AliColour alColour = fColours[i];
    if (name == alColour.GetName()) 
    { return alColour.GetColour(); }
  }
  
  G4String text = "Colour " + name + " is not defined.";
  AliGlobals::Exception(text);
  return 0;
}
    
G4String AliColourStore::GetColoursList() const
{
// Returns the list of all defined colours names.
// ---

  G4String list = "";
  G4int nofCol = fColours.entries();
  for (G4int i=0; i<nofCol; i++)
  {
    list += fColours[i].GetName();
    list += " ";
  };
  return list;
} 
       
G4String AliColourStore::GetColoursListWithCommas() const
{
// Returns the list of all defined colours names
// with commas.
// ---

  G4String list = "";
  G4int nofCol = fColours.entries();
  for (G4int i=0; i<nofCol; i++)
  {
    list += fColours[i].GetName();
    if (i < nofCol-1) list += ", ";
  };
  return list;
} 
