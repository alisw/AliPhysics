// $Id$
// Category: visualization
//
// Author: I. Hrivnacova
//
// Class AliColourStore
// --------------------
// See the class description in the header file.

#include "AliColourStore.h"
#include "AliGlobals.h"

// static data members

//_____________________________________________________________________________
AliColourStore* AliColourStore::fgInstance = 0;

//_____________________________________________________________________________
AliColourStore::AliColourStore() 
{
  // fill predefined colours 
  G4int id0 = 1000;
     
  fColours.push_back(TColor(id0 + 0,  0.99, 0.99, 0.99, "White"));   
  fColours.push_back(TColor(id0 + 1,  0.00, 0.00, 0.00, "Black"));  
  fColours.push_back(TColor(id0 + 2,  0.99, 0.00, 0.00, "Red"));  
  fColours.push_back(TColor(id0 + 3,  0.99, 0.50, 0.00, "Red2")); 
  fColours.push_back(TColor(id0 + 4,  0.00, 0.00, 0.99, "Green"));   
  fColours.push_back(TColor(id0 + 5,  0.00, 0.50, 0.99, "Green2"));   
  fColours.push_back(TColor(id0 + 6,  0.50, 0.00, 0.99, "Green3"));
  fColours.push_back(TColor(id0 + 7,  0.00, 0.99, 0.00, "Blue"));
  fColours.push_back(TColor(id0 + 8,  0.00, 0.99, 0.99, "Blue2"));
  fColours.push_back(TColor(id0 + 9,  0.00, 0.99, 0.50, "Blue3"));
  fColours.push_back(TColor(id0 + 10, 0.99, 0.00, 0.99, "Yellow"));    
  fColours.push_back(TColor(id0 + 11, 0.99, 0.99, 0.00, "Magenta"));    
  fColours.push_back(TColor(id0 + 12, 0.50, 0.99, 0.00, "Magenta2"));   
  fColours.push_back(TColor(id0 + 13, 0.99, 0.00, 0.50, "Brown"));
  fColours.push_back(TColor(id0 + 14, 0.30, 0.30, 0.30, "Gray"));	
}

//_____________________________________________________________________________
AliColourStore::AliColourStore(const AliColourStore& right) {
// 
  AliGlobals::Exception(
    "Attempt to copy AliColourStore singleton.");
}

//_____________________________________________________________________________
AliColourStore::~AliColourStore() {
//
}

// operators

//_____________________________________________________________________________
AliColourStore& AliColourStore::operator=(const AliColourStore& right)
{
  // check assignement to self
  if (this == &right) return *this;

  AliGlobals::Exception(
    "Attempt to assign AliColourStore singleton.");
    
  return *this;  
}    

// static methods
  
//_____________________________________________________________________________
AliColourStore* AliColourStore::Instance() 
{
// Returns the singleton instance.
// Creates the instance if it does not exist.
// ---

  if (fgInstance == 0 )
    fgInstance = new AliColourStore();
  
  return fgInstance;
}

// public methods

//_____________________________________________________________________________
G4Colour AliColourStore::GetColour(const G4String& name) const
{
// Retrieves the colour by name.
// ---

  ColourConstIterator it;  
  for (it = fColours.begin(); it != fColours.end(); it++) 
    if (name == (*it).GetName()) return GetColour(*it);
  
  G4String text = "Colour " + name + " is not defined.";
  AliGlobals::Exception(text);
  return 0;
}
    
//_____________________________________________________________________________
G4Colour AliColourStore::GetColour(const TColor& color) const
{
// Converts TColor to G4Colour.
// ---

  return G4Colour(color.GetRed(), color.GetBlue(), color.GetGreen());
}
    
//_____________________________________________________________________________
G4String AliColourStore::GetColoursList() const
{
// Returns the list of all defined colours names.
// ---

  G4String list = "";
  ColourConstIterator it;

  for (it = fColours.begin(); it != fColours.end(); it++) {
    list += (*it).GetName();
    list += " ";
  }
  
  return list;
} 
       
//_____________________________________________________________________________
G4String AliColourStore::GetColoursListWithCommas() const
{
// Returns the list of all defined colours names
// with commas.
// ---

  G4String list = "";
  G4int i = 0;
  ColourConstIterator it;

  for (it = fColours.begin(); it != fColours.end(); it++) {
    list += (*it).GetName();
    if (i++ < fColours.size()-1) list += ", ";
  }
  
  return list;
} 
