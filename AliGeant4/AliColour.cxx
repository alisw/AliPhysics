// $Id$
// Category: visualization
//
// Author: I. Hrivnacova
//
// Class AliColour
// ---------------
// See the class description in the header file.

#include "AliColour.h"
#include <G4Colour.hh>
  
//_____________________________________________________________________________
AliColour::AliColour()
  : fName (""),
    fRed(0.),
    fBlue(0.),
    fGreen(0.)
{
//
}

//_____________________________________________________________________________
AliColour::AliColour(G4String name, G4double red, G4double blue, G4double green)
  : fName(name),
    fRed(red),
    fBlue(blue),
    fGreen(green)
{
//
}

//_____________________________________________________________________________
AliColour::AliColour(const AliColour& right) {
//
  // copy stuff
  *this = right;
}

//_____________________________________________________________________________
AliColour::~AliColour() {
//
}

// operators

//_____________________________________________________________________________
AliColour& AliColour::operator=(const AliColour& right)
{    
  // check assignement to self
  if (this == &right) return *this;
  
  fName = right.fName;
  fRed = right.fRed;
  fBlue = right.fBlue;
  fGreen = right.fGreen;
  
  return *this;
}

//_____________________________________________________________________________
G4int AliColour::operator==(const AliColour& right) const
{    
//
  G4int returnValue = 0;
  if ( fName == right.fName && 
       fRed == right.fRed && fBlue == right.fBlue && fGreen == right.fGreen) 
  { returnValue = 1; };

  return returnValue;  
}

//_____________________________________________________________________________
G4int AliColour::operator!=(const AliColour& right) const
{
//    
  G4int returnValue = 1;
  if (*this == right) returnValue = 0; 
  
  return returnValue;
}
