// $Id$ //
// Category: geometry
//
// See the class description in the header file.

#include "TG4VSensitiveDetector.h"

#include "TG4StepManager.h"

G4int TG4VSensitiveDetector::fgSDCounter = 0;

TG4VSensitiveDetector::TG4VSensitiveDetector(G4String sdName)
  : G4VSensitiveDetector(sdName)
{
//
  fID = fgSDCounter++;
}

TG4VSensitiveDetector::TG4VSensitiveDetector(G4String sdName, G4int id)
  : G4VSensitiveDetector(sdName),
    fID(id)
{
//
  fgSDCounter++;
}

TG4VSensitiveDetector::TG4VSensitiveDetector(
                                    const TG4VSensitiveDetector& right)
  : G4VSensitiveDetector(right)
{  				    
//
  fID = right.fID;

  fgSDCounter++;;
}

TG4VSensitiveDetector::TG4VSensitiveDetector()
  : G4VSensitiveDetector("") 
{
//
}

TG4VSensitiveDetector::~TG4VSensitiveDetector() {
//
}

// operators

TG4VSensitiveDetector& TG4VSensitiveDetector::operator=(
                                    const TG4VSensitiveDetector &right)
{
  // check assignement to self
  if (this == &right) return *this;

  // base class assignement
  TG4VSensitiveDetector::operator=(right);
  
  fID = right.fID;
  
  return *this;
}
