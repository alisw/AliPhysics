// $Id$ //
// Category: digits+hits
//
// Author: I. Hrivnacova
//
// Class AliLegoSensitiveDetector
// ------------------------------
// See the class description in the header file.

#include "AliLegoSensitiveDetector.h"
#include "AliLego.h"

//_____________________________________________________________________________
AliLegoSensitiveDetector::AliLegoSensitiveDetector(
                             G4String name, AliLego* lego, 
			     G4VSensitiveDetector* standardSD)
  : TG4VSensitiveDetector(name),
    fLego(lego),
    fStandardSD(standardSD)
{
//
}

//_____________________________________________________________________________
AliLegoSensitiveDetector::AliLegoSensitiveDetector(
                                 const AliLegoSensitiveDetector& right)
  : TG4VSensitiveDetector(right)			     
{
//
  // copy stuff
  *this = right;
}

//_____________________________________________________________________________
AliLegoSensitiveDetector::AliLegoSensitiveDetector() {
//
}

//_____________________________________________________________________________
AliLegoSensitiveDetector::~AliLegoSensitiveDetector() {
//
}

// operators

//_____________________________________________________________________________
AliLegoSensitiveDetector& 
AliLegoSensitiveDetector::operator=(const AliLegoSensitiveDetector &right)
{
  // check assignement to self
  if (this == &right) return *this;
  
  // base class assignement
  TG4VSensitiveDetector::operator=(right);

  fLego = right.fLego;
  fStandardSD = right.fStandardSD;
  
  return *this;
}

// public methods

//_____________________________________________________________________________
void AliLegoSensitiveDetector::UserProcessHits(const G4Track* track,
                                               const G4Step* step)
{
// Calls StepManager of associated lego.
// ---

  // let lego process step
  fLego->StepManager();
}

