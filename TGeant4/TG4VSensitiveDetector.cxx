// $Id$ //
// Category: geometry
//
// See the class description in the header file.

#include "TG4VSensitiveDetector.h"
#include "TG4StepManager.h"

G4int TG4VSensitiveDetector::fgSDCounter = 0;

TG4VSensitiveDetector::TG4VSensitiveDetector(G4String sdName)
  : G4VSensitiveDetector(sdName),
    fStepManager(TG4StepManager::Instance())
{
//
  fID = ++fgSDCounter;
}

TG4VSensitiveDetector::TG4VSensitiveDetector(G4String sdName, G4int id)
  : G4VSensitiveDetector(sdName),
    fID(id),
    fStepManager(TG4StepManager::Instance())

{
//
  ++fgSDCounter;
}

TG4VSensitiveDetector::TG4VSensitiveDetector(
                                    const TG4VSensitiveDetector& right)
  : G4VSensitiveDetector(right)
{  				    
//
  fID = right.fID;
  fStepManager = right.fStepManager;

  ++fgSDCounter;;
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
  fStepManager = right.fStepManager;
  
  return *this;
}

// public methods

G4bool TG4VSensitiveDetector::ProcessHits(G4Step* step, G4TouchableHistory*)
{
// Calls StepManager of associated AliModule.
// ---

  // let AliModule process step
  fStepManager->SetStep(step, kNormalStep);
  UserProcessHits(step->GetTrack(), step);

  return true;
}

G4bool TG4VSensitiveDetector::ProcessHitsOnBoundary(G4Step* step)
{
// Calls StepManager of associated AliModule
// when crossing a geometrical boundary.
// ---

  fStepManager->SetStep(step, kBoundary);
  UserProcessHits(step->GetTrack(), step);

  return true;
}

