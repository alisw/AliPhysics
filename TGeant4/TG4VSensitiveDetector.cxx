// $Id$ //
// Category: digits+hits
//
// See the class description in the header file.

#include "TG4VSensitiveDetector.h"
#include "TG4StepManager.h"

G4int TG4VSensitiveDetector::fgSDCounter = 0;

//_____________________________________________________________________________
TG4VSensitiveDetector::TG4VSensitiveDetector(G4String sdName)
  : G4VSensitiveDetector(sdName),
    fStepManager(TG4StepManager::Instance())
{
//
  fID = ++fgSDCounter;
}

//_____________________________________________________________________________
TG4VSensitiveDetector::TG4VSensitiveDetector(G4String sdName, G4int id)
  : G4VSensitiveDetector(sdName),
    fID(id),
    fStepManager(TG4StepManager::Instance())

{
//
  ++fgSDCounter;
}

//_____________________________________________________________________________
TG4VSensitiveDetector::TG4VSensitiveDetector(
                                    const TG4VSensitiveDetector& right)
  : G4VSensitiveDetector(right)
{  				    
//
  // copy stuff
  *this = right;

  ++fgSDCounter;;
}

//_____________________________________________________________________________
TG4VSensitiveDetector::TG4VSensitiveDetector()
  : G4VSensitiveDetector("") 
{
//
}

//_____________________________________________________________________________
TG4VSensitiveDetector::~TG4VSensitiveDetector() {
//
}

// operators

//_____________________________________________________________________________
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

//_____________________________________________________________________________
G4bool TG4VSensitiveDetector::ProcessHits(G4Step* step, G4TouchableHistory*)
{
// Calls StepManager of associated AliModule.
// ---

  // let user sensitive detector process normal step
  fStepManager->SetStep(step, kNormalStep);
  UserProcessHits(step->GetTrack(), step);

  return true;
}

//_____________________________________________________________________________
G4bool TG4VSensitiveDetector::ProcessHitsOnBoundary(G4Step* step)
{
// Calls StepManager of associated AliModule
// when crossing a geometrical boundary.
// ---

  // let user sensitive detector process boundary step
  fStepManager->SetStep(step, kBoundary);
  UserProcessHits(step->GetTrack(), step);

  return true;
}

