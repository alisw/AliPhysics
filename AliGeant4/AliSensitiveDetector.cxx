// $Id$ //
// Category: geometry
//
// See the class description in the header file.

#include "AliSensitiveDetector.h"
#include "AliModule.h" 
#include "AliRun.h"

#include "TG3Units.h"

AliSensitiveDetector::AliSensitiveDetector(G4String sdName, AliModule* module)
  : TG4VSensitiveDetector(sdName),
    fModule(module)
{
//
}

AliSensitiveDetector::AliSensitiveDetector(G4String sdName, AliModule* module, 
                                           G4int id)
  : TG4VSensitiveDetector(sdName, id),
    fModule(module)
{
//
}

AliSensitiveDetector::AliSensitiveDetector(const AliSensitiveDetector& right)
  : TG4VSensitiveDetector(right)
{
//
  fModule = right.fModule;
}  
  
AliSensitiveDetector::AliSensitiveDetector(){
//
}

AliSensitiveDetector::~AliSensitiveDetector() {
//
}

// operators

AliSensitiveDetector& 
AliSensitiveDetector::operator=(const AliSensitiveDetector& right)
{
  // check assignement to self
  if (this == &right) return *this;

  // base class assignement
  TG4VSensitiveDetector::operator=(right);

  fModule = right.fModule;

  return *this;  
}    
          
// public methods

void AliSensitiveDetector::UserProcessHits(const G4Track* track, 
                                           const G4Step* step)
{
// Calls StepManager of associated AliModule.
// ---

  // add energy deposit of the current step
  // directly to AliRun
  if (step) 
    gAlice->AddEnergyDeposit(
      fID, step->GetTotalEnergyDeposit()/TG3Units::Energy());

  // let AliModule process the step
  fModule->StepManager();
}

