// $Id$ //
// Category: geometry
//
// See the class description in the header file.

#include "AliLegoSensitiveDetector.h"
#include "AliLego.h"

#include "TG4StepManager.h"

AliLegoSensitiveDetector::AliLegoSensitiveDetector(
                             G4String name, AliLego* lego, 
			     G4VSensitiveDetector* standardSD)
  : TG4VSensitiveDetector(name),
    fLego(lego),
    fStandardSD(standardSD),
    fStepManager(TG4StepManager::Instance())
{
//
}

AliLegoSensitiveDetector::AliLegoSensitiveDetector(
                                 const AliLegoSensitiveDetector& right)
  : TG4VSensitiveDetector(right)			     
{
//
  fLego = right.fLego;
  fStandardSD = right.fStandardSD;
  fStepManager = right.fStepManager;
}

AliLegoSensitiveDetector::AliLegoSensitiveDetector() {
//
}

AliLegoSensitiveDetector::~AliLegoSensitiveDetector() {
//
}

// operators

AliLegoSensitiveDetector& 
AliLegoSensitiveDetector::operator=(const AliLegoSensitiveDetector &right)
{
  // check assignement to self
  if (this == &right) return *this;
  
  // base class assignement
  TG4VSensitiveDetector::operator=(right);

  fLego = right.fLego;
  fStandardSD = right.fStandardSD;
  fStepManager = right.fStepManager;
  
  return *this;
}

// public methods

void AliLegoSensitiveDetector::Initialize(G4HCofThisEvent* hc)
{
// This method is called at the beginning of event action
// before user defined BeginOfEventAction() method.
}

G4bool AliLegoSensitiveDetector::ProcessHits(G4Step* step, G4TouchableHistory*)
{
// Calls StepManager of associated lego.
// ---

  // let lego process step
  fStepManager->SetStep(step, kPostStepPoint);
  fLego->StepManager();

  return true;
}

void AliLegoSensitiveDetector::EndOfEvent(G4HCofThisEvent* hce){
//
}

//void AliLegoSensitiveDetector::clear()
//{} 

void AliLegoSensitiveDetector::PrintAll() {
//
} 

void AliLegoSensitiveDetector::DrawAll() {
//
} 
