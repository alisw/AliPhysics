// $Id$ //
// Category: geometry
//
// See the class description in the header file.

#include "AliSensitiveDetector.h"
#include "AliModule.h" 
#include "AliRun.h"

#include "TG4StepManager.h"

AliSensitiveDetector::AliSensitiveDetector(G4String sdName, AliModule* module)
  : TG4VSensitiveDetector(sdName),
    fModule(module),
    fStepManager(TG4StepManager::Instance())
{
//
}

AliSensitiveDetector::AliSensitiveDetector(G4String sdName, AliModule* module, 
                                           G4int id)
  : TG4VSensitiveDetector(sdName, id),
    fModule(module),
    fStepManager(TG4StepManager::Instance())
{
//
}

AliSensitiveDetector::AliSensitiveDetector(const AliSensitiveDetector& right)
  : TG4VSensitiveDetector(right)
{
//
  fModule = right.fModule;
  fStepManager = right.fStepManager;
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
  fStepManager = right.fStepManager;

  return *this;  
}    
          
// public methods

void AliSensitiveDetector::Initialize(G4HCofThisEvent* hc)
{
// This method is called at the beginning of event action
// before user defined BeginOfEventAction() method.
}

G4bool AliSensitiveDetector::ProcessHits(G4Step* step, G4TouchableHistory*)
{
// Calls StepManager of associated AliModules.
// ---

  // add energy deposit of the current step
  // directly to AliRun
  G4int copy;
  gAlice->AddEnergyDeposit(fID, step->GetTotalEnergyDeposit());

  // parent ID -> shunt
  G4int parentID
    = step->GetTrack()->GetParentID();
  Int_t shunt = 0;
  if (parentID==0) shunt = 1;
  fModule->SetIshunt(shunt);

  // let AliModule process step
  fStepManager->SetStep(step);
  fModule->StepManager();

  return true;
}

void AliSensitiveDetector::EndOfEvent(G4HCofThisEvent* hce){
//
}

//void AliSensitiveDetector::clear()
//{} 

void AliSensitiveDetector::PrintAll() {
//
} 

void AliSensitiveDetector::DrawAll() {
//
} 
