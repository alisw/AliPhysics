// $Id$
// Category: run
//
// See the class description in the header file.

#include <G4Timer.hh>
   // in order to avoid the odd dependency for the
   // times system function this include must be the first

#include "AliRunAction.h"
#include "AliRunActionMessenger.h"
#include "AliSDManager.h"
#include "AliGlobals.h"
#include "AliRun.h"
#include "AliLego.h"

#include "TG4GeometryManager.h"

#include <G4Run.hh>
#include <G4UImanager.hh>

AliRunAction::AliRunAction()
  : fRunID(-1),
    fVerboseLevel(0)
{
//
  fMessenger = new AliRunActionMessenger(this);
  fTimer = new G4Timer;
}

AliRunAction::AliRunAction(const AliRunAction& right) {
//
  AliGlobals::Exception("AliRunAction is protected from copying.");
}

AliRunAction::~AliRunAction() {
//
  delete fMessenger;
  delete fTimer;
}

// operators

AliRunAction& AliRunAction::operator=(const AliRunAction &right)
{
  // check assignement to self
  if (this == &right) return *this;
  
  AliGlobals::Exception("AliRunAction is protected from assigning.");

  return *this;
}

// public methods

void AliRunAction::BeginOfRunAction(const G4Run* run)
{
// Called by G4 kernel at the beginning of run.
// ---

  fRunID++;
  
  // aliroot
  // store runID in the event header
  gAlice->GetHeader()->SetRun(fRunID);

  // clear remaining G3 tables
  if (fRunID == 0)
    TG4GeometryManager::Instance()->ClearG3TablesFinal();

  // create lego sensitive detectors 
  // if lego is instantiated
  AliLego* lego = gAlice->Lego();
  if (lego) {
    AliSDManager::Instance()->SetLego(lego);
    G4UImanager::GetUIpointer()->ApplyCommand("/aliEvent/verbose 0");
    G4UImanager::GetUIpointer()->ApplyCommand("/aliGenerator/set AliGenerator");
  }  

  G4cout << "### Run " << run->GetRunID() << " start." << G4endl;
  fTimer->Start();
}

void AliRunAction::EndOfRunAction(const G4Run* run)
{
// Called by G4 kernel at the end of run.
// ---

  fTimer->Stop();

  // delete lego sensitive detectors 
  // if lego is instantiated
  AliLego* lego = gAlice->Lego();
  if (lego) {
    AliSDManager::Instance()->UnsetLego();
    G4UImanager::GetUIpointer()->ApplyCommand("/aliEvent/verbose 1");
  }  

  G4cout << "Time of this run:   " << *fTimer << G4endl;
  G4cout << "Number of events processed: " << run->GetNumberOfEvent() << G4endl;
}    
