// $Id$
// Category: run
//
// See the class description in the header file.

#include <G4Timer.hh>
   // in order to avoid the odd dependency for the
   // times system function this include must be the first

#include "AliRunAction.h"
#include "AliRunActionMessenger.h"
#include "AliSDConstruction.h"
#include "AliGlobals.h"
#include "AliRun.h"
#include "AliHeader.h"
#include "AliLego.h"

#include "TG4GeometryManager.h"
#include "TG4SDManager.h"
#include "TG4VSDConstruction.h"

#include <G4Run.hh>
#include <G4UImanager.hh>

//_____________________________________________________________________________
AliRunAction::AliRunAction()
  : fRunID(-1),
    fVerboseLevel(0)
{
//
  fMessenger = new AliRunActionMessenger(this);
  fTimer = new G4Timer;
}

//_____________________________________________________________________________
AliRunAction::AliRunAction(const AliRunAction& right) {
//
  AliGlobals::Exception("AliRunAction is protected from copying.");
}

//_____________________________________________________________________________
AliRunAction::~AliRunAction() {
//
  delete fMessenger;
  delete fTimer;
}

// operators

//_____________________________________________________________________________
AliRunAction& AliRunAction::operator=(const AliRunAction &right)
{
  // check assignement to self
  if (this == &right) return *this;
  
  AliGlobals::Exception("AliRunAction is protected from assigning.");

  return *this;
}

// private methods

//_____________________________________________________________________________
AliSDConstruction* AliRunAction::GetSDConstruction() const
{
// Gets sensitive detectors construction and checks type.
// ---

  TG4VSDConstruction* tg4SDConstruction 
     = TG4SDManager::Instance()->GetSDConstruction();

  AliSDConstruction* aliSDConstruction
     = dynamic_cast<AliSDConstruction*>(tg4SDConstruction);

  if (!aliSDConstruction) {
     G4String text = "AliRunAction::GetSDConstruction:\n";
     text = text + "    Unknown type.";
     AliGlobals::Exception(text);
     return 0;
  }
  
  return aliSDConstruction;
}

// public methods

//_____________________________________________________________________________
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
    GetSDConstruction()->SetLego(lego);
    G4UImanager::GetUIpointer()->ApplyCommand("/aliEvent/verbose 0");
    G4UImanager::GetUIpointer()->ApplyCommand("/aliGenerator/set AliGenerator");
  }  

  G4cout << "### Run " << run->GetRunID() << " start." << G4endl;
  fTimer->Start();
}

//_____________________________________________________________________________
void AliRunAction::EndOfRunAction(const G4Run* run)
{
// Called by G4 kernel at the end of run.
// ---

  fTimer->Stop();

  // delete lego sensitive detectors 
  // if lego is instantiated
  AliLego* lego = gAlice->Lego();
  if (lego) {
    GetSDConstruction()->UnsetLego();
    G4UImanager::GetUIpointer()->ApplyCommand("/aliEvent/verbose 1");
  }  

  G4cout << "Time of this run:   " << *fTimer << G4endl;
  G4cout << "Number of events processed: " << run->GetNumberOfEvent() << G4endl;
}    
