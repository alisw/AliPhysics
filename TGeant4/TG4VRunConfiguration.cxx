// $Id$
// Category: run
//
// See the class description in the header file.

#include "TG4VRunConfiguration.h"
#include "TG4Globals.h"

#include <G4VUserDetectorConstruction.hh>
#include <G4VModularPhysicsList.hh>
#include <G4VUserPrimaryGeneratorAction.hh>
#include <G4UserRunAction.hh>
#include <G4UserEventAction.hh>
#include <G4UserTrackingAction.hh>
#include <G4UserSteppingAction.hh>
#include <G4UserStackingAction.hh>
#include <G4RunManager.hh>

TG4VRunConfiguration::TG4VRunConfiguration()
  : fDetectorConstruction(0),
    fPhysicsList(0),
    fPrimaryGenerator(0),
    fRunAction(0),
    fEventAction(0),
    fTrackingAction(0),
    fSteppingAction(0),
    fStackingAction(0)
{
//
}

TG4VRunConfiguration::TG4VRunConfiguration(const TG4VRunConfiguration& right)
{
//
  TG4Globals::Exception("TG4VRunConfiguration is protected from copying.");
}

TG4VRunConfiguration::~TG4VRunConfiguration(){
//
}

// operators

TG4VRunConfiguration& TG4VRunConfiguration::operator=(
                                const TG4VRunConfiguration& right)
{
  // check assignement to self
  if (this == &right) return *this;
  
  TG4Globals::Exception("TG4VRunConfiguration is protected from assigning.");

  return *this;  
}    
          
// public methods

void TG4VRunConfiguration::ConfigureRunManager(G4RunManager* runManager)
{
// Sets the user action classes to G4RunManager.
// --- 

  //if (!fDetectorConstruction || !fPhysicsList || !fPrimaryGenerator)
  //TG4Globals::Exception("Mandatory user classes are missing.");    
  
  // set mandatory classes
  runManager->SetUserInitialization(fDetectorConstruction);
  runManager->SetUserInitialization(fPhysicsList);
  runManager->SetUserAction(fPrimaryGenerator);      

  // user other action classes 
  if (fRunAction)      runManager->SetUserAction(fRunAction);  
  if (fEventAction)    runManager->SetUserAction(fEventAction);   
  if (fTrackingAction) runManager->SetUserAction(fTrackingAction);   
  if (fSteppingAction) runManager->SetUserAction(fSteppingAction);
  if (fStackingAction) runManager->SetUserAction(fStackingAction);
}

G4VModularPhysicsList* TG4VRunConfiguration::GetPhysicsList() const
{
// Returns the modular physics list.
// ---
  
  return fPhysicsList;
}

