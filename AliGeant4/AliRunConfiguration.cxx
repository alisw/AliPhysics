// $Id$
// Category: run
//
// See the class description in the header file.

#include "AliRunConfiguration.h"
#include "AliRunMessenger.h"

#include "AliDetConstruction.h"
#include "AliPrimaryGeneratorAction.h"
#include "AliRunAction.h"
#include "AliEventAction.h"
#include "AliTrackingAction.h"
#include "AliStackingAction.h"
#include "AliSteppingAction.h"

#ifdef ALICE_EMPTY_PHYSICS_LIST
#include "AliEmptyPhysicsList.h"
#else
#include "TG4PhysicsList.h"
#endif

AliRunConfiguration::AliRunConfiguration(){
//
  fRunMessenger = new AliRunMessenger();
 
  CreateUserConfiguration();
}

AliRunConfiguration::AliRunConfiguration(const AliRunConfiguration& right)
  : TG4VRunConfiguration(right)
{
  // TG4VRunConfiguration is protected from copying
}

AliRunConfiguration::~AliRunConfiguration() {
//
  delete fRunMessenger;

  // all user action data members are deleted 
  // in G4RunManager::~G4RunManager()
}

// operators

AliRunConfiguration& 
AliRunConfiguration::operator=(const AliRunConfiguration& right)
{
  // check assignement to self
  if (this == &right) return *this;
  
  // base class assignement
  // TG4VRunConfiguration is protected from assigning
  TG4VRunConfiguration::operator=(right);

  return *this;  
}    
          
// protected methods

void AliRunConfiguration::CreateUserConfiguration()
{
// Creates the mandatory Geant4 classes and 
// the other user action classes. 
// ---

  // create mandatory Geant4 classes
  fDetectorConstruction = new AliDetConstruction();
#ifndef ALICE_EMPTY_PHYSICS_LIST
  fPhysicsList = new TG4PhysicsList();
#else
  fPhysicsList = new AliEmptyPhysicsList();
#endif
  fPrimaryGenerator = new AliPrimaryGeneratorAction();

  // create the other user action classes
  fRunAction = new AliRunAction();
  fEventAction = new AliEventAction();
  fTrackingAction = new AliTrackingAction();
  fSteppingAction = new AliSteppingAction();
#ifdef ALICE_STACKING
  fStackingAction = new AliStackingAction();
#endif
}
