// $Id$
// Category: run
//
// Author: I. Hrivnacova
//
// Class AliRunConfiguration
// -------------------------
// See the class description in the header file.

#include "AliRunConfiguration.h"
#include "AliRunMessenger.h"
#include "AliDetConstruction.h"
#include "TG4SDConstruction.h"
#include "AliPrimaryGeneratorAction.h"
#include "TG4RunAction.h"
#include "TG4EventAction.h"
#include "TG4TrackingAction.h"
#include "TG4SteppingAction.h"
#include "TG4SpecialStackingAction.h"
#include "AliFiles.h"

#include "TG4ModularPhysicsList.h"

ClassImp(AliRunConfiguration)

Bool_t  AliRunConfiguration::fgIsHoles = true;

//_____________________________________________________________________________
AliRunConfiguration::AliRunConfiguration()
  : TG4VRunConfiguration()
{
//
  CreateUserConfiguration();

  fRunMessenger = new AliRunMessenger();
  fFiles = new AliFiles(); 
}

//_____________________________________________________________________________
AliRunConfiguration::AliRunConfiguration(const AliRunConfiguration& right)
  : TG4VRunConfiguration(right)
{
  // TG4VRunConfiguration is protected from copying
}

//_____________________________________________________________________________
AliRunConfiguration::~AliRunConfiguration() {
//
  delete fRunMessenger;
  delete fFiles;

  // all user action data members are deleted 
  // in G4RunManager::~G4RunManager()
}

// operators

//_____________________________________________________________________________
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

//_____________________________________________________________________________
void AliRunConfiguration::CreateUserConfiguration()
{
// Creates the mandatory Geant4 classes and 
// the other user action classes. 
// ---

  G4cout << "AliRunConfiguration::CreateUserConfiguration()" << G4endl;

  // create mandatory Geant4 classes
  fDetectorConstruction = new AliDetConstruction();
  fSDConstruction = new TG4SDConstruction();
  fPhysicsList = new TG4ModularPhysicsList();
  fPrimaryGenerator = new AliPrimaryGeneratorAction();

  // create the other user action classes
  fRunAction  = new TG4RunAction();
  fEventAction  = new TG4EventAction();
  fTrackingAction = new TG4TrackingAction();
  fSteppingAction = new TG4SteppingAction();
  fStackingAction = new TG4SpecialStackingAction();
}

// public methods

//_____________________________________________________________________________
void AliRunConfiguration::SetConfigName(const char* name)
{
// Sets the configuration macro name 
// ---
  fFiles->SetMacroName(name);
}  

//_____________________________________________________________________________
void AliRunConfiguration::SetG3CallsName(const char* name)
{
// Sets the configuration macro name 
// ---
  fFiles->SetG3CallsName(name);
}  

