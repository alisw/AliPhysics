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
#include "AliFiles.h"

#include "TG4ModularPhysicsList.h"

AliRunConfiguration::AliRunConfiguration(){
//
  fRunMessenger = new AliRunMessenger();
  fFiles = new AliFiles();
 
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
  delete fFiles;

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
  fPhysicsList = new TG4ModularPhysicsList();
  fPrimaryGenerator = new AliPrimaryGeneratorAction();

  // create the other user action classes
  fRunAction = new AliRunAction();
  fEventAction = new AliEventAction();
  fTrackingAction = new AliTrackingAction();
  fSteppingAction = new AliSteppingAction();
  fStackingAction = new AliStackingAction();
}

// public methods

void AliRunConfiguration::SetConfigName(const char* name)
{
// Sets the configuration macro name 
// ---
  fFiles->SetMacroName(name);
}  

void AliRunConfiguration::SetG3CallsName(const char* name)
{
// Sets the configuration macro name 
// ---
  fFiles->SetG3CallsName(name);
}  

