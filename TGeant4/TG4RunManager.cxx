// $Id$
// Category: run
//
// See the class description in the header file.

#include "TG4RunManager.h"
#include "TG4RunMessenger.h"
#include "TG4VRunConfiguration.h"
#include "TG4Globals.h"
#include "TG4GeometryManager.h"
#include "TG4PhysicsManager.h"
#include "TG4G3PhysicsManager.h"

#include <G4RunManager.hh>
#include <G4UIsession.hh>
#include <G4UImanager.hh>
#include <G4UIterminal.hh>
#include <G4UIXm.hh>
#ifdef G4UI_USE_WO
#include <G4UIWo.hh>
#endif
#ifdef G4UI_USE_GAG
#include <G4UIGAG.hh>
#endif

#include <TRint.h>
#include <TROOT.h> 
#include <TCint.h> 

TG4RunManager* TG4RunManager::fgInstance = 0;

TG4RunManager::TG4RunManager(TG4VRunConfiguration* runConfiguration, 
                             int argc, char** argv)		  
  : fRunConfiguration(runConfiguration),
    fGeantUISession(0),
    fRootUISession(0),
    fRootUIOwner(false),
    fARGC(argc),
    fARGV(argv)
  
{
// 
  if (fgInstance) {
    TG4Globals::Exception(
      "TG4RunManager: attempt to create two instances of singleton.");
  };
      
  if (!fRunConfiguration) {
    TG4Globals::Exception(
      "TG4RunManager: attempt to create instance without runConfiguration.");
  };
      
  fgInstance = this;
  
  // create and configure geant4 run manager
  fRunManager =  new G4RunManager();
  fRunConfiguration->ConfigureRunManager(fRunManager);
  // add verbose level
  //G4cout << "G4RunManager has been created." << G4endl;

  // create geant4 UI
  CreateGeantUI();
      // must be created before TG4VisManager::Initialize()
      // (that is invoked in TGeant4 constructor)

  // create root UI
  CreateRootUI();

  // messenger
  fMessenger = new TG4RunMessenger(this);
}

TG4RunManager::TG4RunManager(TG4VRunConfiguration* runConfiguration)
  : fRunConfiguration(runConfiguration),
    fGeantUISession(0),
    fRootUISession(0),
    fRootUIOwner(false),
    fARGC(0),
    fARGV(0)
  
{
//
  if (fgInstance) {
    TG4Globals::Exception(
      "TG4RunManager: attempt to create two instances of singleton.");
  };
      
  if (!fRunConfiguration) {
    TG4Globals::Exception(
      "TG4RunManager: attempt to create instance without runConfiguration.");
  };
      
  fgInstance = this;
  
  // set primary UI
  fRootUISession = gROOT->GetApplication();
  if (fRootUISession) {
    fARGC = fRootUISession->Argc();
    fARGV = fRootUISession->Argv();
  }

  // create and configure geant4 run manager
  fRunManager =  new G4RunManager();
  fRunConfiguration->ConfigureRunManager(fRunManager);
  // add verbose level
  //G4cout << "G4RunManager has been created." << G4endl;

  // create geant4 UI
  CreateGeantUI();
      // must be created before TG4VisManager::Initialize()
      // (that is invoked in TGeant4 constructor)

  // create root UI
  CreateRootUI();

  // messenger
  fMessenger = new TG4RunMessenger(this);
}

TG4RunManager::TG4RunManager() {
//
}

TG4RunManager::TG4RunManager(const TG4RunManager& right) {
// 
  TG4Globals::Exception(
    "Attempt to copy TG4RunManager singleton.");
}

TG4RunManager::~TG4RunManager() {
//  
  delete fRunConfiguration;
  delete fRunManager;
  delete fGeantUISession;
  if (fRootUIOwner) delete fRootUISession;
  delete fMessenger;
}

// operators

TG4RunManager& TG4RunManager::operator=(const TG4RunManager& right)
{
  // check assignement to self
  if (this == &right) return *this;

  TG4Globals::Exception(
    "Attempt to assign TG4RunManager singleton.");
    
  return *this;  
}    

// private methods

void TG4RunManager::CreateGeantUI()
{
// Creates interactive Geant4.
// ---

  if (!fGeantUISession)
  {
    // create geant4 UI
    G4UImanager* pUI = G4UImanager::GetUIpointer();  
    if (fARGC == 1) {
#ifdef G4UI_USE_GAG
      fGeantUISession = new G4UIGAG();
#else      
      fGeantUISession = new G4UIterminal(); 
#endif      
    }  
    else if (strcmp (fARGV[1], "dumb") == 0) {
      fGeantUISession = new G4UIterminal(); 
    }
#ifdef G4UI_USE_WO
    else if (strcmp (fARGV[1], "Wo") == 0) {
      fGeantUISession = new G4UIWo(fARGC, fARGV); 
    }
#endif
#ifdef G4UI_USE_XM
    else if (strcmp (fARGV[1], "Xm") == 0) {
      fGeantUISession = new G4UIXm(fARGC, fARGV); 
    }
#endif
#ifdef G4UI_USE_XAW
    else if (strcmp (fARGV[1], "Xaw") == 0) {
      fGeantUISession = new G4UIXaw(fARGC, fARGV); 
    }
#endif 
#ifdef G4UI_USE_GAG
    else if (strcmp (fARGV[1], "GAG") == 0) {
      fGeantUISession = new G4UIGAG (); 
    }
#endif 
    if (fGeantUISession) {   
      pUI->SetSession(fGeantUISession); 
    }
  }
}

void TG4RunManager::CreateRootUI()
{
// Creates interactive Root.
// ---

  if (!fRootUISession) 
  {
    // create session if it does not exist  
    fRootUISession = new TRint("aliroot", 0, 0, 0, 0);

    // set ownership of Root UI
    fRootUIOwner = true;
  }
}

// public methods

void TG4RunManager::Initialize()
{
// Initializes G4.
// ---

  // create physics constructor
  // (this operation has to precede the "Init" phase)
  TG4PhysicsManager* physicsManager = TG4PhysicsManager::Instance();
  physicsManager->CreatePhysicsConstructors();

  // initialize Geant4 
  fRunManager->Initialize();
  
  // activate/inactivate physics processes
  // (this operation is not allowed in "Init" phase)
  physicsManager->SetProcessActivation();
}

void TG4RunManager::ProcessEvent()
{
// Not yet implemented.
// ---

  TG4Globals::Warning("TG4RunManager::ProcessEvent(): is not yet implemented.");
}
    
void TG4RunManager::ProcessRun(G4int nofEvents)
{
// Processes Geant4 run.
// ---

  fRunManager->BeamOn(nofEvents); 
}
    
void TG4RunManager::StartGeantUI()
{ 
// Starts interactive/batch Geant4.
// ---

  if (!fGeantUISession) CreateGeantUI();
  if (fGeantUISession) {  
    // interactive session
    G4cout << "Welcome back in Geant4" << G4endl;
    fGeantUISession->SessionStart();
    G4cout << "Welcome back in Root" << G4endl;  
  }
  else {
    // execute Geant4 macro if file is specified as an argument 
    G4String fileName = fARGV[1];
    ProcessGeantMacro(fileName);
  }
}

void TG4RunManager::StartRootUI()
{
// Starts interactive Root.
// ---

  if (!fRootUISession) CreateRootUI();
  if (fRootUISession) { 
    G4cout << "Welcome back in Root" << G4endl;
    fRootUISession->Run(kTRUE);
    G4cout << "Welcome back in Geant4" << G4endl;  
  }
}
 
void TG4RunManager::ProcessGeantMacro(G4String macroName)
{
// Processes Geant4 macro.
// ---

  G4String command = "/control/execute " + macroName;
  ProcessGeantCommand(command);
}
 
void TG4RunManager::ProcessRootMacro(G4String macroName)
{
// Processes Root macro.
// ---

  // load macro file
  G4String macroFile = macroName;
  macroFile.append(".C");
  gROOT->LoadMacro(macroFile);

  // execute macro function
  G4String macroFunction = macroName;
  macroFunction.append("()");
  gInterpreter->ProcessLine(macroFunction);
}
 
void TG4RunManager::ProcessGeantCommand(G4String command)
{
// Processes Geant4 command.
// ---

  G4UImanager* pUI = G4UImanager::GetUIpointer();  
  pUI->ApplyCommand(command);
}

void TG4RunManager::ProcessRootCommand(G4String command)
{
// Processes Root command.
// ---

  gInterpreter->ProcessLine(command);
}

void TG4RunManager::UseG3Defaults() 
{
// Controls G3 defaults usage.
// ---

  TG4GeometryManager::Instance()->UseG3TrackingMediaLimits();
  TG4G3PhysicsManager::Instance()->SetG3DefaultCuts();
  TG4G3PhysicsManager::Instance()->SetG3DefaultControls();
}

Int_t TG4RunManager::CurrentEvent() const
{
// Returns the number of the current event.
// ---

  G4int eventID = fRunManager->GetCurrentEvent()->GetEventID();
  return eventID;
}
