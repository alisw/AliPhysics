// $Id$

#include "AliRunConfiguration.h"
#include "AliFiles.h"
#include "AliRun.h"

#include "TGeant4.h"
#include "TG4RunManager.h"

#include <TROOT.h>
#include <TRint.h>

extern void InitGui();

#include <globals.hh>

int main(int argc, char** argv) 
{
  // ROOT  ===================
#ifdef G4VIS_USE_OPACS
  // Root graphics does not work when OPACS graphics is build 
  TROOT aTROOT("Alice","Alice G4 prototype Root I/O");
#else
  VoidFuncPtr_t initfuncs[] = { InitGui, 0 };
  TROOT aTROOT("Alice","Alice G4 prototype Root I/O",initfuncs);
#endif

  // ALICE ======================

  // AliRun
  AliRun* run
    = new AliRun("gAlice","The Alice run manager");
  G4cout << "AliRun has been created." << G4endl;

  // AliRunConfiguration for Geant4
  AliRunConfiguration* runConfiguration 
    = new AliRunConfiguration();
  G4cout << "AliRunConfiguration has been created." << G4endl;
   
  // Geant4 ======================

  // TGeant4
  TGeant4* geant4 
    = new TGeant4("TGeant4", "The Geant4 Monte Carlo",
                   runConfiguration, argc, argv );
  G4cout << "TGeant4 has been created." << G4endl;
  
  // start UI ===================

  TG4RunManager* runManager = TG4RunManager::Instance();

  // Root interactive session
  //runManager->StartRootUI();

  // Geant4 interactive session
  runManager->StartGeantUI();

  delete run;
  //runConfiguration is deleted in TG4RunManager
  //geant4 is deleted in AliRun destructor 

  G4cout << "Everything has been deleted." << G4endl;
  return 0;
}

