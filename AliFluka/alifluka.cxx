// $Id$

#include "TFluka.h"
#include "FGeometryInit.hh"

#include "AliDetConstructionHelper.h"
#include "AliRun.h"

#include <TROOT.h>

#define flukam flukam_

extern "C" void flukam(const G4int & GeoFlag);

int main(int argc, char** argv) 
{
  // ROOT  ===================
  TROOT aTROOT("Alice","Alice Flugg prototype Root I/O");

  // AliRun
  AliRun* run
    = new AliRun("gAlice","The Alice run manager");
  G4cout << "AliRun has been created." << G4endl;
  
  // TFluka
  TFluka* fluka
    = new TFluka("TFluka", "The Fluka Monte Carlo");
  G4cout << "TFluka has been created." << G4endl;

  // Detector construction helper
  AliDetConstructionHelper* detHelper = new AliDetConstructionHelper();
  G4cout << "Detector construction helper has been created." << G4endl;

  // Flugg
  FGeometryInit* theFGeometryInit = FGeometryInit::GetInstance();
  theFGeometryInit->setDetConstruction(detHelper->DetConstruction());
  G4cout << "Detector construction has been set to Flugg." << G4endl;

  //flag for geometry:
  const G4int flag = 1;
        // 1 for GEANT4
        // 0 for FLUKA
        // 2 for Rubia

  //call fluka
  flukam(flag);

  delete run;
  
  return 0;
}
