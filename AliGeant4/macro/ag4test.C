// $Id$
//
// Root macro that runs G4 test macro (alirunN.in) specified by macroName.
// Available tests:
// alirun0.in - standard AliRoot run - all detectors, AliRoot event generator
// alirun1.in - interactive particle gun, event visualization - TPC
// alirun2.in - lego run - TRD
// alirun3.in - geometry test - all detectors (excluded MANY)");
// alirun4.in - geometry browser - FRAME;

void ag4test(const char* macroName)
{
  // Load Geant4 + Geant4 VMC libraries
  //
  if (!gInterpreter->IsLoaded("$ALICE/geant4_vmc/examples/macro/g4libs.C")) 
    gROOT->LoadMacro("$ALICE/geant4_vmc/examples/macro/g4libs.C");
  g4libs();
 
  // Load AliGeant4 library
  //
  gSystem->Load("libAliGeant4.so");

  // Create AliGeant4 VMC
  //  
  if (!gMC) {
  
    // AliRoot run configuration for Geant4
    AliRunConfiguration* runConfiguration 
      = new AliRunConfiguration();
  
    // TGeant4
    TGeant4* geant4
      = new TGeant4("TGeant4", "The Geant4 Monte Carlo", runConfiguration);

    cout << "Geant4 has been created." << endl;
  }      
  else 
    cout << "Monte Carlo has been already created." << endl;

  // Process G4 macro
  // 
  ((TGeant4*)gMC)->ProcessGeantMacro(macroName);
} 
