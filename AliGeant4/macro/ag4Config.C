// $Id$
//
// Configuration macro for running aliroot with AliGeant4.
// All AliRoot specifics are defined in g4ConfigCommon.C
// 
// In order to run aliroot with AliGeant4, you have to Initialize
// AliRun with this ag4Config.C
//
// aliroot
// root [0] gAlice->Init("ag4Config.C");	       
// root [1] gAlice->Run();	
//
// In difference from g4Config.C, TGeant4 is initialized
// with AliRunConfiguration instead of TG4RunConfiguration.


void Config(Bool_t interactiveSetup = true)
{
  // ============================= 
  // Geant4
  // ============================= 

  // Load Geant4 + Geant4 VMC libraries
  //
  if (gClassTable->GetID("TGeant4") == -1) {
    // load dynamically only if executable was not linked
    // with geant4 libraries

    // Load Geant4 and AliRoot steer libraries
    if (!gInterpreter->IsLoaded("$ALICE/geant4_vmc/examples/macro/g4libs.C")) 
      gROOT->LoadMacro("$ALICE/geant4_vmc/examples/macro/g4libs.C");
    gInterpreter->ProcessLine("g4libs()");

    // Load AliGeant4 library
    //
    gSystem->Load("libAliGeant4.so");
  }

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

  // ============================= 
  // AliRoot setup
  // ============================= 

  gROOT->LoadMacro("g4ConfigCommon.C");
  ConfigCommon(interactiveSetup);
}
