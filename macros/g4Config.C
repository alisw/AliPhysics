// $Id$
//
// Configuration macro for running aliroot with Geant4.
// All AliRoot specifics are defined in g4ConfigCommon.C
// 
// In order to run aliroot with Geant4, you have to Initialize
// AliRun with this g4Config.C
//
// aliroot
// root [0] AliSimulation sim("$ALICE_ROOT/macros/g4Config.C");	       
// root [1] .... customize your simulation setting here	       
// root [1] sim.Run(10);
//
// You can also start from a mini GUI - g4menu.C.
// See description in this macro.
//
// By I. Hrivnacova, IPN Orsay
 	

void Config()
{
  cout << "Running g4Config.C ... " << endl;

  // AliRoot setup
  //
  gROOT->LoadMacro("$ALICE_ROOT/macros/g4ConfigCommon1.C");
  ConfigCommon1();

  // Load Geant4 + Geant4 VMC libraries
  //
  if (gClassTable->GetID("TGeant4") == -1) {
    // Load Geant4 libraries 
    if (!gInterpreter->IsLoaded("$ALICE/geant4_vmc/examples/macro/g4libs.C")) {
      gROOT->LoadMacro("$ALICE/geant4_vmc/examples/macro/g4libs.C");
      gInterpreter->ProcessLine("g4libs()");
    }
  }    

  // Create Geant4 VMC
  //  
  TGeant4 *geant4 = 0;
  if ( ! gMC ) {
    TG4RunConfiguration* runConfiguration 
      = new TG4RunConfiguration("geomRoot", 
                                "QGSP_BERT_EMV+optical", 
                                "specialCuts+stackPopper+stepLimiter",
                                 true, false);
//      = new TG4RunConfiguration("geomRootToGeant4",
//                                "emStandard+optical", 
//                                "specialCuts+specialControls+stackPopper+stepLimiter",
//                                 true, false);
      

    geant4 = new TGeant4("TGeant4", "The Geant4 Monte Carlo", runConfiguration);
    cout << "Geant4 has been created." << endl;
  } 
  else {
    cout << "Monte Carlo has been already created." << endl;
  }  

  // Customization of Geant4 VMC
  //

  geant4->ProcessGeantCommand("/mcPhysics/rangeCuts 0.01 mm"); 
  geant4->ProcessGeantCommand("/mcVerbose/all 1");  

  // Uncomment this line to get a detail info from each step 
  //geant4->ProcessGeantCommand("/tracking/verbose 1");  
  
  // More info from the physics list
  // the verbosity level is passed to all contained physics lists and their
  // physics builders
  //geant4->ProcessGeantCommand("/mcVerbose/composedPhysicsList 2");  
  
  // More info from optical processes
  //geant4->ProcessGeantCommand("/mcVerbose/opticalPhysicsList 3");  
  
  // More info from geometry building
  //geant4->ProcessGeantCommand("/mcVerbose/geometryManager 1");  

  // More info from setting geometry properties (in materials and surfaces)
  // for optical physics
  //geant4->ProcessGeantCommand("/mcVerbose/opGeometryManager 1");  
  
  // More info about regions construction 
  // and conversion of VMC cuts in cuts in range per regions 
  //geant4->ProcessGeantCommand("/mcVerbose/regionsManager 2");
  
  // Suppress verbose info from tracks which reached maximum number of steps
  // (default value is 30000)  
  geant4->ProcessGeantCommand("/mcTracking/loopVerbose 0"); 
    
  // AliRoot event generator
  // (it has to be created after MC, as it may use decayer via VMC)
  //
  gROOT->LoadMacro("$ALICE_ROOT/macros/g4ConfigCommon2.C");
  ConfigCommon2();

  cout << "Running g4Config.C finished ... " << endl;
}
