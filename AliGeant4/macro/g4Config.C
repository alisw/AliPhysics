// $Id$
//
// Configuration macro for running aliroot with Geant4.
// All AliRoot specifics are defined in g4ConfigCommon.C
// 
// In order to run aliroot with Geant4, you have to Initialize
// AliRun with this g4Config.C
//
// aliroot
// root [0] gAlice->Init("g4Config.C");	       
// root [1] gAlice->Run();	
//
// You can also start from a mini GUI - g4menu.C.
// See description in this macro.
 	

void Config()
{
  // ============================= 
  // Geant4
  // ============================= 

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
  if (!gMC) {

    TG4RunConfiguration* runConfiguration 
      = new TG4RunConfiguration(true);
        
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
  ConfigCommon(false);
}
