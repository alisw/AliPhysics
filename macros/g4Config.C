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
//
// By I. Hrivnacova, IPN Orsay
 	

void Config()
{
  cout << "Running g4Config.C ... " << endl;

  // AliRoot setup
  //
  gROOT->LoadMacro("$ALICE_ROOT/macros/g4ConfigCommon.C");
  ConfigCommon();

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
  if (!gMC) {
    TG4RunConfiguration* runConfiguration 
      = new TG4RunConfiguration("geomRoot", 
                                "emStandard+optical", 
                                "specialCuts+specialControls+stackPopper",
                                 true);
//      = new TG4RunConfiguration("geomRootToGeant4",
//                                "emStandard+optical", 
//                                "specialCuts+specialControls+stackPopper",
//                                 true);
      

    geant4 = new TGeant4("TGeant4", "The Geant4 Monte Carlo", runConfiguration);
    cout << "Geant4 has been created." << endl;
  } else {
    cout << "Monte Carlo has been already created." << endl;
  }  
    
  //
  // Set External decayer
  TVirtualMCDecayer *decayer = new AliDecayerPythia();

  decayer->SetForceDecay(kAll);
  decayer->Init();
  gMC->SetExternalDecayer(decayer);
  
  //=======================================================================
  // ************* STEERING parameters FOR ALICE SIMULATION **************
  // --- Specify event type to be tracked through the ALICE setup
  // --- All positions are in cm, angles in degrees, and P and E in GeV

  gMC->SetProcess("DCAY",1);
  gMC->SetProcess("PAIR",1);
  gMC->SetProcess("COMP",1);
  gMC->SetProcess("PHOT",1);
  gMC->SetProcess("PFIS",0);
  gMC->SetProcess("DRAY",0);
  gMC->SetProcess("ANNI",1);
  gMC->SetProcess("BREM",1);
  gMC->SetProcess("MUNU",1);
  gMC->SetProcess("CKOV",1);
  gMC->SetProcess("HADR",1);
  gMC->SetProcess("LOSS",2);
  gMC->SetProcess("MULS",1);
  //gMC->SetProcess("RAYL",1);

  Float_t cut = 1.e-3;        // 1MeV cut by default
  Float_t tofmax = 1.e10;

  gMC->SetCut("CUTGAM", cut);
  gMC->SetCut("CUTELE", cut);
  gMC->SetCut("CUTNEU", cut);
  gMC->SetCut("CUTHAD", cut);
  gMC->SetCut("CUTMUO", cut);
  gMC->SetCut("BCUTE",  cut); 
  gMC->SetCut("BCUTM",  cut); 
  gMC->SetCut("DCUTE",  cut); 
  gMC->SetCut("DCUTM",  cut); 
  gMC->SetCut("PPCUTM", cut);
  gMC->SetCut("TOFMAX", tofmax); 

  // Set apply cuts 
  geant4->ProcessGeantCommand("/run/particle/applyCuts");  
  // geant4->ProcessGeantCommand("/mcVerbose/geometryManager 2");  
  geant4->ProcessGeantCommand("/mcVerbose/composedPhysicsList 2");  
  geant4->ProcessGeantCommand("/mcDet/volNameSeparator !");
  geant4->ProcessGeantCommand("/mcControl/g3Defaults");
  geant4->ProcessGeantCommand("/mcPhysics/setStackPopperSelection e+ e- pi+ pi- kaon+ kaon- gamma");
  // geant4->ProcessGeantCommand("/tracking/verbose 1");  

  cout << "Running g4Config.C finished ... " << endl;
}
