// $Id$
//
// Configuration macro for running aliroot with Geant4
// with primary events read from external file.
//
// By I. Hrivnacova, IPN Orsay
 	

void Config(const char * det)
{
  cout << "Running g4Config.C ... " << endl;

  // AliRoot setup
  //
  commonConfig(det);

  // TPC primary ionization 
  AliTPC* TPC = (AliTPC*)gAlice->GetDetector("TPC");
  if ( ! TPC )
    cerr << "Cannot get TPC detector" << endl;
  else  { 
    cerr << "Setting TPC primary ionization" << endl;
    TPC->SetPrimaryIonisation(); // not used with Geant3
  }  

  // Load Geant4 + Geant4 VMC libraries
  //
  // std::string g4libmacro("$G4VMCINSTALL/share/examples/macro/g4libs.C");
  // if (gClassTable->GetID("TGeant4") == -1) {
  //   // Load Geant4 libraries
  //   if (!gInterpreter->IsLoaded(g4libmacro.c_str())) {
  //     gROOT->LoadMacro(g4libmacro.c_str());
  //     gInterpreter->ProcessLine("g4libs()");
  //   }
  // }

  // Create Geant4 VMC
  //  
  TGeant4 *geant4 = 0;
  if ( ! gMC ) {
    TG4RunConfiguration* runConfiguration 
      = new TG4RunConfiguration("geomRoot", 
                                "FTFP_BERT_EMV+optical", 
                                "specialCuts+stackPopper+stepLimiter",
                                 true);
//      = new TG4RunConfiguration("geomRootToGeant4",
//                                "FTFP_BERT_EMV+optical", 
//                                "specialCuts+specialControls+stackPopper+stepLimiter",
//                                 true);
      
    geant4 = new TGeant4("TGeant4", 
                         "The Geant4 Monte Carlo : FTFP_BERT_EMV+optical", 
                         runConfiguration);
             // Repeat physics selection in the title; to be removed
             // with new geant4_vmc tag (1.13)            
                         
    cout << "Geant4 has been created." << endl;
  } 
  else {
    cout << "Monte Carlo has been already created." << endl;
  }  

  // Customization of Geant4 VMC
  //

    geant4->ProcessGeantCommand("/mcVerbose/all 1");  
    geant4->ProcessGeantCommand("/mcVerbose/geometryManager 1");  
    geant4->ProcessGeantCommand("/mcVerbose/opGeometryManager 1");  
    geant4->ProcessGeantCommand("/mcTracking/loopVerbose 1");     
    geant4->ProcessGeantCommand("/mcPhysics/rangeCuts 0.01 mm"); 
    geant4->ProcessGeantCommand("/mcTracking/skipNeutrino true");
    
    // Activate step limit defined in low density materials
    // (the default value is 10 cm)
    //geant4->ProcessGeantCommand("/mcDet/setIsMaxStepInLowDensityMaterials true");
    //geant4->ProcessGeantCommand("/mcDet/setMaxStepInLowDensityMaterials 1 cm");

    // Activate step limit defined in tracking media
    // (Note: this slows down simulation significantly)
    //geant4->ProcessGeantCommand("/mcDet/setIsUserMaxStep true");
   
    // for Geant4 <= 9.4.p03
    //geant4->ProcessGeantCommand("/mcPhysics/selectOpProcess Scintillation");
    //geant4->ProcessGeantCommand("/mcPhysics/setOpProcessActivation false");
    // for Geant4 >= 9.5
    geant4->ProcessGeantCommand("/optics_engine/selectOpProcess Scintillation");
    geant4->ProcessGeantCommand("/optics_engine/setOpProcessUse false");
    geant4->ProcessGeantCommand("/optics_engine/selectOpProcess OpWLS");
    geant4->ProcessGeantCommand("/optics_engine/setOpProcessUse false");
    geant4->ProcessGeantCommand("/optics_engine/selectOpProcess OpMieHG");
    geant4->ProcessGeantCommand("/optics_engine/setOpProcessUse false");
    geant4->ProcessGeantCommand("/optics_engine/selectOpProcess Cerenkov");
    geant4->ProcessGeantCommand("/optics_engine/setTrackSecondariesFirst false");   

    // Activate saving random engine status
    // (the file per event will be re-written with each new event)
    //gAlice->GetMCApp()->SetSaveRndmStatus(kTRUE);
    //geant4->ProcessGeantCommand("/mcRun/saveRandom true");

    // Activate saving random engine status for each event
    // (a new file will be written for each event)
    //gAlice->GetMCApp()->SetSaveRndmStatusPerEvent(kTRUE);
    //geant4->ProcessGeantCommand("/mcRun/saveRandom true");
    //geant4->ProcessGeantCommand("/mcEvent/saveRandom true");
    
    // Activate printing size of used memory per event
    geant4->ProcessGeantCommand("/mcEvent/printMemory true");

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
  // geant4->ProcessGeantCommand("/mcVerbose/regionsManager 2");
  // geant4->ProcessGeantCommand("/mcRegions/print true");
  
  // Suppress verbose info from tracks which reached maximum number of steps
  // (default value is 30000)  
  //geant4->ProcessGeantCommand("/mcTracking/loopVerbose 0"); 
    
  //
  // Set apply cuts 
/*
  geant4->ProcessGeantCommand("/run/particle/applyCuts");  
  // geant4->ProcessGeantCommand("/mcVerbose/geometryManager 2");  

  geant4->ProcessGeantCommand("/mcVerbose/composedPhysicsList 2");  
  geant4->ProcessGeantCommand("/mcDet/volNameSeparator !");
  geant4->ProcessGeantCommand("/mcPhysics/setStackPopperSelection e+ e- pi+ pi- kaon+ kaon- gamma");
  //geant4->ProcessGeantCommand("/tracking/verbose 1");  

  geant4->ProcessGeantCommand("/mcControl/g3Defaults");
!!!!n Generates warnings:
>>> Event 0
G4ProcessTable::Insert : arguments are 0 pointer 
G4ProcessTable::Insert : arguments are 0 pointer 
G4ProcessTable::Insert : arguments are 0 pointer 
G4ProcessTable::Insert : arguments are 0 pointer 
G4ProcessTable::Insert : arguments are 0 pointer 

*/

  cout << "Running g4Config.C finished ... " << endl;
}
