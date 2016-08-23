// Configuration of simulation

void Config()
{
  // The common part contains the description of geometry,
  // magnetic field and generator. The initialization of
  // the random generator seed is also there

  gROOT->Macro("../commonConfig.C");


  //=======================================================================
  // Create Geant4 transport
  TG4RunConfiguration* runConfiguration 
    = new TG4RunConfiguration("geomRoot", 
			      "FTFP_BERT_EMV+optical", 
			      "specialCuts+stackPopper+stepLimiter",
			      true);
  
  TGeant4 * geant4 = new TGeant4("TGeant4", 
				 "The Geant4 Monte Carlo : FTFP_BERT_EMV-EMCAL", 
				 runConfiguration);

  // Customization of Geant4 VMC from the GRID test production
  //
  geant4->ProcessGeantCommand("/control/verbose 2");  
  geant4->ProcessGeantCommand("/mcVerbose/all 1");  
  geant4->ProcessGeantCommand("/mcVerbose/geometryManager 1");  
  geant4->ProcessGeantCommand("/mcVerbose/opGeometryManager 1");  
  geant4->ProcessGeantCommand("/mcTracking/loopVerbose 1");     
  geant4->ProcessGeantCommand("/mcPhysics/rangeCuts 0.01 mm"); 
  
  geant4->ProcessGeantCommand("/mcVerbose/composedPhysicsList 2");  
  geant4->ProcessGeantCommand("/mcTracking/skipNeutrino true");
  geant4->ProcessGeantCommand("/mcDet/setIsMaxStepInLowDensityMaterials true");
  geant4->ProcessGeantCommand("/mcDet/setMaxStepInLowDensityMaterials 10 m");
  geant4->ProcessGeantCommand("/mcMagField/setConstDistance 1 mm");
  //
  // optical
  //
  geant4->ProcessGeantCommand("/process/optical/verbose 0");
  geant4->ProcessGeantCommand("/process/optical/processActivation Scintillation 0");
  geant4->ProcessGeantCommand("/process/optical/processActivation OpWLS 0");
  geant4->ProcessGeantCommand("/process/optical/processActivation OpMieHG 0");
  geant4->ProcessGeantCommand("/process/optical/setTrackSecondariesFirst Cerenkov 0");
  geant4->ProcessGeantCommand("/mcMagField/stepperType NystromRK4");
  
  //
  // PAI for TRD
  //
  // Geant4 VMC >= v3.2
  geant4->ProcessGeantCommand("/mcPhysics/emModel/setEmModel  PAI");
  geant4->ProcessGeantCommand("/mcPhysics/emModel/setRegions  TRD_Gas-mix");
  geant4->ProcessGeantCommand("/mcPhysics/emModel/setParticles  all");
  
  //
  // Precise Msc for EMCAL
  //
  // Geant4 VMC >= v3.2
  geant4->ProcessGeantCommand("/mcPhysics/emModel/setEmModel  SpecialUrbanMsc");
  geant4->ProcessGeantCommand("/mcPhysics/emModel/setRegions  EMCAL_Lead$ EMCAL_Scintillator$");
  geant4->ProcessGeantCommand("/mcPhysics/emModel/setParticles  e- e+");
  
  // Activate printing size of used memory per event
  geant4->ProcessGeantCommand("/mcEvent/printMemory true");

  cout << "End of Geant4 section" << endl;


  //=======================================================================
  // Set External decayer
  AliDecayer *decayer = new AliDecayerPythia();
  decayer->SetForceDecay(kAll);
  decayer->Init();

  //forbid some decays
  AliPythia * py= AliPythia::Instance();
  py->SetMDME(737,1,0); //forbid D*+->D+ + pi0
  py->SetMDME(738,1,0); //forbid D*+->D+ + gamma
  
  for(Int_t d=747; d<=762; d++){ 
    py->SetMDME(d,1,0);
  }
  
  for(Int_t d=764; d<=807; d++){ 
    py->SetMDME(d,1,0);
  }

  TVirtualMC * vmc = TVirtualMC::GetMC();
  vmc->SetExternalDecayer(decayer);
  
  //=======================================================================
  // ************* STEERING parameters FOR ALICE SIMULATION **************
  // --- Specify event type to be tracked through the ALICE setup
  // --- All positions are in cm, angles in degrees, and P and E in GeV

  vmc->SetProcess("DCAY",1);
  vmc->SetProcess("PAIR",1);
  vmc->SetProcess("COMP",1);
  vmc->SetProcess("PHOT",1);
  vmc->SetProcess("PFIS",0);
  vmc->SetProcess("DRAY",0);
  vmc->SetProcess("ANNI",1);
  vmc->SetProcess("BREM",1);
  vmc->SetProcess("MUNU",1);
  vmc->SetProcess("CKOV",1);
  vmc->SetProcess("HADR",1);
  vmc->SetProcess("LOSS",2);
  vmc->SetProcess("MULS",1);
  vmc->SetProcess("RAYL",1);
  
  Float_t cut = 1.e-3;        // 1MeV cut by default
  Float_t tofmax = 1.e10;
  
  vmc->SetCut("CUTGAM", cut);
  vmc->SetCut("CUTELE", cut);
  vmc->SetCut("CUTNEU", cut);
  vmc->SetCut("CUTHAD", cut);
  vmc->SetCut("CUTMUO", cut);
  vmc->SetCut("BCUTE",  cut); 
  vmc->SetCut("BCUTM",  cut); 
  vmc->SetCut("DCUTE",  cut); 
  vmc->SetCut("DCUTM",  cut); 
  vmc->SetCut("PPCUTM", cut);
  vmc->SetCut("TOFMAX", tofmax); 
  
}
