/**
 * @file   PWGLF/FORWARD/analysis2/sim/Config.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Wed Oct 15 13:03:28 2014
 * 
 * @brief  Configuration of the simulation back-end 
 *
 * @note Do not modify this script. 
 *
 * This script depends on the two global variables detCfg and grp
 * already being defined, typically by executing the scripts GRP.C and
 * DetConfig.C from Simulate.C. 
 *
 * New event generator set-ups should be added to the class Setup. 
 * 
 * 
 */

// -------------------------------------------------------------------
/** 
 * Class that defines the set-up.  
 *
 * - The seed of the random number generator is read from the
 *   environment if present.
 *
 * - The event generator type is read from the environment if present.
 *   Otherwise we try to deduce it from the global object "grp".
 * 
 * - The impact parameter range is read from the environment if present. 
 */
struct Setup
{
  TString runType;    // Event generator chosen
  UInt_t  seed;       // Random number seed - (env)
  Float_t minB;       // Least imp.param. (env)
  Float_t maxB;       // Largest imp.param. (env)  
  TString backend;    // The backend to use (geant3 or geant4)
  /** 
   * Get the value of an environment variable as a float 
   * 
   * @param envName Enviroment variable name 
   * @param def     Default value 
   *
   * @return As float, or default
   */
  static Float_t Env2Float(const char* envName, Float_t def) 
  { 
    TString val(gSystem->Getenv(envName));
    if (val.IsNull()) return def;
    return val.Atof();
  }
  /** 
   * Get the value of an environment variable as a unsigned int 
   * 
   * @param envName Enviroment variable name 
   * @param def     Default value 
   *
   * @return As unsigned int, or default
   */
  static UInt_t Env2UInt(const char* envName, UInt_t def)
  {
    TString val(gSystem->Getenv(envName));
    if (val.IsNull()) return def;
    return UInt_t(val.Atoll());
  }
  /** 
   * Get the value of an environment variable as a int
   * 
   * @param envName Enviroment variable name 
   * @param def     Default value 
   *
   * @return As int, or default
   */
  static UInt_t Env2Int(const char* envName, Int_t def)
  {
    TString val(gSystem->Getenv(envName));
    if (val.IsNull()) return def;
    return val.Atoi();
  }

  /** 
   * Constructor - retrieves needed information from CDB manager and
   * environment.
   */
  Setup() 
    : runType(""),
      seed(0),
      minB(0), 
      maxB(100)
  {
    TDatime now;
    runType = gSystem->Getenv("CONFIG_RUN_TYPE");
    backend = gSystem->Getenv("CONFIG_BACKEND");
    seed    = Env2UInt("CONFIG_SEED", now.Get());
    minB    = Env2Float("CONFIG_BMIN", 0);
    maxB    = Env2Float("CONFIG_BMAX", 100);
    if (runType[0] == 'k') runType.Remove(0,1);
    runType.ToLower();

    // gROOT->Macro("GetGRP.C");


    if (runType.IsNull() || runType == "default") DeduceRunType();

    Print();
  }
  /** 
   * Prinf information 
   * 
   */
  void Print()
  {
    Printf("=======================================================\n"
	   "           Set-up of the simulation\n");
    grp->Print();
    Printf("Run type: '%s'", runType.Data());
    Printf("Backend:  '%s'", backend.Data());
    Printf("b range:  [%4.1f,%4.1f]fm", minB, maxB);
    Printf("Seed:     %d", seed);
    Printf("\n"
	   "=======================================================");
  }
  Bool_t IsGeant3() const 
  { 
    return (backend.EqualTo("geant3", TString::kIgnoreCase) ||     
	    backend.EqualTo("g3", TString::kIgnoreCase));
  }
  Bool_t IsGeant4() const 
  { 
    return (backend.EqualTo("geant4", TString::kIgnoreCase) ||     
	    backend.EqualTo("g4", TString::kIgnoreCase));
  }
  /** 
   * Create our simulation back-end (Geant3 or 4)
   */
  void MakeBackend()
  {
    if (!IsGeant4()) {
      gSystem->Load("libgeant321");
      new TGeant3TGeo("C++ Interface to Geant3");
    }
    else { 
      gROOT->Macro("$ALICE/geant4_vmc/examples/macro/g4libs.C");
      // Create G4 VMC
      TGeant4 *g4 = 0;
      
      TG4RunConfiguration* runConfiguration 
	= new TG4RunConfiguration("geomRoot", 
				  "FTFP_BERT_EMV+optical", 
				  "specialCuts+stackPopper+stepLimiter",
				  true);
      TGeant4* g4 = new TGeant4("TGeant4", 
				"The Geant4 Monte Carlo : "
				"FTFP_BERT_EMV+optical", 
				runConfiguration);
      // Customization of Geant4 VMC
      g4->ProcessGeantCommand("/mcVerbose/all 1");  
      g4->ProcessGeantCommand("/mcVerbose/geometryManager 1");  
      g4->ProcessGeantCommand("/mcVerbose/opGeometryManager 1");  
      g4->ProcessGeantCommand("/mcTracking/loopVerbose 1");     
      g4->ProcessGeantCommand("/mcPhysics/rangeCuts 0.01 mm"); 
      g4->ProcessGeantCommand("/mcTracking/skipNeutrino true");
      
      // Activate step limit defined in low density materials
      // (the default value is 10 cm)
      g4->ProcessGeantCommand("/mcDet/setIsMaxStepInLowDensityMaterials true");
      g4->ProcessGeantCommand("/mcDet/setMaxStepInLowDensityMaterials 10 m");
      
      // Activate step limit defined in tracking media
      // (Note: this slows down simulation significantly)
      // g4->ProcessGeantCommand("/mcDet/setIsUserMaxStep true");
      
      // for G4 <= 9.4.p03
      // g4->ProcessGeantCommand("/mcPhysics/selectOpProcess Scintillation");
      // g4->ProcessGeantCommand("/mcPhysics/setOpProcessActivation false");
      // for G4 >= 9.5
      g4->ProcessGeantCommand("/optics_engine/selectOpProcess Scintillation");
      g4->ProcessGeantCommand("/optics_engine/setOpProcessUse false");
      g4->ProcessGeantCommand("/optics_engine/selectOpProcess OpWLS");
      g4->ProcessGeantCommand("/optics_engine/setOpProcessUse false");
      g4->ProcessGeantCommand("/optics_engine/selectOpProcess OpMieHG");
      g4->ProcessGeantCommand("/optics_engine/setOpProcessUse false");
      g4->ProcessGeantCommand("/optics_engine/selectOpProcess Cerenkov");
      
      // Activate saving random engine status
      // (the file per event will be re-written with each new event)
      // gAlice->GetMCApp()->SetSaveRndmStatus(kTRUE);
      // g4->ProcessGeantCommand("/mcRun/saveRandom true");
      
      // Activate saving random engine status for each event
      // (a new file will be written for each event)
      // gAlice->GetMCApp()->SetSaveRndmStatusPerEvent(kTRUE);
      // g4->ProcessGeantCommand("/mcRun/saveRandom true");
      // g4->ProcessGeantCommand("/mcEvent/saveRandom true");
      
      // Activate printing size of used memory per event
      g4->ProcessGeantCommand("/mcEvent/printMemory true");
      
      // Uncomment this line to get a detail info from each step 
      // g4->ProcessGeantCommand("/tracking/verbose 1");  
      
      // More info from the physics list
      // the verbosity level is passed to all contained physics lists and their
      // physics builders
      // g4->ProcessGeantCommand("/mcVerbose/composedPhysicsList 2");  
      
      // More info from optical processes
      // g4->ProcessGeantCommand("/mcVerbose/opticalPhysicsList 3");  
      
      // More info from geometry building
      // g4->ProcessGeantCommand("/mcVerbose/geometryManager 1");  
      
      // More info from setting geometry properties (in materials and surfaces)
      // for optical physics
      // g4->ProcessGeantCommand("/mcVerbose/opGeometryManager 1");  
      
      // More info about regions construction 
      // and conversion of VMC cuts in cuts in range per regions 
      // g4->ProcessGeantCommand("/mcVerbose/regionsManager 2");
      // g4->ProcessGeantCommand("/mcRegions/print true");
      
      // Suppress verbose info from tracks which reached maximum number of steps
      // (default value is 30000)  
      // g4->ProcessGeantCommand("/mcTracking/loopVerbose 0"); 
      
      //
      // Set apply cuts 
      // g4->ProcessGeantCommand("/process/em/applyCuts true");
      // g4->ProcessGeantCommand("/mcVerbose/geometryManager 2");  
      /*
	g4->ProcessGeantCommand("/mcVerbose/composedPhysicsList 2");  
	g4->ProcessGeantCommand("/mcDet/volNameSeparator !");
	g4->ProcessGeantCommand("/mcPhysics/setStackPopperSelection "
	                        "e+ e- pi+ pi- kaon+ kaon- gamma");
	//g4->ProcessGeantCommand("/tracking/verbose 1");  
	
	g4->ProcessGeantCommand("/mcControl/g3Defaults");
	!!!!n Generates warnings:
	>>> Event 0
	G4ProcessTable::Insert : arguments are 0 pointer 
	G4ProcessTable::Insert : arguments are 0 pointer 
	G4ProcessTable::Insert : arguments are 0 pointer 
	G4ProcessTable::Insert : arguments are 0 pointer 
	G4ProcessTable::Insert : arguments are 0 pointer 
	
      */
    } 
  }

  /** 
   * Set the default generator based on the beam type 
   *
   * - p-p PYTHIA
   * - p-A or A-p DPMJet
   * - A-A Hijing 
   */
  void DeduceRunType()
  {
    if      (grp->IsPP())                runType = "pythia";
    else if (grp->IsPA() || grp->IsAP()) runType = "dpmjet";
    else if (grp->IsAA())                runType = "hijing";
  }
};

/** 
 * Configure the simulation backend 
 * 
 */
void Config()
{
  // --- Get settings from environment variables --------------------
  Setup    s;
  gROOT->Macro("EGConfig.C");
  

  // ---- Seed random number generator -------------------------------
  gRandom->SetSeed(s.seed);
  std::cerr << "Seed for random number generation= " << s.seed << std::endl; 

  //------------------------------------------------------------------
  // 
  // Geometry and tracking 
  //
  // --- Libraries required by geant321 ------------------------------
  VirtualEGCfg::LoadGen(s.runType);
  s.MakeBackend();
  // gSystem->Load("libgeant321");
  // new TGeant3TGeo("C++ Interface to Geant3");

  // -----------------------------------------------------------------
  //  Create the output file
  std::cout<< "Config.C: Creating Run Loader ..." << std::endl;
  AliRunLoader* rl = AliRunLoader::Open("galice.root",
					AliConfig::GetDefaultEventFolderName(),
					"recreate");
  if (!rl) Fatal("Config","Can not instatiate the Run Loader");

  rl->SetCompressionLevel(2);
  rl->SetNumberOfEventsPerFile(1000);
  gAlice->SetRunLoader(rl);

  //
  //=======================================================================
  // Steering parameters for ALICE simulation
  // 
  // --- Specify event type to be tracked through the ALICE setup
  // --- All positions are in cm, angles in degrees, and P and E in GeV

  // --- Process switches --------------------------------------------
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
  gMC->SetProcess("RAYL",1);
  

  // --- Tracking cuts -----------------------------------------------
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
  
  // --- Set External decayer ----------------------------------------
  TVirtualMCDecayer* decayer = egCfg->MakeDecayer(s.runType);
  if (decayer) gMC->SetExternalDecayer(decayer);  

  //------------------------------------------------------------------
  // 
  // Generator Configuration 
  //
  // --- Make the generator - this loads libraries
  AliGenerator* gener = egCfg->MakeGenerator(s.runType,
					     s.minB,
					     s.maxB);
  if (!egCfg->IsLego())  {   
    gener->Init();
    if (gener->IsA()->InheritsFrom("AliGenHijing")) {
      Info("", "Setting Hijing debug");
      // static_cast<AliGenHijing*>(gener->GetTHijing()->SetIHPR2(10,1));
    }
  }
      

  // --- Go back to galice.root --------------------------------------
  rl->CdGAFile();
  
  //=================== Alice BODY parameters =============================
  AliBODY *BODY = new AliBODY("BODY", "Alice envelop");
  
  
  if (detCfg->UseMAG())   new AliMAG("MAG", "Magnet");
  if (detCfg->UseABSO())  new AliABSOv3("ABSO", "Muon Absorber");
  if (detCfg->UseDIPO())  new AliDIPOv3("DIPO", "Dipole version 3");
  if (detCfg->UseHALL())  new AliHALLv3("HALL", "Alice Hall");
  if (detCfg->UseFRAME()) (new AliFRAMEv2("FRAME", "Space Frame"))->SetHoles(1);
  if (detCfg->UseSHIL())  new AliSHILv3("SHIL", "Shielding Version 3");
  if (detCfg->UsePIPE())  new AliPIPEv3("PIPE", "Beam Pipe");
  if (detCfg->UseITS())   new AliITSv11("ITS","ITS v11");
  // if (detCfg->UseITS())   new AliITSv11Hybrid("ITS","ITS v11Hybrid");
  if (detCfg->UseTPC()) {
    AliTPC* tpc = new AliTPCv2("TPC", "Default");
    if (s.IsGeant4()) tpc->SetPrimaryIonisation();
  }
  if (detCfg->UseTOF())   new AliTOFv6T0("TOF", "normal TOF");
  if (detCfg->UseHMPID()) new AliHMPIDv3("HMPID", "normal HMPID");
  if (detCfg->UseZDC()) {
    AliZDC *ZDC = 0;
    if (grp->period.EqualTo("LHC10h")) {
      // Need to use older ZDC for PbPb 
      ZDC = new AliZDCv3("ZDC", "normal ZDC");
      ZDC->SetSpectatorsTrack();
    }
    else 
      ZDC = new AliZDCv4("ZDC", "normal ZDC");
    if (grp->Year() < 2011) { //?
      // What are these? Do they need to be set properly? 
      //Collimators aperture
      ZDC->SetVCollSideCAperture(0.85);
      ZDC->SetVCollSideCCentre(0.);
      ZDC->SetVCollSideAAperture(0.75);
      ZDC->SetVCollSideACentre(0.);
      //Detector position
      ZDC->SetYZNC(1.6);
      ZDC->SetYZNA(1.6);
      ZDC->SetYZPC(1.6);
      ZDC->SetYZPA(1.6);
    }
    ZDC->SetLumiLength(0.);
    if (grp->IsPA() || grp->IsAP()) { 
      ZDC->SetpAsystem();
      ZDC->SetBeamEnergy(82.*grp->beamEnergy/208.);
    }
  }
  if (detCfg->UseTRD()) {
    AliTRD *TRD = new AliTRDv1("TRD", "TRD slow simulator");
    AliTRDgeometry *geoTRD = TRD->GetGeometry();
    // Total of 18 super modules. We turn them all off by default 
    for (Int_t i = 0; i < 18; i++) geoTRD->SetSMstatus(i, 0);

    // '09-'10 had 7 super modules 
    geoTRD->SetSMstatus( 0,1);
    geoTRD->SetSMstatus( 1,1);
    geoTRD->SetSMstatus( 7,1);
    geoTRD->SetSMstatus( 8,1);
    geoTRD->SetSMstatus( 9,1);
    geoTRD->SetSMstatus(14,1);//?
    geoTRD->SetSMstatus(17,1);

    // In jan '11 3 more were added 
    if (grp->Year() > 2010) { 
      geoTRD->SetSMstatus(11, 1);
      geoTRD->SetSMstatus(15, 1);
      geoTRD->SetSMstatus(16, 1);//?
    }

    // In the 2012 shutdow 3 more were added 
    if (grp->Year() > 2012) { 
      geoTRD->SetSMstatus( 2,1);
      geoTRD->SetSMstatus( 3,1);
      geoTRD->SetSMstatus( 6,1);
    }
    if (grp->Year() > 2014) {
      geoTRD->SetSMstatus( 4,1);
      geoTRD->SetSMstatus( 5,1);
      geoTRD->SetSMstatus(10,1);
      geoTRD->SetSMstatus(12,1);
      geoTRD->SetSMstatus(13,1);
    }      
  }
  if (detCfg->UseFMD())    new AliFMDv1("FMD", "normal FMD");
  if (detCfg->UseMUON()) {
    AliMUON *MUON = new AliMUONv1("MUON", "default");
    MUON->SetTriggerEffCells(1); // not needed if raw masks
    MUON->SetTriggerResponseV1(2);
  }
  if (detCfg->UsePHOS())   new AliPHOSv1("PHOS", "noCPV_Modules123");
  if (detCfg->UsePMD())    new AliPMDv1("PMD", "normal PMD");
  if (detCfg->UseT0())     new AliT0v1("T0", "T0 Detector");
  if (detCfg->UseEMCAL())  {
    TString var;
    if       (grp->run <= 140000) var="EMCAL_FIRSTYEARV1";
    else if  (grp->run <= 170593) var="COMPLETEV1";
    else if  (grp->run <= 197692) var="EMCAL_COMPLETE12SMV1";
    else                          var="EMCAL_COMPLETE12SMV1_DCAL_8SM";
    new AliEMCALv2("EMCAL", var.Data());
  }
  if (detCfg->UseACORDE()) new AliACORDEv1("ACORDE", "normal ACORDE");
  if (detCfg->UseVZERO())  new AliVZEROv7("VZERO", "normal VZERO");
}




// 
// EOF
// 
