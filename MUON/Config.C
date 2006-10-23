// Config file test for MUON spectormeter
// Remember to define the directory and option
// gAlice->SetConfigFunction("Config('$HOME','box');");

void Config(char directory[100]="", char option[6]="param")
{
  //
  // Config file for MUON test
  //

  //=====================================================================
  //  Libraries required by geant321
  gSystem->Load("libgeant321.so");
  new TGeant3TGeo("C++ Interface to Geant3");
  //=======================================================================
  //  Create the output file    
  Text_t filename[100];
  sprintf(filename,"%sgalice.root",directory);
  cout << ">>> Output file is " << filename << endl;   
  cout << ">>> Config.C: Creating Run Loader ..."<<endl;
  AliRunLoader* rl=0x0;
  rl = AliRunLoader::Open(
	filename, AliConfig::GetDefaultEventFolderName(), "recreate");
  if (rl == 0x0) {
    gAlice->Fatal("Config.C","Can not instatiate the Run Loader");
    return;
  }
  rl->SetCompressionLevel(2);
  rl->SetNumberOfEventsPerFile(100);
  gAlice->SetRunLoader(rl);

  //AliLog::SetModuleDebugLevel("MUON", 1);
  
  //=======================================================================
  // Set External decayer
  TVirtualMCDecayer *decayer = new AliDecayerPythia();
  decayer->SetForceDecay(kAll);
  decayer->Init();
  gMC->SetExternalDecayer(decayer);

  //
  //=======================================================================
  // ******* GEANT STEERING parameters FOR ALICE SIMULATION *******
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
  //
  //=======================================================================
  // ************* STEERING parameters FOR ALICE SIMULATION **************
  // Chamber positions
  // From AliMUONConstants class we get :
  //   Position Z (along beam) of the chambers (in cm) 
  //        (from AliMUONConstants class):  
  //    533.5,  546.5,  678.5, 693.5,  964.0, 986.0, 1251.5, 1278.5, 
  //   1416.5, 1443.5,  1610, 1625.,  1710., 1725. 
  //   Internal Radius (in cm)   
  //     36.4,  46.2,  66.0,  80.,  80., 100., 100.    
  //   External Radius (in cm)
  //    183.,  245.,  395.,  560., 563., 850., 900.  
  //=======================================================================
  if (!strcmp(option,"box")) {
    AliGenBox * gener = new AliGenBox(1);
    gener->SetMomentumRange(20.,20.1);
    gener->SetPhiRange(0., 360.);         
    gener->SetThetaRange(171.000,178.001);
    gener->SetPart(13);           // Muons
    gener->SetOrigin(0.,0., 0.);  //vertex position
    gener->SetSigma(0.0, 0.0, 0.0);         //Sigma in (X,Y,Z) (cm) on IP position
  }
  if (!strcmp(option,"gun")) {
    //*********************************************
    // Example for Fixed Particle Gun             *
    //*********************************************
    AliGenFixed *gener = new AliGenFixed(ntracks);
    gener->SetMomentum(10);
    gener->SetPhiRange(0.);
    gener->SetThetaRange(0.);
    gener->SetOrigin(30,30,1200);//vertex position
    gener->SetPart(13);          //GEANT particle type  13 is muons
  }
  if (!strcmp(option,"scan")) {
    AliGenScan *gener = new AliGenScan(-1);
    gener->SetMomentumRange(10,10);
    gener->SetPhiRange(0, 0);
    gener->SetThetaRange(-180, -180);
    //vertex position
    //gener->SetSigma(1,1,0);           //Sigma in (X,Y,Z) (cm) on IP position
    gener->SetPart(kRootino); 
    gener->SetRange(100, -300., 300., 100, -300., 300., 1, 2000, 2000);
  }  
  if (!strcmp(option,"param")) {
    //*******************************************************
    // Example for J/psi or Upsilon Production from  Parameterisation *
    //*******************************************************
    AliGenParam *gener = new AliGenParam(1, AliGenMUONlib::kUpsilon);
    gener->SetMomentumRange(0,999);
    gener->SetPtRange(0,100.);
    gener->SetPhiRange(0., 360.);
    gener->SetCutOnChild(1);
    gener->SetChildPhiRange(0.,360.);
    gener->SetChildThetaRange(171.0,178.0);
    gener->SetOrigin(0,0,0);          //vertex position    gener->SetSigma(0,0,0);           //Sigma in (X,Y,Z) (cm) on IP position
    gener->SetForceDecay(kDiMuon);
    gener->SetTrackingFlag(1);
    gener->Init();
  }
  //============================================================= 
  // Field (L3 0.4 T)
  AliMagFMaps* field = new AliMagFMaps("Maps","Maps", 1, 1., 10., AliMagFMaps::k4kG);
  gAlice->SetField(field);

  //=================== Alice BODY parameters =============================
  AliBODY *BODY = new AliBODY("BODY","Alice envelop");
  //=================== ABSO parameters ============================
  AliABSO *ABSO = new AliABSOv0("ABSO", "Muon Absorber");
  //=================== DIPO parameters ============================
  AliDIPO *DIPO = new AliDIPOv2("DIPO", "Dipole version 2");
  //================== HALL parameters ============================
  AliHALL *HALL = new AliHALL("HALL", "Alice Hall");
  //=================== PIPE parameters ============================
  AliPIPE *PIPE = new AliPIPEv0("PIPE", "Beam Pipe");
  //=================== SHIL parameters ============================
  AliSHIL *SHIL = new AliSHILv2("SHIL", "Shielding Version 2");

  //=================== MUON Subsystem ===========================
  cout << ">>> Config.C: Creating AliMUONv1 ..."<<endl;

  AliMUON *MUON = new AliMUONv1("MUON", "default");

  // The 3 switches below are to be used for the trigger code
  // their default value is set in AliMUON.h
  // activate trigger cluster-size (0=default, 1=cluster-size according to AliMUONResponseTriggerV1
  //  MUON->SetTriggerResponseV1(0);
  // activate 4/4 trigger coincidence (0=default (coinc 3/4), 1=coinc 4/4)
  //  MUON->SetTriggerCoinc44(0);
  // activate trigger chamber efficiency by cells (0=default, 1=trigger efficiency according to AliMUONTriggerEfficiencyCells
  //  MUON->SetTriggerEffCells(0);

  // To get same as above w/o noise-only digits for the tracker do  :
  // MUON->SetDigitizerWithNoise(kKALSE);

  //
  // If SetAlign, the detection elements transformations
  // are taken from the input file and not from the code
  // MUON->SetAlign("transform.dat");

  // To generate and read scaler trigger events in rawdata
  // MUON->SetTriggerScalerEvent();

  // If you want to play with builders, first reset the geometry builder,
  // and then add yours.
  //  MUON->ResetGeometryBuilder();
  //  MUON->AddGeometryBuilder(new AliMUONSt1GeometryBuilderV2(MUON));
  //  MUON->AddGeometryBuilder(new AliMUONSt2GeometryBuilderV2(MUON));
  //  MUON->AddGeometryBuilder(new AliMUONSlatGeometryBuilder(MUON));
  //  MUON->AddGeometryBuilder(new AliMUONTriggerGeometryBuilder(MUON));
}

Float_t EtaToTheta(Float_t arg){
  return (180./TMath::Pi())*2.*atan(exp(-arg));
}
