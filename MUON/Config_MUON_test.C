// Config file test for MUON spectormeter
// Remember to define the directory and option
// gAlice->SetConfigFunction("Config('$HOME','box');");

void Config(char directory[100]="", char option[6]="box")
{
  //
  // Config file for MUON test
  // Gines MARITNEZ, Subatech, mai 2003, august 2003
  // 

  //=====================================================================
  //  Libraries required by geant321
  gSystem->Load("libgeant321.so");
  new TGeant3("C++ Interface to Geant3");
  //=======================================================================
  //  Create the output file    
  Text_t filename[100];
  sprintf(filename,"%sgalice.root",directory);
  cout << ">>> Output file is " << filename << endl;   
  cout << ">>> Config_MUON_test.C: Creating Run Loader ..."<<endl;
  AliRunLoader* rl=0x0;
  rl = AliRunLoader::Open(
	filename, AliConfig::fgkDefaultEventFolderName, "recreate");
  if (rl == 0x0) {
    gAlice->Fatal("Config_MUON_test.C","Can not instatiate the Run Loader");
    return;
  }
  rl->SetCompressionLevel(2);
  rl->SetNumberOfEventsPerFile(100);
  gAlice->SetRunLoader(rl);

  
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
    gener->SetMomentumRange(7.,7.1);
    gener->SetPhiRange(-180., 180.);         
    gener->SetThetaRange(2.000,9.000);
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
  if (!strcmp(option,"param")) {
    //*******************************************************
    // Example for J/psi or Upsilon Production from  Parameterisation *
    //*******************************************************
    AliGenParam *gener = new AliGenParam(1, AliGenMUONlib::kUpsilon);
    gener->SetMomentumRange(0,999);
    gener->SetPtRange(0,100.);
    gener->SetPhiRange(-180, 180);
    gener->SetYRange(2.5,4);
    gener->SetCutOnChild(1);
    gener->SetChildThetaRange(2.0,9);
    gener->SetOrigin(0,0,0);          //vertex position    gener->SetSigma(0,0,0);           //Sigma in (X,Y,Z) (cm) on IP position
    gener->SetForceDecay(kDiMuon);
    gener->SetTrackingFlag(1);
  }
     
  //=============================================================
  //Specify maximum magnetic field in Tesla (neg. ==> default field)   
  AliMagFMaps* field = new AliMagFMaps("Maps","Maps", 2, 1., 10., 1);
  gAlice->SetField(field); 

  //=================== Alice BODY parameters =============================
  AliBODY *BODY = new AliBODY("BODY","Alice envelop");


  //=================== MUON Subsystem ===========================
  AliMUONv1 *MUON  = new AliMUONv1("MUON","default");
  
}


Float_t EtaToTheta(Float_t arg){
  return (180./TMath::Pi())*2.*atan(exp(-arg));
}
