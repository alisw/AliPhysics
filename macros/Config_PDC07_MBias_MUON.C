// Config file MUON + ITS (for vertex) + VZERO + FMD for PDC07
// Tuned for p+p min biais and quarkonia production, heavy quark production 
// from AliGenCorrHF(AliGenMUONCocktailpp)
// Remember to define the directory and option
// gAlice->SetConfigFunction("Config('$HOME','box');");
// april 3rd: added L3 magnet 

void Config(char directory[100]="", char option[6]="trgAll")
{
    
    static Int_t sseed = 0; // Set 0 to use the current time

  // Load Pythia libraries
    LoadPythia();
    
  //=====================================================================
  //  Libraries required by geant321
    gSystem->Load("libgeant321.so");
    
    new TGeant3TGeo("C++ Interface to Geant3");
  //=======================================================================

    if(!AliCDBManager::Instance()->IsDefaultStorageSet()){
	AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT");
	AliCDBManager::Instance()->SetRun(0);
    }

  // Set Random Number seed
    gRandom->SetSeed(sseed);

  //  Create the output file    
    Text_t filename[100];
    sprintf(filename,"%sgalice.root",directory);

    AliRunLoader* rl=0x0;
    rl = AliRunLoader::Open(
	filename, AliConfig::GetDefaultEventFolderName(), "recreate");
    if (rl == 0x0) {
	gAlice->Fatal("Config.C","Can not instatiate the Run Loader");
	return;
    }
    rl->SetCompressionLevel(2);
    rl->SetNumberOfEventsPerFile(1000);
    gAlice->SetRunLoader(rl);

  //=======================================================================
  // Set the trigger configuration - not needed here for MUON in
  // stand-alone mode
  // gAlice->SetTriggerDescriptor("p-p");
  // cout<<"Trigger configuration is set to  p-p"<<endl;

  //=======================================================================
  // Set External decayer
    TVirtualMCDecayer *decayer = new AliDecayerPythia();
    decayer->SetForceDecay(kAll);
    decayer->Init();
    gMC->SetExternalDecayer(decayer);

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

    AliGenPythia* PythiaForMUONCocktail(Decay_t dt)
	{
	    AliGenPythia *pythia = new AliGenPythia(1);
	    pythia->SetProcess(kPyMbMSEL1);
	    pythia->SetStrucFunc(kCTEQ5L);
	    pythia->SetEnergyCMS(14000.);
	    pythia->SetForceDecay(dt);
	    pythia->SetPtRange(0.,100.);
	    pythia->SetYRange(-8.,8.);
	    pythia->SetPhiRange(0.,360.);
	    pythia->SetPtHard(2.76,-1.0);
	    pythia->SwitchHFOff();
	    return pythia;
	}

    if (!strcmp(option,"trg2mu")) {
	AliGenMUONCocktailpp *gener = new AliGenMUONCocktailpp();
	gener->SetPtRange(0.,100.);
	gener->SetYRange(-4.,-2.4);
	gener->SetPhiRange(0., 360.);
	gener->SetMuonMultiplicity(2);  
	gener->SetMuonPtCut(0.5);
	gener->SetMuonThetaRange(171.,178.);      
	gener->SetOrigin(0.,0.,0.); 
	gener->SetSigma(0.,0.,5.);
	gener->SetVertexSmear(kPerEvent);
	Decay_t dt = gener->GetDecayModePythia();
	AliGenPythia* pythia = PythiaForMUONCocktail(dt);
	pythia->Init();     
	gener->AddGenerator(pythia,"Pythia",1);
	gener->Init(); 
    }
 
    if (!strcmp(option,"trg1mu")) {
	AliGenMUONCocktailpp *gener = new AliGenMUONCocktailpp();
	gener->SetPtRange(0.,100.);
	gener->SetYRange(-4.,-2.4);
	gener->SetPhiRange(0., 360.);
	gener->SetMuonMultiplicity(1);  
	gener->SetMuonPtCut(0.5);
	gener->SetMuonThetaRange(171.,178.);      
	gener->SetOrigin(0.,0.,0.); 
	gener->SetSigma(0.,0.,5.);
	gener->SetVertexSmear(kPerEvent);
	Decay_t dt = gener->GetDecayModePythia();
	AliGenPythia* pythia = PythiaForMUONCocktail(dt);
	pythia->Init();     
	gener->AddGenerator(pythia,"Pythia",1);
	gener->Init(); 
    }
    if (!strcmp(option,"trgAll")) {
	AliGenMUONCocktailpp *gener = new AliGenMUONCocktailpp();
	gener->SetPtRange(0.,100.);
	gener->SetYRange(-4.,-2.4);
	gener->SetPhiRange(0., 360.);
	gener->SetMuonMultiplicity(0);  
	gener->SetMuonPtCut(0.);    
	gener->SetOrigin(0.,0.,0.); 
	gener->SetSigma(0.,0.,5.3);
	gener->SetVertexSmear(kPerEvent);
	Decay_t dt = gener->GetDecayModePythia();
	AliGenPythia* pythia = PythiaForMUONCocktail(dt);
	pythia->Init();     
	gener->AddGenerator(pythia,"Pythia",1);
	gener->Init();
    }
  //============================================================= 
  // Field (L3 0.5 T) outside dimuon spectrometer
    AliMagF* field = new AliMagF("Maps","Maps", 2, 1., 10., AliMagF::k5kG);
    field->SetL3ConstField(0); // Using const. field in the barrel 
    TGeoGlobalMagField::Instance()->SetField(field);

    Int_t   iITS = 1;
    Int_t   iFMD = 1;
    Int_t   iVZERO = 1;

    rl->CdGAFile();

  //=================== Alice BODY parameters =============================
    AliBODY *BODY = new AliBODY("BODY","Alice envelop");
  //=================== ABSO parameters ============================
    AliABSO *ABSO = new AliABSOv3("ABSO", "Muon Absorber");
  //=================== DIPO parameters ============================
    AliDIPO *DIPO = new AliDIPOv3("DIPO", "Dipole version 3");
  //================== HALL parameters ============================
    AliHALL *HALL = new AliHALLv3("HALL", "Alice Hall");
  //================== The L3 Magnet ==============================
    AliMAG *MAG = new AliMAG("MAG", "L3 Magnet");
  //=================== PIPE parameters ============================
    AliPIPE *PIPE = new AliPIPEv3("PIPE", "Beam Pipe");
  //=================== SHIL parameters ============================
    AliSHIL *SHIL = new AliSHILv3("SHIL", "Shielding Version 3");
  //=================== ITS parameters =============================
    if(iITS) {
	AliITS *ITS  = new AliITSv11Hybrid("ITS","ITS v11Hybrid");
    }      
 //=================== FMD parameters =============================
    if(iFMD) {
	AliFMD *FMD = new AliFMDv1("FMD", "normal FMD");
    }
 //=================== VZERO parameters =============================
    if (iVZERO) {
	AliVZERO *VZERO = new AliVZEROv7("VZERO", "normal VZERO");
    }

  //=================== MUON Subsystem ===========================
    cout << ">>> Config.C: Creating AliMUONv1 ..."<<endl;

    AliMUON *MUON = new AliMUONv1("MUON");
}

Float_t EtaToTheta(Float_t arg){
    return (180./TMath::Pi())*2.*atan(exp(-arg));
}

  

void LoadPythia()
{
    // Load Pythia related libraries
    gSystem->Load("liblhapdf.so");      // Parton density functions
    gSystem->Load("libEGPythia6.so");   // TGenerator interface
    gSystem->Load("libpythia6.so");     // Pythia
    gSystem->Load("libAliPythia6.so");  // ALICE specific implementations
}
