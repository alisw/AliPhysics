// $Id$
//
// AliRoot Configuration for running aliroot with Monte Carlo.
// Called from g4Config.C

Float_t EtaToTheta(Float_t arg);
static Int_t    eventsPerRun = 50;

void ConfigCommon()
{
  // ============================= 
  // Root file
  // ============================= 

  // Create the output file
  AliRunLoader* rl = 0;
  rl = AliRunLoader::Open("galice.root",
			   AliConfig::GetDefaultEventFolderName(),
			   "recreate");
  if (!rl) {
    gAlice->Fatal("Config.C","Can not instatiate the Run Loader");
    return;
  }
  rl->SetCompressionLevel(2);
  rl->SetNumberOfEventsPerFile(3);
  gAlice->SetRunLoader(rl);
//  gAlice->SetRootGeometry();
//  gAlice->SetGeometryFileName("geometry.root");

  // Set the trigger configuration
  gAlice->SetTriggerDescriptor("Pb-Pb");
  cout<<"Trigger configuration is set to  Pb-Pb"<<endl;

    // Set Random Number seed
    gRandom->SetSeed(123456); // Set 0 to use the currecnt time
    AliLog::Message(AliLog::kInfo, Form("Seed for random number generation = %d",gRandom->GetSeed()), "Config.C", "Config.C", "Config()","Config.C", __LINE__);
    int     nParticles = 500;
    if (gSystem->Getenv("CONFIG_NPARTICLES"))
    {
        nParticles = atoi(gSystem->Getenv("CONFIG_NPARTICLES"));
    }

///*
    AliGenCocktail *gener = new AliGenCocktail();
    gener->SetPhiRange(0, 360);
    // Set pseudorapidity range from -8 to 8.
    Float_t thmin = EtaToTheta(8);   // theta min. <---> eta max
    Float_t thmax = EtaToTheta(-8);  // theta max. <---> eta min 
    gener->SetThetaRange(thmin,thmax);
    gener->SetOrigin(0, 0, 0);  //vertex position
    gener->SetSigma(0, 0, 0);   //Sigma in (X,Y,Z) (cm) on IP position

    AliGenHIJINGpara *hijingparam = new AliGenHIJINGpara(nParticles);
    hijingparam->SetMomentumRange(0.2, 999);
    gener->AddGenerator(hijingparam,"HIJING PARAM",1);
    gener->Init();

    // 
    // Activate this line if you want the vertex smearing to happen
    // track by track
    //
    //gener->SetVertexSmear(perTrack); 

  // ============================= 
  // Magnetic field
  // ============================= 

    // Field (L3 0.4 T)
    AliMagFMaps* field = new AliMagFMaps("Maps","Maps", 2, 1., 10., 1);
    gAlice->SetField(field);    

    Int_t   iABSO  =  1;
    Int_t   iDIPO  =  1;
    Int_t   iFMD   =  1;
    Int_t   iFRAME =  1;
    Int_t   iHALL  =  1;
    Int_t   iITS   =  1;
    Int_t   iMAG   =  1;
    Int_t   iMUON  =  1;
    Int_t   iPHOS  =  1;
    Int_t   iPIPE  =  1;
    Int_t   iPMD   =  1;
    Int_t   iHMPID =  1;
    Int_t   iSHIL  =  1;
    Int_t   iT0    =  1;
    Int_t   iTOF   =  1;
    Int_t   iTPC   =  1;
    Int_t   iTRD   =  1;
    Int_t   iZDC   =  1;
    Int_t   iEMCAL =  1;
    Int_t   iACORDE = 0;
    Int_t   iVZERO =  1;

    rl->CdGAFile();
    //=================== Alice BODY parameters =============================
    AliBODY *BODY = new AliBODY("BODY", "Alice envelop");

    if (iMAG)
    {
        //=================== MAG parameters ============================
        // --- Start with Magnet since detector layouts may be depending ---
        // --- on the selected Magnet dimensions ---
        AliMAG *MAG = new AliMAG("MAG", "Magnet");
    }


    if (iABSO)
    {
        //=================== ABSO parameters ============================
        AliABSO *ABSO = new AliABSOv3("ABSO", "Muon Absorber");
    }

    if (iDIPO)
    {
        //=================== DIPO parameters ============================

        AliDIPO *DIPO = new AliDIPOv3("DIPO", "Dipole version 3");
    }

    if (iHALL)
    {
        //=================== HALL parameters ============================

        AliHALL *HALL = new AliHALLv3("HALL", "Alice Hall");
    }


    if (iFRAME)
    {
        //=================== FRAME parameters ============================

        AliFRAMEv2 *FRAME = new AliFRAMEv2("FRAME", "Space Frame");
    }

    if (iSHIL)
    {
        //=================== SHIL parameters ============================

        AliSHIL *SHIL = new AliSHILv3("SHIL", "Shielding Version 3");
    }


    if (iPIPE)
    {
        //=================== PIPE parameters ============================

        AliPIPE *PIPE = new AliPIPEv3("PIPE", "Beam Pipe");
    }
 
    if(iITS) {

    //=================== ITS parameters ============================
    //
    // As the innermost detector in ALICE, the Inner Tracking System "impacts" on
    // almost all other detectors. This involves the fact that the ITS geometry
    // still has several options to be followed in parallel in order to determine
    // the best set-up which minimizes the induced background. All the geometries
    // available to date are described in the following. Read carefully the comments
    // and use the default version (the only one uncommented) unless you are making
    // comparisons and you know what you are doing. In this case just uncomment the
    // ITS geometry you want to use and run Aliroot.
    //
    // Detailed geometries:         
    //
    //
    //AliITS *ITS  = new AliITSv5symm("ITS","Updated ITS TDR detailed version with symmetric services");
    //
    //AliITS *ITS  = new AliITSv5asymm("ITS","Updates ITS TDR detailed version with asymmetric services");
    //
	AliITSvPPRasymmFMD *ITS  = new AliITSvPPRasymmFMD("ITS","New ITS PPR detailed version with asymmetric services");
	ITS->SetMinorVersion(2);  // don't touch this parameter if you're not an ITS developer
	ITS->SetReadDet(kFALSE);	  // don't touch this parameter if you're not an ITS developer
    //    ITS->SetWriteDet("$ALICE_ROOT/ITS/ITSgeometry_vPPRasymm2.det");  // don't touch this parameter if you're not an ITS developer
	ITS->SetThicknessDet1(200.);   // detector thickness on layer 1 must be in the range [100,300]
	ITS->SetThicknessDet2(200.);   // detector thickness on layer 2 must be in the range [100,300]
	ITS->SetThicknessChip1(150.);  // chip thickness on layer 1 must be in the range [150,300]
	ITS->SetThicknessChip2(150.);  // chip thickness on layer 2 must be in the range [150,300]
	ITS->SetRails(0);	     // 1 --> rails in ; 0 --> rails out
	ITS->SetCoolingFluid(1);   // 1 --> water ; 0 --> freon

    // Coarse geometries (warning: no hits are produced with these coarse geometries and they unuseful 
    // for reconstruction !):
    //                                                     
    //
    //AliITSvPPRcoarseasymm *ITS  = new AliITSvPPRcoarseasymm("ITS","New ITS PPR coarse version with asymmetric services");
    //ITS->SetRails(0);                // 1 --> rails in ; 0 --> rails out
    //ITS->SetSupportMaterial(0);      // 0 --> Copper ; 1 --> Aluminum ; 2 --> Carbon
    //
    //AliITS *ITS  = new AliITSvPPRcoarsesymm("ITS","New ITS PPR coarse version with symmetric services");
    //ITS->SetRails(0);                // 1 --> rails in ; 0 --> rails out
    //ITS->SetSupportMaterial(0);      // 0 --> Copper ; 1 --> Aluminum ; 2 --> Carbon
    //                      
    //
    //
    // Geant3 <-> EUCLID conversion
    // ============================
    //
    // SetEUCLID is a flag to output (=1) or not to output (=0) both geometry and
    // media to two ASCII files (called by default ITSgeometry.euc and
    // ITSgeometry.tme) in a format understandable to the CAD system EUCLID.
    // The default (=0) means that you dont want to use this facility.
    //
	ITS->SetEUCLID(0);  
    }

    if (iTPC)
    {
        //============================ TPC parameters ===================
        AliTPC *TPC = new AliTPCv2("TPC", "Default");
    }


    if (iTOF) {
        //=================== TOF parameters ============================
	AliTOF *TOF = new AliTOFv6T0("TOF", "normal TOF");
	// Partial geometry: modules at 2,3,4,6,7,11,12,14,15,16
	// starting at 6h in positive direction
	//	Int_t TOFSectors[18]={-1,-1,0,0,0,-1,0,0,-1,-1,-1,0,0,-1,0,0,0,0};
	// Partial geometry: modules at 1,2,6,7,9,10,11,12,15,16,17
	// (ALICE numbering convention)
       	Int_t TOFSectors[18]={-1,0,0,-1,-1,-1,0,0,-1,0,0,0,0,-1,-1,0,0,0};
	TOF->SetTOFSectors(TOFSectors);
    }


    if (iHMPID)
    {
        //=================== HMPID parameters ===========================
        AliHMPID *HMPID = new AliHMPIDv2("HMPID", "normal HMPID");

    }


    if (iZDC)
    {
        //=================== ZDC parameters ============================

        AliZDC *ZDC = new AliZDCv2("ZDC", "normal ZDC");
    }

    if (iTRD)
    {
        //=================== TRD parameters ============================

        AliTRD *TRD = new AliTRDv1("TRD", "TRD slow simulator");
    }

    if (iFMD)
    {
        //=================== FMD parameters ============================
	AliFMD *FMD = new AliFMDv1("FMD", "normal FMD");
   }

    if (iMUON)
    {
        //=================== MUON parameters ===========================
        // New MUONv1 version (geometry defined via builders)
        AliMUON *MUON = new AliMUONv1("MUON", "default");
    }
    //=================== PHOS parameters ===========================

    if (iPHOS)
    {
        AliPHOS *PHOS = new AliPHOSv1("PHOS", "IHEP");
    }


    if (iPMD)
    {
        //=================== PMD parameters ============================
        AliPMD *PMD = new AliPMDv1("PMD", "normal PMD");
    }

    if (iT0)
    {
        //=================== T0 parameters ============================
        AliT0 *T0 = new AliT0v1("T0", "T0 Detector");
    }

    if (iEMCAL)
    {
        //=================== EMCAL parameters ============================
        AliEMCAL *EMCAL = new AliEMCALv2("EMCAL", "SHISH_77_TRD1_2X2_FINAL_110DEG");
    }

     if (iACORDE)
    {
        //=================== ACORDE parameters ============================
        AliACORDE *ACORDE = new AliACORDEv0("ACORDE", "normal ACORDE");
    }

     if (iVZERO)
    {
        //=================== ACORDE parameters ============================
        AliVZERO *VZERO = new AliVZEROv7("VZERO", "normal VZERO");
    }

     AliLog::Message(AliLog::kInfo, "End of Config", "Config.C", "Config.C", "Config()"," Config.C", __LINE__);

}

Float_t EtaToTheta(Float_t arg){
  return (180./TMath::Pi())*2.*atan(exp(-arg));
}
