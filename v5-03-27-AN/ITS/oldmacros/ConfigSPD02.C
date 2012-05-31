#include <TPDGCode.h>

void Config(){
    // Set Random Number seed
    // gRandom->SetSeed(12345);
    // libraries required by geant321
    gSystem->Load("libgeant321");
    new TGeant3("C++ Interface to Geant3");
    if (!gSystem->Getenv("CONFIG_FILE")){
        cout<<"Config.C: Creating Run Loader ..."<<endl;
        AliRunLoader *rl = AliRunLoader::Open("galice.root",
                      AliConfig::fgkDefaultEventFolderName,"recreate");
        if (rl == 0x0){
	    gAlice->Fatal("Config.C","Can not instatiate the Run Loader");
	    return;
	} // end if rl==0x0
        rl->SetCompressionLevel(2);
        rl->SetNumberOfEventsPerFile(1000);
        gAlice->SetRunLoader(rl);
    } // end if !gSystem
    TGeant3 *geant3 = (TGeant3 *) gMC;
    // Set External decayer
    AliDecayer *decayer = new AliDecayerPythia();
    decayer->SetForceDecay(kAll);
    decayer->Init();
    gMC->SetExternalDecayer(decayer);
    //=======================================================================
    // ******* GEANT STEERING parameters FOR ALICE SIMULATION *******
    geant3->SetTRIG(1);         //Number of events to be processed 
    geant3->SetSWIT(4, 10);
    geant3->SetDEBU(0, 0, 1);
    //geant3->SetSWIT(2,2);
    geant3->SetDCAY(1);
    geant3->SetPAIR(1);
    geant3->SetCOMP(1);
    geant3->SetPHOT(1);
    geant3->SetPFIS(0);
    geant3->SetDRAY(0);
    geant3->SetANNI(1);
    geant3->SetBREM(1);
    geant3->SetMUNU(1);
    geant3->SetCKOV(1);
    geant3->SetHADR(1);//Select pure GEANH (HADR 1) or GEANH/NUCRIN (HADR 3)
    geant3->SetLOSS(2);
    geant3->SetMULS(1);
    geant3->SetRAYL(1);
    geant3->SetAUTO(1);//Select automatic STMIN etc... calc. (AUTO 1) or manual (AUTO 0)
    geant3->SetABAN(0);//Restore 3.16 behaviour for abandoned tracks
    geant3->SetOPTI(2);//Select optimisation level for GEANT geometry searches (0,1,2)
    geant3->SetERAN(5.e-7);
    Float_t cut = 1.e-3;        // 1MeV cut by default
    Float_t tofmax = 1.e10;
    //             GAM ELEC NHAD CHAD MUON EBREM MUHAB EDEL MUDEL MUPA TOFMAX
    geant3->SetCUTS(cut, cut, cut, cut, cut, cut, cut, cut, cut, cut,
                    tofmax);
    //=======================================================================
    // ************* STEERING parameters FOR ALICE SIMULATION **************
    // --- Specify event type to be tracked through the ALICE setup
    // --- All positions are in cm, angles in degrees, and P and E in GeV
    if (gSystem->Getenv("CONFIG_NPARTICLES")){
        int     nParticles = atoi(gSystem->Getenv("CONFIG_NPARTICLES"));
    }else{
        int     nParticles = 1;
    } // end if
    //*********************************************
    // Example for Moving Particle Gun            *
    //*********************************************
    AliGenBox *gener = new AliGenBox(nParticles);
    gener->SetMomentumRange(100.,300.);
    gener->SetPhiRange(0,0.1);
    gener->SetThetaRange(0.0, .1);
    gener->SetOrigin(0.,0.,-50.);
    //vertex position
    gener->SetSigma(0.1,0.1,0.0); //Sigma in (X,Y,Z) (cm) on IP position
    gener->SetPart(kPiPlus);
    gener->Init();
    // Activate this line if you want the vertex smearing to happen
    // track by track
    //
    //gener->SetVertexSmear(perTrack); 
    // Field (L3 0.4 T)
    rootfile->cd();
    //gAlice->SetField(field);

    Int_t   iHALL  =  0;
    Int_t   iITS   =  1;
    rl->CdGAFile();
    //=================== Alice BODY parameters =============================
    AliBODY *BODY = new AliBODY("BODY", "Alice envelop");

    if (iHALL){
        //=================== HALL parameters ============================
        AliHALL *HALL = new AliHALL("HALL", "Alice Hall");
    } // end if
    if(iITS) {
	//=================== ITS parameters ============================
	AliITSvSPD02 *ITS  = new AliITSvSPD02("SPD test beam 2002");
    }
    return;
}

Float_t EtaToTheta(Float_t arg){
  return (180./TMath::Pi())*2.*atan(exp(-arg));
}
