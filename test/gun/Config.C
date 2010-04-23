// One can use the configuration macro in compiled mode by
// root [0] gSystem->Load("libgeant321");
// root [0] gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include\
//                   -I$ALICE_ROOT -I$ALICE/geant3/TGeant3");
// root [0] .x grun.C(1,"Config.C++")

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TPDGCode.h>
#include <TRandom.h>
#include <TSystem.h>
#include <TVirtualMC.h>
#include <TGeant3TGeo.h>
#include "STEER/AliRunLoader.h"
#include "STEER/AliRun.h"
#include "STEER/AliConfig.h"
#include "PYTHIA6/AliDecayerPythia.h"
#include "EVGEN/AliGenCocktail.h"
#include "EVGEN/AliGenHIJINGpara.h"
#include "EVGEN/AliGenFixed.h"
#include "EVGEN/AliGenBox.h"
#include "STEER/AliMagF.h"
#include "STEER/AliPID.h"
#include "STRUCT/AliBODY.h"
#include "STRUCT/AliMAG.h"
#include "STRUCT/AliABSOv3.h"
#include "STRUCT/AliDIPOv3.h"
#include "STRUCT/AliHALLv3.h"
#include "STRUCT/AliFRAMEv2.h"
#include "STRUCT/AliSHILv3.h"
#include "STRUCT/AliPIPEv3.h"
#include "ITS/AliITSv11Hybrid.h"
#include "TPC/AliTPCv2.h"
#include "TOF/AliTOFv6T0.h"
#include "HMPID/AliHMPIDv3.h"
#include "ZDC/AliZDCv3.h"
#include "TRD/AliTRDv1.h"
#include "TRD/AliTRDgeometry.h"
#include "FMD/AliFMDv1.h"
#include "MUON/AliMUONv1.h"
#include "PHOS/AliPHOSv1.h"
#include "PMD/AliPMDv1.h"
#include "T0/AliT0v1.h"
#include "EMCAL/AliEMCALv2.h"
#include "ACORDE/AliACORDEv1.h"
#include "VZERO/AliVZEROv7.h"
#endif

enum PprTrigConf_t
{
    kDefaultPPTrig, kDefaultPbPbTrig
};

const char * pprTrigConfName[] = {
    "p-p","Pb-Pb"
};

Float_t EtaToTheta(Float_t arg);

static PprTrigConf_t strig = kDefaultPPTrig;// default PP trigger configuration

void Config()
{
    // ThetaRange is (0., 180.). It was (0.28,179.72) 7/12/00 09:00
    // Theta range given through pseudorapidity limits 22/6/2001

    // Set Random Number seed
  gRandom->SetSeed(123456); // Set 0 to use the currecnt time


   // libraries required by geant321
#if defined(__CINT__)
    gSystem->Load("liblhapdf");
    gSystem->Load("libEGPythia6");
    gSystem->Load("libpythia6");
    gSystem->Load("libAliPythia6");
    gSystem->Load("libgeant321");
#endif

    new     TGeant3TGeo("C++ Interface to Geant3");

    AliRunLoader* rl=0x0;


    rl = AliRunLoader::Open("galice.root",
			    AliConfig::GetDefaultEventFolderName(),
			    "recreate");
    if (rl == 0x0)
      {
	gAlice->Fatal("Config.C","Can not instatiate the Run Loader");
	return;
      }
    rl->SetCompressionLevel(2);
    rl->SetNumberOfEventsPerFile(3);
    gAlice->SetRunLoader(rl);

    // Set the trigger configuration
    AliSimulation::Instance()->SetTriggerConfig(pprTrigConfName[strig]);
    cout<<"Trigger configuration is set to  "<<pprTrigConfName[strig]<<endl;

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

    // Special generation for Valgrind tests
    // Each detector is fired by few particles selected 
    // to cover specific cases


    // The cocktail itself

    AliGenCocktail *gener = new AliGenCocktail();
    gener->SetPhiRange(0, 360);
    // Set pseudorapidity range from -8 to 8.
    Float_t thmin = EtaToTheta(8);   // theta min. <---> eta max
    Float_t thmax = EtaToTheta(-8);  // theta max. <---> eta min 
    gener->SetThetaRange(thmin,thmax);
    gener->SetOrigin(0, 0, 0);  //vertex position
    gener->SetSigma(0, 0, 0);   //Sigma in (X,Y,Z) (cm) on IP position


    // Particle guns for the barrel part (taken from RichConfig)

    AliGenFixed *pG1=new AliGenFixed(1);
    pG1->SetPart(kProton);
    pG1->SetMomentum(2.5);
    pG1->SetTheta(109.5-3);
    pG1->SetPhi(10);
    gener->AddGenerator(pG1,"g1",1);
    
    AliGenFixed *pG2=new AliGenFixed(1);
    pG2->SetPart(kPiPlus);
    pG2->SetMomentum(1.0);
    pG2->SetTheta( 90.0-3);
    pG2->SetPhi(10);
    gener->AddGenerator(pG2,"g2",1);

    AliGenFixed *pG3=new AliGenFixed(1);
    pG3->SetPart(kPiMinus);
    pG3->SetMomentum(1.5);
    pG3->SetTheta(109.5-3);
    pG3->SetPhi(30);
    gener->AddGenerator(pG3,"g3",1);
    
    AliGenFixed *pG4=new AliGenFixed(1);
    pG4->SetPart(kKPlus);
    pG4->SetMomentum(0.7);
    pG4->SetTheta( 90.0-3);
    pG4->SetPhi(30);
    gener->AddGenerator(pG4,"g4",1);
    
    AliGenFixed *pG5=new AliGenFixed(1);
    pG5->SetPart(kKMinus);
    pG5->SetMomentum(1.0);
    pG5->SetTheta( 70.0-3);
    pG5->SetPhi(30);
    gener->AddGenerator(pG5,"g5",1);
    
    AliGenFixed *pG6=new AliGenFixed(1);
    pG6->SetPart(kProtonBar);
    pG6->SetMomentum(2.5);
    pG6->SetTheta( 90.0-3);
    pG6->SetPhi(50);
    gener->AddGenerator(pG6,"g6",1);
    
    AliGenFixed *pG7=new AliGenFixed(1);
    pG7->SetPart(kPiMinus);
    pG7->SetMomentum(0.7);
    pG7->SetTheta( 70.0-3);
    pG7->SetPhi(50);
    gener->AddGenerator(pG7,"g7",1);

    // Electrons for TRD

    AliGenFixed *pG8=new AliGenFixed(1);
    pG8->SetPart(kElectron);
    pG8->SetMomentum(1.2);
    pG8->SetTheta( 95.0);
    pG8->SetPhi(190);
    gener->AddGenerator(pG8,"g8",1);

    AliGenFixed *pG9=new AliGenFixed(1);
    pG9->SetPart(kPositron);
    pG9->SetMomentum(1.2);
    pG9->SetTheta( 85.0);
    pG9->SetPhi(190);
    gener->AddGenerator(pG9,"g9",1);

    // PHOS

    AliGenBox *gphos = new AliGenBox(1);
    gphos->SetMomentumRange(10,11.);
    gphos->SetPhiRange(270.5,270.7);
    gphos->SetThetaRange(90.5,90.7);
    gphos->SetPart(kGamma);
    gener->AddGenerator(gphos,"GENBOX GAMMA for PHOS",1);

    // EMCAL

    AliGenBox *gemcal = new AliGenBox(1);
    gemcal->SetMomentumRange(10,11.);
    gemcal->SetPhiRange(90.5,199.5);
    gemcal->SetThetaRange(90.5,90.7);
    gemcal->SetPart(kGamma);
    gener->AddGenerator(gemcal,"GENBOX GAMMA for EMCAL",1);

    // MUON
    AliGenBox * gmuon1 = new AliGenBox(1);
    gmuon1->SetMomentumRange(20.,20.1);
    gmuon1->SetPhiRange(0., 360.);         
    gmuon1->SetThetaRange(171.000,178.001);
    gmuon1->SetPart(kMuonMinus);           // Muons
    gener->AddGenerator(gmuon1,"GENBOX MUON1",1);

    AliGenBox * gmuon2 = new AliGenBox(1);
    gmuon2->SetMomentumRange(20.,20.1);
    gmuon2->SetPhiRange(0., 360.);         
    gmuon2->SetThetaRange(171.000,178.001);
    gmuon2->SetPart(kMuonPlus);           // Muons
    gener->AddGenerator(gmuon2,"GENBOX MUON1",1);

    //TOF
    AliGenFixed *gtof=new AliGenFixed(1);
    gtof->SetPart(kProton);
    gtof->SetMomentum(2.5);
    gtof->SetTheta(95);
    gtof->SetPhi(340);
    gener->AddGenerator(gtof,"Proton for TOF",1);

    //FMD1
    AliGenFixed *gfmd1=new AliGenFixed(1);
    gfmd1->SetPart(kGamma);
    gfmd1->SetMomentum(25);
    gfmd1->SetTheta(1.8);
    gfmd1->SetPhi(10);
    gener->AddGenerator(gfmd1,"Gamma for FMD1",1);
    
    //FMD2i
    AliGenFixed *gfmd2i=new AliGenFixed(1);
    gfmd2i->SetPart(kPiPlus);
    gfmd2i->SetMomentum(1.5);
    gfmd2i->SetTheta(7.3);
    gfmd2i->SetPhi(20);
    gener->AddGenerator(gfmd2i,"Pi+ for FMD2i",1);
    
    //FMD2o
    AliGenFixed *gfmd2o=new AliGenFixed(1);
    gfmd2o->SetPart(kPiMinus);
    gfmd2o->SetMomentum(1.5);
    gfmd2o->SetTheta(16.1);
    gfmd2o->SetPhi(30);
    gener->AddGenerator(gfmd2o,"Pi- for FMD2o",1);
    
    //FMD3o
    AliGenFixed *gfmd3o=new AliGenFixed(1);
    gfmd3o->SetPart(kPiPlus);
    gfmd3o->SetMomentum(1.5);
    gfmd3o->SetTheta(163.9);
    gfmd3o->SetPhi(40);
    gener->AddGenerator(gfmd3o,"Pi+ for FMD3o",1);
    
    //FMD3i
    AliGenFixed *gfmd3i=new AliGenFixed(1);
    gfmd3i->SetPart(kPiMinus);
    gfmd3i->SetMomentum(1.5);
    gfmd3i->SetTheta(170.5);
    gfmd3i->SetPhi(50);
    gener->AddGenerator(gfmd3i,"Pi- for FMD3i",1);
    
    //VZERO C
    AliGenFixed *gv0c=new AliGenFixed(1);
    gv0c->SetPart(kPiPlus);
    gv0c->SetMomentum(1.5);
    gv0c->SetTheta(170);
    gv0c->SetPhi(50);
    gener->AddGenerator(gv0c,"Pi+ for V0C",1);
    
    //VZERO A
    AliGenFixed *gv0a=new AliGenFixed(1);
    gv0a->SetPart(kPiMinus);
    gv0a->SetMomentum(1.5);
    gv0a->SetTheta(1.5);
    gv0a->SetPhi(70);
    gener->AddGenerator(gv0a,"Pi- for V0A",1);


    //PMD
    AliGenFixed *gpmd=new AliGenFixed(1);
    gpmd->SetPart(kGamma);
    gpmd->SetMomentum(2);
    gpmd->SetTheta(12.6);
    gpmd->SetPhi(60);
    gener->AddGenerator(gpmd,"Gamma for PMD",1);

    //ZDC
    AliGenFixed *gzdc1=new AliGenFixed(1);
    gzdc1->SetPart(kProton);
    gzdc1->SetMomentum(700);
    gzdc1->SetTheta(0.6);
    gzdc1->SetPhi(60);
    gener->AddGenerator(gzdc1,"Proton for ZDC",1);

    AliGenFixed *gzdc2=new AliGenFixed(1);
    gzdc2->SetPart(kNeutron);
    gzdc2->SetMomentum(500);
    gzdc2->SetTheta(0.6);
    gzdc2->SetPhi(60);
    gener->AddGenerator(gzdc2,"Neutron for ZDC",1);

    //T0
    AliGenFixed *gt0=new AliGenFixed(1);
    gt0->SetPart(kPiPlus);
    gt0->SetMomentum(2);
    gt0->SetTheta(5.1);
    gt0->SetPhi(60);
    gener->AddGenerator(gt0,"Pi+ for T0",1);

    AliGenFixed *gt01=new AliGenFixed(1);
    gt01->SetPart(kPiMinus);
    gt01->SetMomentum(2);
    gt01->SetTheta(5.1);
    gt01->SetPhi(60);
    gener->AddGenerator(gt01,"Pi- for T0",1);


    //ACORDE
    AliGenFixed *gacorde=new AliGenFixed(1);
    gacorde->SetPart(kMuonPlus);
    gacorde->SetMomentum(20);
    gacorde->SetTheta(90.);
    gacorde->SetPhi(90);
    gener->AddGenerator(gacorde,"Muon+ for ACORDE",1);

    AliGenFixed *gacorde1=new AliGenFixed(1);
    gacorde1->SetPart(kMuonMinus);
    gacorde1->SetMomentum(20);
    gacorde1->SetTheta(90.);
    gacorde1->SetPhi(90);
    gener->AddGenerator(gacorde1,"Muon- for ACORDE",1);

    gener->Init();


    // 
    // Activate this line if you want the vertex smearing to happen
    // track by track
    //
    //gener->SetVertexSmear(perTrack); 
    // Field (L3 0.5 T)
    TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", -1., -1., AliMagF::k5kG));

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
    Int_t   iACORDE=  1;
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
	FRAME->SetHoles(1);
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
 
    if (iITS)
    {
        //=================== ITS parameters ============================

	AliITS *ITS  = new AliITSv11Hybrid("ITS","ITS v11Hybrid");
    }

    if (iTPC)
    {
        //============================ TPC parameters ===================
        AliTPC *TPC = new AliTPCv2("TPC", "Default");
    }


    if (iTOF) {
        //=================== TOF parameters ============================
	AliTOF *TOF = new AliTOFv6T0("TOF", "normal TOF");
     }


    if (iHMPID)
    {
        //=================== HMPID parameters ===========================
        AliHMPID *HMPID = new AliHMPIDv3("HMPID", "normal HMPID");

    }


    if (iZDC)
    {
        //=================== ZDC parameters ============================

        AliZDC *ZDC = new AliZDCv3("ZDC", "normal ZDC");
    }

    if (iTRD)
    {
        //=================== TRD parameters ============================

        AliTRD *TRD = new AliTRDv1("TRD", "TRD slow simulator");
        AliTRDgeometry *geoTRD = TRD->GetGeometry();
	// Partial geometry: modules at 0,1,7,8,9,10,17
	// starting at 3h in positive direction
	geoTRD->SetSMstatus(2,0);
	geoTRD->SetSMstatus(3,0);
	geoTRD->SetSMstatus(4,0);
        geoTRD->SetSMstatus(5,0);
	geoTRD->SetSMstatus(6,0);
        geoTRD->SetSMstatus(11,0);
        geoTRD->SetSMstatus(12,0);
        geoTRD->SetSMstatus(13,0);
        geoTRD->SetSMstatus(14,0);
        geoTRD->SetSMstatus(15,0);
        geoTRD->SetSMstatus(16,0);
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
        AliMUON *MUON = new AliMUONv1("MUON","default");
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
        AliEMCAL *EMCAL = new AliEMCALv2("EMCAL", "EMCAL_COMPLETE");
    }

     if (iACORDE)
    {
        //=================== ACORDE parameters ============================
        AliACORDE *ACORDE = new AliACORDEv1("ACORDE", "normal ACORDE");
    }

     if (iVZERO)
    {
        //=================== VZERO parameters ============================
        AliVZERO *VZERO = new AliVZEROv7("VZERO", "normal VZERO");
    }


}

Float_t EtaToTheta(Float_t arg){
  return (180./TMath::Pi())*2.*atan(exp(-arg));
}
