// One can use the configuration macro in compiled mode by
// root [0] gSystem->Load("libgeant321");
// root [0] gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include\
//                   -I$ALICE_ROOT -I$ALICE/geant3/TGeant3");
// root [0] .x grun.C(1,"Config.C++")

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
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
#include "STEER/AliMagFMaps.h"
#include "STRUCT/AliBODY.h"
#include "STRUCT/AliMAG.h"
#include "STRUCT/AliABSOv0.h"
#include "STRUCT/AliDIPOv2.h"
#include "STRUCT/AliHALL.h"
#include "STRUCT/AliFRAMEv2.h"
#include "STRUCT/AliSHILv2.h"
#include "STRUCT/AliPIPEv0.h"
#include "ITS/AliITSvPPRasymmFMD.h"
#include "TPC/AliTPCv2.h"
#include "TOF/AliTOFv4T0.h"
#include "HMPID/AliHMPIDv1.h"
#include "ZDC/AliZDCv2.h"
#include "TRD/AliTRDv1.h"
#include "FMD/AliFMDv1.h"
#include "MUON/AliMUONv1.h"
#include "MUON/AliMUONSt1GeometryBuilderV2.h"
#include "MUON/AliMUONSt2GeometryBuilderV2.h"
#include "MUON/AliMUONSlatGeometryBuilder.h"
#include "MUON/AliMUONTriggerGeometryBuilder.h"
#include "PHOS/AliPHOSv1.h"
#include "PMD/AliPMDv1.h"
#include "T0/AliT0v1.h"
#include "EMCAL/AliEMCALv2.h"
#include "ACORDE/AliACORDEv0.h"
#include "VZERO/AliVZEROv5.h"
#endif

Float_t EtaToTheta(Float_t arg);
enum PprGeo_t 
{
    kHoles, kNoHoles
};
static PprGeo_t geo = kHoles;

void Config()
{
    // ThetaRange is (0., 180.). It was (0.28,179.72) 7/12/00 09:00
    // Theta range given through pseudorapidity limits 22/6/2001

    // Set Random Number seed
  gRandom->SetSeed(123456); // Set 0 to use the currecnt time
  AliLog::Message(AliLog::kInfo, Form("Seed for random number generation = %d",gRandom->GetSeed()), "Config.C", "Config.C", "Config()","Config.C", __LINE__);


   // libraries required by geant321
#if defined(__CINT__)
    gSystem->Load("libgeant321");
#endif

    new     TGeant3TGeo("C++ Interface to Geant3");

    AliRunLoader* rl=0x0;

    AliLog::Message(AliLog::kInfo, "Creating Run Loader", "Config.C", "Config.C", "Config()"," Config.C", __LINE__);

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


    int     nParticles = 100;
    if (gSystem->Getenv("CONFIG_NPARTICLES"))
    {
        nParticles = atoi(gSystem->Getenv("CONFIG_NPARTICLES"));
    }

    AliGenCocktail *gener = new AliGenCocktail();
    gener->SetPhiRange(220, 320);
    // Set pseudorapidity range from -0.13 to 0.13.
    Float_t thmin = EtaToTheta(+0.13);  // theta min. <---> eta max
    Float_t thmax = EtaToTheta(-0.13);  // theta max. <---> eta min 
    gener->SetThetaRange(thmin,thmax);
    gener->SetOrigin(0, 0, 0);  //vertex position
    gener->SetSigma(0, 0, 0);   //Sigma in (X,Y,Z) (cm) on IP position

    AliGenBox *genbox = new AliGenBox(nParticles);
    genbox->SetPart(11);
    genbox->SetPtRange(0.3, 10.00);
    gener->AddGenerator(genbox,"GENBOX GAMMA for PHOS",1);
    gener->Init();


    // 
    // Activate this line if you want the vertex smearing to happen
    // track by track
    //
    //gener->SetVertexSmear(perTrack); 
    // Field (L3 0.4 T)
    AliMagFMaps* field = new AliMagFMaps("Maps","Maps", 2, 1., 10., 1);
    gAlice->SetField(field);    


    Int_t   iABSO  =  0;
    Int_t   iDIPO  =  0;
    Int_t   iFMD   =  0;
    Int_t   iFRAME =  0;
    Int_t   iHALL  =  0;
    Int_t   iITS   =  0;
    Int_t   iMAG   =  0;
    Int_t   iMUON  =  0;
    Int_t   iPHOS  =  1;
    Int_t   iPIPE  =  1;
    Int_t   iPMD   =  0;
    Int_t   iHMPID  =  0;
    Int_t   iSHIL  =  0;
    Int_t   iT0 =  0;
    Int_t   iTOF   =  0;
    Int_t   iTPC   =  0;
    Int_t   iTRD   =  0;
    Int_t   iZDC   =  0;
    Int_t   iEMCAL =  0;
    Int_t   iACORDE   =  0;
    Int_t   iVZERO =  0;
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
        AliABSO *ABSO = new AliABSOv0("ABSO", "Muon Absorber");
    }

    if (iDIPO)
    {
        //=================== DIPO parameters ============================

        AliDIPO *DIPO = new AliDIPOv2("DIPO", "Dipole version 2");
    }

    if (iHALL)
    {
        //=================== HALL parameters ============================

        AliHALL *HALL = new AliHALL("HALL", "Alice Hall");
    }


    if (iFRAME)
    {
        //=================== FRAME parameters ============================

        AliFRAMEv2 *FRAME = new AliFRAMEv2("FRAME", "Space Frame");
	if (geo == kHoles) {
	    FRAME->SetHoles(1);
	} else {
	    FRAME->SetHoles(0);
	}
    }

    if (iSHIL)
    {
        //=================== SHIL parameters ============================

        AliSHIL *SHIL = new AliSHILv2("SHIL", "Shielding Version 2");
    }


    if (iPIPE)
    {
        //=================== PIPE parameters ============================

        AliPIPE *PIPE = new AliPIPEv0("PIPE", "Beam Pipe");
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
    //
	AliITSvPPRasymmFMD *ITS  = new AliITSvPPRasymmFMD("ITS","ITS PPR detailed version with asymmetric services");
	ITS->SetMinorVersion(2);  // don't touch this parameter if you're not an ITS developer
	ITS->SetReadDet(kTRUE);	  // don't touch this parameter if you're not an ITS developer
	//    ITS->SetWriteDet("$ALICE_ROOT/ITS/ITSgeometry_vPPRasymm2.det");  // don't touch this parameter if you're not an ITS developer
	ITS->SetThicknessDet1(200.);   // detector thickness on layer 1 must be in the range [100,300]
	ITS->SetThicknessDet2(200.);   // detector thickness on layer 2 must be in the range [100,300]
	ITS->SetThicknessChip1(200.);  // chip thickness on layer 1 must be in the range [150,300]
	ITS->SetThicknessChip2(200.);  // chip thickness on layer 2 must be in the range [150,300]
	ITS->SetRails(0);	       // 1 --> rails in ; 0 --> rails out
	ITS->SetCoolingFluid(1);       // 1 --> water ; 0 --> freon

 
    //
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
	AliTOF *TOF = new AliTOFv4T0("TOF", "normal TOF");
    }


    if (iHMPID)
    {
        //=================== HMPID parameters ===========================
        AliHMPID *HMPID = new AliHMPIDv1("HMPID", "normal HMPID");

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

        // Select the gas mixture (0: 97% Xe + 3% isobutane, 1: 90% Xe + 10% CO2)
        TRD->SetGasMix(1);
	if (geo == kHoles) {
	    // With hole in front of PHOS
	    TRD->SetPHOShole();
	    // With hole in front of HMPID
	    TRD->SetHMPIDhole();
	}
	    // Switch on TR
	    AliTRDsim *TRDsim = TRD->CreateTR();
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
	((AliMUONv1*)MUON)->SetStepManagerVersionDE(true);
	MUON->AddGeometryBuilder(new AliMUONSt1GeometryBuilderV2(MUON));
	MUON->AddGeometryBuilder(new AliMUONSt2GeometryBuilderV2(MUON));
	MUON->AddGeometryBuilder(new AliMUONSlatGeometryBuilder(MUON));
	MUON->AddGeometryBuilder(new AliMUONTriggerGeometryBuilder(MUON));
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
        AliEMCAL *EMCAL = new AliEMCALv2("EMCAL", "SHISH");
    }

     if (iACORDE)
    {
        //=================== ACORDE parameters ============================
        AliACORDE *ACORDE = new AliACORDEv0("ACORDE", "normal ACORDE");
    }

     if (iVZERO)
    {
        //=================== ACORDE parameters ============================
        AliVZERO *VZERO = new AliVZEROv5("VZERO", "normal VZERO");
    }

     AliLog::Message(AliLog::kInfo, "End of Config", "Config.C", "Config.C", "Config()"," Config.C", __LINE__);

}

Float_t EtaToTheta(Float_t arg){
  return (180./TMath::Pi())*2.*atan(exp(-arg));
}
