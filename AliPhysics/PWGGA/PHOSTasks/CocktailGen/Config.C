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
#include "STEER/AliMagF.h"
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
#include "TOF/AliTOFv6T0.h"
#include "RICH/AliRICHv1.h"
#include "ZDC/AliZDCv2.h"
#include "TRD/AliTRDv1.h"
#include "TRD/AliTRDgeometry.h"
#include "FMD/AliFMDv1.h"
#include "MUON/AliMUONv1.h"
#include "MUON/AliMUONSt1GeometryBuilderV2.h"
#include "MUON/AliMUONSt2GeometryBuilderV2.h"
#include "MUON/AliMUONSlatGeometryBuilder.h"
#include "MUON/AliMUONTriggerGeometryBuilder.h"
#include "PHOS/AliPHOSv1.h"
#include "PMD/AliPMDv1.h"
#include "START/AliSTARTv1.h"
#include "EMCAL/AliEMCALv2.h"
#include "CRT/AliCRTv0.h"
#include "VZERO/AliVZEROv5.h"
#endif

static Int_t    sseed = 12345; //Set 0 to use the current time
static Int_t mesonPDG =   111; // PDG code of a neutral meson to generate
static TString beams  = "pp" ; // colliding system

Float_t EtaToTheta(Float_t arg);
enum PprGeo_t 
{
    kHoles, kNoHoles
};
static PprGeo_t geo = kHoles;

void ProcessEnvironmentVars();

void Config()
{
    // ThetaRange is (0., 180.). It was (0.28,179.72) 7/12/00 09:00
    // Theta range given through pseudorapidity limits 22/6/2001

  ProcessEnvironmentVars();

    // Set Random Number seed
  gRandom->SetSeed(sseed);
  AliLog::Message(AliLog::kInfo, Form("Seed for random number generation = %d",gRandom->GetSeed()), "Config.C", "Config.C", "Config()","Config.C", __LINE__);


   // libraries required by geant321
#if defined(__CINT__)
  gSystem->Load("liblhapdf");      // Parton density functions
  gSystem->Load("libEGPythia6");   // TGenerator interface
  gSystem->Load("libpythia6");     // Pythia
  gSystem->Load("libAliPythia6");  // ALICE specific implementations

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
    rl->SetNumberOfEventsPerFile(2000);
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
    
    Double_t yMin = -0.30;
    Double_t yMax =  0.30;
    
    Double_t phiMin = 240;
    Double_t phiMax = 340;
    
    AliGenCocktail *gener  = new AliGenCocktail();

    TString mesonParamFile = "_2.76TeV_SpectraParam.root";
    mesonParamFile.Prepend(beams);

    TString ptSpectrumName = "pi0Spectrum";
    if (mesonPDG == 111)
      ptSpectrumName = "pi0Spectrum";
    else if (mesonPDG == 221)
      ptSpectrumName = "etaSpectrum";
    else if (mesonPDG == 223)
      ptSpectrumName = "omegaSpectrum";
    else if (mesonPDG == 331)
      ptSpectrumName = "etaprimeSpectrum";
    else
      Fatal("",Form("Unknown meson PDF code to generate: %d",mesonPDG));
    
    TFile* fpp = TFile::Open(mesonParamFile);
    if (!fpp)
      Fatal("",Form("Cannot open file %s",mesonParamFile));

    TF1*   ptSpectrum = (TF1*)fpp->Get(ptSpectrumName);     
    Info("",Form("\n\tGenerator %s is initialized\n",ptSpectrum->GetTitle()));
    
    AliGenPHOSlibPlus * myGener = new AliGenPHOSlibPlus(mesonPDG,ptSpectrum) ;
    AliGenParam *genMeson = new AliGenParam(1,AliGenPHOSlib::kPion,
					    myGener->GetPt(AliGenPHOSlib::kPion, ""),
					    myGener->GetY (AliGenPHOSlib::kPion, ""),
					    myGener->GetIp(1, ""));
    
    genMeson->SetPhiRange(phiMin,phiMax);
    genMeson->SetYRange  (yMin,yMax) ;
    genMeson->SetPtRange (0.5,25.) ; 
    genMeson->SetOrigin  (0, 0, 0);     
    genMeson->SetSigma   (0., 0., 0.);   
    
    gener->AddGenerator(genMeson,ptSpectrumName,1) ;
    gener->Init();

    TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", -1., -1., AliMagF::k5kG));

    Int_t   iABSO  =  0;
    Int_t   iDIPO  =  0;
    Int_t   iFMD   =  0;
    Int_t   iFRAME =  1;
    Int_t   iHALL  =  0;
    Int_t   iITS   =  1;
    Int_t   iMAG   =  1;
    Int_t   iMUON  =  0;
    Int_t   iPHOS  =  1;
    Int_t   iPIPE  =  1;
    Int_t   iPMD   =  0;
    Int_t   iRICH  =  0;
    Int_t   iSHIL  =  0;
    Int_t   iSTART =  0;
    Int_t   iTOF   =  1;
    Int_t   iTPC   =  1;
    Int_t   iTRD   =  1;
    Int_t   iZDC   =  0;
    Int_t   iEMCAL =  0;
    Int_t   iCRT   =  0;
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
	FRAME->SetHoles(1);
    }

    if (iSHIL)
    {
        //=================== SHIL parameters ============================

        AliSHIL *SHIL = new AliSHILv2("SHIL", "Shielding Version 2");
    }


    if (iPIPE)
    {
        //=================== PIPE parameters ============================

      AliPIPE *PIPE = new AliPIPEv3("PIPE", "Beam Pipe");
    }
 
    if(iITS) {
      //=================== ITS parameters ============================
      
      AliITS *ITS  = new AliITSv11("ITS","ITS v11");
    }
    
    if (iTPC)
    {
      //============================ TPC parameters ================================
      
      AliTPC *TPC = new AliTPCv2("TPC", "Default");
    }


    if (iTOF) {
      //=================== TOF parameters ============================

      AliTOF *TOF = new AliTOFv6T0("TOF", "normal TOF");
    }


    if (iRICH)
    {
        //=================== RICH parameters ===========================
        AliRICH *RICH = new AliRICHv1("RICH", "normal RICH");

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
      //AliPHOS *PHOS = new AliPHOSv1("PHOS", "noCPV");
    }


    if (iPMD)
    {
        //=================== PMD parameters ============================
        AliPMD *PMD = new AliPMDv1("PMD", "normal PMD");
    }

    if (iSTART)
    {
        //=================== START parameters ============================
        AliSTART *START = new AliSTARTv1("START", "START Detector");
    }

    if (iEMCAL)
    {
        //=================== EMCAL parameters ============================
        AliEMCAL *EMCAL = new AliEMCALv2("EMCAL", "SHISH");
    }

     if (iCRT)
    {
        //=================== CRT parameters ============================
        AliCRT *CRT = new AliCRTv0("CRT", "normal ACORDE");
    }

     if (iVZERO)
    {
        //=================== CRT parameters ============================
        AliVZERO *VZERO = new AliVZEROv5("VZERO", "normal VZERO");
    }

     AliLog::Message(AliLog::kInfo, "End of Config", "Config.C", "Config.C", "Config()"," Config.C", __LINE__);

}

Float_t EtaToTheta(Float_t arg){
  return (180./TMath::Pi())*2.*atan(exp(-arg));
}

void ProcessEnvironmentVars()
{
  // Colliding system
  if (gSystem->Getenv("CONFIG_BEAMS")) {
    beams = gSystem->Getenv("CONFIG_BEAMS");
  }
  
  // PDG code of neutral meson
  if (gSystem->Getenv("CONFIG_MESON_PDG")) {
    mesonPDG = atoi(gSystem->Getenv("CONFIG_MESON_PDG"));
  }
  
  // Random Number seed
  if (gSystem->Getenv("CONFIG_SEED")) {
    sseed = atoi(gSystem->Getenv("CONFIG_SEED"));
  }


}
