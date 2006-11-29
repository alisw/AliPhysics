#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TRandom.h>
#include <TSystem.h>
#include <TVirtualMC.h>
#include <TGeant3.h>
#include "STEER/AliRunLoader.h"
#include "STEER/AliRun.h"
#include "STEER/AliConfig.h"
#include "STEER/AliGenerator.h"
#include "PYTHIA6/AliDecayerPythia.h"
#include "EVGEN/AliGenHIJINGpara.h"
#include "THijing/AliGenHijing.h"
#include "EVGEN/AliGenCocktail.h"
#include "EVGEN/AliGenSlowNucleons.h"
#include "EVGEN/AliSlowNucleonModelExp.h"
#include "PYTHIA6/AliGenPythia.h"
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
#include "PHOS/AliPHOSv1.h"
#include "PMD/AliPMDv1.h"
#include "START/AliSTARTv1.h"
#include "EMCAL/AliEMCALv1.h"
#include "CRT/AliCRTv0.h"
#include "VZERO/AliVZEROv3.h"
#endif

enum PprRun_t 
{
    test50,
    kParam_8000,   kParam_4000,  kParam_2000, 
    kHijing_cent1, kHijing_cent2, 
    kHijing_per1,  kHijing_per2, kHijing_per3, kHijing_per4,  kHijing_per5,
    kHijing_jj25,  kHijing_jj50, kHijing_jj75, kHijing_jj100, kHijing_jj200, 
    kHijing_gj25,  kHijing_gj50, kHijing_gj75, kHijing_gj100, kHijing_gj200,
    kHijing_pA, kPythia6, kPythia6Jets, kD0PbPb5500
};

enum PprGeo_t 
{
    kHoles, kNoHoles
};

enum PprRad_t
{
    kGluonRadiation, kNoGluonRadiation
};

enum PprMag_t
{
    k2kG, k4kG, k5kG
};


// This part for configuration    
static PprRun_t srun = test50;
//static PprRun_t srun = kPythia6;
static PprGeo_t sgeo = kNoHoles;
static PprRad_t srad = kGluonRadiation;
static PprMag_t smag = k5kG;

// Comment line 
static TString  comment;

// Functions
Float_t EtaToTheta(Float_t arg);
AliGenerator* GeneratorFactory(PprRun_t srun);
AliGenHijing* HijingStandard();

void Config()
{
    // ThetaRange is (0., 180.). It was (0.28,179.72) 7/12/00 09:00
    // Theta range given through pseudorapidity limits 22/6/2001

    // Set Random Number seed
  gRandom->SetSeed(12345); //Set 0 to use the current time
    cout<<"Seed for random number generation= "<<gRandom->GetSeed()<<endl; 


   // libraries required by geant321
#if defined(__CINT__)
    gSystem->Load("libgeant321");
#endif

    new     TGeant3("C++ Interface to Geant3");

    AliRunLoader* rl=0x0;

    cout<<"Config.C: Creating Run Loader ..."<<endl;
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
    AliDecayer *decayer = new AliDecayerPythia();

    decayer->SetForceDecay(kAll);
    decayer->Init();
    gMC->SetExternalDecayer(decayer);
    //
    //
    //=======================================================================
    //
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


// Generator Configuration
    gAlice->SetDebug(1);
    AliGenerator* gener = GeneratorFactory(srun);
    gener->SetOrigin(0, 0, 0);    // vertex position
    gener->SetSigma(0, 0, 5.3);   // Sigma in (X,Y,Z) (cm) on IP position
    gener->SetCutVertexZ(1.);     // Truncate at 1 sigma
    gener->SetVertexSmear(kPerEvent); 
    gener->SetTrackingFlag(1);
    gener->Init();
    
    if (smag == k2kG) {
	comment = comment.Append(" | L3 field 0.2 T");
    } else if (smag == k4kG) {
	comment = comment.Append(" | L3 field 0.4 T");
    } else if (smag == k5kG) {
	comment = comment.Append(" | L3 field 0.5 T");
    }
    
    
    if (srad == kGluonRadiation)
    {
	comment = comment.Append(" | Gluon Radiation On");
	
    } else {
	comment = comment.Append(" | Gluon Radiation Off");
    }

    if (sgeo == kHoles)
    {
	comment = comment.Append(" | Holes for PHOS/HMPID");
	
    } else {
	comment = comment.Append(" | No holes for PHOS/HMPID");
    }

    printf("\n \n Comment: %s \n \n", comment.Data());
    
    
// Field (L3 0.4 T)
    AliMagFMaps* field = new AliMagFMaps("Maps","Maps", 2, 1., 10., smag);
    field->SetL3ConstField(0); //Using const. field in the barrel
    rl->CdGAFile();
    gAlice->SetField(field);    
//
    Int_t   iABSO   = 0;
    Int_t   iDIPO   = 0;
    Int_t   iFMD    = 0;
    Int_t   iFRAME  = 1;
    Int_t   iHALL   = 0;
    Int_t   iITS    = 1;
    Int_t   iMAG    = 1;
    Int_t   iMUON   = 0;
    Int_t   iPHOS   = 0;
    Int_t   iPIPE   = 1;
    Int_t   iPMD    = 0;
    Int_t   iHMPID   = 0;
    Int_t   iSHIL   = 0;
    Int_t   iSTART  = 0;
    Int_t   iTOF    = 1;
    Int_t   iTPC    = 1;
    Int_t   iTRD    = 1;
    Int_t   iZDC    = 0;
    Int_t   iEMCAL  = 0;
    Int_t   iVZERO  = 0;
    Int_t   iCRT    = 0;

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
	if (sgeo == kHoles) {
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
    //AliITS *ITS  = new AliITSv5symm("ITS","Updated ITS TDR detailed version with symmetric services");
    //
    //AliITS *ITS  = new AliITSv5asymm("ITS","Updates ITS TDR detailed version with asymmetric services");
    //
	AliITSvPPRasymmFMD *ITS  = new AliITSvPPRasymmFMD("ITS","New ITS PPR detailed version with asymmetric services");
	ITS->SetMinorVersion(2);  // don't touch this parameter if you're not an ITS developer
	ITS->SetReadDet(kTRUE);	  // don't touch this parameter if you're not an ITS developer
    //    ITS->SetWriteDet("$ALICE_ROOT/ITS/ITSgeometry_vPPRasymm2.det");  // don't touch this parameter if you're not an ITS developer
	ITS->SetThicknessDet1(200.);   // detector thickness on layer 1 must be in the range [100,300]
	ITS->SetThicknessDet2(200.);   // detector thickness on layer 2 must be in the range [100,300]
	ITS->SetThicknessChip1(200.);  // chip thickness on layer 1 must be in the range [150,300]
	ITS->SetThicknessChip2(200.);  // chip thickness on layer 2 must be in the range [150,300]
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
        //============================ TPC parameters ================================
        // --- This allows the user to specify sectors for the SLOW (TPC geometry 2)
        // --- Simulator. SecAL (SecAU) <0 means that ALL lower (upper)
        // --- sectors are specified, any value other than that requires at least one 
        // --- sector (lower or upper)to be specified!
        // --- Reminder: sectors 1-24 are lower sectors (1-12 -> z>0, 13-24 -> z<0)
        // ---           sectors 25-72 are the upper ones (25-48 -> z>0, 49-72 -> z<0)
        // --- SecLows - number of lower sectors specified (up to 6)
        // --- SecUps - number of upper sectors specified (up to 12)
        // --- Sens - sensitive strips for the Slow Simulator !!!
        // --- This does NOT work if all S or L-sectors are specified, i.e.
        // --- if SecAL or SecAU < 0
        //
        //
        //-----------------------------------------------------------------------------

        //  gROOT->LoadMacro("SetTPCParam.C");
        //  AliTPCParam *param = SetTPCParam();
        AliTPC *TPC = new AliTPCv2("TPC", "Default");

        // All sectors included 
        TPC->SetSecAL(-1);
        TPC->SetSecAU(-1);

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
	if (sgeo == kHoles) {
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

    if (iSTART)
    {
        //=================== START parameters ============================
        AliSTART *START = new AliSTARTv1("START", "START Detector");
    }

    if (iEMCAL)
    {
        //=================== EMCAL parameters ============================
        AliEMCAL *EMCAL = new AliEMCALv1("EMCAL", "G56_2_55_19");
    }

     if (iCRT)
    {
        //=================== CRT parameters ============================
        AliCRT *CRT = new AliCRTv0("CRT", "normal ACORDE");
    }

     if (iVZERO)
    {
        //=================== CRT parameters ============================
        AliVZERO *VZERO = new AliVZEROv3("VZERO", "normal VZERO");
    }
 
             
}

Float_t EtaToTheta(Float_t arg){
  return (180./TMath::Pi())*2.*atan(exp(-arg));
}



AliGenerator* GeneratorFactory(PprRun_t srun) {
    Int_t isw = 3;
    if (srad == kNoGluonRadiation) isw = 0;
    

    AliGenerator * gGener = 0x0;
    switch (srun) {
    case test50:
      {
	comment = comment.Append(":HIJINGparam test 50 particles");
	AliGenHIJINGpara *gener = new AliGenHIJINGpara(50);
	gener->SetMomentumRange(0, 999999.);
	gener->SetPhiRange(0., 360.);
	// Set pseudorapidity range from -8 to 8.
	Float_t thmin = EtaToTheta(1);   // theta min. <---> eta max
	Float_t thmax = EtaToTheta(-1);  // theta max. <---> eta min 
	gener->SetThetaRange(thmin,thmax);
	gGener=gener;
      }
      break;
    case kParam_8000:
      {
	comment = comment.Append(":HIJINGparam N=8000");
	AliGenHIJINGpara *gener = new AliGenHIJINGpara(86030);
	gener->SetMomentumRange(0, 999999.);
	gener->SetPhiRange(0., 360.);
	// Set pseudorapidity range from -8 to 8.
	Float_t thmin = EtaToTheta(8);   // theta min. <---> eta max
	Float_t thmax = EtaToTheta(-8);  // theta max. <---> eta min 
	gener->SetThetaRange(thmin,thmax);
	gGener=gener;
      }
      break;
    case kParam_4000:
      {
	comment = comment.Append("HIJINGparam N=4000");
	AliGenHIJINGpara *gener = new AliGenHIJINGpara(43015);
	gener->SetMomentumRange(0, 999999.);
	gener->SetPhiRange(0., 360.);
	// Set pseudorapidity range from -8 to 8.
	Float_t thmin = EtaToTheta(8);   // theta min. <---> eta max
	Float_t thmax = EtaToTheta(-8);  // theta max. <---> eta min 
	gener->SetThetaRange(thmin,thmax);
	gGener=gener;
      }
	break;
    case kParam_2000:
      {
	comment = comment.Append("HIJINGparam N=2000");
	AliGenHIJINGpara *gener = new AliGenHIJINGpara(21507);
	gener->SetMomentumRange(0, 999999.);
	gener->SetPhiRange(0., 360.);
	// Set pseudorapidity range from -8 to 8.
	Float_t thmin = EtaToTheta(8);   // theta min. <---> eta max
	Float_t thmax = EtaToTheta(-8);  // theta max. <---> eta min 
	gener->SetThetaRange(thmin,thmax);
	gGener=gener;
      }
      break;
//
//  Hijing Central
//
    case kHijing_cent1:
      {
	comment = comment.Append("HIJING cent1");
	AliGenHijing *gener = HijingStandard();
// impact parameter range
	gener->SetImpactParameterRange(0., 5.);
	gGener=gener;
      }
      break;
    case kHijing_cent2:
      {
	comment = comment.Append("HIJING cent2");
	AliGenHijing *gener = HijingStandard();
// impact parameter range
	gener->SetImpactParameterRange(0., 2.);
	gGener=gener;
      }
      break;
//
// Hijing Peripheral 
//
    case kHijing_per1:
      {
	comment = comment.Append("HIJING per1");
	AliGenHijing *gener = HijingStandard();
// impact parameter range
	gener->SetImpactParameterRange(5., 8.6);
	gGener=gener;
      }
      break;
    case kHijing_per2:
      {
	comment = comment.Append("HIJING per2");
	AliGenHijing *gener = HijingStandard();
// impact parameter range
	gener->SetImpactParameterRange(8.6, 11.2);
	gGener=gener;
      }
      break;
    case kHijing_per3:
      {
	comment = comment.Append("HIJING per3");
	AliGenHijing *gener = HijingStandard();
// impact parameter range
	gener->SetImpactParameterRange(11.2, 13.2);
	gGener=gener;
      }
      break;
    case kHijing_per4:
      {
	comment = comment.Append("HIJING per4");
	AliGenHijing *gener = HijingStandard();
// impact parameter range
	gener->SetImpactParameterRange(13.2, 15.);
	gGener=gener;
      }
      break;
    case kHijing_per5:
      {
	comment = comment.Append("HIJING per5");
	AliGenHijing *gener = HijingStandard();
// impact parameter range
	gener->SetImpactParameterRange(15., 100.);
	gGener=gener;
      }
      break;
//
//  Jet-Jet
//
    case kHijing_jj25:
      {
	comment = comment.Append("HIJING Jet 25 GeV");
	AliGenHijing *gener = HijingStandard();
// impact parameter range
	gener->SetImpactParameterRange(0., 5.);
	// trigger
	gener->SetTrigger(1);
	gener->SetPtJet(25.);
	gener->SetRadiation(isw);
	gener->SetSimpleJets(!isw);
	gener->SetJetEtaRange(-0.3,0.3);
	gener->SetJetPhiRange(75., 165.);   
	gGener=gener;
      }
      break;

    case kHijing_jj50:
      {
	comment = comment.Append("HIJING Jet 50 GeV");
	AliGenHijing *gener = HijingStandard();
// impact parameter range
	gener->SetImpactParameterRange(0., 5.);
	// trigger
	gener->SetTrigger(1);
	gener->SetPtJet(50.);
	gener->SetRadiation(isw);
	gener->SetSimpleJets(!isw);
	gener->SetJetEtaRange(-0.3,0.3);
	gener->SetJetPhiRange(75., 165.);   
	gGener=gener;
      }
	break;

    case kHijing_jj75:
      {
	comment = comment.Append("HIJING Jet 75 GeV");
	AliGenHijing *gener = HijingStandard();
// impact parameter range
	gener->SetImpactParameterRange(0., 5.);
	// trigger
	gener->SetTrigger(1);
	gener->SetPtJet(75.);
	gener->SetRadiation(isw);
	gener->SetSimpleJets(!isw);
	gener->SetJetEtaRange(-0.3,0.3);
	gener->SetJetPhiRange(75., 165.);   
	gGener=gener;
      }
      break;

    case kHijing_jj100:
      {
	comment = comment.Append("HIJING Jet 100 GeV");
	AliGenHijing *gener = HijingStandard();
// impact parameter range
	gener->SetImpactParameterRange(0., 5.);
	// trigger
	gener->SetTrigger(1);
	gener->SetPtJet(100.);
	gener->SetRadiation(isw);
	gener->SetSimpleJets(!isw);
	gener->SetJetEtaRange(-0.3,0.3);
	gener->SetJetPhiRange(75., 165.);   
	gGener=gener;
      }
      break;

    case kHijing_jj200:
      {
	comment = comment.Append("HIJING Jet 200 GeV");
	AliGenHijing *gener = HijingStandard();
// impact parameter range
	gener->SetImpactParameterRange(0., 5.);
	// trigger
	gener->SetTrigger(1);
	gener->SetPtJet(200.);
	gener->SetRadiation(isw);
	gener->SetSimpleJets(!isw);
	gener->SetJetEtaRange(-0.3,0.3);
	gener->SetJetPhiRange(75., 165.);   
	gGener=gener;
      }
      break;
//
// Gamma-Jet
//
    case kHijing_gj25:
      {
	comment = comment.Append("HIJING Gamma 25 GeV");
	AliGenHijing *gener = HijingStandard();
// impact parameter range
	gener->SetImpactParameterRange(0., 5.);
	// trigger
	gener->SetTrigger(2);
	gener->SetPtJet(25.);
	gener->SetRadiation(isw);
	gener->SetSimpleJets(!isw);
	gener->SetJetEtaRange(-0.12, 0.12);
        gener->SetJetPhiRange(220., 320.);
	gGener=gener;
      }
      break;

    case kHijing_gj50:
      {
	comment = comment.Append("HIJING Gamma 50 GeV");
	AliGenHijing *gener = HijingStandard();
// impact parameter range
	gener->SetImpactParameterRange(0., 5.);
	// trigger
	gener->SetTrigger(2);
	gener->SetPtJet(50.);
	gener->SetRadiation(isw);
	gener->SetSimpleJets(!isw);
	gener->SetJetEtaRange(-0.12, 0.12);
        gener->SetJetPhiRange(220., 320.);
	gGener=gener;
      }
      break;

    case kHijing_gj75:
      {
	comment = comment.Append("HIJING Gamma 75 GeV");
	AliGenHijing *gener = HijingStandard();
// impact parameter range
	gener->SetImpactParameterRange(0., 5.);
	// trigger
	gener->SetTrigger(2);
	gener->SetPtJet(75.);
	gener->SetRadiation(isw);
	gener->SetSimpleJets(!isw);
	gener->SetJetEtaRange(-0.12, 0.12);
        gener->SetJetPhiRange(220., 320.);
	gGener=gener;
      }
      break;

    case kHijing_gj100:
      {
	comment = comment.Append("HIJING Gamma 100 GeV");
	AliGenHijing *gener = HijingStandard();
// impact parameter range
	gener->SetImpactParameterRange(0., 5.);
	// trigger
	gener->SetTrigger(2);
	gener->SetPtJet(100.);
	gener->SetRadiation(isw);
	gener->SetSimpleJets(!isw);
	gener->SetJetEtaRange(-0.12, 0.12);
        gener->SetJetPhiRange(220., 320.);
	gGener=gener;
      }
      break;

    case kHijing_gj200:
      {
	comment = comment.Append("HIJING Gamma 200 GeV");
	AliGenHijing *gener = HijingStandard();
// impact parameter range
	gener->SetImpactParameterRange(0., 5.);
	// trigger
	gener->SetTrigger(2);
	gener->SetPtJet(200.);
	gener->SetRadiation(isw);
	gener->SetSimpleJets(!isw);
	gener->SetJetEtaRange(-0.12, 0.12);
        gener->SetJetPhiRange(220., 320.);
	gGener=gener;
      }
      break;
    case kHijing_pA:
      {
	comment = comment.Append("HIJING pA");

	AliGenCocktail *gener  = new AliGenCocktail();
//	gener->SetTrackingFlag(0);
	gener->SetOrigin(0, 0, 0);    // vertex position
	gener->SetSigma(0, 0, 5.3);   // Sigma in (X,Y,Z) (cm) on IP position
	gener->SetCutVertexZ(1.);     // Truncate at 1 sigma
	gener->SetVertexSmear(kPerEvent); 
	AliGenHijing   *hijing = new AliGenHijing(-1);
// centre of mass energy 
	hijing->SetEnergyCMS(TMath::Sqrt(82./208.) * 14000.);
// impact parameter range
	hijing->SetImpactParameterRange(0., 15.);
// reference frame
	hijing->SetReferenceFrame("CMS");
	hijing->SetBoostLHC(1);
// projectile
	hijing->SetProjectile("P", 1, 1);
	hijing->SetTarget    ("A", 208, 82);
// tell hijing to keep the full parent child chain
	hijing->KeepFullEvent();
// enable jet quenching
	hijing->SetJetQuenching(0);
// enable shadowing
	hijing->SetShadowing(1);
// Don't track spectators
	hijing->SetSpectators(0);
// kinematic selection
	hijing->SetSelectAll(0);
//
	AliGenSlowNucleons*  gray    = new AliGenSlowNucleons(1);
	AliSlowNucleonModel* model   = new AliSlowNucleonModelExp();
	gray->SetSlowNucleonModel(model);
	gray->SetDebug(1);
	gener->AddGenerator(hijing,"Hijing pPb", 1);
	gener->AddGenerator(gray,  "Gray Particles",1);
	gGener=gener;
      }
      break;
    case kPythia6:
      {
        comment = comment.Append(":Pythia p-p @ 14 TeV");
        AliGenPythia *gener = new AliGenPythia(-1); 
        gener->SetMomentumRange(0,999999);
        gener->SetPhiRange(-180,180);
        gener->SetThetaRange(0., 180.);
        gener->SetYRange(-12,12);
        gener->SetPtRange(0,1000);
        gener->SetStrucFunc(kCTEQ4L);   
        gener->SetProcess(kPyMb);
        gener->SetEnergyCMS(14000.);
	gGener=gener;
      }
      break;
    case kPythia6Jets:
      {
        comment = comment.Append(":Pythia jets @ 5.5 TeV");
        AliGenPythia * gener = new AliGenPythia(-1);
	gener->SetEnergyCMS(5500.);//        Centre of mass energy
	gener->SetProcess(kPyJets);//        Process type
	gener->SetJetEtaRange(-0.5, 0.5);//  Final state kinematic cuts
	gener->SetJetPhiRange(0., 360.);
	gener->SetJetEtRange(10., 1000.);
	gener->SetGluonRadiation(1,1);
	//    gener->SetPtKick(0.);
	//   Structure function
	gener->SetStrucFunc(kCTEQ4L);
	gener->SetPtHard(86., 104.);// Pt transfer of the hard scattering
	gener->SetPycellParameters(2., 274, 432, 0., 4., 5., 1.0);
	gener->SetForceDecay(kAll);//  Decay type (semielectronic, etc.)
	gGener=gener;
      }
      break;
    case kD0PbPb5500:
      {
	comment = comment.Append(" D0 in Pb-Pb at 5.5 TeV");
	AliGenPythia * gener = new AliGenPythia(10);
	gener->SetProcess(kPyD0PbPbMNR);
	gener->SetStrucFunc(kCTEQ4L);
	gener->SetPtHard(2.1,-1.0);
	gener->SetEnergyCMS(5500.);
	gener->SetNuclei(208,208);
	gGener=gener;
      }
      break;
    }
    return gGener;
}

AliGenHijing* HijingStandard()
{
    AliGenHijing *gener = new AliGenHijing(-1);
// centre of mass energy 
    gener->SetEnergyCMS(5500.);
// reference frame
    gener->SetReferenceFrame("CMS");
// projectile
     gener->SetProjectile("A", 208, 82);
     gener->SetTarget    ("A", 208, 82);
// tell hijing to keep the full parent child chain
     gener->KeepFullEvent();
// enable jet quenching
     gener->SetJetQuenching(4);
// enable shadowing
     gener->SetShadowing(1);
// neutral pion and heavy particle decays switched off
     gener->SetDecaysOff(1);
// Don't track spectators
     gener->SetSpectators(0);
// kinematic selection
     gener->SetSelectAll(0);
     return gener;
}



