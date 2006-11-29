// One can use the configuration macro in compiled mode by
// root [0] gSystem->Load("libgeant321");
// root [0] gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include\
//                   -I$ALICE_ROOT -I$ALICE/geant3/TGeant3");
// root [0] .x grun.C(1,"ConfigHBT.C++")

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TRandom.h>
#include <TSystem.h>
#include <TVirtualMC.h>
#include <TGeant3.h>
#include "STEER/AliRunLoader.h"
#include "STEER/AliRun.h"
#include "STEER/AliConfig.h"
#include "PYTHIA6/AliDecayerPythia.h"
#include "EVGEN/AliGenCocktailAfterBurner.h"
#include "TMEVSIM/AliMevSimConfig.h"
#include "TMEVSIM/AliMevSimParticle.h"
#include "TMEVSIM/AliGenMevSim.h"
#include "THbtp/AliGenHBTprocessor.h"
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
#include "TOF/AliTOFv5T0.h"
#include "HMPID/AliHMPIDv1.h"
#include "ZDC/AliZDCv1.h"
#include "TRD/AliTRDv1.h"
#include "FMD/AliFMDv1.h"
#include "MUON/AliMUONv1.h"
#include "PHOS/AliPHOSv1.h"
#include "PMD/AliPMDv1.h"
#include "T0/AliT0v1.h"
#include "EMCAL/AliEMCALv2.h"
#include "CRT/AliCRTv1.h"
#endif

void Config()
{
    // Set Random Number seed
    // gRandom->SetSeed(12345);

   // libraries required by geant321
#if defined(__CINT__)
    gSystem->Load("libgeant321");
#endif

    new     TGeant3TGeo("C++ Interface to Geant3");

    cout<<"Config.C: Creating Run Loader ..."<<endl;
    AliRunLoader* rl = AliRunLoader::Open("galice.root",AliConfig::GetDefaultEventFolderName(),
                                              "recreate");
    if (rl == 0x0)
      {
	gAlice->Fatal("Config.C","Can not instatiate the Run Loader");
	return;
      }
    rl->SetCompressionLevel(2);
    rl->SetNumberOfEventsPerFile(6);        
    gAlice->SetRunLoader(rl);

    // Set the trigger configuration
    gAlice->SetTriggerDescriptor("Pb-Pb");
    cout<<"Trigger configuration is set to  Pb-Pb"<<endl;

    //
    // Set External decayer
    AliDecayer *decayer = new AliDecayerPythia();

    decayer->SetForceDecay(kAll);
    decayer->Init();
    gMC->SetExternalDecayer(decayer);
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


    if (gSystem->Getenv("CONFIG_NPARTICLES"))
    {
        int     nParticles = atoi(gSystem->Getenv("CONFIG_NPARTICLES"));
    } else
    {
        int     nParticles = 10000;
    }
    
    /*********************************************************************************/
    /*********************************************************************************/
    /*****************     G E N E R A T O R S       *********************************/    

    /**  C O C K T A I L   W I T H   A F T E R B U R N E R S   A F T E R W A R D S  **/    
    /*****************        (WRAPPER GENERATOR)    *********************************/    


    
    AliGenCocktailAfterBurner* gener = new AliGenCocktailAfterBurner();
    gener->SetPhiRange(0, 360);
    gener->SetThetaRange(40., 140.);
    gener->SetOrigin(0, 0, 0);  //vertex position
    gener->SetSigma(0, 0, 0);   //Sigma in (X,Y,Z) (cm) on IP position
    
    /*****************     HIJING PARAMETRIZATION    ********************************/    
    
    //AliGenHIJINGpara *g2 = new AliGenHIJINGpara(5);
    //g2->SetMomentumRange(0, 999);
    
    /*****************        G E N   B O X       ***********************************/    
    
   // AliGenBox *g1 = new AliGenBox(5);
   // g1->SetMomentumRange(0, 3);
   // g1->SetOrigin(0,0,0);        //vertex position
   // g1->SetSigma(0,0,0);         //Sigma in (X,Y,Z) (cm) on IP position
   // g1->SetPart(310);        //K0short(310), Lambda(3122)
   

    /*****************     M E V S I M    ********************************/    
    
    AliMevSimConfig *c = new AliMevSimConfig(1);
    c->SetRectPlane(1);                                 // reaction plane control, model 4
    c->SetGrid(80,80);
    
    AliGenMevSim *g3 = new AliGenMevSim(c);
    g3->SetPtRange(0.001, 3); 
    g3->SetMomentumRange(0.1, 3); 
    
    g3->SetTrackingFlag(1);
    g3->SetOrigin(0.0, 0.0, 0.0);
    g3->SetSigma(0.0, 0.0, 0.0);
    //g3->SetEtaCutRange(-4,4);       
 
    AliMevSimParticle *piPlus = new AliMevSimParticle(kPiPlus, 100, 3, 0.2, 0.01, 3, 0.1, 0.0, 0.0);
    AliMevSimParticle *piMinus = new AliMevSimParticle(kPiMinus, 100, 3, 0.2, 0.01, 3, 0.1, 0.0, 0.0);
    //AliMevSimParticle *KPlus = new AliMevSimParticle(kKPlus, 10, 0, 0.25, 0.0, 2, 0.15, 0.0, 0.0 );
    //AliMevSimParticle *KMinus = new AliMevSimParticle(kKMinus, 10, 0, 0.25, 0.0, 2, 0.15, 0.0, 0.0 );
    //AliMevSimParticle *protonPlus = new AliMevSimParticle(kProton, 3, 0,  0.4, 0.0, 2, 0.15, 0.0, 0.0);
    //AliMevSimParticle *protonMinus = new AliMevSimParticle(kProtonBar, 2, 0, 0.4, 0.0, 2, 0.15, 0.0, 0.0);
 
    g3->AddParticleType(piPlus);
    g3->AddParticleType(piMinus);
    //g3->AddParticleType(KPlus);
    //g3->AddParticleType(KMinus);
    //g3->AddParticleType(protonPlus);
    //g3->AddParticleType(protonMinus);
    
    
    /*****************     H B T   P R O C E S S O R    ********************************/    
    
    AliGenHBTprocessor *hbtp = new AliGenHBTprocessor();
    hbtp->SetRefControl(2);
    hbtp->SetSwitch1D(1);
    hbtp->SetSwitch3D(1);
    hbtp->SetSwitchCoulomb(1);
    hbtp->SetSwitchCoherence(0);
    hbtp->SetSwitchFermiBose(-1);
    hbtp->SetDeltap(0.1);
    hbtp->SetDelChi(0.1);
    hbtp->SetLambda(0.5);
    hbtp->SetQ0(8.0);
    hbtp->SetR0(8.0);
    hbtp->SetRParallel(8.0);
    hbtp->SetR1d(8.0);
    hbtp->SetRSide(8.0);
    hbtp->SetRLong(8.0);
    hbtp->SetROut(8.0);
    hbtp->SetPtRange(0.1,0.98);
    hbtp->SetPIDs(211,-211); //pi+ pi-
    hbtp->SetSwitchType(3);  // fit both the like and unlike pair correl
    hbtp->SetMaxIterations(300);
    /***********************************************************************************/
    
   // gener->AddGenerator(g2,"HIJING PARAMETRIZATION",1);
   // gener->AddGenerator(g1,"BOX (K0short)",1);
    gener->AddGenerator(g3,"MEVSIM",1);
    
    gener->AddAfterBurner(hbtp,"HBT PROCESSOR",1);
    
    gener->Init();
    // 
    // Activate this line if you want the vertex smearing to happen
    // track by track
    //
    //gener->SetVertexSmear(perTrack); 

    gAlice->SetField(-999, 2);  //Specify maximum magnetic field in Tesla (neg. ==> default field)
//    gAlice->SetField(-999, 2, 2.);  //Specify maximum magnetic field in Tesla (neg. ==> default field)
     //Last number indicates the scale factor 

    Int_t   iABSO = 1;
    Int_t   iCRT = 0;
    Int_t   iDIPO = 1;
    Int_t   iFMD = 0;
    Int_t   iFRAME = 1;
    Int_t   iHALL = 1;
    Int_t   iITS = 1;
    Int_t   iMAG = 1;
    Int_t   iMUON = 0;
    Int_t   iPHOS = 0;
    Int_t   iPIPE = 1;
    Int_t   iPMD = 0;
    Int_t   iHMPID = 0;
    Int_t   iSHIL = 1;
    Int_t   iT0 = 0;
    Int_t   iTOF = 0;
    Int_t   iTPC = 1;
    Int_t   iTRD = 0;
    Int_t   iZDC = 0;
    Int_t   iEMCAL = 0;

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

        AliFRAME *FRAME = new AliFRAMEv2("FRAME", "Space Frame");

    }

    if (iSHIL)
    {
        //=================== SHIL parameters ============================

        AliSHIL *SHIL = new AliSHILv2("SHIL", "Shielding");
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
    ITS->SetMinorVersion(2); // don't touch this parameter if you're not an ITS developer
    ITS->SetReadDet(kFALSE);  // don't touch this parameter if you're not an ITS developer
    //    ITS->SetWriteDet("$ALICE_ROOT/ITS/ITSgeometry_vPPRasymm2.det");  // don't touch this parameter if you're not an ITS developer
    ITS->SetThicknessDet1(200.);   // detector thickness on layer 1 must be in the range [100,300]
    ITS->SetThicknessDet2(200.);   // detector thickness on layer 2 must be in the range [100,300]
    ITS->SetThicknessChip1(150.);  // chip thickness on layer 1 must be in the range [150,300]
    ITS->SetThicknessChip2(150.);  // chip thickness on layer 2 must be in the range [150,300]
    ITS->SetRails(1);	           // 1 --> rails in ; 0 --> rails out
    ITS->SetCoolingFluid(1);       // 1 --> water ; 0 --> freon
    //
    //
    // Coarse geometries (warning: no hits are produced with these coarse geometries and they unuseful 
    // for reconstruction !):
    //                                                     
    //
    //AliITSvPPRcoarseasymm *ITS  = new AliITSvPPRcoarseasymm("ITS","New ITS PPR coarse version with asymmetric services");
    //ITS->SetRails(1);                // 1 --> rails in ; 0 --> rails out
    //ITS->SetSupportMaterial(0);      // 0 --> Copper ; 1 --> Aluminum ; 2 --> Carbon
    //
    //AliITS *ITS  = new AliITSvPPRcoarsesymm("ITS","New ITS PPR coarse version with symmetric services");
    //ITS->SetRails(1);                // 1 --> rails in ; 0 --> rails out
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

    if (iTOF)
    {
        //=================== TOF parameters ============================
        AliTOF *TOF = new AliTOFv5T0("TOF", "normal TOF");
    }

    if (iHMPID)
    {
        //=================== HMPID parameters ===========================
        AliHMPID *HMPID = new AliHMPIDv1("HMPID", "normal HMPID");

    }


    if (iZDC)
    {
        //=================== ZDC parameters ============================

        AliZDC *ZDC = new AliZDCv1("ZDC", "normal ZDC");
    }

    if (iCRT)
    {
        //=================== CRT parameters ============================

        AliCRT *CRT = new AliCRTv1("CRT", "normal CRT");
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

        AliMUON *MUON = new AliMUONv1("MUON", "default");
    }
    //=================== PHOS parameters ===========================

    if (iPHOS)
    {
        AliPHOS *PHOS = new AliPHOSv1("PHOS", "GPS2");
    }


    if (iPMD)
    {
        //=================== PMD parameters ============================

        AliPMD *PMD = new AliPMDv1("PMD", "normal PMD");

        PMD->SetPAR(1., 1., 0.8, 0.02);
        PMD->SetIN(6., 18., -580., 27., 27.);
        PMD->SetGEO(0.0, 0.2, 4.);
        PMD->SetPadSize(0.8, 1.0, 1.0, 1.5);

    }
    if (iEMCAL && !iHMPID)
    {
        //=================== EMCAL parameters ============================
        AliEMCAL *EMCAL = new AliEMCALv2("EMCAL", "SHISH_77_TRD1_2X2_FINAL_110DEG");
    }

    if (iT0)
    {
        //=================== T0 parameters ============================
        AliT0 *T0 = new AliT0v1("T0", "T0 Detector");
    }


}
