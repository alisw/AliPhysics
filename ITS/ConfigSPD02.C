//
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TPDGCode.h>
#include <TRandom.h>
#include <TSystem.h>
#include <TVirtualMC.h>
#include <TGeant3.h>
#include "STEER/AliRunLoader.h"
#include "STEER/AliRun.h"
#include "STEER/AliConfig.h"
#include "PYTHIA6/AliDecayerPythia.h"
#include "EVGEN/AliGenCocktail.h"
#include "EVGEN/AliGenHIJINGpara.h"
#include "STRUCT/AliBODY.h"
#include "STRUCT/AliABSOv0.h"
#include "STRUCT/AliDIPOv2.h"
#include "STRUCT/AliHALL.h"
#include "STRUCT/AliFRAMEv2.h"
#include "STRUCT/AliSHILv2.h"
#include "STRUCT/AliPIPEv0.h"
#include "ITS/AliITSvSPD02.h"
#endif

Float_t EtaToTheta(Float_t arg){
  return (180./TMath::Pi())*2.*atan(exp(-arg));
}
//------------------------------------------------------------------
void Config(){
    // Set Random Number seed
    // gRandom->SetSeed(12345);
    // libraries required by geant321
#if defined(__CINT__)
    gSystem->Load("libgeant321");
#endif
    new TGeant3("C++ Interface to Geant3");
    AliRunLoader *rl = 0;
    rl = AliRunLoader::Open("galice.root",
			    AliConfig::GetDefaultEventFolderName(),"recreate");
    if (rl == 0x0){
      gAlice->Fatal("Config.C","Can not instatiate the Run Loader");
      return;
    } // end if rl==0x0
    rl->SetCompressionLevel(2);
    rl->SetNumberOfEventsPerFile(1000);
    gAlice->SetRunLoader(rl);
    //
    // Set External decayer
    TVirtualMCDecayer *decayer = new AliDecayerPythia();
    decayer->SetForceDecay(kAll);
    decayer->Init();
    gMC->SetExternalDecayer(decayer);
    //=======================================================================
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
    int     nParticles = 1;
    if (gSystem->Getenv("CONFIG_NPARTICLES")){
      nParticles = atoi(gSystem->Getenv("CONFIG_NPARTICLES"));
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
	AliITSvSPD02 *ITS  = new AliITSvSPD02("SPD test beam 2002",2002);
    }
    return;
}
