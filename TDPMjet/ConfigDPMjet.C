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
#include "EVGEN/AliGenCocktail.h"
#include "EVGEN/AliGenFixed.h"
#include "EVGEN/AliGenBox.h"
#include "EVGEN/AliGenScan.h"
#include "EVGEN/AliGenHIJINGpara.h"
#include "THijing/AliGenHijing.h"
#include "PYTHIA6/AliGenPythia.h"
#include "THerwig/AliGenHerwig.h"
#include "TIsajet/AliGenIsajet.h"
#include "TDPMjet/AliGenDPMjet.h"
#include "EVGEN/AliGenParam.h"
#include "EVGEN/AliGenMUONlib.h"
#include "EVGEN/AliGenPHOSlib.h"
#include "EVGEN/AliGenGSIlib.h"
#include "EVGEN/AliGenFLUKAsource.h"
#include "EVGEN/AliGenExtFile.h"
#include "EVGEN/AliGenHalo.h"
#include "EVGEN/AliGenReaderTreeK.h"
#include "STRUCT/AliBODY.h"
#include "STRUCT/AliMAG.h"

#endif

enum gentype_t {hijing, hijingParam, gun, box, pythia, herwig, 
		param1, param2, param3, param4, 
		cocktail, fluka, halo, ntuple, scan, isajet, dpmjet};

gentype_t gentype = dpmjet;

Int_t ntracks=1;

void Config()
{

  // Set Random Number seed
  gRandom->SetSeed(12345); //Set 0 to use the current time
  cout<<"Seed for random number generation= "<<gRandom->GetSeed()<<endl; 


  // libraries required by geant321
#if defined(__CINT__)
  gSystem->Load("libgeant321");
#endif
gSystem->Load("libherwig.so");
gSystem->Load("libTHerwig.so");
gSystem->Load("libisajet.so");
gSystem->Load("libTIsajet.so");
gSystem->Load("libdpmjet.so");
gSystem->Load("libTDPMjet.so");

  new TGeant3TGeo("C++ Interface to Geant3");

//=======================================================================
//  Create the output file
   
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
  rl->SetNumberOfEventsPerFile(100);
  gAlice->SetRunLoader(rl);

//
// Set External decayer
  AliDecayer* decayer = new AliDecayerPythia();
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
  

  AliGenerator * gGener = 0x0;
  switch(gentype)
    {
    case gun:
//*********************************************
// Example for Fixed Particle Gun             
//*********************************************
      {
	AliGenFixed *gener = new AliGenFixed(ntracks);
	gener->SetMomentum(50);
	gener->SetPhi(180.);
	gener->SetTheta(5.);
	gener->SetOrigin(0,0,0);        //vertex position
	gener->SetPart(13);             //GEANT particle type
	gGener = gener;
      }
      break;
    case box:  
//*********************************************
// Example for Moving Particle Gun            *
//*********************************************
      {
	AliGenBox *gener = new AliGenBox(ntracks);
	gener->SetMomentumRange(3,4);
	gener->SetPhiRange(0,360);
	gener->SetThetaRange(90, 180. );
	gener->SetOrigin(0,0,0);   
	//vertex position
	gener->SetSigma(0,0,0);         //Sigma in (X,Y,Z) (cm) on IP position
	gener->SetPart(5);              //GEANT particle type
	gGener = gener;
      }
      break;
    case scan:  
//*********************************************
// Scanning on a grid                         *
//*********************************************
      {
	AliGenScan *gener = new AliGenScan(-1);
	gener->SetMomentumRange(4,4);
	gener->SetPhiRange(0,360);
	gener->SetThetaRange(9,9);
	//vertex position
	gener->SetSigma(6,6,0);         //Sigma in (X,Y,Z) (cm) on IP position
	gener->SetPart(5); 
	gener->SetRange(20, -100, 100, 20, -100, 100, 1, 500, 500);
	gGener = gener;
      }
      break;
     
    case hijingParam:
      {
	AliGenHIJINGpara *gener = new AliGenHIJINGpara(ntracks);
	gener->SetMomentumRange(0,999);
	gener->SetPhiRange(0,360);
	gener->SetThetaRange(2,10);
	gener->SetOrigin(0,0,0);        //vertex position
	gener->SetSigma(0,0,0);         //Sigma in (X,Y,Z) (cm) on IP position
	gGener = gener;
      }
      break;
    case hijing:
      {
	AliGenHijing *gener = new AliGenHijing(-1);
// centre of mass energy 
	gener->SetEnergyCMS(5500);
// reference frame
	gener->SetReferenceFrame("CMS     ");
// projectile
	gener->SetProjectile("A       ", 208, 82);
	gener->SetTarget    ("A       ", 208, 82);
// impact parameter range
	gener->SetImpactParameterRange(0, 3.);
// evaluate cross section before run
	gener->SetEvaluate(0);
// tell hijing to keep the full parent child chain
	gener->KeepFullEvent();
// enable jet quenching
	gener->SetJetQuenching(1);
// enable shadowing
	gener->SetShadowing(1);
// neutral pion and heavy particle decays switched off
	gener->SetDecaysOff(1);
// trigger
	gener->SetTrigger(0);
// kinematic selection
	gener->SetSelectAll(0);
// momentum range
	gener->SetMomentumRange(0,999);
// phi range
	gener->SetPhiRange(0.,360.);
// theta range 
	gener->SetThetaRange(0,180.);
// select flavor (0: no, 4: charm+beauty, 5:beauty)
	gener->SetFlavor(0);
//     
	gener->SetOrigin(0., 0.0 ,0);
	gener->SetSigma(0,0,5.3);
	gener->SetVertexSmear(kPerEvent); 
// no tracking
	gener->SetTrackingFlag(0);
	gGener = gener;
      }
      break;
     
    case pythia:
//********************************************
// Example for Charm  Production with Pythia *
//********************************************
      {
	AliGenPythia *gener = new AliGenPythia(-1);
//   final state kinematic cuts
	gener->SetMomentumRange(0,999);
	gener->SetPhiRange(0. ,360.);
	gener->SetThetaRange(0., 180.);
	gener->SetYRange(-10,10);
	gener->SetPtRange(0,100);
//   vertex position and smearing 
	gener->SetOrigin(0,0,0);       // vertex position
	gener->SetVertexSmear(kPerEvent);
	gener->SetSigma(0,0,5.6);      // Sigma in (X,Y,Z) (cm) on IP position
//   Structure function. See the list in EVGEN/AliStructFuncType.h
	gener->SetStrucFunc(kGRVHO);
// Select corection for nuclear structure functions
//     gener->SetNuclei(208,208);
//
//   Process type. See the list in PYTHIA6/AliPythia.h
	gener->SetProcess(kPyBeauty);
//   
//   Pt transfer of the hard scattering
	gener->SetPtHard(0.,5.);
//   Decay type (semielectronic, semimuonic, nodecay)
	gener->SetForceDecay(kSemiElectronic);
//   Centre of mass energy 
	gener->SetEnergyCMS(5500.);
//   No Tracking 
	gener->SetTrackingFlag(0);
	gGener = gener;
      }
      break;              

    case herwig:
//********************************************
// Example for Charm  Production with Pythia *
//********************************************
      {
	AliGenHerwig *gener = new AliGenHerwig(-1);
//   final state kinematic cuts
	gener->SetMomentumRange(0,7000);
	gener->SetPhiRange(0. ,360.);
	gener->SetThetaRange(0., 180.);
	gener->SetYRange(-10,10);
	gener->SetPtRange(0,7000);
//   vertex position and smearing 
	gener->SetOrigin(0,0,0);       // vertex position
	gener->SetVertexSmear(kPerEvent);
	gener->SetSigma(0,0,5.6);      // Sigma in (X,Y,Z) (cm) on IP position
//   Beam momenta
	gener->SetBeamMomenta(7000,7000);
//   Beams
	gener->SetProjectile("P");
	gener->SetTarget("P");
//   Structure function
	gener->SetStrucFunc(kGRVHO);
//   Hard scatering
	gener->SetPtHardMin(200);
	gener->SetPtRMS(20);
//   Min bias
	gener->SetProcess(8000);
//   No Tracking 
	gener->SetTrackingFlag(0);
	gGener = gener;
      }
      break;              
    case isajet:
	AliGenIsajet *gener = new AliGenIsajet(-1);
        gGener = gener;
      break;
   case dpmjet:
	AliGenDPMjet *gener = new AliGenDPMjet(-1);
	//kDpmSingleDiffr, kDpmDoubleDiffr, kDpmDiffr, kDpmNonDiffr
		gener->SetProcess(kDpmMb);
	gGener = gener;
      break;



    case param1:
//*******************************************************
// Example for J/psi  Production from  Parameterisation 
// using default library (AliMUONlib)                                       
//*******************************************************
      {
	AliGenParam *gener =
	  new AliGenParam(ntracks, AliGenMUONlib::kUpsilon);
	gener->SetMomentumRange(0,999);
	gener->SetPtRange(0,999);     
	gener->SetPhiRange(0. , 360.);
	gener->SetYRange(2.5,4);
	gener->SetCutOnChild(1);
	gener->SetChildThetaRange(2,9);
	gener->SetOrigin(0,0,0);        //vertex position
	gener->SetSigma(0,0,5.3);       //Sigma in (X,Y,Z) (cm) on IP position
	gener->SetForceDecay(kDiMuon);
	gener->SetTrackingFlag(0);
	gGener = gener;
      }
      break;

    case param2:
//*******************************************************
// Example for Omega  Production from  Parameterisation 
// specifying library.                                       
//*******************************************************
      {
	AliGenParam *gener = new AliGenParam(1000,new AliGenPHOSlib(), 
					     AliGenPHOSlib::kOmega);
	gener->SetWeighting(kNonAnalog);
	gener->SetForceDecay(kNoDecay);
	gener->SetPtRange(0,100);
	gener->SetThetaRange(45,135);
	gener->SetTrackingFlag(0);
	gGener = gener;
      }
      break;

    case param3:
//*******************************************************
// Example for Upsilon  Production from  Parameterisation 
// specifying library.                                       
// GSI style
//*******************************************************
      {
	AliGenParam *gener = new AliGenParam(1000,new AliGenGSIlib(), 
					     AliGenGSIlib::kUpsilon, "MUON");
	gener->SetMomentumRange(0,999);
	gener->SetPtRange(0,999);     
	gener->SetPhiRange(0., 360.);
	gener->SetYRange(2.5,4);
	gener->SetCutOnChild(1);
	gener->SetChildThetaRange(2,9);
	gener->SetOrigin(0,0,0);        //vertex position
	gener->SetSigma(0,0,5.3);       //Sigma in (X,Y,Z) (cm) on IP position
	gener->SetForceDecay(kDiMuon);
	gener->SetTrackingFlag(0);
	gGener = gener;
      }
      break;
     
    case param4:
//*******************************************************
// Example for Omega  Production from  Parameterisation 
// specifying library.
// The alternative way.                                       
//*******************************************************
      {
	AliGenLib* Lib=new AliGenPHOSlib();
	Int_t iOmega = AliGenPHOSlib::kOmega;
	AliGenParam *gener = new AliGenParam(50, iOmega,            
					     Lib->GetPt(iOmega, ""),
					     Lib->GetY (iOmega, ""),
					     Lib->GetIp(iOmega, ""));
	gener->SetPtRange(0,999);     
	gener->SetWeighting(kNonAnalog);
	gener->SetForceDecay(kNoDecay);
	gener->SetTrackingFlag(0);
	gGener = gener;
      }
      break;
      
    case fluka:
//*******************************************************
// Example for a FLUKA Boundary Source                  *
//*******************************************************
      {
	AliGenFLUKAsource *gener = new AliGenFLUKAsource(-1);
	gener->SetFileName("$(ALICE_ROOT)/data/all32.root"); 
	gener->SetPartFlag(9);
	gener->SetAgeMax(1.e-5);
//  31.7 events     
	gener->SetFraction(0.0315);     
//     gener->SetFraction(0.75*0.0315);     
	rl->CdGAFile();
//     gener->SetPartFlag(10);
	gener->SetMomentumRange(0,999);
	gener->SetPhiRange(0.,360.);
	gener->SetThetaRange(0., 180.); 
	gener->SetAgeMax(1.e-5);
     
//  31.7 events     
//     gener->SetFraction(0.0315);     
	gGener = gener;
      }
      break;

    case ntuple:
//*******************************************************
// Example for reading from a external file                  *
//*******************************************************
      {
	AliGenExtFile *gener = new AliGenExtFile(-1); 
	gener->SetVertexSmear(kPerEvent); 
	gener->SetTrackingFlag(1);
	
	AliGenReaderTreeK * reader = new AliGenReaderTreeK();
	reader->SetFileName("$(ALICE_ROOT)/data/dtujet93.root");
	gener->SetReader(reader);
	gGener = gener;
      }
      break;

    case halo:
//*******************************************************
// Example for Tunnel Halo Source                       *
//*******************************************************
      {
	AliGenHalo *gener = new AliGenHalo(ntracks); 
	gener->SetFileName("/h1/morsch/marsip/marsip5.mu");
	gGener = gener;
      }
      break;
      
    case cocktail:
//*******************************************************
// Example for a Cocktail                               *
//*******************************************************
      {
	AliGenCocktail *gener = new AliGenCocktail(); 

	gener->SetPhiRange(0,360);
	gener->SetYRange(2.5,4);
	gener->SetThetaRange(2,9);
	gener->SetPtRange(0,10);
	gener->SetOrigin(0,0,0);        //vertex position
	gener->SetSigma(0,0,0);         //Sigma in (X,Y,Z) (cm) on IP position
	gener->SetMomentumRange(0,999);

	AliGenParam *jpsi = new AliGenParam(1,AliGenMUONlib::kJpsi);
	jpsi->SetForceDecay(kDiMuon);
	jpsi->SetCutOnChild(1);

     
	AliGenFLUKAsource *bg = new AliGenFLUKAsource(-1);
	bg->AddFile("$(ALICE_ROOT)/data/all32.root"); 
	rl->CdGAFile();
	bg->SetPartFlag(9);
	bg->SetAgeMax(1.e-5);
//  31.7 events     
//     gener->SetFraction(0.0315);     
	bg->SetFraction(0.01*0.0315);     
      
	gener->AddGenerator(jpsi,"J/Psi", 1);
	gener->AddGenerator(bg,"Background",1);

	gGener = gener;
      }
      break;
    }
 
// Activate this line if you want the vertex smearing to happen
// track by track
//
// gener->SetVertexSmear(kPerTrack); 

  gGener->Init();

  gAlice->SetField(-999,2);    //Specify maximum magnetic field in Tesla (neg. ==> default field)

  Int_t iMAG=1;
  rl->CdGAFile();

//=================== Alice BODY parameters =============================
  AliBODY *BODY = new AliBODY("BODY","Alice envelop");


  if(iMAG) {
//=================== MAG parameters ============================
// --- Start with Magnet since detector layouts may be depending ---
// --- on the selected Magnet dimensions ---
    AliMAG *MAG  = new AliMAG("MAG","Magnet");
  }
}
