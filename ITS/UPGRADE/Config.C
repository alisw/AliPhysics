/// One can use the configuration macro in compiled mode by
// root [0] gSystem->Load("libgeant321");
// root [0] gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include\
//                   -I$ALICE_ROOT -I$ALICE/geant3/TGeant3");
// root [0] .x grun.C(1,"Config.C++"

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
#include "STEER/AliMagWrapCheb.h"
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
#include <TVirtualMagField.h>
#endif


Int_t generatorFlag = 2;
Int_t ITSflag = 3;



/* $Id$ */
enum PprTrigConf_t
{
  kDefaultPPTrig, kDefaultPbPbTrig
};

const char * pprTrigConfName[] = {
  "p-p","Pb-Pb"
};

Float_t EtaToTheta(Float_t arg);
AliGenerator *Hijing();

static PprTrigConf_t strig = kDefaultPPTrig;// default PP trigger configuration

void Config()
{
  // ThetaRange is (0., 180.). It was (0.28,179.72) 7/12/00 09:00
  // Theta range given through pseudorapidity limits 22/6/2001

  // Set Random Number seed
  gRandom->SetSeed(0); // Set 0 to use the currecnt time


  // libraries required by geant321
#if defined(__CINT__)
  gSystem->Load("liblhapdf");
  gSystem->Load("libEGPythia6");
  gSystem->Load("libpythia6");
  gSystem->Load("libAliPythia6");
  gSystem->Load("libgeant321");
  gSystem->Load("libhijing");	
  gSystem->Load("libTHijing");
 
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
  rl->SetNumberOfEventsPerFile(2000);
  gAlice->SetRunLoader(rl);

  // Set the trigger configuration
  // gAlice->SetTriggerDescriptor(pprTrigConfName[strig]);
  //cout<<"Trigger configuration is set to  "<<pprTrigConfName[strig]<<endl;
  AliSimulation::Instance()->SetTriggerConfig(pprTrigConfName[strig]);
  cout<<"Trigger configuration is set to  pprTrigConfName[strig] "<<endl;

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
  gMC->SetProcess("HADR",0);
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

  if (generatorFlag==0) {
    
    AliGenCocktail *gener = new AliGenCocktail();
    gener->SetPhiRange(0, 360);
    // Set pseudorapidity range from -8 to 8.
    Float_t thmin = EtaToTheta(8);   // theta min. <---> eta max
    Float_t thmax = EtaToTheta(-8);  // theta max. <---> eta min 
    gener->SetThetaRange(thmin,thmax);
    gener->SetOrigin(0., 0., 0);  //vertex position
    gener->SetSigma(0., 0., 0);   //Sigma in (X,Y,Z) (cm) on IP position
    
  /* 
     AliGenBox *gbox2 = new AliGenBox(2300*1.8);
     gbox2->SetPtRange(0.1,2.0);  //
     gbox2->SetPhiRange(0,360);
     gbox2->SetThetaRange(44.,135.); // +/- 0.9 eta
     //gbox2->SetPhiRange(250,250.5);
     //  gbox2->SetThetaRange(50.,50.1);
     gbox2->SetPart(kPiMinus);
     gener->AddGenerator(gbox2,"GENBOX PIONS for ITS",1);
  */


  
    AliGenHIJINGpara *gbox2=new AliGenHIJINGpara(2300*1.8); // 2300 dNdy in 2*0.9eta
    gbox2->SetPtRange(0.06,9999999); 
    gbox2->SetPhiRange(0,360.);
    gbox2->SetThetaRange(44.,135.); // +/- 0.9 eta
    gener->AddGenerator(gbox2,"GenHIJINGpara PIONS for ITS",1);
    
    gener->Init();


  } else if (generatorFlag==1) {


    AliGenHijing *generHijing = new AliGenHijing(-1);
    // centre of mass energy 
    generHijing->SetEnergyCMS(5500.);
    generHijing->SetImpactParameterRange(0,1);	
    // reference frame
    generHijing->SetReferenceFrame("CMS");
    // projectile
    generHijing->SetProjectile("A", 208, 82);
    generHijing->SetTarget    ("A", 208, 82);
    // tell hijing to keep the full parent child chain
    generHijing->KeepFullEvent();
    // enable jet quenching
    generHijing->SetJetQuenching(1);
    // enable shadowing
    generHijing->SetShadowing(1);
    // neutral pion and heavy particle decays switched off
    //  generHijing->SetDecaysOff(3); // 3 requested by Ana Marin
    // Don't track spectators
    generHijing->SetSpectators(0);
    // kinematic selection
    generHijing->SetSelectAll(0);
 
    AliGenerator*  gener = generHijing;
    gener->SetSigma(0, 0, 0);      // Sigma in (X,Y,Z) (cm) on IP position
    gener->SetVertexSmear(kPerEvent);
    gener->Init();


  } else if (generatorFlag==2) {


    // Francesco ...
    int     nParticles = 14022;
    AliGenHIJINGpara *gener = new AliGenHIJINGpara(nParticles);
    gener->SetMomentumRange(0.1, 10.);
    gener->SetPhiRange(0., 360.);
    Float_t thmin = EtaToTheta(2.5);   // theta min. <---> eta max
    Float_t thmax = EtaToTheta(-2.5);  // theta max. <---> eta min
    gener->SetThetaRange(thmin,thmax);
    gener->SetOrigin(0, 0, 0);  //vertex position
    gener->SetSigma(0, 0, 0);   //Sigma in (X,Y,Z) (cm) on IP position
    gener->Init();


  } else if (generatorFlag==3) {

    // Annalisa ...
    AliGenHijing *generHijing = new AliGenHijing(-1);
    generHijing->SetEnergyCMS(5500.);
    generHijing->SetImpactParameterRange(0,5);
    generHijing->SetReferenceFrame("CMS");
    generHijing->SetProjectile("A", 208, 82);
    generHijing->SetTarget    ("A", 208, 82);
    generHijing->KeepFullEvent();
    generHijing->SetJetQuenching(1);
    generHijing->SetShadowing(1);
    generHijing->SetSpectators(0);
    generHijing->SetSelectAll(0);
    generHijing->SetPtHardMin(4.5);

    AliGenerator*  gener = generHijing;
    gener->SetSigma(0, 0, 0);      // Sigma in (X,Y,Z) (cm) on IP position
    gener->SetVertexSmear(kPerEvent);
    gener->Init();


  }


  // 
  // Activate this line if you want the vertex smearing to happen
  // track by track
  //
  //VertexSmear_t perTrack;
  //gener->SetVertexSmear(perTrack); 
  // Field (L3 0.5 T)
  //AliMagF* field = new AliMagF("map","map",2, -1.,1., 15, AliMagF::k5kGUniform);
  //TGeoGlobalMagField::Instance()->SetField(field);
  TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", -1., -1., AliMagF::k5kG));

  Int_t   iABSO  =  0;
  Int_t   iDIPO  =  0;
  Int_t   iFMD   =  0;
  Int_t   iFRAME =  0;
  Int_t   iHALL  =  0;
  Int_t   iITS   =  1;
  Int_t  isUpgrade = ITSflag;
  Int_t   iMAG   =  0;
  Int_t   iMUON  =  0;
  Int_t   iPHOS  =  0;
  Int_t   iPIPE  =  0;
  Int_t   iPMD   =  0;
  Int_t   iHMPID =  0;
  Int_t   iSHIL  =  0;
  Int_t   iT0    =  0;
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

  
      if (isUpgrade==1) { // current ITS geometry ...


	AliITSupgrade *ITSu  = new AliITSupgrade("ITS","ITS upgrade",kTRUE);
	  

      } else if(isUpgrade==2) { // "New SPDs" Configuration 

	// BeamPipe
	Bool_t bp=kTRUE;
	Double_t widthBP = 0.08;
	Double_t radiusBP = 2.0-widthBP;
	Double_t halflengthBP = 200.;


	// BeamPipe
	Int_t  nlayers =7;
	TArrayD xsize(nlayers); 	TArrayD zsize(nlayers);
	TArrayD width(nlayers); 	TArrayD radii(nlayers);
       
	TArrayD widthCu(nlayers); 	TArrayD radiiCu(nlayers);
	TArrayS copper(nlayers);
	TArrayD halfLengths(nlayers);


	for(Int_t i=0; i<nlayers; i++){

	  Double_t xsz[7]={14*1e-04, 14*1e-04, 14*1e-04, 121.2e-04,121.4*1e-04,69.3*1e-04,69.3*1e-04};
	  Double_t zsz[7]={14*1e-04, 14*1e-04, 14*1e-04,  97*1e-04,97*1e-04,2875*1e-04,2875*1e-04};

	  Double_t halfL[7]={21./2, 25./2., 32./2., 22.2,29.7,43.1,48.9};


	  Double_t r[7]={2.2,3.8,6.8,14.9,23.8,39.1,43.6};
	  radii.AddAt(r[i],i);

	  Double_t thick[7]={50.*1e-04,50.*1e-04,50.*1e-04,150.*1e-04,150.*1e-04,150.*1e-04,150.*1e-04};
	  width.AddAt(thick[i],i);

	  // Copper radius & thickness
	  radiiCu.AddAt(r[i]+thick[i],i);
	  if (i<3) 
	    widthCu.AddAt(0.0036,i);
	  else
	    widthCu.AddAt(0.0150,i);//micron
	  
	  // no idea ????????
	  Int_t c[7]={1,1,1,1,1,1,1};
	  copper.AddAt(c[i],i);


	  // Silicon pixel sizes and length in Z

	  Int_t npixHalf[7];
	  npixHalf[i]=(Int_t)(halfL[i]/zsz[i]); 
	  Double_t HalfL[7];
	  HalfL[i]=npixHalf[i]*zsz[i];
	
	  Int_t npixR[7];
	  npixR[i] = (Int_t)(2*TMath::Pi()*r[i]/xsz[i]);
	  Double_t xszInt[7];
	  xszInt[i]= 2*TMath::Pi()*r[i]/npixR[i];
	  xsize.AddAt(xszInt[i],i);
	  zsize.AddAt(zsz[i],i);

	  halfLengths.AddAt(HalfL[i],i);


	}
	AliITSupgrade *ITSu  = 
	  new AliITSupgrade("ITS","ITS upgrade",
			    width,radii,halfLengths,
			    radiiCu,widthCu,copper,
			    bp,radiusBP, widthBP, halflengthBP );
     

	ITSu->SetFullSegmentation(xsize,zsize);

      } else if (isUpgrade==3) { // "All New" Configuration 

	// BeamPipe
	Bool_t bp=kTRUE;
	Double_t widthBP = 0.08;
	Double_t radiusBP = 2.0-widthBP;
	Double_t halflengthBP = 200.;


	// BeamPipe
	Int_t  nlayers =7;
	TArrayD xsize(nlayers); 	TArrayD zsize(nlayers);
	TArrayD width(nlayers); 	TArrayD radii(nlayers);
       
	TArrayD widthCu(nlayers); 	TArrayD radiiCu(nlayers);
	TArrayS copper(nlayers);
	TArrayD halfLengths(nlayers);


	for(Int_t i=0; i<nlayers; i++){

	  Double_t xsz[7]={14*1e-04, 14*1e-04, 14*1e-04, 14*1e-04, 14*1e-04, 14*1e-04, 14*1e-04};
	  Double_t zsz[7]={14*1e-04, 14*1e-04, 14*1e-04, 14*1e-04, 14*1e-04, 14*1e-04, 14*1e-04};

	  Double_t halfL[7]={21./2, 25./2., 32./2., 45./2, 67./2, 107./2, 114./2};


	  Double_t r[7]={2.2,3.8,6.8,12.4,23.5,39.6,43.0};
	  radii.AddAt(r[i],i);

	  Double_t thick[7]={50.*1e-04,50.*1e-04,50.*1e-04,50.*1e-04,50.*1e-04,50.*1e-04,50.*1e-04};
	  width.AddAt(thick[i],i);

	  // Copper radius & thickness
	  radiiCu.AddAt(r[i]+thick[i],i);
	  widthCu.AddAt(0.0036,i);
	  
	  // no idea ????????
	  Int_t c[7]={1,1,1,1,1,1,1};
	  copper.AddAt(c[i],i);


	  // Silicon pixel sizes and length in Z

	  Int_t npixHalf[7];
	  npixHalf[i]=(Int_t)(halfL[i]/zsz[i]); 
	  Double_t HalfL[7];
	  HalfL[i]=npixHalf[i]*zsz[i];
	
	  Int_t npixR[7];
	  npixR[i] = (Int_t)(2*TMath::Pi()*r[i]/xsz[i]);
	  Double_t xszInt[7];
	  xszInt[i]= 2*TMath::Pi()*r[i]/npixR[i];
	  xsize.AddAt(xszInt[i],i);
	  zsize.AddAt(zsz[i],i);

	  halfLengths.AddAt(HalfL[i],i);


	}
	AliITSupgrade *ITSu  = 
	  new AliITSupgrade("ITS","ITS upgrade",
			    width,radii,halfLengths,
			    radiiCu,widthCu,copper,
			    bp,radiusBP, widthBP, halflengthBP );
     

	ITSu->SetFullSegmentation(xsize,zsize);

      } else if (isUpgrade==4) { // "All New" Configuration with 8 layers

	// BeamPipe
	Bool_t bp=kTRUE;
	Double_t widthBP = 0.08;
	Double_t radiusBP = 2.0-widthBP;
	Double_t halflengthBP = 200.;


	// BeamPipe
	Int_t  nlayers =8;
	TArrayD xsize(nlayers); 	TArrayD zsize(nlayers);
	TArrayD width(nlayers); 	TArrayD radii(nlayers);
       
	TArrayD widthCu(nlayers); 	TArrayD radiiCu(nlayers);
	TArrayS copper(nlayers);
	TArrayD halfLengths(nlayers);


	for(Int_t i=0; i<nlayers; i++){

	  Double_t xsz[8]={14*1e-04, 14*1e-04, 14*1e-04, 14*1e-04, 14*1e-04, 14*1e-04, 14*1e-04, 14*1e-04};
	  Double_t zsz[8]={14*1e-04, 14*1e-04, 14*1e-04, 14*1e-04, 14*1e-04, 14*1e-04, 14*1e-04, 14*1e-04};

	  Double_t halfL[8]={21./2, 25./2., 32./2., 45./2, 67./2, 114./2, 114./2, 114./2};


	  Double_t r[8]={2.2,3.8,6.8,12.4,23.5, 42.6,43.0,43.4};
	  radii.AddAt(r[i],i);

	  Double_t thick[8]={50.*1e-04,50.*1e-04,50.*1e-04,50.*1e-04,50.*1e-04,50.*1e-04,50.*1e-04,50.*1e-04};
	  width.AddAt(thick[i],i);

	  // Copper radius & thickness
	  radiiCu.AddAt(r[i]+thick[i],i);
	  widthCu.AddAt(0.0036,i);
	  
	  // no idea ????????
	  Int_t c[8]={1,1,1,1,1,1,1};
	  copper.AddAt(c[i],i);


	  // Silicon pixel sizes and length in Z

	  Int_t npixHalf[8];
	  npixHalf[i]=(Int_t)(halfL[i]/zsz[i]); 
	  Double_t HalfL[8];
	  HalfL[i]=npixHalf[i]*zsz[i];
	
	  Int_t npixR[8];
	  npixR[i] = (Int_t)(2*TMath::Pi()*r[i]/xsz[i]);
	  Double_t xszInt[8];
	  xszInt[i]= 2*TMath::Pi()*r[i]/npixR[i];
	  xsize.AddAt(xszInt[i],i);
	  zsize.AddAt(zsz[i],i);

	  halfLengths.AddAt(HalfL[i],i);


	}
	AliITSupgrade *ITSu  = 
	  new AliITSupgrade("ITS","ITS upgrade",
			    width,radii,halfLengths,
			    radiiCu,widthCu,copper,
			    bp,radiusBP, widthBP, halflengthBP );
     

	ITSu->SetFullSegmentation(xsize,zsize);

      } else if (isUpgrade==5) { // "All New" Configuration with 8 layers

	// BeamPipe
	Bool_t bp=kTRUE;
	Double_t widthBP = 0.08;
	Double_t radiusBP = 2.0-widthBP;
	Double_t halflengthBP = 200.;


	// BeamPipe
	Int_t  nlayers =9;
	TArrayD xsize(nlayers); 	TArrayD zsize(nlayers);
	TArrayD width(nlayers); 	TArrayD radii(nlayers);
       
	TArrayD widthCu(nlayers); 	TArrayD radiiCu(nlayers);
	TArrayS copper(nlayers);
	TArrayD halfLengths(nlayers);


	for(Int_t i=0; i<nlayers; i++){

	  Double_t xsz[9]={14*1e-04, 14*1e-04, 14*1e-04, 14*1e-04, 14*1e-04, 14*1e-04, 14*1e-04, 14*1e-04, 14*1e-04};
	  Double_t zsz[9]={14*1e-04, 14*1e-04, 14*1e-04, 14*1e-04, 14*1e-04, 14*1e-04, 14*1e-04, 14*1e-04, 14*1e-04};

	  Double_t halfL[9]={21./2, 25./2., 32./2., 45./2, 67./2, 107./2, 107./2, 114./2, 114./2};


	  Double_t r[9]={2.2,3.8,6.8,12.4,23.5, 39.6,40.0,43.0,43.4};
	  radii.AddAt(r[i],i);

	  Double_t thick[9]={50.*1e-04,50.*1e-04,50.*1e-04,50.*1e-04,50.*1e-04,50.*1e-04,50.*1e-04,50.*1e-04,50.*1e-04};
	  width.AddAt(thick[i],i);

	  // Copper radius & thickness
	  radiiCu.AddAt(r[i]+thick[i],i);
	  widthCu.AddAt(0.0036,i);
	  
	  // no idea ????????
	  Int_t c[9]={1,1,1,1,1,1,1};
	  copper.AddAt(c[i],i);


	  // Silicon pixel sizes and length in Z

	  Int_t npixHalf[9];
	  npixHalf[i]=(Int_t)(halfL[i]/zsz[i]); 
	  Double_t HalfL[9];
	  HalfL[i]=npixHalf[i]*zsz[i];
	
	  Int_t npixR[9];
	  npixR[i] = (Int_t)(2*TMath::Pi()*r[i]/xsz[i]);
	  Double_t xszInt[9];
	  xszInt[i]= 2*TMath::Pi()*r[i]/npixR[i];
	  xsize.AddAt(xszInt[i],i);
	  zsize.AddAt(zsz[i],i);

	  halfLengths.AddAt(HalfL[i],i);


	}
	AliITSupgrade *ITSu  = 
	  new AliITSupgrade("ITS","ITS upgrade",
			    width,radii,halfLengths,
			    radiiCu,widthCu,copper,
			    bp,radiusBP, widthBP, halflengthBP );
     

	ITSu->SetFullSegmentation(xsize,zsize);


      } else {

      } else {
	AliITS *ITS =  new AliITSv11Hybrid("ITS","ITS v11Hybrid");
      }
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
      // Partial geometry: modules at 2,3,4,6,11,12,14,15
      // starting at 6h in positive direction
      geoTRD->SetSMstatus(0,0);
      geoTRD->SetSMstatus(1,0);
      geoTRD->SetSMstatus(5,0);
      geoTRD->SetSMstatus(7,0);
      geoTRD->SetSMstatus(8,0);
      geoTRD->SetSMstatus(9,0);
      geoTRD->SetSMstatus(10,0);
      geoTRD->SetSMstatus(13,0);
      geoTRD->SetSMstatus(16,0);
      geoTRD->SetSMstatus(17,0);
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
      AliEMCAL *EMCAL = new AliEMCALv2("EMCAL", "EMCAL_COMPLETEV1");
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


AliGenerator* Hijing()
{
     return gener;
}

