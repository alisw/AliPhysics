/***************************************************************************
 *  Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                         *
 *  Author: ALICE OFFLICE.                                                 *
 *  Contributors are mentioned in the code where appropriate.              *
 *                                                                         *
 *  Permission to use, copy, modify and distribute this software and its   *
 *  documentation strictly for non-commercial purposes is hereby granted   *
 *  without fee, provided that the above copyright notice appears in all   *
 *  copies and that both the copyright notice and this permission notice   *
 *  appear in the supporting documentation. The authors make no claims     *
 *  about the suitability of this software for any purpose. It is          *
 *  provided "as is" without express or implied warranty.                  *
 *                                                                         *
 ***************************************************************************
     $Satyajit Jena || alien:sjena Sun Apr 21 14:05:19 CEST 2013$
                                                                          
       Sterring macro for fast production of MC events - Followed from 
       the original macro of fastGenAmpt.C     
     
     Implemented Generators: (Version 1.0: Sun Apr 21 14:05:19 CEST 2013)
     ----------------------------------------------------------------
      kPythia6,            kPythia8,               kPythia6D6T,       
      kPythiaPerugia0,     kPythia6ATLAS,          
      kPythiaJets,         
      kPhojet,                
      kDPMjet,             kDPMjet_pA, 
      kAmptDefault,        kAmptStringMelting,          kAmptStringMeltingNoART,
      kAmptpA,                kAmptFlow,
      kAmptReducedFlow,

     FIXME: 
     kPythia6ATLAS_Flat, 
     kHijing,             kHijing2000,            kHijing_pA,             

                                                                
 ***************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TStopwatch.h>
#include <TDatime.h>
#include <TRandom.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TParticle.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>
#include <TRandom3.h>
#include <TChain.h>
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliHeader.h"
#include "AliStack.h"
#include "AliPDG.h"
#include "AliGenAmpt.h"
#include "TAmpt.h"
#include "PYTHIA6/AliDecayerPythia.h"
#include <TVirtualMC.h>
#include <TF1.h>
#include "STEER/AliConfig.h"
#include "STEER/AliGenerator.h"
#include "STEER/AliLog.h"
#include "EVGEN/AliGenHIJINGpara.h"
#include "THijing/AliGenHijing.h"
#include "EVGEN/AliGenCocktail.h"
#include "EVGEN/AliGenSlowNucleons.h"
#include "EVGEN/AliSlowNucleonModelExp.h"
#include "EVGEN/AliGenParam.h"
#include "EVGEN/AliGenMUONlib.h"
#include "EVGEN/AliGenSTRANGElib.h"
#include "EVGEN/AliGenMUONCocktail.h"
#include "EVGEN/AliGenCocktail.h"
#include "EVGEN/AliGenGeVSim.h"
#include "EVGEN/AliGeVSimParticle.h"
#include "PYTHIA6/AliGenPythia.h"
#endif

//___________________________________________________________________________
enum PDC06Proc_t { 
  kPythia6,
  kPythia8, 
  kPythia6D6T, 
  kPythiaPerugia0, 
  kPythia6ATLAS, 
  kPythia6ATLAS_Flat, 
  kPythiaJets,
  kPhojet, 
  kHijing, 
  kHijing2000, 
  kHijing_pA,
  kDPMjet, 
  kDPMjet_pA,
  kAmptDefault,
  kAmpt, 
  kAmptpA,
  kAmptFlowStringMelting,
  kAmptStringMeltingNoART,
  kAmptFlow,
  kAmptReducedFlow,
  kRunMax
};

//___________________________________________________________________________
const char * pprRunName[] = {
  "kPythia6", 
  "kPythia8",
  "kPythia6D6T", 
  "kPythiaPerugia0", 
  "kPythia6ATLAS", 
  "kPythia6ATLAS_Flat", 
  "kPythiaJets", 
  "kPhojet",  
  "kHijing", 
  "kHijing2000", 
  "kHijing_pA",
  "kDPMjet", 
  "kDPMjet_pA", 
  "kAmpt",
  "kAmptpA", 
  "kAmptFlow",
  "kAmptReducedFlow"  
};

enum PprTrigConf_t {kDefaultPPTrig, kDefaultPbPbTrig };
const char * pprTrigConfName[] = {"p-p","Pb-Pb"};

//___________________________________________________________________________
void ProcessEnvironmentVars();
class AliGenPythia;
AliGenerator *MbPythia();
AliGenerator *Pythia8();
AliGenerator *MbPythiaTuneD6T();
AliGenerator *MbPythiaTunePerugia0();
AliGenerator *MbPythiaTuneATLAS();
AliGenerator *MbPythiaTuneATLAS_Flat();
AliGenerator *PythiaJets();
AliGenerator *MbPhojet();
AliGenerator *Hijing();
AliGenerator *Hijing2000();
AliGenerator *Hijing_pA(Bool_t kSlowN);
AliGenerator *DPMjet();
AliGenerator *DPMjet_pA(Bool_t fragments);
AliGenerator *Ampt();
AliGenerator *AmptpA();
AliGenerator* AmptFlow();
AliGenerator *AmptReducedFlow();


//_________________________________________________________________________
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// Geterator, field, beam energy

//static Double_t    pBeamEnergy = 4000.0;  // Used during pA runs
//static Double_t  energy  = 2.*pBeamEnergy*TMath::Sqrt(82./208.); //energy in CMS 

static Double_t      energy    = 7000.;
static PDC06Proc_t   proc      = kPythia6;
static Float_t       bMin      = 0.;
static Float_t       bMax      = 100.;
static PprTrigConf_t strig     = kDefaultPPTrig; // default pp trigger configuration


static Double_t  JpsiPol      = 0; // Jpsi polarisation
static Bool_t    JpsiHarderPt = kFALSE; // Jpsi harder pt spectrum (8.8 TeV)
static TString comment;
//static PprTrigConf_t strig    = kDefaultPbPbTrig; // default pp trigger configuration
TDatime dt; 
static UInt_t seed    = dt.Get();

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//_________________________________________________________________________
void fastMcProduction(Int_t nev = 300) {

  ProcessEnvironmentVars();

  gRandom->SetSeed(seed);
  cerr<<" Seed for random number generation = "<<seed<<endl; 
  
#if defined(__CINT__)
  gSystem->Load("liblhapdf");  
  gSystem->Load("libEGPythia6"); 
  
  if (proc == kPythia6 || proc == kPhojet || proc == kDPMjet || proc==kDPMjet_pA) {
    gSystem->Load("libpythia6");        // Pythia 6.2
    gSystem->Load("libAliPythia6");     // ALICE specific implementations
  }
    
  if (proc == kHijing || proc == kHijing2000 || proc == kHijing_pA ) {
    gSystem->Load("libHIJING");	
    gSystem->Load("libTHijing");
  } 
  
  else if ( proc == kDPMjet || proc== kDPMjet_pA ) {
    gSystem->Load("libDPMJET"); 
    gSystem->Load("libTDPMjet");
  } 
  
  else if (proc == kAmptDefault || kAmptFlowStringMelting || proc ==  kAmptStringMeltingNoART || proc == kAmptpA || proc == kAmptReducedFlow) {
    gSystem->Load("libampt");  
    gSystem->Load("libTAmpt");
    gSystem->Load("libpythia6");
    gSystem->Load("libAliPythia6");
  } 

  if (proc == kPythia8) {
    gSystem->Load("libpythia8");
    gSystem->Load("libAliPythia8");
    gSystem->Setenv("PYTHIA8DATA", gSystem->ExpandPathName("$ALICE_ROOT/PYTHIA8/pythia8145/xmldoc"));
    gSystem->Setenv("LHAPDF",      gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF"));
    gSystem->Setenv("LHAPATH",     gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF/PDFsets"));
  }
#endif


 AliGenerator* gener = 0x0;
 
 cout<<"Run type set to ------------- "<<pprRunName[proc]<<"   " << proc << "    " << kDPMjet_pA<< endl;

 if (proc == kPythia6) {
   Printf("<<<<<<<<<<<<<<<<<<< Processing Mb. PYTHIA >>>>>>>>>>>>>>>>>>>>"); 
   gener = MbPythia();
 } 
 
 else if (proc == kPythia8) {
   Printf("<<<<<<<<<<<<<<<<<<< Processing Mb. Pythia 8 >>>>>>>>>>>>>>>>>>>>"); 
   gener = Pythia8();
 }  
 
 else if (proc == kPythia6D6T) {
   Printf("<<<<<<<<<<<<<<<<<<< Processing Mb. PYTHIA D6T >>>>>>>>>>>>>>>>>>>>"); 
   gener = MbPythiaTuneD6T();
 } 

 else if (proc == kPythiaPerugia0) {
   Printf("<<<<<<<<<<<<<<<<<<< Processing Mb. PYTHIA Perugia0 >>>>>>>>>>>>>>>>>>>>"); 
   gener = MbPythiaTunePerugia0();
 } 
 
 else if (proc == kPythia6ATLAS) {
   Printf("<<<<<<<<<<<<<<<<<<< Processing Mb. PYTHIA ATLAS>>>>>>>>>>>>>>>>>>>>"); 
   gener = MbPythiaTuneATLAS();
 } 
 
 else if (proc == kPythia6ATLAS_Flat) {
   Printf("<<<<<<<<<<<<<<<<<<< Processing Mb. PYTHIA ATLAS_FLAT >>>>>>>>>>>>>>>>>>>>"); 
   gener = MbPythiaTuneATLAS_Flat();
 } 
 
 else if (proc == kPythiaJets ) {
   Printf("<<<<<<<<<<<<<<<<<<< Processing Pythia Jets >>>>>>>>>>>>>>>>>>>>"); 
   gener = PythiaJets();
 } 
 
 else if (proc == kPhojet) {
   Printf("<<<<<<<<<<<<<<<<<<< Processing Mb. PHOJET >>>>>>>>>>>>>>>>>>>>"); 
   gener = MbPhojet();
 } 

 else if (proc == kHijing) {
   Printf("<<<<<<<<<<<<<<<<<<< Processing Mb. HIJING >>>>>>>>>>>>>>>>>>>>"); 
   gener = Hijing();	
 } 
 
 else if (proc == kHijing2000) {
   Printf("<<<<<<<<<<<<<<<<<<< Processing Mb. HIJING 2000 >>>>>>>>>>>>>>>>>>>>"); 
   gener = Hijing2000();	
 }
 
 else if (proc ==kHijing_pA) {
   Printf("<<<<<<<<<<<<<<<<<<< Processing Mb. pA Hijing >>>>>>>>>>>>>>>>>>>>"); 
   gener = Hijing_pA(kTRUE);
 } 
 
 else if (proc == kDPMjet) {
   Printf("<<<<<<<<<<<<<<<<<<< Processing  DMPJet  >>>>>>>>>>>>>>>>>>>>");
   gener = DPMjet();	
 } 

 else if (proc == kDPMjet_pA) {
   Printf("<<<<<<<<<<<<<<<<<<< Processing  DMPJet pA >>>>>>>>>>>>>>>>>>>>");
   gener = DPMjet_pA(kFALSE);	
 } 
 
 else if (proc == kAmptDefault) {
   Printf("<<<<<<<<<<<<<<<<<<< Processing  AMPT Default >>>>>>>>>>>>>>>>>>>>");
   gener = AmptDefault();
 } 

 else if (proc == kAmptStringMelting) {
     Printf("<<<<<<<<<<<<<<<<<<< Processing  AMPT With Flow  >>>>>>>>>>>>>>>>>>>>");
     gener = AmptStringMelting();
 }

 else if (proc == kAmptStringMeltingNoART) {
     Printf("<<<<<<<<<<<<<<<<<<< Processing  AMPT With Flow  >>>>>>>>>>>>>>>>>>>>");
     gener = AmptStringMeltingNoART();
 }

 else if (proc == kAmptpA) {
   Printf("<<<<<<<<<<<<<<<<<<< Processing  AMPT pA  >>>>>>>>>>>>>>>>>>>>");
   gener = AmptpA();
 } 
 

 else if (proc == kAmptReducedFlow) {
   // Specific Fastgen
   Printf("<<<<<<<<<<<<<<<<<<< Processing  AMPT With Reduced Flow >>>>>>>>>>>>>>>>>>>>");
   gener = AmptReducedFlow();
 } 

 else {
   cout << "ERROR : Wrong Procss Selcted !!!" << endl;
   return;
 }


 AliPDG::AddParticlesToPdgDataBase();
 TDatabasePDG::Instance();
 
 const char* filename = "galice.root";
 AliRunLoader* rl = AliRunLoader::Open(filename,"FASTRUN","recreate");
 
 rl->SetCompressionLevel(2);
 rl->SetNumberOfEventsPerFile(nev);
 rl->LoadKinematics("RECREATE");
 rl->MakeTree("E");
 gAlice->SetRunLoader(rl);
 rl->MakeStack();
 AliStack* stack = rl->Stack();
 
 AliHeader* header = rl->GetHeader();
 
 /*
   Float_t sigmaz  = 5.4 / TMath::Sqrt(2.); // [cm]
   Float_t betast  = 3.5;                      // beta* [m]
   Float_t eps     = 3.75e-6;                   // emittance [m]
   Float_t gamma   = energy / 2.0 / 0.938272;  // relativistic gamma [1]
   Float_t sigmaxy = TMath::Sqrt(eps * betast / gamma) / TMath::Sqrt(2.) * 100.;  // [cm]
   printf("\n \n Diamond size x-y: %10.3e z: %10.3e\n \n", sigmaxy, sigmaz);
   gener->SetSigma(sigmaxy, sigmaxy, sigmaz);      // Sigma in (X,Y,Z) (cm) on IP position
   gener->SetVertexSmear(kPerEvent);
 */
 

 gener->Init();
 gener->SetStack(stack);
 
 rl->CdGAFile();
 
 TStopwatch timer;
 timer.Start();
 for (Int_t iev = 0; iev < nev; iev++) {
   cout <<"============================================= Event number "<< iev << endl;
   //  Initialize event
   header->Reset(0,iev);
   rl->SetEventNumber(iev);
   stack->Reset();
   rl->MakeTree("K");
   Int_t nprim = 0;
   Int_t ntrial = 0;
   //  Int_t ndstar = 0;
   stack->Reset();
   stack->ConnectTree(rl->TreeK());
   gener->Generate();
   ntrial++;
   nprim = stack->GetNprimary();
   cout << "Number of particles " << nprim << endl;
   cout << "Number of trials " << ntrial << endl;
   header->SetNprimary(stack->GetNprimary());
   header->SetNtrack(stack->GetNtrack());  
   stack->FinishEvent();
    header->SetStack(stack);
    rl->TreeE()->Fill();
    rl->WriteKinematics("OVERWRITE");
 } // event loop
  timer.Stop();
  timer.Print();
  gener->FinishRun();
  rl->WriteHeader("OVERWRITE");
  gener->Write();
  rl->Write();
}

//___________________________________________________//
void ProcessEnvironmentVars() {
    // Run type
    if (gSystem->Getenv("CONFIG_RUN_TYPE")) {
      for (Int_t iRun = 0; iRun < kRunMax; iRun++) {
        if (strcmp(gSystem->Getenv("CONFIG_RUN_TYPE"), pprRunName[iRun]) == 0) {
          proc = (PDC06Proc_t)iRun;
          cout<<"Run type set to "<<pprRunName[iRun]<<endl;
        }
      }
    }

    
    // Energy
    if (gSystem->Getenv("CONFIG_ENERGY")) {
      energy = atoi(gSystem->Getenv("CONFIG_ENERGY"));
      cout<<"Energy set to "<<energy<<" GeV"<<endl;
    }

    // Random Number seed
    if (gSystem->Getenv("CONFIG_SEED")) {
      seed = atoi(gSystem->Getenv("CONFIG_SEED"));
    }

    // Impact param
    if (gSystem->Getenv("CONFIG_BMIN")) {
      bMin = atof(gSystem->Getenv("CONFIG_BMIN"));
    }

    if (gSystem->Getenv("CONFIG_BMAX")) {
      bMax = atof(gSystem->Getenv("CONFIG_BMAX"));
    }
    cout<<"Impact parameter in ["<<bMin<<","<<bMax<<"]"<<endl;
}



//______________________________________________________________________
AliGenerator* MbPythia() // Mb Pythia
{
      comment = comment.Append(" pp: Pythia low-pt");
      AliGenPythia* pythia = new AliGenPythia(-1); 
      /* pythia->SetMomentumRange(0, 999999.);
      pythia->SetThetaRange(0., 180.);
      pythia->SetYRange(-12.,12.);
      pythia->SetPtRange(0,1000.);*/
      pythia->SetProcess(kPyMb);
      pythia->SetEnergyCMS(energy);
      
      return pythia;
}


//______________________________________________________________________
AliGenerator* Pythia8()
{
  AliGenPythiaPlus* gener = new AliGenPythiaPlus(AliPythia8::Instance());
  gener->SetProcess(kPyMbDefault);
  gener->SetEnergyCMS(energy);
  gener->SetEventListRange(-1, 2);
  return gener;
}



//______________________________________________________________________
AliGenerator* MbPythiaTuneD6T()
{
      comment = comment.Append(" pp: Pythia low-pt");
      AliGenPythia* pythia = new AliGenPythia(-1); 
      pythia->SetMomentumRange(0, 999999.);
      pythia->SetThetaRange(0., 180.);
      pythia->SetYRange(-12.,12.);
      pythia->SetPtRange(0,1000.);
      pythia->SetProcess(kPyMb);
      pythia->SetEnergyCMS(energy);
//    Tune
//    109     D6T : Rick Field's CDF Tune D6T (NB: needs CTEQ6L pdfs externally)
      pythia->SetTune(109); // F I X 
      pythia->SetStrucFunc(kCTEQ6l);
//
      return pythia;
}

//______________________________________________________________________
AliGenerator* MbPythiaTunePerugia0()
{
      comment = comment.Append(" pp: Pythia low-pt (Perugia0)");
//
//    Pythia
      AliGenPythia* pythia = new AliGenPythia(-1); 
      /* pythia->SetMomentumRange(0, 999999.);
      pythia->SetThetaRange(0., 180.);
      pythia->SetYRange(-12.,12.);
      pythia->SetPtRange(0,1000.);*/
      pythia->SetProcess(kPyMb);
      pythia->SetEnergyCMS(energy);
//    Tune
//    320     Perugia 0
      pythia->SetTune(320); 
      pythia->UseNewMultipleInteractionsScenario();
//
      return pythia;
}

//______________________________________________________________________
AliGenerator* MbPythiaTuneATLAS()
{
      comment = comment.Append(" pp: Pythia low-pt");
//
//    Pythia
      AliGenPythia* pythia = new AliGenPythia(-1); 
      /*   pythia->SetMomentumRange(0, 999999.);
      pythia->SetThetaRange(0., 180.);
      pythia->SetYRange(-12.,12.);
      pythia->SetPtRange(0,1000.);*/
      pythia->SetProcess(kPyMb);
      pythia->SetEnergyCMS(energy);
//    Tune
//    C   306 ATLAS-CSC: Arthur Moraes' (new) ATLAS tune (needs CTEQ6L externally)
      pythia->SetTune(306);
      pythia->SetStrucFunc(kCTEQ6l);
//
      return pythia;
}


//______________________________________________________________________
AliGenerator* MbPythiaTuneATLAS_Flat()
{
  AliGenPythia* pythia = MbPythiaTuneATLAS();
  
  comment = comment.Append("; flat multiplicity distribution");
  
  // set high multiplicity trigger
  // this weight achieves a flat multiplicity distribution
  TH1 *weight = new TH1D("weight","weight",201,-0.5,200.5);
  weight->SetBinContent(1,5.49443);
  weight->SetBinContent(2,8.770816);
  weight->SetBinContent(6,0.4568624);
  weight->SetBinContent(7,0.2919915);
  weight->SetBinContent(8,0.6674189);
  weight->SetBinContent(9,0.364737);
  weight->SetBinContent(10,0.8818444);
  weight->SetBinContent(11,0.531885);
  weight->SetBinContent(12,1.035197);
  weight->SetBinContent(13,0.9394057);
  weight->SetBinContent(14,0.9643193);
  weight->SetBinContent(15,0.94543);
  weight->SetBinContent(16,0.9426507);
  weight->SetBinContent(17,0.9423649);
  weight->SetBinContent(18,0.789456);
  weight->SetBinContent(19,1.149026);
  weight->SetBinContent(20,1.100491);
  weight->SetBinContent(21,0.6350525);
  weight->SetBinContent(22,1.351941);
  weight->SetBinContent(23,0.03233504);
  weight->SetBinContent(24,0.9574557);
  weight->SetBinContent(25,0.868133);
  weight->SetBinContent(26,1.030998);
  weight->SetBinContent(27,1.08897);
  weight->SetBinContent(28,1.251382);
  weight->SetBinContent(29,0.1391099);
  weight->SetBinContent(30,1.192876);
  weight->SetBinContent(31,0.448944);
  weight->SetBinContent(32,1);
  weight->SetBinContent(33,1);
  weight->SetBinContent(34,1);
  weight->SetBinContent(35,1);
  weight->SetBinContent(36,0.9999997);
  weight->SetBinContent(37,0.9999997);
  weight->SetBinContent(38,0.9999996);
  weight->SetBinContent(39,0.9999996);
  weight->SetBinContent(40,0.9999995);
  weight->SetBinContent(41,0.9999993);
  weight->SetBinContent(42,1);
  weight->SetBinContent(43,1);
  weight->SetBinContent(44,1);
  weight->SetBinContent(45,1);
  weight->SetBinContent(46,1);
  weight->SetBinContent(47,0.9999999);
  weight->SetBinContent(48,0.9999998);
  weight->SetBinContent(49,0.9999998);
  weight->SetBinContent(50,0.9999999);
  weight->SetBinContent(51,0.9999999);
  weight->SetBinContent(52,0.9999999);
  weight->SetBinContent(53,0.9999999);
  weight->SetBinContent(54,0.9999998);
  weight->SetBinContent(55,0.9999998);
  weight->SetBinContent(56,0.9999998);
  weight->SetBinContent(57,0.9999997);
  weight->SetBinContent(58,0.9999996);
  weight->SetBinContent(59,0.9999995);
  weight->SetBinContent(60,1);
  weight->SetBinContent(61,1);
  weight->SetBinContent(62,1);
  weight->SetBinContent(63,1);
  weight->SetBinContent(64,1);
  weight->SetBinContent(65,0.9999999);
  weight->SetBinContent(66,0.9999998);
  weight->SetBinContent(67,0.9999998);
  weight->SetBinContent(68,0.9999999);
  weight->SetBinContent(69,1);
  weight->SetBinContent(70,1);
  weight->SetBinContent(71,0.9999997);
  weight->SetBinContent(72,0.9999995);
  weight->SetBinContent(73,0.9999994);
  weight->SetBinContent(74,1);
  weight->SetBinContent(75,1);
  weight->SetBinContent(76,1);
  weight->SetBinContent(77,1);
  weight->SetBinContent(78,0.9999999);
  weight->SetBinContent(79,1);
  weight->SetBinContent(80,1);
  weight->SetEntries(526);
  Int_t limit = weight->GetRandom();
  pythia->SetTriggerChargedMultiplicity(limit, 1.4);
  comment = comment.Append(Form("; multiplicity threshold set to %d in |eta| < 1.4", limit));
  return pythia;
}

//______________________________________________________________________
AliGenerator* PythiaJets()
{
      comment = comment.Append(" pp: Pythia low-pt");
//
//    Pythia
      AliGenPythia* pythia = new AliGenPythia(-1); 
      /*   pythia->SetMomentumRange(0, 999999.);
      pythia->SetThetaRange(0., 180.);
      pythia->SetYRange(-12.,12.);
      pythia->SetPtRange(0,1000.);*/
      pythia->SetProcess(kPyJets);
      pythia->SetEnergyCMS(energy);
      pythia->SetStrucFunc(kCTEQ6l);
      //  pythia->SetPtHard(50, 1000);
//
      return pythia;
}



//___________________________________________________________________
AliGenerator* MbPhojet()
{
  comment = comment.Append(" pp: Pythia low-pt");
#if defined(__CINT__)
  gSystem->Load("libDPMJET");      // Parton density functions
  gSystem->Load("libTDPMjet");      // Parton density functions
#endif
  AliGenDPMjet* dpmjet = new AliGenDPMjet(-1); 

  /*  dpmjet->SetMomentumRange(0, 999999.);
  dpmjet->SetThetaRange(0., 180.);
  dpmjet->SetYRange(-12.,12.);
  dpmjet->SetPtRange(0,1000.);*/

  dpmjet->SetProcess(kDpmMb);
  dpmjet->SetEnergyCMS(energy);
  return dpmjet;
}

//__________________________________________________________________
AliGenerator* Hijing()
{
    AliGenHijing *gener = new AliGenHijing(-1);
    // centre of mass energy 
    gener->SetEnergyCMS(energy);
    gener->SetImpactParameterRange(bMin, bMax);	
    // reference frame
    gener->SetReferenceFrame("CMS");
    // projectile
    gener->SetProjectile("A", 208, 82);
    gener->SetTarget    ("A", 208, 82);
    // tell hijing to keep the full parent child chain
    gener->KeepFullEvent();
    // enable jet quenching
    gener->SetJetQuenching(1);
    // enable shadowing
    gener->SetDecaysOff(0);
    gener->SetShadowing(1);
    // Don't track spectators
    gener->SetSpectators(0);
    // kinematic selection
    gener->SetSelectAll(0);
    return gener;
}


//__________________________________________________________________
AliGenerator* Hijing2000()
{
  AliGenHijing *gener = (AliGenHijing*) Hijing();
  gener->SetJetQuenching(0);	
  gener->SetPtHardMin (2.3);
  return gener;
}



//_____________________________________________________________
AliGenerator* Hijing_pA(Bool_t kSlowN) {
  AliGenCocktail *gener = 0x0;
  if (kSlowN) {
    gener  = new AliGenCocktail();
    gener->SetProjectile("A", 208, 82);
    gener->SetTarget    ("P",   1,  1);
    gener->SetEnergyCMS(energy);
  }
  
  AliGenHijing   *hijing = new AliGenHijing(-1);
  // centre of mass energy 
  hijing->SetEnergyCMS(energy);
  // impact parameter range
  hijing->SetImpactParameterRange(0., 100.);
  // reference frame
  hijing->SetReferenceFrame("CMS");
  hijing->SetBoostLHC(1);
  // projectile
  hijing->SetProjectile("A", 208, 82);
  hijing->SetTarget    ("P", 1, 1);
  // tell hijing to keep the full parent child chain
  hijing->KeepFullEvent();
  // enable jet quenching
  hijing->SetJetQuenching(0);
  // enable shadowing
  hijing->SetShadowing(1);
  // kinematic selection
  hijing->SetSelectAll(0);
  
  if (!kSlowN) {
    // DO track spectators
    hijing->SetSpectators(1);
    return hijing;
  }
  else {
    // Cocktail with slow nucleons generator
    // DO NOT track spectators
    hijing->SetSpectators(0);
    AliGenSlowNucleons* gray   = new AliGenSlowNucleons(1);
    AliCollisionGeometry* coll = hijing->CollisionGeometry();
    AliSlowNucleonModelExp* model = new AliSlowNucleonModelExp();
    //  Not yet in the release...
    //	  model->SetSaturation(kTRUE);
    gray->SetSlowNucleonModel(model);
    gray->SetTarget(208,82);
    gray->SetThetaDist(1);
    gray->SetProtonDirection(1);
    //	  gray->SetDebug(1);
    gray->SetNominalCmsEnergy(2.*pBeamEnergy);
    gray->NeedsCollisionGeometry();
    gray->SetCollisionGeometry(coll);
    
    gener->AddGenerator(hijing, "Hijing pPb", 1);
    gener->AddGenerator(gray, "Gray Particles", 1);
    
    return gener;
  }
}



//__________________________________________________________________
AliGenerator* DPMjet()
{
  AliGenDPMjet* dpmjet = new AliGenDPMjet(-1); 
  dpmjet->SetEnergyCMS(energy);
  dpmjet->SetProjectile("A", 208, 82);
  dpmjet->SetTarget    ("A", 208, 82);
  dpmjet->SetImpactParameterRange(bMin, bMax);
  dpmjet->SetPi0Decay(0);
  return dpmjet;
}


//____________________________________________
AliGenerator* DPMjet_pA(Bool_t fragments)
{
  AliGenDPMjet *gener = new AliGenDPMjet(-1);
  //  A-p
  gener->SetProjectile("A", 208, 82);
  gener->SetTarget("P", 1, 1);
  //
  gener->SetEnergyCMS(energy);
  gener->SetProjectileBeamEnergy(82.*pBeamEnergy/208.);

  //
  gener->SetProcess(kDpmMb);
  gener->SetImpactParameterRange(0., 100.);
  // DO NOT BOOST !
  //  gener->SetBoostLHC(1);
  //
  gener->SetFragmentProd(fragments);
  
  return gener;

}

//__________________________________________________________________
AliGenerator* AmptDefault()
{
  AliGenAmpt *genHi = new AliGenAmpt(-1);
  //=========================================================================
  // THE DECAYER
  AliDecayer *decayer = new AliDecayerPythia();
  cout << "*****************************************" << endl;
  genHi->SetForceDecay( kHadronicD );
  genHi->SetDecayer( decayer );
  //=========================================================================
  genHi->SetEnergyCMS(energy);
  genHi->SetReferenceFrame("CMS");
  genHi->SetProjectile("A",208,82);
  genHi->SetTarget("A",208,82);
    
  genHi->SetIsoft(1); //1: defaul - 4: string melting
  genHi->SetStringFrag(0.5,0.9); //Lund string frgmentation parameters
  genHi->SetPtHardMin (3);
  //genHi->SetImpactParameterRange(9.,9.5);
  genHi->SetImpactParameterRange(bMin,bMax);
    
  // Xmu = 3.2 fm^-1 and as = 0.33 ==> sigma_{partonic} = 1.5mb
  // Ntmax = 150
  // v_2{2} = 0.102105 +/- 0.000266894
  // v_2{4} = 0.0829477 +/- 0.00106158
    
  genHi->SetNtMax(150);        // NTMAX: number of timesteps (D=150)
  genHi->SetXmu(3.2264);        // parton screening mass in fm^(-1) (D=3.2264d0)
    
  genHi->SetJetQuenching(0);  // enable jet quenching
  genHi->SetShadowing(1);     // enable shadowing
  genHi->SetDecaysOff(1);     // neutral pion and heavy particle decays switched off
  genHi->SetSpectators(0);    // track spectators
  //Boost into LHC lab frame
  genHi->SetBoostLHC(1);
  //  genHi->Init();
  genHi->SetRandomReactionPlane(kTRUE);
 
  return genHi;
}

//__________________________________________________________________
AliGenerator* AmptStringMelting()
{
  AliGenAmpt *genHi = new AliGenAmpt(-1);
  //=========================================================================
  // THE DECAYER
  AliDecayer *decayer = new AliDecayerPythia();
  cout << "*****************************************" << endl;
  genHi->SetForceDecay( kHadronicD );
  genHi->SetDecayer( decayer );
  //=========================================================================
  genHi->SetEnergyCMS(energy);
  genHi->SetReferenceFrame("CMS");
  genHi->SetProjectile("A",208,82);
  genHi->SetTarget("A",208,82);
    
  genHi->SetIsoft(4); //1: defaul - 4: string melting
  genHi->SetStringFrag(0.5,0.9); //Lund string frgmentation parameters
  genHi->SetPtHardMin (3);
  //genHi->SetImpactParameterRange(9.,9.5);
  genHi->SetImpactParameterRange(bMin,bMax);

  // Xmu = 3.2 fm^-1 and as = 0.33 ==> sigma_{partonic} = 1.5mb
  // Ntmax = 150
  // v_2{2} = 0.102105 +/- 0.000266894
  // v_2{4} = 0.0829477 +/- 0.00106158

  genHi->SetNtMax(150);        // NTMAX: number of timesteps (D=150)
  genHi->SetXmu(3.2264);        // parton screening mass in fm^(-1) (D=3.2264d0)

  genHi->SetJetQuenching(0);  // enable jet quenching
  genHi->SetShadowing(1);     // enable shadowing
  genHi->SetDecaysOff(1);     // neutral pion and heavy particle decays switched off
  genHi->SetSpectators(0);    // track spectators
  //Boost into LHC lab frame
  genHi->SetBoostLHC(1);
//  genHi->Init();
  genHi->SetRandomReactionPlane(kTRUE);
  return genHi;

}

//__________________________________________________________________
AliGenerator* AmptStringMeltingNoART()
{
    AliGenAmpt *genHi = new AliGenAmpt(-1);
    //=========================================================================
    // THE DECAYER
    AliDecayer *decayer = new AliDecayerPythia();
    cout << "*****************************************" << endl;
    genHi->SetForceDecay( kHadronicD );
    genHi->SetDecayer( decayer );
    //=========================================================================
    genHi->SetEnergyCMS(energy);
    genHi->SetReferenceFrame("CMS");
    genHi->SetProjectile("A",208,82);
    genHi->SetTarget("A",208,82);
    
    genHi->SetIsoft(4); //1: defaul - 4: string melting
    genHi->SetStringFrag(0.5,0.9); //Lund string frgmentation parameters
    genHi->SetPtHardMin (3);
    //genHi->SetImpactParameterRange(9.,9.5);
    genHi->SetImpactParameterRange(bMin,bMax);
    
    // Xmu = 3.2 fm^-1 and as = 0.33 ==> sigma_{partonic} = 1.5mb
    // Ntmax = 150
    // v_2{2} = 0.102105 +/- 0.000266894
    // v_2{4} = 0.0829477 +/- 0.00106158
    
    genHi->SetNtMax(3);        // NTMAX: number of timesteps (D=150)
    genHi->SetXmu(3.2264);        // parton screening mass in fm^(-1) (D=3.2264d0)
    
    genHi->SetJetQuenching(0);  // enable jet quenching
    genHi->SetShadowing(1);     // enable shadowing
    genHi->SetDecaysOff(1);     // neutral pion and heavy particle decays switched off
    genHi->SetSpectators(0);    // track spectators
    //Boost into LHC lab frame
    genHi->SetBoostLHC(1);
    //  genHi->Init();
    genHi->SetRandomReactionPlane(kTRUE);
    return genHi;
    
}


//__________________________________________________________________
AliGenerator* AmptReducedFlow()
{
  AliGenAmpt *genHi = new AliGenAmpt(-1);
  //=========================================================================
  // THE DECAYER
  AliDecayer *decayer = new AliDecayerPythia();
  cout << "*****************************************" << endl;
  genHi->SetForceDecay( kHadronicD );
  genHi->SetDecayer( decayer );
  //=========================================================================
  genHi->SetEnergyCMS(energy);
  genHi->SetReferenceFrame("CMS");
  genHi->SetProjectile("A",208,82);
  genHi->SetTarget("A",208,82);
    
  genHi->SetIsoft(4); //1: defaul - 4: string melting
  genHi->SetStringFrag(0.5,0.9); //Lund string frgmentation parameters
  genHi->SetPtHardMin (3);
  //genHi->SetImpactParameterRange(9.,9.5);
  genHi->SetImpactParameterRange(bMin,bMax);

  // Xmu = 12.4 fm^-1 and as = 0.33 ==> sigma_{partonic} = 0.1mb
  // Ntmax = 20
  // flow estimates from Q-cumulants
  // (POI, without weights)
  // v_2{2} = 0.0549735 +/- 0.000270249
  // v_2{4} = 0.0421905 +/- 0.00189449

  genHi->SetNtMax(20);        // NTMAX: number of timesteps (D=150)
  genHi->SetXmu(12.4);        // parton screening mass in fm^(-1) (D=3.2264d0)

  genHi->SetJetQuenching(0);  // enable jet quenching
  genHi->SetShadowing(1);     // enable shadowing
  genHi->SetDecaysOff(1);     // neutral pion and heavy particle decays switched off
  genHi->SetSpectators(0);    // track spectators
  //Boost into LHC lab frame
  genHi->SetBoostLHC(1);
 // genHi->Init();
  genHi->SetRandomReactionPlane(kTRUE);
  return genHi;

}

//__________________________________________________________________
AliGenerator* AmptpA()
{
    AliGenAmpt *genHi = new AliGenAmpt(-1);
    //=========================================================================
    // THE DECAYER
    AliDecayer *decayer = new AliDecayerPythia();
    cout << "*****************************************" << endl;
    genHi->SetForceDecay( kHadronicD );
    genHi->SetDecayer( decayer );
    //=========================================================================
    genHi->SetEnergyCMS(energy);
    genHi->SetReferenceFrame("CMS");
    genHi->SetProjectile("A", 208, 82);
    genHi->SetTarget    ("P", 1, 1);
    genHi->SetIsoft(4); //1: defaul - 4: string melting
    genHi->SetStringFrag(0.5,0.9); //Lund string frgmentation parameters
    genHi->SetPtHardMin (3);
    //genHi->SetImpactParameterRange(9.,9.5);
    genHi->SetImpactParameterRange(bMin,bMax);
    genHi->SetNtMax(1500); //NTMAX: number of timesteps (D=150)
    genHi->SetXmu(3.2264); //parton screening mass in fm^(-1) (D=3.2264d0)
    genHi->SetJetQuenching(0); // enable jet quenching
    genHi->SetShadowing(1);    // enable shadowing
    genHi->SetDecaysOff(1);    // neutral pion and heavy particle decays switched off
    genHi->SetSpectators(0);   // track spectators
    //Boost into LHC lab frame
    genHi->SetBoostLHC(1);
    //  genHi->Init();
    genHi->SetRandomReactionPlane(kTRUE);
    return genHi;
    
}
