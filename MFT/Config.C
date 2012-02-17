#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TRandom.h>
#include <TDatime.h>
#include <TSystem.h>
#include <TVirtualMC.h>
#include <TGeant3TGeo.h>
#include "AliRunLoader.h"
#include "AliRun.h"
#include "AliConfig.h"
#include "AliDecayerPythia.h"
#include "AliGenPythia.h"
#include "AliGenDPMjet.h"
#include "AliMagFCheb.h"
#include "AliBODY.h"
#include "AliMAG.h"
#include "AliABSOv3.h"
#include "AliDIPOv3.h"
#include "AliHALLv3.h"
#include "AliFRAMEv2.h"
#include "AliSHILv3.h"
#include "AliPIPEv3.h"
#include "AliPIPEv4.h"
#include "AliITSv11.h"
#include "AliTPCv2.h"
#include "AliTOFv6T0.h"
#include "AliHMPIDv3.h"
#include "AliZDCv3.h"
#include "AliTRDv1.h"
#include "AliTRDgeometry.h"
#include "AliFMDv1.h"
#include "AliMUONv1.h"
#include "AliPHOSv1.h"
#include "AliPHOSSimParam.h"
#include "AliPMDv1.h"
#include "AliT0v1.h"
#include "AliEMCALv2.h"
#include "AliACORDEv1.h"
#include "AliVZEROv7.h"
#include "AliMFT.h"
#endif

enum PDCProc_t {kGenBox,
		kGenMuonLMR,
		kGenPionKaon,
		kGenCorrHF,
                kPythia6,
		kPythiaPerugia0, 
		kPythiaPerugia0Jpsi2e, 
		kPythiaPerugia0BtoJpsi2e, 
		kHijing, 
		kHijing2500,
		kHijing2500Cocktail};

const Char_t* pprRunName[] = {"kGenBox",
			      "kGenMuonLMR",
			      "kGenPionKaon",
			      "kGenCorrHF",
			      "kPythia6",
			      "kPythiaPerugia0", 
			      "kPythiaPerugia0Jpsi2e", 
			      "kPythiaPerugia0BtoJpsi2e", 
			      "kHijing", 
			      "kHijing2500", 
			      "kHijing2500Cocktail"};

enum Mag_t { kNoField, k5kG, kFieldMax };

const Char_t* pprField[] = { "kNoField", "k5kG", "kFieldMax" };

void LoadLibs();

// ----------------------- Generator, field, beam energy,... ------------------------------------------------------------
static PDCProc_t     proc     = kGenBox;
static PDCProc_t     signal   = kGenBox;    // only in case kHijing2500Cocktail is the proc
static Mag_t         mag      = k5kG;
static Float_t       energy   = 14000.; // energy in CMS
static Float_t       bMin     = 0.;
static Float_t       bMax =   = 5.; // 0-5 fm corresponds to around 0-10% (see https://twiki.cern.ch/twiki/bin/viewauth/ALICE/CentStudies#Tables_with_centrality_bins_for)
static Double_t      JpsiPol  = 0; // Jpsi polarisation
static Bool_t        JpsiHarderPt = kFALSE; // Jpsi harder pt spectrum (8.8 TeV)
// ----------------------------------------------------------------------------------------------------------------------

static TString comment;

//====================================================================================================================================================

void Config() {

  //  AliLog::SetClassDebugLevel("AliMFT", 1);

  LoadLibs();

  new TGeant3TGeo("C++ Interface to Geant3");

  // Create the output file

  AliRunLoader* rl=0x0;

  cout<<"Config.C: Creating Run Loader ..."<<endl;
  rl = AliRunLoader::Open("galice.root", AliConfig::GetDefaultEventFolderName(), "recreate");
  if (rl == 0x0) {
    gAlice->Fatal("Config.C","Can not instatiate the Run Loader");
    return;
  }
  rl->SetCompressionLevel(2);
  rl->SetNumberOfEventsPerFile(1000);
  gAlice->SetRunLoader(rl);

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
  
  TVirtualMCDecayer *decayer = new AliDecayerPythia();
  decayer->SetForceDecay(kAll);
  decayer->Init();
  gMC->SetExternalDecayer(decayer);
  
  // Generator
  AliGenerator* gener = 0x0;
  if (proc == kPythia6)                   gener = MbPythia(); 
  else if (proc == kPythiaPerugia0)       gener = MbPythiaTunePerugia0();
  else if (proc == kHijing)               gener = Hijing();	
  else if (proc == kHijing2500)           gener = Hijing2500();	
  else if (proc == kHijing2500Cocktail)   gener = Hijing2500Cocktail();
  else if (proc == kGenBox)               gener = GenBox();
  else if (proc == kGenMuonLMR)           gener = GenMuonLMR();
  else if (proc == kGenCorrHF)            gener = GenCorrHF();
  else if (proc == kGenPionKaon)          gener = GenParamPionKaon();

  // Size of the interaction diamond
  Float_t sigmaz  = 5.4 / TMath::Sqrt(2.);     // [cm]
  Float_t betast  = 3.5;                       // beta* [m]
  Float_t eps     = 3.75e-6;                   // emittance [m]
  Float_t gamma   = energy / 2.0 / 0.938272;   // relativistic gamma [1]
  Float_t sigmaxy = TMath::Sqrt(eps * betast / gamma) / TMath::Sqrt(2.) * 100.;  // [cm]

  printf("\n \n Diamond size x-y: %10.3e z: %10.3e\n \n", sigmaxy, sigmaz);
    
//   gener->SetSigma(sigmaxy, sigmaxy, sigmaz);      // Sigma in (X,Y,Z) (cm) on IP position
//   gener->SetVertexSmear(kPerEvent);
  gener->Init();

  printf("\n \n Comment: %s \n \n", comment.Data());

  //  TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", -1., -1., AliMagF::k5kG, AliMagF::kBeamTypeAA, 2750.));
  TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", -1., -1., AliMagF::k5kG, AliMagF::kBeamTypepp, 7000.));

  rl->CdGAFile();
  
  // Detector Setup

  Int_t iABSO  = 1;
  Int_t iDIPO  = 1;
  Int_t iHALL  = 1;
  Int_t iMUON  = 1;
  Int_t iPIPE  = 1;
  Int_t iSHIL  = 1;
  Int_t iT0    = 0;
  Int_t iVZERO = 1;
  Int_t iMFT   = 1;
  Int_t iACORDE= 0;
  Int_t iEMCAL = 0;
  Int_t iFMD   = 0;
  Int_t iFRAME = 0;
  Int_t iITS   = 0;
  Int_t iMAG   = 0;
  Int_t iPHOS  = 0;
  Int_t iPMD   = 0;
  Int_t iHMPID = 0;
  Int_t iTOF   = 0;
  Int_t iTPC   = 0;
  Int_t iTRD   = 0;
  Int_t iZDC   = 0;
  

  AliBODY *BODY = new AliBODY("BODY", "Alice envelop");

  if (iMAG)       AliMAG    *MAG    = new AliMAG("MAG", "Magnet");
  if (iABSO)      AliABSO   *ABSO   = new AliABSOv3("ABSO", "Muon Absorber");
  if (iDIPO)      AliDIPO   *DIPO   = new AliDIPOv3("DIPO", "Dipole version 3");
  if (iHALL)      AliHALL   *HALL   = new AliHALLv3("HALL", "Alice Hall");
  if (iSHIL)      AliSHIL   *SHIL   = new AliSHILv3("SHIL", "Shielding Version 3");
  if (iITS)       AliITS    *ITS    = new AliITSv11("ITS","ITS v11");
  if (iTPC)       AliTPC    *TPC    = new AliTPCv2("TPC", "Default");
  if (iTOF)       AliTOF    *TOF    = new AliTOFv6T0("TOF", "normal TOF");
  if (iHMPID)     AliHMPID  *HMPID  = new AliHMPIDv3("HMPID", "normal HMPID");
  if (iFMD)       AliFMD    *FMD    = new AliFMDv1("FMD", "normal FMD");
  if (iPHOS)      AliPHOS   *PHOS   = new AliPHOSv1("PHOS", "noCPV_Modules123");
  if (iPMD)       AliPMD    *PMD    = new AliPMDv1("PMD", "normal PMD");
  if (iT0)        AliT0     *T0     = new AliT0v1("T0", "T0 Detector");
  if (iEMCAL)     AliEMCAL  *EMCAL  = new AliEMCALv2("EMCAL", "EMCAL_FIRSTYEARV1");
  if (iACORDE)    AliACORDE *ACORDE = new AliACORDEv1("ACORDE", "normal ACORDE");
  if (iVZERO)     AliVZERO  *VZERO  = new AliVZEROv7("VZERO", "normal VZERO");
  if (iFRAME) {
    AliFRAMEv2 *FRAME = new AliFRAMEv2("FRAME", "Space Frame");
    FRAME->SetHoles(1);
  }
  if (iPIPE) {
    //    AliPIPE *PIPE = new AliPIPEv3("PIPE", "Beam Pipe");
    AliPIPE *PIPE = new AliPIPEv4("PIPE", "Beam Pipe", 1.98, 0.08);
  }
  if (iZDC) {
    AliZDC *ZDC = new AliZDCv3("ZDC", "normal ZDC");
    ZDC->SetSpectatorsTrack();	
    ZDC->SetLumiLength(0.);
  }
  if (iTRD) {
    AliTRD *TRD = new AliTRDv1("TRD", "TRD slow simulator");
  }
  if (iMUON) {
    AliMUON *MUON = new AliMUONv1("MUON", "default");
    MUON->SetTriggerEffCells(1); // not needed if raw masks 
    Char_t* digitstore="AliMUONDigitStoreV2S";    
    MUON->SetDigitStoreClassName(digitstore);
  }
  if (iMFT) {
    AliMFT *MFT = new AliMFT("MFT", "normal MFT");
    MFT->SetNSlices(1);
    MFT->SetChargeDispersion(20.e-4);
    MFT->SetNStepForChargeDispersion(4);
  }

  TIter next(gAlice->Modules());
  AliModule *detector;
  printf("gAlice->Modules:\n");
  while((detector = (AliModule*)next())) printf("%s\n",detector->GetName());

}

//====================================================================================================================================================

AliGenerator* GenBox() {

  AliGenBox *gener = new AliGenBox(2);
  gener->SetMomentumRange(10.0, 10.1);
  gener->SetPhiRange(0., 360.);         
  gener->SetThetaRange(171.0,178.0);
  gener->SetPart(kMuonMinus);           // Muons
  gener->SetOrigin(0., 0., 0.);         // vertex position
  
  return gener;
  
}

//====================================================================================================================================================

AliGenerator* GenMuonLMR() {
  
  AliGenMUONLMR *gener = new AliGenMUONLMR();
  gener->SetMomentumRange(0,999);
  gener->SetPtRange(0,100.);
  gener->SetPhiRange(0., 360.);
  gener->SetThetaRange(0., 180.);
  gener->SetChildThetaRange(170.0,179.0);
  gener->SetYRange(-4.5, -2.0);
  gener->SetOrigin(0.0, 0.0, 0.0);             // vertex position
  gener->SetSigma(0.0, 0.0, 0.0);              // vertex position smearing
  gener->SetVertexSmear(kPerEvent);
  enum {kEta2Body, kEtaDalitz, kRho2Body, kOmega2Body, kOmegaDalitz, kPhi2Body, kEtaPrimeDalitz, kPionLMR, kKaonLMR}; 
  gener->GenerateSingleProcess(kPhi2Body, 500);
  gener->SetCutOnChild(1);

  return gener;

}

//====================================================================================================================================================

AliGenerator* GenParamPionKaon() {
  
  AliGenParamPionsKaons *gener = new AliGenParamPionsKaons(100,"$ALICE_ROOT/MFT/PionKaonKinematics.root");
  gener->SetPtRange(0, 5.);
  gener->SetPhiRange(0., 360.);
  gener->SetYRange(-10., 0.);
  gener->SetOrigin(0.0, 0.0, 0.0);  // vertex position
  gener->SetSigma(0.0, 0.0, 0.0);   // vertex position smearing
  //  gener->SetCutOnChild(1);

  return gener;

}

//====================================================================================================================================================

AliGenerator* GenCorrHF() {
  
  AliGenCorrHF *gener = new AliGenCorrHF(1, 4, 14);  // for charm, 1 pair per event
  // AliGenCorrHF *gener = new AliGenCorrHF(1, 5, 14);  // for beauty, 1 pair per event
  
  gener->SetMomentumRange(0,9999);
  gener->SetCutOnChild(1);          // 1/0 means cuts on children enable/disable
  gener->SetChildThetaRange(170.0,179.0);
  gener->SetOrigin(0,0,0);          //vertex position    
  gener->SetForceDecay(kSemiMuonic);
  gener->SetTrackingFlag(1);
  gener->Init();
  
  return gener;
  
}

//====================================================================================================================================================

AliGenerator* MbPythia() {
  
  comment = comment.Append(" pp: Pythia low-pt");
  
  //    Pythia
  AliGenPythia* pythia = new AliGenPythia(-1); 
  pythia->SetMomentumRange(0, 999999.);
  //  pythia->SetThetaRange(0., 180.);
  //  pythia->SetChildYRange(-12.,0.);
  //  pythia->SetPtRange(0,1000.);
  //  pythia->SetCutOnChild(1);
  pythia->SetProcess(kPyMb);
  pythia->SetEnergyCMS(energy);
  pythia->SetForceDecay(kSemiMuonic);
  
  return pythia;
}

//====================================================================================================================================================

AliGenerator* MbPythiaTunePerugia0() {
  
  comment = comment.Append(" pp: Pythia low-pt (Perugia0)");
  
  //    Pythia
  AliGenPythia* pythia = new AliGenPythia(-1); 
  pythia->SetMomentumRange(0, 999999.);
  pythia->SetThetaRange(0., 180.);
  pythia->SetYRange(-12.,12.);
  pythia->SetPtRange(0,1000.);
  pythia->SetProcess(kPyMb);
  pythia->SetEnergyCMS(energy);
  //    Tune
  //    320     Perugia 0
  pythia->SetTune(320); 
  pythia->UseNewMultipleInteractionsScenario();
  
  return pythia;
}

//====================================================================================================================================================

AliGenerator* MbPythiaTunePerugia0Jpsi2e() {
  
  comment = comment.Append("Jpsi forced to dielectrons");
  AliGenParam *jpsi=0x0;
  if(JpsiHarderPt) jpsi = new AliGenParam(1, AliGenMUONlib::kJpsi, "CDF pp 8.8", "Jpsi");  // 8.8 TeV
  else jpsi = new AliGenParam(1, AliGenMUONlib::kJpsi, "CDF pp 7", "Jpsi");  // 7 TeV
  jpsi->SetPtRange(0.,999.);
  jpsi->SetYRange(-1.0, 1.0);
  jpsi->SetPhiRange(0.,360.);
  jpsi->SetForceDecay(kDiElectron);
  return jpsi;

}

//====================================================================================================================================================

AliGenerator* MbPythiaTunePerugia0BtoJpsi2e() {

  comment = comment.Append(" pp: Pythia (Perugia0) BtoJpsi (1 bbbar per event, 1 b-hadron in |y|<2, 1 J/psi in |y|<2");
  
  //    Pythia
  AliGenPythia* pythia = new AliGenPythia(-1);
  pythia->SetMomentumRange(0, 999999.);
  pythia->SetThetaRange(0., 180.);
  pythia->SetYRange(-2.,2.);
  pythia->SetPtRange(0,1000.);
  pythia->SetProcess(kPyBeautyppMNRwmi);
  pythia->SetEnergyCMS(energy);
  //    Tune
  //    320     Perugia 0
  pythia->SetTune(320);
  pythia->UseNewMultipleInteractionsScenario();
  //
  //    decays
  pythia->SetCutOnChild(1);
  pythia->SetPdgCodeParticleforAcceptanceCut(443);
  pythia->SetChildYRange(-2,2);
  pythia->SetChildPtRange(0,10000.);
  //
  //    decays
  pythia->SetForceDecay(kBJpsiDiElectron);
  
  return pythia;
}

//====================================================================================================================================================

AliGenerator* HijingParam() {

  AliGenHIJINGpara *gener = new AliGenHIJINGpara(2000);
  gener->SetMomentumRange(0,999);
  gener->SetPhiRange(0,360);
  gener->SetThetaRange(171, 179);
  gener->SetOrigin(0,0,0);        // vertex position
  gener->SetSigma(0,0,0);         // Sigma in (X,Y,Z) (cm) on IP position
  gGener = gener;
  
  return gener;

}

//====================================================================================================================================================

AliGenerator* Hijing() {
  
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
  gener->SetShadowing(1);
  // Don't track spectators
  gener->SetSpectators(0);
  // kinematic selection
  gener->SetSelectAll(0);
  return gener;
}

//====================================================================================================================================================

AliGenerator* Hijing2500() {
  
  AliGenHijing *gener = (AliGenHijing*) Hijing();
  gener->SetJetQuenching(0);	
  gener->SetPtHardMin (4.5);
  return gener;
  
}

//====================================================================================================================================================

AliGenerator* Hijing2500Cocktail() {
  
  comment = comment.Append(" PbPb: Hjing2500 at 5.5 + muon signals");

  AliGenCocktail *cocktail = new AliGenCocktail();
  cocktail->SetProjectile("A", 208, 82);
  cocktail->SetTarget    ("A", 208, 82);
  cocktail->SetEnergyCMS(energy);

  // 1 Hijing event: provides underlying event and collision geometry 
  AliGenHijing *hijing = Hijing2500();
  cocktail->AddGenerator(hijing,"hijing",1);

  // generator for the custom signal
  AliGenerator* signalGen = 0x0;      
  switch (signal) {
  case 0:
    signalGen = GenBox();
    break;
  case 1:
    signalGen = GenMuonLMR(); 
    break;
  case 2:
    signalGen = GenCorrHF();
    break;
  default:
    signalGen = GenBox();
    break;
  }
  cocktail->AddGenerator(signalGen, "signal", 1);

  cocktail->SetTrackingFlag(1);

  return cocktail;

}

//====================================================================================================================================================

void LoadLibs() {

#if defined(__CINT__)

  gSystem->Load("liblhapdf");      // Parton density functions
  gSystem->Load("libEGPythia6");   // TGenerator interface
  if (proc == kPythia6) {
    gSystem->Load("libpythia6");        // Pythia 6.2
    gSystem->Load("libAliPythia6");     // ALICE specific implementations
  } 
  else {
    gSystem->Load("libpythia6.4.21");   // Pythia 6.4
    gSystem->Load("libAliPythia6");     // ALICE specific implementations	
  }
  
  if (proc == kHijing || proc == kHijing2500 || proc == kHijing2500Cocktail) {
    gSystem->Load("libhijing");	
    gSystem->Load("libTHijing");
  } 
  
  gSystem->Load("libgeant321");
  
#endif

}

//====================================================================================================================================================

