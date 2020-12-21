//
// Template for a Config.C file
// to be used by AliMuonAccEffSubmitter, which will
// replace all instances of VAR_ with something
// relevant and useable...

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
#include "AliBODY.h"
#include "AliMAG.h"
#include "AliABSOv3.h"
#include "AliDIPOv3.h"
#include "AliHALLv3.h"
#include "AliFRAMEv2.h"
#include "AliSHILv3.h"
#include "AliPIPEv3.h"
#include "AliITSv11Hybrid.h"
#include "AliZDCv3.h"
#include "AliFMDv1.h"
#include "AliMUONv1.h"
#include "AliT0v1.h"
#include "AliVZEROv7.h"
#endif


//--- Functions ---
class AliGenPythia;

//_________________________________________________________
void LoadPhotos( void )
{
  // old evtgen libraries
    std::cout << "LoadPhotos" << std::endl;
     gSystem->Setenv("PYTHIA8DATA", gSystem->ExpandPathName("$ALICE_ROOT/PYTHIA8/pythia8210/xmldoc"));
  gSystem->Load("libPhotos" );
  gSystem->Load("libEvtGen");
  gSystem->Load("libTEvtGen");
}

void Config()
{
  // Libraries required by geant321
#if defined(__CINT__)
  gSystem->AddIncludePath("-I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
  gSystem->Load("libgeant321");
  LoadPhotos();
#endif

  
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
  rl->SetNumberOfEventsPerFile(5000);
  gAlice->SetRunLoader(rl);
  
  if ( TString("VAR_TRIGGER_CONFIGURATION").Length() > 0 )
  {
    AliSimulation::Instance()->SetTriggerConfig("VAR_TRIGGER_CONFIGURATION");
    cout<<"Trigger configuration is set to VAR_TRIGGER_CONFIGURATION" << std::endl;
  }
  
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

  //======================//
  // Set External decayer //
  //======================//
  TVirtualMCDecayer* decayer = new AliDecayerPythia;
  
  /* FIXME: put back polarization switch ?
  if (polar == kNO_pol){
    decayer = new AliDecayerPythia();
  } else if (polar == kTH_pol){
    decayer = new AliDecayerPolarized(1.,AliDecayerPolarized::kHelicity,AliDecayerPolarized::kMuon);
  } else if (polar == kT_pol){
    decayer = new AliDecayerPolarized(1.,AliDecayerPolarized::kColSop,AliDecayerPolarized::kMuon);
  } else if (polar == kLH_pol){
    decayer = new AliDecayerPolarized(-1.,AliDecayerPolarized::kHelicity,AliDecayerPolarized::kMuon);
  } else if (polar == kL_pol){
    decayer = new AliDecayerPolarized(-1.,AliDecayerPolarized::kColSop,AliDecayerPolarized::kMuon);
  }
   */
  
  decayer->SetForceDecay(kAll);
  decayer->Init();
  gMC->SetExternalDecayer(decayer);

  //=========================//
  // Generator Configuration //
  //=========================//

  //AliGenerator* gener = CreateGenerator();

  std::cout << "VAR_GENERATOR settings " << std::endl;
  gROOT->LoadMacro("VAR_GENERATOR.C+");
  AliGenerator* gener = VAR_GENERATOR();
  
  TString slibs = gSystem->GetLibraries();
  TObjArray* olibs = slibs.Tokenize(" ");
  TObjString* s;
  TIter next(olibs);
  std::cout << "List of libraries=" << std::endl;
  while ( ( s = static_cast<TObjString*>(next())) )
  {
    std::cout << s->String().Data() << std::endl;
  }
  

  gener->SetOrigin(0., 0., 0.); // Taken from OCDB
  
  gener->SetSigma(VAR_VERTEX_SIGMA_X, VAR_VERTEX_SIGMA_Y, 0.);      // Sigma in (X,Y,Z) (cm) on IP position, sigmaz taken from OCDB
  gener->SetVertexSmear(kPerEvent);
  gener->Init();
  
  gener->Print();
    
  rl->CdGAFile();
  
  Int_t iABSO  = 1;
  Int_t iDIPO  = 1;
  Int_t iFMD   = 1;
  Int_t iFRAME = 1;
  Int_t iHALL  = 1;
  Int_t iITS   = 1;
  Int_t iMAG   = 1;
  Int_t iMUON  = 1;
  Int_t iPIPE  = 1;
  Int_t iSHIL  = 1;
  Int_t iT0    = 1;
  Int_t iVZERO = 1;
  Int_t iZDC   = 0;
  Int_t iAD    = 1;

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

        AliFRAMEv3 *FRAME = new AliFRAMEv3("FRAME", "Space Frame");
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

	AliITS *ITS  = new AliITSv11("ITS","ITS v11");
    }

    if (iZDC)
    {
        //=================== ZDC parameters ============================
	
      AliZDC *ZDC = new AliZDCv4("ZDC", "normal ZDC");
      ZDC->SetLumiLength(0.);
      ZDC->SetVCollSideCAperture(2.8);
      ZDC->SetVCollSideCApertureNeg(2.8);
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
      MUON->SetTriggerEffCells(1);
      MUON->SetTriggerResponseV1(2);
    }

    if (iT0)
    {
        //=================== T0 parameters ============================
        AliT0 *T0 = new AliT0v1("T0", "T0 Detector");
    }

    if (iVZERO)
    {
        //=================== ACORDE parameters ============================

        AliVZERO *VZERO = new AliVZEROv7("VZERO", "normal VZERO");
    }

    if (iAD)
    {
        //=================== AD parameters ============================

        AliAD *AD = new AliADv1("AD", "normal AD");
    }

}
