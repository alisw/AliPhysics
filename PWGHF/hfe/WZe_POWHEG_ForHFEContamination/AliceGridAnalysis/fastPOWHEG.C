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
#include "STEER/AliRunLoader.h"
#include "STEER/AliRun.h"
#include "STEER/AliConfig.h"
#include "PYTHIA6/AliDecayerPythia.h"
#include "PYTHIA6/AliGenPythia.h"
#include "TDPMjet/AliGenDPMjet.h"
#include "STEER/AliMagFCheb.h"
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
#include "PHOS/AliPHOSSimParam.h"
#include "PMD/AliPMDv1.h"
#include "T0/AliT0v1.h"
#include "EMCAL/AliEMCALv2.h"
#include "ACORDE/AliACORDEv1.h"
#include "VZERO/AliVZEROv7.h"
#endif


//--- Functions ---
class AliGenPythia;

void fastPOWHEG()
{
  // Libraries required by geant321
#if defined(__CINT__)
  gSystem->Load("liblhapdf5_9_1");      // Parton density functions
  gSystem->Load("libEGPythia6");   // TGenerator interface
  if ( TString("GenW_Pythia6_POWHEG").Contains("pythia6",TString::kIgnoreCase) )
  {
    std::cout << "Setting up Pythia6 required env. variables" << std::endl;
    gSystem->AddIncludePath("-I$ALICE_ROOT/PYTHIA6 -I$ALICE_ROOT/STEER/STEER -I$ALICE_ROOT/STEER/STEERBase -I$ALICE_ROOT/LHAPDF -I$ALICE_ROOT/FASTSIM");
    gSystem->Load("libpythia6.4.25");
  }
  else  gSystem->Load("libpythia6");     // Pythia 6.2 (for decayer)
  gSystem->Load("libAliPythia6");  // ALICE specific implementations
  gSystem->Load("libgeant321");

  if ( TString("GenW_Pythia6_POWHEG").Contains("pythia8",TString::kIgnoreCase) )
  {
    std::cout << "Setting up Pythia8 required libraries and env. variables" << std::endl;
    gSystem->Load("libpythia8");
    gSystem->Load("libAliPythia8");
    
    
  }

#endif


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
  rl->SetNumberOfEventsPerFile(100000);
  rl->LoadKinematics("RECREATE");
  rl->MakeTree("E");
  gAlice->SetRunLoader(rl);
  
  if ( TString("p-p").Length() > 0 )
  {
    AliSimulation::Instance()->SetTriggerConfig("p-p");
    cout<<"Trigger configuration is set to p-p" << std::endl;
  }
 

//  Create stack
    rl->MakeStack();
    AliStack* stack      = rl->Stack();
 
//  Header
    AliHeader* header = rl->GetHeader();
//
  //=========================//
  // Generator Configuration //
  //=========================//

  //AliGenerator* gener = CreateGenerator();

  std::cout << "GenW_Pythia6_POWHEG settings " << std::endl;
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gSystem->AddIncludePath("-I$ALICE_ROOT/EVGEN");
  gSystem->AddIncludePath("-I$ALICE_ROOT/STEER/STEER");
  gROOT->LoadMacro("GenW_Pythia6_POWHEG.C+");
  AliGenerator* gener = GenW_Pythia6_POWHEG();
  
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

  Float_t sigmax = 0.0025;
  Float_t sigmay = 0.0029;
  
  gener->SetSigma(sigmax, sigmay, 0.);      // Sigma in (X,Y,Z) (cm) on IP position, sigmaz taken from OCDB
  gener->SetVertexSmear(kPerEvent);
  gener->Init();
  gener->SetStack(stack);
  
  gener->Print();
    
  //rl->CdGAFile();
 
    Int_t iev = 0;
    //Int_t nev = 80;  // ok
    Int_t nev = 3500;  // ok
    //Int_t nev = 1000;  // ok
     
    for (iev = 0; iev < nev; iev++) {

	printf("\n \n Event number %d \n \n", iev);
	
//  Initialize event
	header->Reset(0,iev);
	rl->SetEventNumber(iev);
	stack->Reset();
	rl->MakeTree("K");
//	stack->ConnectTree();
    
//  Generate event
	gener->Generate();
//  Analysis
	Int_t npart = stack->GetNprimary();
	printf("Analyse %d Particles\n", npart);
	for (Int_t part=0; part<npart; part++) {
	    TParticle *MPart = stack->Particle(part);
	    Int_t mpart  = MPart->GetPdgCode();
	    printf("Particle %d\n", mpart);
	}
	
//  Finish event
	header->SetNprimary(stack->GetNprimary());
	header->SetNtrack(stack->GetNtrack());  
//      I/O
//	
	stack->FinishEvent();
	header->SetStack(stack);
	rl->TreeE()->Fill();
	rl->WriteKinematics("OVERWRITE");

    } // event loop
//
//                         Termination
//  Generator
    gener->FinishRun();
//  Write file
    rl->WriteHeader("OVERWRITE");
    gener->Write();
    rl->Write();

 
}
