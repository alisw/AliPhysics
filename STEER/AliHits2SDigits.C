////////////////////////////////////////////////////////////////////////
//
// name: AliHits2SDigits
// date: 4.4.2002
// last update: 4.4.2002
// author: Jiri Chudoba
// version: 1.0
//
// description: 
//       creates sdigits for several detectors
//
// input:
//       TString fileName ... galice input file
//       Int_t nEvents  ... how many events to proceed
//       Int_t firstEvent  ... first event number
//       Int_t ITS, TPC, ...   many flags for diff. detectors
//
// History:
//
// 21.07.03 - changes for NewIO
//
// 04.04.02 - first version
// 
////////////////////////////////////////////////////////////////////////

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "iostream.h"
#include "TTree.h"
#include "TBranch.h"
#include "TDirectory.h"
#include "TFile.h"
#include "AliRun.h"
#include "TParticle.h"
#include "TPC/AliTPCDigitsArray.h"
#include "AliHeader.h"
#include "TGeometry.h"
#include "TObjArray.h"
#include "TString.h"
#include "ITS/AliITS.h"
#include "TPC/AliTPC.h"
#include "PHOS/AliPHOSSDigitizer.h"
#include "TRD/AliTRDdigitizer.h"
#include "TStopwatch.h"
#include "TRD/AliTRDparameter.h"
#endif

AliRunLoader* Init(TString fileName);

AliTRDdigitizer *InitTRDdigitizer();

// global variables

TFile *gFileHits = 0;
Bool_t gSameFiles = kFALSE;
Int_t gDEBUG = 1;


Int_t AliHits2SDigits(TString fileName="galice.root", 
                  Int_t nEvents = 1, Int_t firstEvent = 0, Int_t iITS = 0,
                  Int_t iTPC = 0, Int_t iTRD = 0,Int_t iPHOS = 0)
{
//
// Initialization
//
  AliRunLoader* rl = Init(fileName);
  if (!rl) return 1;

// ITS
  AliITS *ITS = NULL;
  if (iITS) {
    ITS  = (AliITS*) gAlice->GetModule("ITS");
    if (!ITS) {
      iITS = 0;
      cerr<<"AliITS object not found on file." << endl;
    } else if (!ITS->GetITSgeom()) {
      cerr<<"AliITSgeom not found." << endl;
      iITS = 0;
    }
  }

// TPC
  AliTPC *TPC = NULL;
  if (iTPC) {
    TPC = (AliTPC*)gAlice->GetDetector("TPC");
    if (!TPC) {
      iTPC = 0;
      cerr<<"AliTPC object not found"<<endl;
    }
  }

// TRD
  AliTRDdigitizer *sdTRD = NULL;
  if (iTRD) {
    sdTRD = InitTRDdigitizer();
  }


// PHOS
  AliPHOSSDigitizer *sdPHOS = NULL;
  if (iPHOS) {
    sdPHOS = new AliPHOSSDigitizer(fileName.Data());
  }



//
// loop over events
//
  TStopwatch timer;
  timer.Start();
  for (Int_t iEvent = firstEvent;iEvent<firstEvent+nEvents;iEvent++){
    rl->GetEvent(iEvent);
//    gAlice->MakeTree("S",fileSDigits);
    
// ITS
    if (iITS) {
      if (gDEBUG) {cout<<"  Create ITS sdigits: ";}
      AliLoader* loader = rl->GetLoader("ITSLoader");
      if (loader)
       { 
        loader->LoadHits("read");
        loader->LoadSDigits("recreate");
	if(!loader->TreeS()) loader->MakeTree("S");
	ITS->MakeBranch("S");
        ITS->SetTreeAddress();
        ITS->Hits2SDigits();
        loader->UnloadHits();
        loader->UnloadSDigits();
        if (gDEBUG) {cout<<"done"<<endl;}
       }
      else if (gDEBUG) {cout<<"Did not get loader"<<endl;}
    }

// TPC
    if (iTPC) {
      if (gDEBUG) {cout<<"  Create TPC sdigits: ";}
      AliLoader* loader = rl->GetLoader("TPCLoader");
      if (loader)
       { 
        loader->LoadHits("read");
        loader->LoadSDigits("recreate");
      
	TPC->SetTreeAddress();
        TPC->SetActiveSectors(1);
        TPC->Hits2SDigits2(iEvent);
        loader->UnloadHits();
        loader->UnloadSDigits();
        if (gDEBUG) {cout<<"done"<<endl;}
       }
      else if (gDEBUG) {cout<<"Did not get loader"<<endl;}
    }

// TRD
    if (iTRD) {
      if (gDEBUG) {cout<<"  Create TRD sdigits: ";}
      AliLoader* loader = rl->GetLoader("TRDLoader");
      if (loader)
       { 
        loader->LoadHits("read");
        loader->LoadSDigits("recreate");
        sdTRD->MakeDigits();
        sdTRD->WriteDigits();
        loader->UnloadHits();
        loader->UnloadSDigits();
        if (gDEBUG) {cout<<"done"<<endl;}
       }
      else if (gDEBUG) {cout<<"Did not get loader"<<endl;}
    }
    
  } // end of loop over events

// PHOS processes always all events
  if (iPHOS) {
    sdPHOS->ExecuteTask("deb all");
  }

//
// finish 
//
  timer.Stop(); 
  timer.Print();

  delete rl;
  return 0;
}
 

////////////////////////////////////////////////////////////////////////
AliRunLoader* Init(TString fileName) 
 {
// open input file, read in gAlice, prepare output file
  if (gAlice) delete gAlice;
  gAlice = 0;
  AliRunLoader*rl = AliRunLoader::Open(fileName);
  if (rl == 0x0) return 0x0;
  rl->LoadgAlice();
  gAlice = rl->GetAliRun();
  return rl;
 
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
AliTRDdigitizer *InitTRDdigitizer() {
// initialization of TRD digitizer
  AliTRDdigitizer *sdTRD = new AliTRDdigitizer("TRDdigitizer"
                                     ,"TRD digitizer class");
  sdTRD->SetDebug(0);
  sdTRD->SetSDigits(kTRUE);
  AliTRDparameter *TRDparam = new AliTRDparameter("TRDparameter"
                                      ,"TRD parameter class");

  sdTRD->SetParameter(TRDparam);
  sdTRD->InitDetector();
  return sdTRD;
}
