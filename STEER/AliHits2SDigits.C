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
//       stores sdigits in separate file (or in the source file
//       with hits). Stores gAlice object and copies TE to the
//       file with sdigits
//
// input:
//       TString fileNameSDigits ... output file with sdigits
//       TString fileNameHits ... input file with hits
//       Int_t nEvents  ... how many events to proceed
//       Int_t firstEvent  ... first event number
//       Int_t ITS, TPC, ...   many flags for diff. detectors
//
// History:
//
// 04.04.02 - first version
// 
////////////////////////////////////////////////////////////////////////

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "Riostream.h"
#include "TTree.h"
#include "TBranch.h"
#include "TDirectory.h"
#include "TFile.h"
#include "AliRun.h"
#include "TParticle.h"
#include "AliTPCDigitsArray.h"
#include "AliHeader.h"
#include "TGeometry.h"
#include "TObjArray.h"
#include "TString.h"
#include "AliITS.h"
#include "AliTPC.h"
#include "AliPHOSSDigitizer.h"
#include "AliTRDdigitizer.h"
#include "TStopwatch.h"
#include "AliTRDparameter.h"
#include "AliTOFSDigitizer.h"
#endif

TFile* Init(TString fileNameSDigits, TString fileNameHits);
TFile* OpenFile(TString fileName);
Bool_t ImportgAlice(TFile *file);
AliTRDdigitizer *InitTRDdigitizer();
void AliCopy(TFile *inputFile, TFile *outputFile);

// global variables

TFile *gFileHits = 0;
Bool_t gSameFiles = kFALSE;
Int_t gDEBUG = 1;


Int_t AliHits2SDigits(TString fileNameSDigits="sdigits.root", 
		      TString fileNameHits="rfio:galice.root", 
		      Int_t nEvents = 1, Int_t firstEvent = 0, Int_t iITS = 0,
		      Int_t iTPC = 0, Int_t iTRD = 0,Int_t iPHOS = 0, 
		      Int_t iTOF = 0, Int_t iCopy = 1)
{
//
// Initialization
//
  TFile *fileSDigits;
  fileSDigits = Init(fileNameSDigits, fileNameHits);
  if (!fileSDigits) return 1;
  if (iCopy) {
    AliCopy(gFileHits,fileSDigits);
    gFileHits->cd();
  }  

// ITS
  AliITS *ITS;
  if (iITS) {
    ITS  = (AliITS*) gAlice->GetModule("ITS");
    if (!ITS) {
      iITS = 0;
      cerr<<"AliITS object not found on file." << endl;
    } else if (!ITS->GetITSgeom()) {
      iITS = 0;
      cerr<<"AliITSgeom not found." << endl;
    }
  }

// TPC
  AliTPC *TPC;
  if (iTPC) {
    TPC = (AliTPC*)gAlice->GetDetector("TPC");
    if (!TPC) {
      iTPC = 0;
      cerr<<"AliTPC object not found"<<endl;
    }
  }

// TRD
  AliTRDdigitizer *sdTRD;
  if (iTRD) {
    sdTRD = InitTRDdigitizer();
  }


// PHOS
  AliPHOSSDigitizer *sdPHOS;
  if (iPHOS) {
    sdPHOS = new AliPHOSSDigitizer(fileNameHits.Data());
  }

// TOF
  AliTOFSDigitizer *sdTOF;
  if (iTOF) {
    sdTOF = new AliTOFSDigitizer(fileNameHits.Data(),firstEvent,nEvents);
  }

//
// loop over events
//
  TStopwatch timer;
  timer.Start();
  for (Int_t iEvent = firstEvent;iEvent<firstEvent+nEvents;iEvent++){
    gAlice->GetEvent(iEvent);
    if (!gAlice->TreeS()) gAlice->MakeTree("S",fileSDigits);
    
// ITS
    if (iITS) {
      if (gDEBUG) {cout<<"  Create ITS sdigits: ";}
      ITS->MakeBranch("S");
      ITS->SetTreeAddress();
      ITS->Hits2SDigits();
      if (gDEBUG) {cout<<"done"<<endl;}
    }

// TPC
    if (iTPC) {
      if (gDEBUG) {cout<<"  Create TPC sdigits: ";}
      if (iTPC == 1) {
// process all sectors
	TPC->SetActiveSectors(1);
	if (gDEBUG) {cout<<"All TPC sectors set active."<<endl;}
      } else if (iTPC == 2) {
// process only sectors with hits
	TPC->SetActiveSectors(0);
	if (gDEBUG) {
	  printf("\nActive sectors\n");
	  Int_t iActive = 0;
	  for (Int_t i=0;i<72;i++) {
	    if (TPC->IsSectorActive(i)) {
	      printf("%2d ",i);
	      iActive++;
	      if (iActive%10 == 0) printf("\n");
	    }
	  }
	  printf("\n");
	}
      }    
      TPC->Hits2SDigits2(iEvent);
      if (gDEBUG) {cout<<"done"<<endl;}
    }

// TRD
    if (iTRD) {
      if (gDEBUG) {cout<<"  Create TRD sdigits: ";}
      sdTRD->InitOutput(fileSDigits, iEvent);
      sdTRD->MakeDigits();
      sdTRD->WriteDigits();
      if (gDEBUG) {cout<<"done"<<endl;}
    }
    
  } // end of loop over events

// PHOS processes always all events
  if (iPHOS) {
    sdPHOS->ExecuteTask("deb all");
  }

// TOF does its own loop
  if (iTOF) {
    sdTOF->Exec("");
  }


//
// finish 
//
  timer.Stop(); 
  timer.Print();

  if (iTRD) { 
    fileSDigits->cd();
    sdTRD->GetParameter()->Write();
  }

  gFileHits->cd();
  delete gAlice;
  gAlice = 0;
  if (!gSameFiles) {
    gFileHits->Close();
    delete gFileHits;
  }

}
 

////////////////////////////////////////////////////////////////////////
TFile* Init(TString fileNameSDigits, TString fileNameHits) {
// open input file, read in gAlice, prepare output file
  if (gAlice) delete gAlice;
  gAlice = 0;

  Bool_t gSameFiles = kFALSE;
  if (fileNameSDigits == fileNameHits || fileNameSDigits == "") gSameFiles = kTRUE;

  TString fileMode = "read";
  if (gSameFiles) fileMode = "update";

  gFileHits =  TFile::Open(fileNameHits.Data(),fileMode.Data());
  if (!gFileHits->IsOpen()) {
    cerr<<"Can't open "<<fileNameHits.Data()<<" !\n";
    return 0;
  }
  if (!ImportgAlice(gFileHits)) return 0;
  if (!gSameFiles) return gAlice->InitTreeFile("S",fileNameSDigits.Data());
  return gFileHits;

}

////////////////////////////////////////////////////////////////////////
TFile* OpenFile(TString fileName) {
// open file fileName
  TFile *file = TFile::Open(fileName.Data());
  if (!file->IsOpen()) {
    cerr<<"Can't open "<<fileName.Data()<<" !\n";
    return 0;
  }
  return file;
}

////////////////////////////////////////////////////////////////////////
Bool_t ImportgAlice(TFile *file) {
// read in gAlice object from the file
  gAlice = (AliRun*)file->Get("gAlice");
  if (!gAlice)  return kFALSE;
  return kTRUE;
}
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
  if (!sdTRD->MakeBranch()) {
    cerr<<"Problems with TRD digitizer initialization."<<endl;
  }
  return sdTRD;
}
////////////////////////////////////////////////////////////////////////
void AliCopy(TFile *inputFile, TFile *outputFile) {
// copy some objects

// copy gAlice
  if (gDEBUG) cout<<"Copy gAlice: ";
  outputFile->cd();
  gAlice->Write();
  if (gDEBUG) cout<<"done"<<endl;

  TTree *treeE  = gAlice->TreeE();
  if (!treeE) {
    cerr<<"No TreeE found "<<endl;
    return;
  }      

// copy TreeE
  if (gDEBUG) cout<<"Copy TreeE: ";
  AliHeader *header = new AliHeader();
  treeE->SetBranchAddress("Header", &header);
  treeE->SetBranchStatus("*",1);
  TTree *treeENew =  treeE->CloneTree();
  treeENew->Write();
  if (gDEBUG) cout<<"done"<<endl;

// copy AliceGeom
  if (gDEBUG) cout<<"Copy AliceGeom: ";
  TGeometry *AliceGeom = static_cast<TGeometry*>(inputFile->Get("AliceGeom"));
  if (!AliceGeom) {
    cerr<<"AliceGeom was not found in the input file "<<endl;
    return;
  }
  AliceGeom->Write();
  if (gDEBUG) cout<<"done"<<endl;

}
