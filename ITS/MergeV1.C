////////////////////////////////////////////////////////////////////////
//
// name: MergeV1
// date: 11.4.2002
// last update: 11.4.2002
// author: Jiri Chudoba
// version: 1.0
//
// description: 
//       creates digits from sdigits for several detectors
//       stores sdigits in separate file (or in the source file
//       with sdigits). Stores gAlice object and copies TE to the
//       file with digits
//       ITS region of Interest is set
//       test
//
// input:
//       TString fileNameSDigits ... input file with sdigits
//       TString fileNameDigits ... output file with digits
//       Int_t nEvents  ... how many events to process
//       Int_t ITS, TPC, ...   many flags for diff. detectors
//
// History:
//
// 04.04.02 - first version
// 
////////////////////////////////////////////////////////////////////////

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "iostream.h"
#include "TDatetime.h"
#include "STEER/AliRun.h"
#include "STEER/AliRunDigitizer.h"
#include "ITS/AliITSDigitizer.h"
#include "ITS/AliITS.h"
#include "ITS/AliITSDetType.h"
#include "ITS/AliITSresponseSDD.h"
#include "TPC/AliTPCDigitizer.h"
#include "TRD/AliTRDdigitizer.h"
#include "PHOS/AliPHOSDigitizer.h"
#include "MUON/AliMUONDigitizer.h"
#include "RICH/AliRICHDigitizer.h"
#include "TStopwatch.h"
#endif

// #include "AliHits2SDigits.C"

// void AliCopyN(TString inputFile, TString outputFile);

Int_t MergeV1(TString fileNameDigits="digits.root", 
	    TString fileNameSDigitsSig="sig.sdigits.root", 
	    TString fileNameSDigitsBgr="bgr.sdigits.root", 
	    Int_t nEvents = 1, Int_t iITS = 2, Int_t iTPC = 0,
	    Int_t iTRD = 0,  Int_t iPHOS = 0, Int_t iMUON = 0,
	    Int_t iRICH = 0, Int_t iCopy = 1)
{
// delete the current gAlice object, the one from input file
//  will be used

  if(gAlice){
    delete gAlice;
    gAlice = 0;
  } // end if gAlice

  // Connect the Root Galice file containing Geometry, Kine and Hits
  TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(fileNameSDigitsSig.Data());
  if(!file) file = new TFile(fileNameSDigitsSig.Data());
  TDatime *ct0 = new TDatime(2002,04,26,00,00,00), ct = file->GetCreationDate();
  
 
  // Get AliRun object from file or create it if not on file
  if(!gAlice) {
      gAlice = (AliRun*)file->Get("gAlice");
      if(gAlice) printf("AliRun object found on file\n");
      if(!gAlice) gAlice = new AliRun("gAlice","Alice test program");
  } // end if !gAlice

  AliRunDigitizer * manager = new AliRunDigitizer(2,1);
  manager->SetInputStream(0,fileNameSDigitsSig.Data());
  manager->SetInputStream(1,fileNameSDigitsBgr.Data());
  if (fileNameDigits != "") {
//    if (iCopy) {
//      AliCopyN(fileNameSDigitsSig,fileNameDigits);
//    }
    manager->SetOutputFile(fileNameDigits);
  }
  manager->SetNrOfEventsToWrite(nEvents);
  
  if (iITS) {
    AliITSDigitizer *dITS  = new AliITSDigitizer(manager);
    if (iITS == 2) dITS->SetByRegionOfInterestFlag(1);
    if(ct0->GetDate()>ct.GetDate()){
	// For old files, must change SDD noise.
	AliITS *ITS = (AliITS*) gAlice->GetDetector("ITS");
	AliITSresponseSDD *resp1 = ITS->DetType(1)->GetResponseModel();
	resp1->SetNoiseParam();
	resp1->SetNoiseAfterElectronics();
	Float_t n,b;
	Int_t cPar[8];
	resp1->GetNoiseParam(n,b);
	n = resp1->GetNoiseAfterElectronics();
	cPar[0]=0;
	cPar[1]=0;
	cPar[2]=(Int_t)(b + 2.*n + 0.5);
	cPar[3]=(Int_t)(b + 2.*n + 0.5);
	cPar[4]=0;
	cPar[5]=0;
	cPar[6]=0;
	cPar[7]=0;
	resp1->SetCompressParam(cPar);
    } // end if
  }
  if (iTPC) AliTPCDigitizer *dTPC  = new AliTPCDigitizer(manager);
  if (iTRD) AliTRDdigitizer *dTRD  = new AliTRDdigitizer(manager);
  if (iPHOS) AliPHOSDigitizer *dPHOS  = new AliPHOSDigitizer(manager);
  if (iMUON) AliMUONDigitizer *dMUON  = new AliMUONDigitizer(manager);
  if (iRICH) AliRICHDigitizer *dRICH  = new AliRICHDigitizer(manager);
  TStopwatch timer;
  timer.Start();
  manager->Exec("deb all");
  timer.Stop(); 
  timer.Print();
//  delete gAlice;  
//  gAlice = 0;
  delete manager;
}


/*
////////////////////////////////////////////////////////////////////////
void AliCopyN(TString inputFileName, TString outputFileName) {
// copy some objects

  TFile *inputFile = OpenFile(inputFileName);
  if (!inputFile) return;

  TFile *outputFile = TFile::Open(outputFileName.Data(),"update");
  if (!outputFile->IsOpen()) {
    cerr<<"Can't open "<<outputFileName.Data()<<" !\n";
    return;
  }
  if (!ImportgAlice(inputFile)) return;
  AliCopy(inputFile, outputFile);
  delete gAlice;
  gAlice=0;
  inputFile->Close();
  delete inputFile;
  outputFile->Close();
  delete outputFile;
}
////////////////////////////////////////////////////////////////////////
*/
