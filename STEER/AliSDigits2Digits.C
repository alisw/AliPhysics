////////////////////////////////////////////////////////////////////////
//
// name: AliSDigits2Digits
// date: 4.4.2002
// last update: 4.4.2002
// author: Jiri Chudoba
// version: 1.0
//
// description: 
//       creates digits from sdigits for several detectors
//       stores sdigits in separate file (or in the source file
//       with sdigits). Stores gAlice object and copies TE to the
//       file with digits
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
#include "STEER/AliRun.h"
#include "STEER/AliRunDigitizer.h"
#include "ITS/AliITSDigitizer.h"
#include "TPC/AliTPCDigitizer.h"
#include "TRD/AliTRDdigitizer.h"
#include "PHOS/AliPHOSDigitizer.h"
#include "MUON/AliMUONDigitizer.h"
#include "RICH/AliRICHDigitizer.h"
#include "TStopwatch.h"
#endif

#include "AliHits2SDigits.C"

void AliCopyN(TString inputFile, TString outputFile);

Int_t AliSDigits2Digits(TString fileNameDigits="digits.root", 
			TString fileNameSDigits="rfio:sdigits.root", 
			Int_t nEvents = 1, Int_t iITS = 0, Int_t iTPC = 0,
			Int_t iTRD = 0,  Int_t iPHOS = 0, Int_t iMUON = 0,
			Int_t iRICH = 0, Int_t iCopy = 1)
{
// delete the current gAlice object, the one from input file
//  will be used

  if(gAlice){
    delete gAlice;
    gAlice = 0;
  } // end if gAlice
  AliRunDigitizer * manager = new AliRunDigitizer(1,1);
  manager->SetInputStream(0,fileNameSDigits.Data());
  if (fileNameDigits != "") {
    if (iCopy) {
      AliCopyN(fileNameSDigits,fileNameDigits);
    }
    manager->SetOutputFile(fileNameDigits);
  }
  manager->SetNrOfEventsToWrite(nEvents);
  if (iITS == 1) AliITSDigitizer *dITS  = new AliITSDigitizer(manager);
  if (iITS == 2) AliITSFDigitizer *dITS  = new AliITSFDigitizer(manager);
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
  delete manager;
}


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
  inputFile->Close();
  delete inputFile;
  outputFile->Close();
  delete outputFile;
}
////////////////////////////////////////////////////////////////////////
