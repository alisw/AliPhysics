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
//
// input:
//       TString input ... galice input file
//       Int_t nEvents  ... how many events to process
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
#include "STEER/AliRun.h"
#include "STEER/AliRunDigitizer.h"
#include "ITS/AliITSDigitizer.h"
#include "ITS/AliITSFDigitizer.h"
#include "TPC/AliTPCDigitizer.h"
#include "TRD/AliTRDdigitizer.h"
#include "PHOS/AliPHOSDigitizer.h"
#include "MUON/AliMUONDigitizer.h"
#include "RICH/AliRICHDigitizer.h"
#include "TStopwatch.h"
#endif

Int_t AliSDigits2Digits(TString input="galice.root", 
                        Int_t nEvents = 1, Int_t iITS = 0, Int_t iTPC = 0,
                        Int_t iTRD = 0,  Int_t iPHOS = 0, Int_t iMUON = 0,
                        Int_t iRICH = 0)
{
// delete the current gAlice object, the one from input file
//  will be used

  if(gAlice){
    delete gAlice;
    gAlice = 0;
  } // end if gAlice
  AliRunDigitizer * manager = new AliRunDigitizer(1,1);
  manager->SetDebug(1000);
  manager->SetInputStream(0,input);
  
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

