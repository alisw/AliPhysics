////////////////////////////////////////////////////////////////////////
//
// name: AliMerge
// date: 21.07.2003
// last update: 21.07.2003
// author: Thomas Kuhr
// version: 1.0
//
// description: 
//    merges sdigits from a signal and a background event
//    to digits for several detectors
//    (for advanced features like more input streams see 
//     AliRunDigitizer.cxx)
//
// input:
//    const char* fileNameSignal : galice file of the signal event
//    const char* fileNameBkgrd  : galice file of the background event
//    Int_t nEvents              : how many events to process (<0 means all)
//    Int_t signalPerBkgrd       : number of signal events merged with 
//                                 the same background event
//    Int_t ITS, TPC, ...        : flags for diff. detectors
//
// History:
//
// 21.07.03 - first version
// 
////////////////////////////////////////////////////////////////////////

#if !defined(__CINT__) || defined(__MAKECINT__)
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

Bool_t AliMerge(const char* fileNameSignal = "signal/galice.root",
		const char* fileNameBkgrd = "bkgrd/galice.root",
		Int_t nEvents = -1, Int_t signalPerBkgrd = 1,
		Int_t iITS = 0, Int_t iTPC = 0, Int_t iTRD = 0,
		Int_t iPHOS = 0, Int_t iMUON = 0, Int_t iRICH = 0)
{
  if (gAlice) delete gAlice;
  gAlice = NULL;

  AliRunDigitizer* manager = new AliRunDigitizer(2, signalPerBkgrd);
  manager->SetInputStream(0, fileNameSignal);
  manager->SetInputStream(1, fileNameBkgrd);
  if (nEvents >= 0) manager->SetNrOfEventsToWrite(nEvents);

  if (iITS == 1) new AliITSDigitizer(manager);
  if (iITS == 2) new AliITSFDigitizer(manager);
  if (iTPC) new AliTPCDigitizer(manager);
  if (iTRD) new AliTRDdigitizer(manager);
  if (iPHOS) new AliPHOSDigitizer(manager);
  if (iMUON) new AliMUONDigitizer(manager);
  if (iRICH) new AliRICHDigitizer(manager);

  TStopwatch timer;
  timer.Start();
  manager->Exec("deb all");
  timer.Stop(); 
  timer.Print();

  delete manager;
  return kTRUE;
}

