////////////////////////////////////////////////////////////////////////
//
// name: AliTOFtestDigitizer
// date: 11-VI-2002
// last update: 11-VI-2002
// author: F. Pierella | pierella@bo.infn.it
// version: 1.0
//
// description: 
//       creates digits from sdigits for TOF detector
//       stores sdigits in separate file (or in the source file
//       with sdigits). Stores gAlice object and copies TE to the
//       file with digits
//
// input:
//       char* fileNameSignal ... input file with sdigits
//       TString fileNameDigits ... output file with digits
//       Int_t nEvents  ... how many events to process
//
// Updated to the new I/O: C. Zampolli
//
/////////////////////////////////////////////////////////////////////

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "iostream.h"
#include "AliTOFDigitizer.h"
#include "../STEER/AliRunDigitizer.h"
#include "../STEER/AliDigitizer.h"
#include "TStopwatch.h"
#endif

Int_t AliTOFtestDigitizer(const char* fileNameSignal = "galice.root",
			  /*const char* fileNameSignal = "signal/galice.root",
			    const char* fileNameBkgrd = "bkgrd/galice.root",*/
			  Int_t nEvents = -1, Int_t signalPerBkgrd = 1,
			  Int_t iTOF = 1)
{
  
// delete the current gAlice object, the one from input file will be used
  if(gAlice)
    {
      delete AliRunLoader::GetRunLoader();
      delete gAlice;
      gAlice=0x0;
    }

  AliRunDigitizer * manager = new AliRunDigitizer(1/*2*/,signalPerBkgrd);
  manager->SetInputStream(0, fileNameSignal);
  //manager->SetInputStream(1, fileNameBkgrd);
  //manager->SetOutputFile(fileNameSignal);
  if (nEvents >= 0) manager->SetNrOfEventsToWrite(nEvents);

  if (iTOF) AliTOFDigitizer *dTOF = new AliTOFDigitizer(manager);

  TStopwatch timer;
  timer.Start();
  manager->Exec("deb all");
  timer.Stop(); 
  timer.Print();

//  delete manager;
//  manager = 0x0;
  
  if(gAlice)
    {
      delete AliRunLoader::GetRunLoader();
      delete gAlice;
      gAlice=0x0;
    }

  return 0;

}
