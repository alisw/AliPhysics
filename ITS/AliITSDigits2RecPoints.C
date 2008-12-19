#if !defined(__CINT__) || defined(__MAKECINT__)

#include <Riostream.h>
#include <TClassTable.h>
#include <TDatime.h>
#include <TGeoManager.h>
#include <TROOT.h>
#include <TStopwatch.h>
#include "AliRun.h"
#include "AliRunDigitizer.h"
#include "AliITSDetTypeRec.h"
#include "AliITSDigitizer.h"
#include "AliITS.h"
#include "AliITSresponseSDD.h"
#include "AliITSreconstruction.h"

#endif
#define DEBUG

Int_t AliITSDigits2RecPoints(TString filename="galice.root",TString fileRP=""){
  // Standard ITS Digits to RecPoints.

  // Get geometry
  TGeoManager::Import("geometry.root");

  // Dynamically link some shared libs
  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->ProcessLine(".x $(ALICE_ROOT)/macros/loadlibs.C");
  }else if (gAlice){
    delete AliRunLoader::GetRunLoader();
    delete gAlice;
    gAlice=0;
  } // end if

  TStopwatch timer;
#ifdef DEBUG
  cout << "Creating reconstructed points from digits for the ITS..." << endl;
#endif
  AliITSreconstruction *itsr = new AliITSreconstruction(filename);

  timer.Start();
  if(!(fileRP.IsNull()))itsr->SetOutputFile(fileRP);
  itsr->Init();
  itsr->Exec(); 
  timer.Stop(); 
  timer.Print();
  delete itsr;
  return 0;
}
