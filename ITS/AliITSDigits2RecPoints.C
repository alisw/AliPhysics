#if !defined(__CINT__) || defined(__MAKECINT__)

#include "iostream.h"
#include "TDatetime.h"
#include "STEER/AliRun.h"
#include "STEER/AliRunDigitizer.h"
#include "ITS/AliITSDigitizer.h"
#include "ITS/AliITS.h"
#include "ITS/AliITSDetType.h"
#include "ITS/AliITSresponseSDD.h"
#include "TStopwatch.h"

#endif
#define DEBUG
Int_t AliITSDigits2RecPoints(TString filename="galice.root"){
    // Standard ITS Digits to RecPoints.

    // Dynamically link some shared libs
    if (gClassTable->GetID("AliRun") < 0) {
     gROOT->LoadMacro("loadlibs.C");
     loadlibs();
    }else if (gAlice){
      delete gAlice->GetRunLoader();
      delete gAlice;
      gAlice=0;
     } // end if

    TStopwatch timer;
#ifdef DEBUG
    cout << "Creating reconstructed points from digits for the ITS..." << endl;
#endif
    AliITSreconstruction *itsr = new AliITSreconstruction(filename);

    timer.Start();
    itsr->Init();
    itsr->Exec(); 
    timer.Stop(); 
    timer.Print();
    delete itsr;
    return 0;
}
