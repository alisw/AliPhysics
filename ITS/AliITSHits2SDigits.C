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

Int_t AliITSHits2SDigits(TString  filename = "galice.root")
 {
    // Standeard ITS Hits to SDigits.

    // Dynamically link some shared libs
    if (gClassTable->GetID("AliRun") < 0) {
     gROOT->LoadMacro("loadlibs.C");
     loadlibs();
    } else if (gAlice){
       delete gAlice->GetRunLoader();
       delete gAlice;
       gAlice=0;
     }

    // Connect the Root Galice file containing Geometry, Kine and Hits

    AliRunLoader* rl = AliRunLoader::Open(filename);
    if (rl == 0x0)
     {
      cerr<<"AliITSHits2DigitsDefault.C : Can not open session RL=NULL"
           << endl;
       return 3;
     }
     
    Int_t retval = rl->LoadgAlice();
    if (retval)
     {
      cerr<<"AliITSHits2DigitsDefault.C : LoadgAlice returned error"
           << endl;
       return 3;
     }
    gAlice=rl->GetAliRun();
    AliLoader* gime = rl->GetLoader("ITSLoader");
    if (gime == 0x0)
     {
      cerr<<"AliITSHits2DigitsDefault.C : can not get ITS loader"
           << endl;
     }
    AliITS *ITS = (AliITS*)gAlice->GetDetector("ITS");      
    if (!ITS) {
	cerr<<"AliITSHits2DigitsDefault.C : AliITS object not found on file"
	    << endl;
	return 3;
    }  // end if !ITS
    if(!(ITS->GetITSgeom())){
       cerr << " AliITSgeom not found. Can't digitize with out it." << endl;
       return 4;
    } // end if

    TStopwatch timer;
    Int_t evNumber1 = 0;
    Int_t evNumber2 = gAlice->GetEventsPerRun();
    timer.Start();
    retval = gime->LoadHits();
    if (retval)
     {
      cerr<<"AliITSHits2DigitsDefault.C : ITSLoader::LoadHits returned error"
           << endl;
       return 3;
     }
    retval = gime->LoadSDigits("recreate");
    if (retval)
     {
      cerr<<"AliITSHits2DigitsDefault.C : ITSLoader::LoadSDigits returned error"
           << endl;
       return 3;
     }
    for(Int_t event = evNumber1; event < evNumber2; event++){
       rl->GetEvent(event);
       if(!gime->TreeS()){ 
           cout << "Having to create the SDigits Tree." << endl;
           gime->MakeTree("S");
       } // end 

       ITS->MakeBranch("S");
       ITS->SetTreeAddress();
       cout<<"Making ITS SDigits for event "<<event<<endl;
       ITS->Hits2SDigits();
    } // end for event
    timer.Stop();
    timer.Print();

    delete rl; // sdigfile is closed by deleting gAlice if != hitfile.
}
