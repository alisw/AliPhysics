#if !defined(__CINT__) || defined(__MAKECINT__)

#include "Riostream.h"
#include "TDatime.h"
#include "TClassTable.h"
#include "AliRun.h"
#include "AliRunDigitizer.h"
#include "AliITSDigitizer.h"
#include "AliITS.h"
#include "AliITSDetType.h"
#include "AliITSLoader.h"
#include "AliITSresponseSDD.h"
#include "TStopwatch.h"

#endif

//#define DEBUG
Int_t AliITSHits2Digits(TString inFile = "galice.root"){
    // Standard ITS Hits to Digits, excluding creation of SDigits.

    // Dynamically link some shared libs
   if (gClassTable->GetID("AliRun") < 0) {
      gROOT->ProcessLine(".x $(ALICE_ROOT)/macros/loadlibs.C");
    }else if (gAlice){
      delete AliRunLoader::GetRunLoader();
      delete gAlice;
      gAlice=0;
     } 

   AliRunLoader* rl = AliRunLoader::Open(inFile.Data());
    if (rl == 0x0)
     {
      cerr<<"AliITSHits2Digits.C : Can not open session RL=NULL"
           << endl;
       return 3;
     }
     
    Int_t retval = rl->LoadgAlice();
    if (retval)
     {
      cerr<<"AliITSHits2Digits.C : LoadgAlice returned error"
           << endl;
       return 3;
     }
    gAlice=rl->GetAliRun();
    AliITSLoader* gime = (AliITSLoader*)rl->GetLoader("ITSLoader");
    if (gime == 0x0)
     {
      cerr<<"AliITSHits2Digits.C : can not get ITS loader"
           << endl;
     }
    AliITS *ITS = (AliITS*)gAlice->GetDetector("ITS");      
    if (!ITS) {
	cerr<<"AliITSHits2Digit.C : AliITS object not found on file"
	    << endl;
	return 3;
    }  // end if !ITS

    if(!(ITS->GetITSgeom())){
	cerr << " AliITSgeom not found. Can't digitize without it." << endl;
	return 4;
    } // end if

    

    TStopwatch timer;
    Int_t evNumber1 = 0;
    Int_t evNumber2 = AliRunLoader::GetNumberOfEvents();
    timer.Start();
    retval = gime->LoadHits();
    if (retval)
     {
      cerr<<"AliITSHits2Digits.C : ITSLoader::LoadHits returned error"
           << endl;
       return 3;
     }

    retval = gime->LoadDigits("recreate");
    if (retval)
     {
      cerr<<"AliITSHits2Digits.C : ITSLoader::LoadDigits returned error"
           << endl;
       return 3;
     }
    for(Int_t nevent = evNumber1; nevent < evNumber2; nevent++){
	// cout<<"Producing Digits for event n."<<nevent<<endl;

	rl->GetEvent(nevent);
	if(!gime->TreeD()){ 
	    cout << "Having to create the Digits Tree." << endl;
	    gime->MakeTree("D");
	} // end if creating digits tree
	ITS->MakeBranch("D");
	ITS->SetTreeAddress();   
	ITS->Hits2Digits();
    } // end for nevent
    timer.Stop();
    timer.Print();

    delete rl;
    return 0;
}
