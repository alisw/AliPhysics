#if !defined(__CINT__) || defined(__MAKECINT__)

#include "iostream.h"
#include "TDatime.h"
#include "TDatime.h"
#include "TFile.h"
#include "TString.h"
#include "../STEER/AliRun.h"
#include "../STEER/AliRunDigitizer.h"
//#include "ITS/AliITSDigitizer.h"
#include "ITS/AliITS.h"
#include "ITS/AliITSDetType.h"
#include "ITS/AliITSresponseSDD.h"
#include "TStopwatch.h"

Bool_t GaliceITSok();
TFile* AccessFile(TString inFile="galice.root", TString acctype="R");
void writeAR(TFile * fin, TFile *fou);
Int_t ChangeITSDefaults(TFile *hitfile,AliITS *ITS,TString opt="");
void loadlibs();

#endif

//#define DEBUG

Int_t AliITSh2d(TString hf="galice.root",TString df="",TString opt=""){
    // Produce ITS SDigits from Hits.

    // Dynamically link some shared libs
    if (gClassTable->GetID("AliRun") < 0) {
	gROOT->LoadMacro("loadlibs.C");
	loadlibs();
    } // end if
    gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSstandard.C");

    TFile *hfp = 0, *dfp = 0;
    if(!GaliceITSok()){
	// gAlice not define. Must open a file and read it in.
	if(hf.CompareTo(df) == 0 || df.CompareTo("") == 0) {
	    // Input file = output file
	    hfp = AccessFile(hf,"U");  // input file open for update.
	}else{ // Input file different from output file.
	    hfp = AccessFile(hf,"R"); // input file open as read only
	    // open output file and create TreeR on it
	    dfp = gAlice->InitTreeFile("D",df);
	} // end if
    } // end if !GALICEITSOK()
    AliITS *ITS = (AliITS*)gAlice->GetDetector("ITS");

    ChangeITSDefaults(hfp,ITS,opt);
    // write the AliRun object to the output file if different from input file.
    if(dfp) writeAR(hfp,dfp);

    TStopwatch timer;
    Int_t evNumber1 = 0;
    Int_t evNumber2 = gAlice->GetEventsPerRun();
    timer.Start();
    for(Int_t event = evNumber1; event < evNumber2; event++){
	gAlice->GetEvent(event);
	if(!gAlice->TreeD() && dfp == 0){ 
	    cout << "Having to create the Digits Tree." << endl;
	    gAlice->MakeTree("S");
	} // end if !gAlice->TreeS()
	if(dfp) gAlice->MakeTree("D",dfp);
	//    make branch
	ITS->MakeBranch("D");
	ITS->SetTreeAddress();
#ifdef DEBUG
	cout<<"Making ITS Digits for event "<<event<<endl;
#endif
	ITS->Hits2Digits();
    } // end for event
    timer.Stop();
    timer.Print();
    if(dfp!=0){
	cout <<"After Digitization "<<df<< " size =" << dfp->GetSize() << endl;
    }else if(hfp!=0){
	cout <<"After digitization "<<hf<< " size =" << hfp->GetSize() << endl;
    } // end if dfp!=0
    return 0;
}
