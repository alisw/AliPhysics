#if !defined(__CINT__) || defined(__MAKECINT__)

#include "iostream.h"
#include "TDatime.h"
#include "TFile.h"
#include "TString.h"
#include "../STEER/AliRun.h"
#include "../STEER/AliRunDigitizer.h"
#include "ITS/AliITSDigitizer.h"
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

Int_t AliITSd2r(TString df="galice.root",TString rf="",TString opt=""){
    // Produce ITS RecPoints from Digits.

    // Dynamically link some shared libs
    if (gClassTable->GetID("AliRun") < 0) {
	gROOT->LoadMacro("loadlibs.C");
	loadlibs();
    } // end if
    gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSstandard.C");

    TFile *dfp = 0, *rfp = 0;
    if(!GaliceITSok()){
	// gAlice not define. Must open a file and read it in.
	if(rf.CompareTo(df) == 0 || rf.CompareTo("") == 0) {
	    // Input file = output file
	    dfp = AccessFile(df,"U");  // input file open for update.
	}else{ // Input file different from output file.
	    dfp = AccessFile(df,"R"); // input file open as read only
	    // open output file and create TreeR on it
	    rfp = new TFile(rp,"NEW");
	} // end if
    } // end if !GALICEITSOK()
    AliITS *ITS = (AliITS*)gAlice->GetDetector("ITS"); 

    ChangeITSDefaults(dfp,ITS,opt);
    // write the AliRun object to the output file if different from input file.
    if(rfp) {
	writeAR(dfp,rfp);
	rfp->Close(); // Manager will open in update mode.
    } // end if

    TStopwatch timer;

#ifdef DEBUG
    cout << "Creating reconstructed points from digits for the ITS..." << endl;
#endif
    if(dfp!=0)
	AliITSreconstruction *itsr = new AliITSreconstruction(df.Data());
    else
	AliITSreconstruction *itsr = new AliITSreconstruction(0);
    if(df.CompareTo(rf)!=0 && rf.CompareTo("")!=0) {
	itsr->SetOutputFile(rf.Data());
    } // end if
    timer.Start();
    itsr->Init();
    itsr->Exec(); 
    timer.Stop(); 
    timer.Print();
    delete itsr;

    if(rfp!=0){
	cout <<"After Reconstruction "<<rf<<" size ="<< rfp->GetSize()<< endl;
    }else if(dfp!=0){
	cout <<"After Reconstruction "<<df<<" size ="<< dfp->GetSize()<< endl;
    } // end if sf1!=0
    return 0;
}
