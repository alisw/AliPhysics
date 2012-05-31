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

#endif

//#define DEBUG

Int_t AliITSsd2d(TString df="galice.root",TString sf1="galice.root",
		 TString sf2="",TString opt=""){
    // Produce ITS Digits from SDigits, with psible merging.

    // Dynamically link some shared libs
    if (gClassTable->GetID("AliRun") < 0) {
	gROOT->LoadMacro("loadlibs.C");
	loadlibs();
    } // end if
    gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSstandard.C");

    TFile *sfp1 = 0, *sfp2 = 0, *dfp = 0;
    if(!GaliceITSok()){
	// gAlice not define. Must open a file and read it in.
	if(df.CompareTo(sf1) == 0 || df.CompareTo("") == 0) {
	    // Input file = output file
	    sfp1 = AccessFile(sf1,"U");  // input file open for update.
	}else{ // Input file different from output file.
	    sfp1 = AccessFile(sf1,"R"); // input file open as read only
	    // open output file and create TreeR on it
	    dfp = gAlice->InitTreeFile("S",df);
	} // end if
    } // end if !GALICEITSOK()
    AliITS *ITS = (AliITS*)gAlice->GetDetector("ITS"); 

    ChangeITSDefaults(sfp1,ITS,opt);
    // write the AliRun object to the output file if different from input file.
    if(dfp) {
	writeAR(sfp1,dfp);
	dfp->Close(); // Manager will open in update mode.
    } // end if

    AliRunDigitizer *manager;
    if(sf2.CompareTo("")==0) { // do not merge.
	manager = new AliRunDigitizer(1,1);
	manager->SetInputStream(0,sf1);
    }else{
	manager = new AliRunDigitizer(2,1);
	manager->SetInputStream(0,sf1.Data());
	manager->SetInputStream(1,sf2.Data());
    } // end if
    if (df.CompareTo(sf1) !=0) {
	manager->SetOutputFile(df);
    } // end if
    AliITSDigitizer *dITS = new AliITSDigitizer(manager);
    if(opt.Contains("ROI")==0) dITS->SetByRegionOfInterestFlag(1);

    TStopwatch timer;
    timer.Start();
    manager->Exec("all");
    timer.Stop(); 
    timer.Print();
    delete manager;

    if(dfp!=0){
	cout << df << " size =" << df->GetSize() << endl;
    }else if(sfp1!=0){
	cout << sf1 << " size =" << sf1->GetSize() << endl;
    } // end if sf1!=0
    return 0;
}
