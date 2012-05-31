#if !defined(__CINT__) || defined(__MAKECINT__)

#include "iostream.h"
#include "TDatime.h"
#include "TFile.h"
#include "TString.h"
#include "../STEER/AliRun.h"
#include "../STEER/AliRunDigitizer.h"
//#include "AliITSDigitizer.h"
#include "AliITS.h"
#include "AliITSDetType.h"
#include "AliITSresponseSDD.h"
#include "TStopwatch.h"

Bool_t GaliceITSok();
TFile* AccessFile(TString inFile="galice.root", TString acctype="R");
void writeAR(TFile * fin, TFile *fou);
Int_t ChangeITSDefaults(TFile *hitfile,AliITS *ITS,TString opt="");
void loadlibs();

#endif


//#define DEBUG

Int_t AliITSh2sd(TString hf="galice.root",TString sf="",TString opt=""){
    // Produce ITS SDigits from Hits.

    // Dynamically link some shared libs
    if (gClassTable->GetID("AliRun") < 0) {
	gROOT->LoadMacro("loadlibs.C");
	loadlibs();
    } // end if
    gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSstandard.C");

    TFile *hfp = 0, *sfp = 0;
    if(!GaliceITSok()){
	// gAlice not define. Must open a file and read it in.
	if(hf.CompareTo(sf) == 0 || sf.CompareTo("") == 0) {
	    // Input file = output file
	    hfp = AccessFile(hf,"U");  // input file open for update.
	}else{ // Input file different from output file.
	    hfp = AccessFile(hf,"R"); // input file open as read only
	    // open output file and create TreeR on it
	    sfp = gAlice->InitTreeFile("S",sf);
	} // end if
    } // end if !GALICEITSOK()
    AliITS *ITS = (AliITS*) (gAlice->GetDetector("ITS"));

    ChangeITSDefaults(hfp,ITS,opt);
    // write the AliRun object to the output file if different from input file.
    if(sfp) writeAR(hfp,sfp);

    TStopwatch timer;
    Int_t evNumber1 = 0;
    Int_t evNumber2 = gAlice->GetEventsPerRun();
    timer.Start();
    for(Int_t event = evNumber1; event < evNumber2; event++){
	gAlice->GetEvent(event);
	if(!gAlice->TreeS() && sfp == 0){ 
	    cout << "Having to create the SDigits Tree." << endl;
	    gAlice->MakeTree("S");
	} // end if !gAlice->TreeS()
	if(sfp) gAlice->MakeTree("S",sfp);
	//    make branch
	ITS->MakeBranch("S");
	ITS->SetTreeAddress();
#ifdef DEBUG
	cout<<"Making ITS SDigits for event "<<event<<endl;
#endif
	ITS->Hits2SDigits();
    } // end for event
    timer.Stop();
    timer.Print();
    if(sfp!=0){
	cout << sf << " size =" << sfp->GetSize() << endl;
    }else if(hfp!=0){
	cout << hfp << " size =" << hfp->GetSize() << endl;
    } // end if sfp!=0
    return 0;
}
