#include "iostream.h"
#include "TFile.h"
#include "TString.h"
#include "TClonesArray.h"
/*
#include "$(ALICE_ROOT)/STEER/AliRun.h"
#include "$(ALICE_ROOT)/ITS/AliITS.h"
#include "$(ALICE_ROOT)/ITS/AliITSgeom.h"
#include "$(ALICE_ROOT)/ITS/AliITSpList.h"
*/
void AliITSPrintSDigits(TString rfn="galice.root",Int_t mod=-1,
			  Int_t evnt=-1){
    // Macro to print out the recpoints for all or a specific module

    // Dynamically link some shared libs
    if (gClassTable->GetID("AliRun") < 0) {
        gROOT->LoadMacro("loadlibs.C");
        loadlibs();
    } // end if
    gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSstandard.C");

    TFile *rf=0;
    rf = AccessFile(rfn,"R"); // Set up to read in Data
    AliITS *ITS = (AliITS*)gAlice->GetDetector("ITS");
    if(!ITS){
	cout << "Error: no ITS found. Aborting"<<endl;
	return;
    } // end if !ITS

    Int_t evNumber1 = 0;
    Int_t evNumber2 = gAlice->GetEventsPerRun();
    if(evnt>=0){
	evNumber1 = evnt;
	evNumber2 = evnt+1;
    } // end if evnt>=0
    Int_t mod1 = 0;
    Int_t mod2 = ITS->GetITSgeom()->GetIndexMax();
    if(mod>=0){
	mod1 = mod;
	mod2 = mode+1;
    } // end if mod>=0
    TClonesArray *sdpa = ITS->GetSDigits();
    AliITSpListItem *sdp = 0;

    Int_t event,m,i,i2;
    for(event = evNumber1; event < evNumber2; event++){
        gAlice->GetEvent(event);
	for(m=mod1;m<mod2;m++){
	    ITS->ResetSDigits();
	    gAlice->TreeS()->GetEvent(m);
	    i2 = sdpa->GetEntriesFast();
	    cout <<  "Event=" << event << " module=" << m <<
		" Number of SDigits=" << i2 <<endl;
	    for(i=0;i<i2;i++){
		sdp = (AliITSpListItem*)(sdpa->At(i));
		cout << i << " ";
		sdp->Print((ostream*)cout);
		cout << endl;
	    } // end for i
	} // end for m
    } // end for event

}
