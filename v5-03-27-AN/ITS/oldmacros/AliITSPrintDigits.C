#include "iostream.h"
#include "TFile.h"
#include "TString.h"
#include "TClonesArray.h"
/*
#include "$(ALICE_ROOT)/STEER/AliRun.h"
#include "$(ALICE_ROOT)/ITS/AliITS.h"
#include "$(ALICE_ROOT)/ITS/AliITSgeom.h"
#include "$(ALICE_ROOT)/ITS/AliITSdigit.h"
*/
void AliITSPrintDigits(TString rfn="galice.root",Int_t mod=-1,
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
    TClonesArray   *pda = 0;
    AliITSdigitSPD *dp0 = 0;
    AliITSdigitSDD *dp1 = 0;
    AliITSdigitSSD *dp2 = 0;

    Int_t event,m,i,i2,id;
    for(event = evNumber1; event < evNumber2; event++){
        gAlice->GetEvent(event);
	for(m=mod1;m<mod2;m++){
	    id = ITS->GetITSgeom()->GetModuleType(m);
	    dpa = ITS->DigitsAddress(id);
	    ITS->ResetDigits();
	    gAlice->TreeD()->GetEvent(m);
	    i2 = dpa->GetEntriesFast();
	    switch (id) {
	    case 0:
		cout <<  "Event=" << event << " module=" << m <<
		    " Number of SPD Digits=" << i2 <<endl;
		for(i=0;i<i2;i++){
		    dp0 = (AliITSdigitSPD*)(dpa->At(i));
		    cout << i << " ";
		    dp0->Print((ostream*)cout);
		    cout << endl;
		} // end for i
		break;
	    case 1:
		cout <<  "Event=" << event << " module=" << m <<
		    " Number of SDD Digits=" << i2 <<endl;
		for(i=0;i<i2;i++){
		    dp1 = (AliITSdigitSDD*)(dpa->At(i));
		    cout << i << " ";
		    dp1->Print((ostream*)cout);
		    cout << endl;
		} // end for i
		break;
	    case 2:
		cout <<  "Event=" << event << " module=" << m <<
		    " Number of SSD Digits=" << i2 <<endl;
		for(i=0;i<i2;i++){
		    dp2 = (AliITSdigitSSD*)(dpa->At(i));
		    cout << i << " ";
		    dp2->Print((ostream*)cout);
		    cout << endl;
		} // end for i
		break;
	    default:
		break;
	    } // end switch
	} // end for m
    } // end for event
}
