#include "iostream.h"
#include "TFile.h"
#include "TString.h"
#include "TClonesArray.h"
/*
#include "$(ALICE_ROOT)/STEER/AliRun.h"
#include "$(ALICE_ROOT)/ITS/AliITS.h"
#include "$(ALICE_ROOT)/ITS/AliITSgeom.h"
#include "$(ALICE_ROOT)/ITS/AliITSHit.h"
*/
void AliITSPrintHits(TString hfn="galice.root",Int_t mod=-1,
			  Int_t evnt=-1){
    // Macro to print out the recpoints for all or a specific module

    // Dynamically link some shared libs
    if (gClassTable->GetID("AliRun") < 0) {
        gROOT->LoadMacro("loadlibs.C");
        loadlibs();
    } // end if
    gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSstandard.C");

    TFile *hf=0;
    hf = AccessFile(hfn,"R"); // Set up to read in Data
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
    TClonesArray *hpa = ITS->RecPoints();
    AliITShit *hp = 0;

    Int_t nmodules,size=-1;
    ITS->InitModules(size,nmodules);
    Int_t event,m,i,i2,hit,trk;
    for(event = evNumber1; event < evNumber2; event++){
        gAlice->GetEvent(event);
	ITS->FillModules(event,0,-1," "," ");
	for(m=mod1;m<mod2;m++){
	    i2 = (ITS->GetModule(m))->GetNhits();
	    cout <<  "Event=" << event << " module=" << m <<
		" Number of Hits=" << i2 <<endl;
	    for(i=0;i<i2;i++){
		trk = (ITS->GetModule(m))->GetHitTrackIndex(i);
		hit = (ITS->GetModule(m))->GetHitHitIndex(i);
		hp  = (ITS->GetModule(m))->GetHit(i);
		cout << i << " trk#="<<trk<<" hit#="<< hit << " ";
		hp->Print((ostream*)cout);
		cout << endl;
	    } // end for i
	} // end for m
	ITS->ClearModules();
    } // end for event
}
