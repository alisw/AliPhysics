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
void AliITSPrintGeom(TString hfn="galice.root",Int_t mod=-1){
    // Macro to print out the information kept in the AliITSgeom class, for 
    //all or a specific module

    // Dynamically link some shared libs
    if (gClassTable->GetID("AliRun") < 0) {
        gROOT->LoadMacro("loadlibs.C");
        loadlibs();
    } // end if

    TFile *hf = (TFile*)gROOT->GetListOfFiles()->FindObject(hfn);
    if(hf) {
        hf->Close();
        delete hf;
        hf = 0;
    } // end if file
    hf = new TFile(hfn,"READ");
    // Get AliRun object from file or return if not on file
    if (gAlice) {delete gAlice; gAlice = 0;}
    gAlice = (AliRun*)hf->Get("gAlice");
    if (!gAlice) {
        cerr << "AliRun object not found on file "<< FileName << "!" << endl;
        file->Close();  // close file and return error.
        return;
    } // end if !gAlice

/*
    gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSstandard.C");
    TFile *hf=0;
    hf = AccessFile(hfn,"R"); // Set up to read in Data
*/
    AliITS *ITS = (AliITS*)gAlice->GetDetector("ITS");
    if(!ITS){
        cout << "Error: no ITS found. Aborting"<<endl;
        return;
    } // end if !ITS
    AliITSgeom *gm = ITS->GetITSgeom();
    Int_t mod1 = 0;
    Int_t mod2 = gm->GetIndexMax();
    if(mod>=0){
        mod1 = mod;
        mod2 = mode+1;
    } // end if mod>=0
    AliITSgeomMatrix *gmm = gm->GetGeomMatrix(0);
    Int_t m;
    gmm->PrintComment(&cout); cout << endl;
    for(m=mod1;m<mod2;m++){
	gmm = gm->GetGeomMatrix(m);
	gmm->Print(&cout); cout << endl;
    } // end for m
}
