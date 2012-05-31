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

Int_t AliITSh2sd(TString hf="galice.root",TString sf="",TString opt="");
Int_t AliITSh2d(TString hf="galice.root",TString df="",TString opt="");
Int_t AliITSsd2d(TString df="galice.root",TString sf1="galice.root",
                 TString sf2="",TString opt="")
Int_t AliITSd2r(TString df="galice.root",TString rf="",TString opt="");
void grun(Int_t nevent=1, const char *config="Config.C");

#endif
//#define DEBUG


Int_t AliITSh2r(TString hf="galice.root",TString opt=""){
    // Runs grun to create the file galice.root and then runs the 
    // nessesary files to creat recpoints form that file.

    gROOT->LoadMacro("$(ALICE_ROOT)/macros/grun.C");
    if(!opt.Contains("Nogrun")) grun(1,"Config.C");
    if(opt.Contains("SDigit")){
	gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSh2sd.C");
	if(AliITSh2sd(hf,hf,opt)!=0){
	    cout << "Error AliITSh2sd failed" << endl;
	    return 1;
	} // end if
	gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSsd2d.C");
	if(AliITSsd2d(hf,hf,"",opt)!=0){
	    cout << "Error AliITSsd2d failed" << endl;
	    return 2;
	} // end if
    }else{
	gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSh2d.C");
	if(AliITSh2d(hf,hf,opt)!=0){
	    cout << "Error AliITSh2d failed" << endl;
	    return 3;
	} // end if
    } // end if
    gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSd2r.C");
    if(AliITSd2r(hf,hf,opt)!=0){
	cout << "Error AliITSd2r failed" << endl;
	return 4;
    } // end if
    return 0;
}
