#if !defined(__CINT__) || defined(__MAKECINT__)

#include "iostream.h"
#include "TDatetime.h"
#include "STEER/AliRun.h"
#include "STEER/AliRunDigitizer.h"
#include "ITS/AliITSDigitizer.h"
#include "ITS/AliITS.h"
#include "ITS/AliITSDetType.h"
#include "ITS/AliITSresponseSDD.h"
#include "TStopwatch.h"

#endif

TFile* AccessFile(TString inFile="galice.root", TString acctype="R");
void writeAR(TFile * fin, TFile *fou);
Int_t ChangeITSDefaults(TFile *hitfile,AliITS *ITS,TString opt="");
#define DEBUG
Int_t AliITSDigits2RecPoints(TString digFile="galiceD.root", 
			    TString recFile="galiceR.root",TString opt=""){
    // Standard ITS Digits to RecPoints.

    // Dynamically link some shared libs
    if (gClassTable->GetID("AliRun") < 0) {
	gROOT->LoadMacro("loadlibs.C");
	loadlibs();
    } // end if

    // Connect the Root Galice file containing Geometry, Kine and Hits

    TFile *digfile = 0;   // pointer to input file.
    TFile *recfile = 0;  // possible output file for TreeD
    if(recFile.CompareTo(digFile) == 0) { //write output to same file as input.
	digfile = AccessFile(digFile,"U");  // input file open for update.
    }else{ // different output file then input file.
	digfile = AccessFile(digFile,"R"); // input file open as read only
	recfile = new TFile(recFile,"NEW");
    } // end if sdigFile == hitFile.
/*
    AliITS *ITS = (AliITS*)gAlice->GetDetector("ITS");      
    if (!ITS) {
	cerr<<"AliITSHits2DigitsDefault.C : AliITS object not found on file"
	    << endl;
	return 3;
    }  // end if !ITS
    if(!(ITS->GetITSgeom())){
	cerr << " AliITSgeom not found. Can't digitize with out it." << endl;
	return 4;
    } // end if

    ChangeITSDefaults(digfile,ITS,opt);
*/
    // write the AliRun object to the output file if different from input file.
    if(recfile) writeAR(digfile,recfile);
    if(recfile){ recfile->Close(); recfile = 0;}
//    if(digfile){ digfile->Close(); digfile = 0;}
//    delete gAlice;
//    gAlice = 0;

    TStopwatch timer;

#ifdef DEBUG
    cout << "Creating reconstructed points from digits for the ITS..." << endl;
#endif
//    const char *nulptr=0;
    AliITSreconstruction *itsr = new AliITSreconstruction(gAlice);
    if(digFile.CompareTo(recFile)!=0) {
#ifdef DEBUG
	cout << itsr << " filename="  << recFile << " compare="<<
	    digFile.CompareTo(recFile) << endl;
#endif
	itsr->SetOutputFile(recFile);
    } // end if
    timer.Start();
    itsr->Init();
    itsr->Exec(); 
    timer.Stop(); 
    timer.Print();
    delete itsr;
//    delete gAlice;
//    gAlice = 0;
    if(digfile){ digfile->Close(); digfile = 0;}
    if(recfile){ recfile->Close(); recfile = 0;}
    return 0;
}
//______________________________________________________________________
TFile * AccessFile(TString FileName, TString acctype){
    // Function used to open the input file and fetch the AliRun object

    TFile *retfil = 0;
    TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(FileName);
    if(file) {
	file->Close();
	delete file;
	file = 0;
    } // end if file
    if(acctype.Contains("U")){
	file = new TFile(FileName,"UPDATE");
    } // end if open for update
    if(acctype.Contains("N") && !file){
	file = new TFile(FileName,"RECREATE");
    } // end if open a new file
    if(!file) file = new TFile(FileName,"READ");   // default readonly
    if (!file->IsOpen()) {
	cerr << "Can't open " << FileName << " !" << endl;
	return retfil;
    } // end if error opeing file

    // Get AliRun object from file or return if not on file
    if (gAlice) {delete gAlice; gAlice = 0;}
    gAlice = (AliRun*)file->Get("gAlice");
    if (!gAlice) {
	cerr << "AliRun object not found on file "<< FileName << "!" << endl;
	file->Close();  // close file and return error.
	return retfil;
    } // end if !gAlice
    return file;
}
//______________________________________________________________________
void writeAR(TFile * fin, TFile *fou) {
    TDirectory *current = gDirectory;
/*
    TTree *TeOld;
    TTree *TeNew;
    AliHeader *alhe = new AliHeader();
    TeOld = (TTree*)fin->Get("TE");
    TeOld->SetBranchAddress("Header",&alhe);
    TeOld->SetBranchStatus("*",1);
    fou->cd();
    TeNew = TeOld->CloneTree();
    TeNew->Write(0,TObject::kOverwrite);
*/
    fou->cd();
    gAlice->TreeE()->SetBranchStatus("*",1);
    gAlice->TreeE()->Write(0,TObject::kOverwrite);
    gAlice->Write(0,TObject::kOverwrite);
    current->cd();
//    delete alhe;
#ifdef DEBUG
    cout << "AliRun object written to file" << endl;
#endif
}
//______________________________________________________________________
Int_t ChangeITSDefaults(TFile *hitfile,AliITS *ITS,TString opt){

    TDatime *ct0 = new TDatime(2002,04,26,00,00,00);
    TDatime ct = hitfile->GetCreationDate();

    if(ct0->GetDate()>ct.GetDate()){
	// For old files, must change SDD noise.
	AliITSresponseSDD *resp1 = (AliITSresponseSDD*)ITS->DetType(1)->
	    GetResponseModel();
	resp1 = new AliITSresponseSDD();
	ITS->SetResponseModel(1,resp1);
	cout << "Changed response class for SDD:" << endl;
	resp1->Print();
    } // end if

    if(opt.Contains("Dubna")){
	AliITSresponseSPDdubna *resp0 = new AliITSresponseSPDdubna();
	if(ITS->DetType(0)->GetResponseModel() !=0){
	    delete ((AliITSresponse*)ITS->DetType(0)->GetResponseModel());
	    ITS->DetType(0)->ResponseModel(0);
	} // end if
	ITS->DetType(0)->ResponseModel(resp0);
	AliITSsegmentationSPD *seg0 = (AliITSsegmentationSPD*)ITS->
	    DetType(0)->GetSegmentationModel();
	AliITSsimulationSPDdubna *sim0 = new AliITSsimulationSPDdubna(seg0,
								      resp0);
	if(ITS->DetType(0)->GetSimulationModel() !=0){
	    delete ((AliITSsimulation*)ITS->DetType(0)->GetSimulationModel());
	    ITS->DetType(0)->SimulationModel(0);
	} // end if
	ITS->DetType(0)->SimulationModel(sim0);
    } // end if Dubna
}
