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

//#define DEBUG

//______________________________________________________________________
Bool_t GaliceITSok(){
    // Checks gAlice to see that ITS and the ITS geometry are properly
    // defined. If not return kFALSE and deletes gAlice and set it to zero.

    if(!gAlice){
	return kFALSE;
    } // end if !gAlice
    // gAlice defined check to see if ITS is properly defined.
    AliITS *ITS = (AliITS*)gAlice->GetDetector("ITS"); 
    if(!ITS){ // ITS not defined, delete and reload gAlice
	delete gAlice;
	gAlice = 0;
	return kFALSE;
    } // end if !ITS
    // ITS defined
    if(!(ITS->GetITSgeom())){
	delete gAlice;
	gAlice = 0;
	return kFALSE;
    } // end if !(ITS->GetITSgeom())
    // ITS and ITS geometry properly defined defined.
    return kTRUE;
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
    fou->cd();
    gAlice->TreeE()->SetBranchStatus("*",1);
    gAlice->TreeE()->Write(0,TObject::kOverwrite);
    gAlice->Write(0,TObject::kOverwrite);
    current->cd();
#ifdef DEBUG
    cout << "AliRun object from file "<<fin->GetName() 
	 << " written to file " << fou->GetName() <<"." << endl;
#endif
}
//______________________________________________________________________
Int_t ChangeITSDefaults(TFile *hitfile,AliITS *ITS,TString opt){

    TDatime *ct0 = new TDatime(2002,04,26,00,00,00);
    TDatime ct = *ct0;
    if(hitfile) ct = hitfile->GetCreationDate();

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
