TFile* AccessFile(TString inFile="galice.root", TString acctype="R");
void writeAR(TFile * fin, TFile *fou);
void AliITSSD2D(TString inFile, TString outFile);

void AliITSSDigits2DigitsDubna(TString inFile= "galiceS.root",
			  TString outFile = "galiceD.root"){
    // This macro takes SDigits and produces Digits. No merging is done
    // and only one galice.root file is used. 
    // Dynamically link some shared libs 
    TStopwatch timer;

    if(gAlice){
	delete gAlice;
	gAlice = 0;
    } // end if gAlice
    cout << "Creating digits from summable digits for the ITS..." << endl;
    AliITSSD2D(inFile,outFile);
    timer.Stop(); 
    timer.Print();
}
//______________________________________________________________________
void AliITSSD2D(TString inFile, TString outFile){
    AliRunDigitizer * manager = new AliRunDigitizer(1,1);
    char ftmp[50];
    sprintf(ftmp,"%s",inFile.Data());
    TFile *file0 = AccessFile(ftmp);
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

    // For old files, must change SPD noise.
    AliITSresponseSPDdubna *resp0 = new AliITSresponseSPDdubna();
    if(ITS->DetType(0)->GetResponseModel() !=0){
	delete ((AliITSresponse*)ITS->DetType(0)->GetResponseModel());
	ITS->DetType(0)->ResponseModel(0);
    } // end if
    ITS->DetType(0)->ResponseModel(resp0);
    AliITSsegmentationSPD *seg0 = (AliITSsegmentationSPD*)ITS->DetType(0)->
	GetSegmentationModel();
    AliITSsimulationSPDdubna *sim0 = new AliITSsimulationSPDdubna(seg0,resp0);
    if(ITS->DetType(0)->GetSimulationModel() !=0){
	delete ((AliITSsimulation*)ITS->DetType(0)->GetSimulationModel());
	ITS->DetType(0)->SimulationModel(0);
    } // end if
    ITS->DetType(0)->SimulationModel(sim0);
    manager->SetInputStream(0,ftmp);
    if(outFile != "")manager->SetOutputFile(outFile);
    AliITSDigitizer *dITS  = new AliITSDigitizer(manager);
    manager->Exec("");
    TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(inFile);
    TFile *file2 = 0;
    if(outFile != ""){ 
	file2 = new TFile(outFile,"UPDATE");
	writeAR(file,file2);
    } // end if outFile!=""
    delete manager;
    if(file){
	file->Write();
    } // end if file
    if(file2){
	file2->Close();
	delete file2;
    } // end if file2
}
//______________________________________________________________________
TFile * AccessFile(TString FileName, TString acctype){

  // Function used to open the input file and fetch the AliRun object

  TFile *retfil = 0;
  TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(FileName);
  if (file) {file->Close(); delete file; file = 0;} 
  if(acctype.Contains("U")){
    file = new TFile(FileName,"update");
  }
  if(acctype.Contains("N") && !file){
    file = new TFile(FileName,"recreate");
  }
  if(!file) file = new TFile(FileName);   // default readonly
  if (!file->IsOpen()) {
        cerr<<"Can't open "<<FileName<<" !" << endl;
        return retfil;
  } 

  // Get AliRun object from file or return if not on file
  if (gAlice) {delete gAlice; gAlice = 0;}
  gAlice = (AliRun*)file->Get("gAlice");
  if (!gAlice) {
        cerr << "AliRun object not found on file"<< endl;
        return retfil;
  } 
  return file;
}
//______________________________________________________________________
void writeAR(TFile * fin, TFile *fou) {
    TDirectory *current = gDirectory;
    TTree *Te;
    TTree *TeNew;
    AliHeader *alhe = new AliHeader();
    Te = (TTree*)fin->Get("TE");
    Te->SetBranchAddress("Header",&alhe);
    Te->SetBranchStatus("*",1);
    fou->cd();
    TeNew = Te->CloneTree();
    TeNew->Write(0,TObject::kOverwrite);
    gAlice->Write(0,TObject::kOverwrite);
    current->cd();
    delete alhe;
    cout<<"AliRun object written to file\n";
}
