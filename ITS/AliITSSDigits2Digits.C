void writeAR(TFile * fin, TFile *fou);
void AliITSSD2D(TString inFile, TString outFile);

void AliITSSDigits2Digits(TString inFile= "galice.root", TString outFile = ""){
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

void AliITSSD2D(TString inFile, TString outFile){
  AliRunDigitizer * manager = new AliRunDigitizer(1,1);
  char ftmp[50];
  sprintf(ftmp,"%s",inFile.Data());
  manager->SetInputStream(0,ftmp);
  if(outFile != "")manager->SetOutputFile(outFile);
  AliITSDigitizer *dITS  = new AliITSDigitizer(manager);
  manager->Exec("");
  TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(inFile);
  TFile *file2 = 0;
  if(outFile != ""){ 
    file2 = new TFile(outFile,"UPDATE");
    writeAR(file,file2);
  }
  delete manager;
  if(file){
    file->Write();
  }
  if(file2){
    file2->Close();
    delete file2;
  }
}

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






