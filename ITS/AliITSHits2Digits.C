TFile* AccessFile(TString inFile="galice.root", TString acctype="R");
void writeAR(TFile * fin, TFile *fou);

Int_t AliITSHits2Digits(Int_t evNumber1=0,Int_t evNumber2=0, TString inFile = "galice.root", TString outFile="galice.root"){

  // Dynamically link some shared libs
  if (gClassTable->GetID("AliRun") < 0) {
	gROOT->LoadMacro("loadlibs.C");
	loadlibs();
  } // end if

  // Connect the Root Galice file containing Geometry, Kine and Hits

  TFile *file;
  if(outFile.Data() == inFile.Data()){
    file = AccessFile(inFile,"U");
  }
  else {
    file = AccessFile(inFile);
  }
  
  TFile *file2 = 0;  // possible output file for TreeD

  if(!(outFile.Data() == inFile.Data())){
    // open output file and create TreeR on it
    file2 = gAlice->InitTreeFile("D",outFile);
  }

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

  TDatime *ct0 = new TDatime(2002,04,26,00,00,00);
  TDatime ct = file->GetCreationDate();

  if(ct0->GetDate()>ct.GetDate()){
	// For old files, must change SDD noise.
    AliITSresponseSDD *resp1 = (AliITSresponseSDD*)ITS->DetType(1)->GetResponseModel();
    resp1 = new AliITSresponseSDD();
    ITS->SetResponseModel(1,resp1);
    cout << "Changed response class for SDD: \n";
    resp1->Print();
  } // end if
  TStopwatch timer;
  timer.Start();

  for(Int_t nevent = evNumber1; nevent <= evNumber2; nevent++){
    cout<<"Producing Digits for event n."<<nevent<<endl;
    gAlice->GetEvent(nevent);
    if(!gAlice->TreeD() && file2 == 0){ 
      cout << "Having to create the Digits Tree." << endl;
      gAlice->MakeTree("D");
    } 
    if(file2)gAlice->MakeTree("D",file2);
    ITS->MakeBranch("D");
    ITS->SetTreeAddress();   
    ITS->Hits2Digits();
  }
  timer.Stop();
  timer.Print();

  // write the AliRun object to the output file
  if(file2)writeAR(file,file2);

  delete gAlice;   
  gAlice=0;
  file->Close();
}

//-------------------------------------------------------------------
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

//-------------------------------------------------------------------
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





