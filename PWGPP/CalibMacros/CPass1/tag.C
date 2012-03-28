void tag() {
  const char* turl = gSystem->Getenv("ALIEN_JDL_OUTPUTDIR");

  gSystem->Load("libNet.so");
  gSystem->Load("libMonaLisa.so");
  new TMonaLisaWriter(0, "GridAliRoot-tag.C", 0, 0, "global");

  TString fESDFileName = "alien://";
  fESDFileName += turl;
  fESDFileName += "/AliESDs.root";  

  TString fGUID = 0;
  GetGUID(fGUID);

  gEnv->Print();

  TString fAliroot, fRoot, fGeant;
  GetVersions(fAliroot,fRoot,fGeant);

  UpdateTag(fAliroot,fRoot,fGeant,fESDFileName,fGUID);
}
  
//_____________________________________//
GetVersions(TString &fAliroot, TString &froot, TString &fgeant) {
  const char* fver = gSystem->Getenv("ALIEN_JDL_PACKAGES");
  TString fS = fver;
  Int_t fFirst = fS.First("#");

  while(fFirst != -1) {
    Int_t fTotalLength = fS.Length();
    TString tmp = fS;
    TString fS1 = fS(0,fFirst);
    tmp = fS(fFirst+2,fTotalLength);
    fS = tmp;

    if(fS1.Contains("Root")) fAliroot = fS1;
    if(fS1.Contains("ROOT")) froot = fS1;
    if(fS1.Contains("GEANT")) fgeant = fS1;

    if(tmp.Contains("Root")) fAliroot = tmp;
    if(tmp.Contains("ROOT")) froot = tmp;
    if(tmp.Contains("GEANT")) fgeant = tmp;
    
    fFirst = tmp.First("#");
  }
}

//_____________________________________//
GetGUID(TString &guid) {
  ofstream myfile ("guid.txt");
  if (myfile.is_open()) {
    TFile *f = TFile::Open("AliESDs.root","read");
    if(f && !f->IsZombie() && f->IsOpen()) {
      guid = f->GetUUID().AsString();
      myfile << "AliESDs.root \t"<<f->GetUUID().AsString();
      cout<<guid.Data()<<endl;
      myfile.close();
    }
    else cout<<"Input file not found"<<endl;
  }
  else cout<<"Output file can't be created..."<<endl;
}


//_____________________________________//
Bool_t UpdateTag(TString faliroot, TString froot, TString fgeant, TString turl, TString guid) {
  cout<<"> Updating tags...."<<endl;

  const char * tagPattern = "tag.root";
  // Open the working directory
  void * dirp = gSystem->OpenDirectory(gSystem->pwd());
  const char * name = 0x0;
  // Add all files matching *pattern* to the chain
  while((name = gSystem->GetDirEntry(dirp))) {
    cout<<">>> Adding to chain file " << name << "...." << endl;
    if (strstr(name,tagPattern)) {
      TFile *f = TFile::Open(name,"read") ;
 
      AliRunTag *tag = 0x0;
      AliEventTag *evTag = 0x0;
      TTree *fTree = (TTree *)f->Get("T");
      if (!fTree) { f->Close(); continue; }
      fTree->SetBranchAddress("AliTAG",&tag);
   
      //Defining new tag objects
      AliRunTag *newTag = 0x0;
      TTree ttag("T","A Tree with event tags");
      TBranch * btag = ttag.Branch("AliTAG", &newTag);
      btag->SetCompressionLevel(9);
      
      cout<<">>>>> Found " << fTree->GetEntries() << " entries...." << endl;
      for(Int_t iTagFiles = 0; iTagFiles < fTree->GetEntries(); iTagFiles++) {
	fTree->GetEntry(iTagFiles);
	newTag = new AliRunTag(*tag);
	newTag->SetAlirootVersion(faliroot);
	newTag->SetRootVersion(froot);
	newTag->SetGeant3Version(fgeant);
 	TClonesArray *tagList = newTag->GetEventTags();
 	cout << "Found " << tagList->GetEntries() << " tags" << endl;
 	for(Int_t j = 0; j < tagList->GetEntries(); j++) {
 	  evTag = (AliEventTag *) tagList->At(j);
 	  evTag->SetTURL(turl);
 	  evTag->SetGUID(guid);
 	}
	ttag.Fill();

	delete tag;
	delete newTag;
      }//tag file loop 
      
      TFile* ftag = TFile::Open(name, "recreate");
      ftag->cd();
      ttag.Write();
      ftag->Close();

    }//pattern check
  }//directory loop
  return kTRUE;
}

