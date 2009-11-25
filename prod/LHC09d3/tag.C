void tag() {
  const char* turl = gSystem->Getenv("ALIEN_JDL_OUTPUTDIR");
  TString fESDFileName = "alien://";
  fESDFileName += turl;
  fESDFileName += "/AliESDs.root";

  TString fGUID = 0;
  GetGUID(fGUID);

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
    if(f->IsOpen()) {
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
  cout<<"Updating tags....."<<endl;

  const char * tagPattern = "tag.root";
  // Open the working directory
  void * dirp = gSystem->OpenDirectory(gSystem->pwd());
  const char * name = 0x0;
  // Add all files matching *pattern* to the chain
  while((name = gSystem->GetDirEntry(dirp))) {
    if (strstr(name,tagPattern)) {
      TFile *f = TFile::Open(name,"read") ;

      AliRunTag *tag = new AliRunTag;
      AliEventTag *evTag = new AliEventTag;
      TTree *fTree = (TTree *)f->Get("T");
      fTree->SetBranchAddress("AliTAG",&tag);

      //Defining new tag objects
      AliRunTag *newTag = new AliRunTag();
      TTree ttag("T","A Tree with event tags");
      TBranch * btag = ttag.Branch("AliTAG", &newTag);
      btag->SetCompressionLevel(9);
      for(Int_t iTagFiles = 0; iTagFiles < fTree->GetEntries(); iTagFiles++) {
	fTree->GetEntry(iTagFiles);
	newTag->SetRunId(tag->GetRunId());
	newTag->SetAlirootVersion(faliroot);
	newTag->SetRootVersion(froot);
	newTag->SetGeant3Version(fgeant);
	const TClonesArray *tagList = tag->GetEventTags();
	for(Int_t j = 0; j < tagList->GetEntries(); j++) {
	  evTag = (AliEventTag *) tagList->At(j);
	  evTag->SetTURL(turl);
	  evTag->SetGUID(guid);
	  newTag->AddEventTag(*evTag);
	}
	ttag.Fill();
	newTag->Clear();
      }//tag file loop

      TFile* ftag = TFile::Open(name, "recreate");
      ftag->cd();
      ttag.Write();
      ftag->Close();

      delete tag;
      delete newTag;
    }//pattern check
  }//directory loop
  return kTRUE;
}
