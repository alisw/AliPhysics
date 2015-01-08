/**
 * @file   Tag.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Wed Oct 15 13:28:47 2014
 * 
 * @brief  Update tags 
 */
//____________________________________________________________________
/** 
 * Tag files 
 * 
 */
void Tag() 
{
  const char* turl = gSystem->Getenv("ALIEN_JDL_OUTPUTDIR");

  gSystem->Load("libNet");
  //  gSystem->Load("libMonaLisa");
  //  new TMonaLisaWriter(0, "GridAliRoot-tag.C", 0, 0, "global");
  
  TString fESDFileName = "alien://";
  fESDFileName += turl;
  fESDFileName += "/AliESDs.root";  

  TString fGUID = 0;
  GetGUID(fGUID);

  // gEnv->Print();

  TString fAliroot, fRoot, fGeant;
  GetVersions(fAliroot,fRoot,fGeant);

  TString fPeriod, fPass, fName;
  GetProductionInfo(fPeriod, fPass, fName);

  UpdateTag(fAliroot,fRoot,fGeant,fESDFileName,fGUID,fPeriod,fPass,fName);
}

//____________________________________________________________________
/** 
 * Extract production information from path 
 * 
 * @param fPeriod  On return, the period
 * @param fPass    On return, the possible pass number 
 * @param fName    On return, the full name 
 */
void GetProductionInfo(TString &fPeriod, TString &fPass, TString &fName) 
{
  const char* turl = gSystem->Getenv("ALIEN_JDL_OUTPUTDIR");
  
  TString fS = turl;
  TObjArray *fDirs = fS.Tokenize("/");
  
  for (int iter=0; iter<fDirs->GetEntries(); iter++) {
    TString fDir = ((TObjString *) fDirs->At(iter))->String();

    if (fDir.Contains("LHC")) fPeriod = fDir;
    if (fDir.Contains("pass")) fPass = fDir;
  }
  fName = fPeriod+"."+fPass;
}
  
//____________________________________________________________________
/** 
 * Get the AliROOT, ROOT, and GEANT3 versions from the JDL packages
 * 
 * @param fAliroot On return, the AliROOT version 
 * @param froot    On return, the ROOT version 
 * @param fgeant   On return, the GEANT3 version 
 */
void GetVersions(TString &fAliroot, TString &froot, TString &fgeant) 
{
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

//____________________________________________________________________
/** 
 * Get the global univeral identifier of the ESD file 
 * 
 * @param guid 
 */
void GetGUID(TString &guid) 
{
  ofstream myfile ("guid.txt");
  if (!myfile.is_open()) {
    Warning("GetGUID", "Couldn't open guid.txt for writing");
    return;
  }

  TFile *f = TFile::Open("AliESDs.root","read");
  if (!f || f->IsZombie() || !f->IsOpen()) {
    Warning("GetGUID", "Input file AliESDs.root not found");
    return;
  }

  guid = f->GetUUID().AsString();
  f->Close();
  Info("", "Global Unique IDentifier: %s", guid.Data());

  myfile << "AliESDs.root \t"<< guid << std::endl;
  myfile.close();
}


//____________________________________________________________________
/** 
 * Update the tags 
 * 
 * @param faliroot  AliROOT version 
 * @param froot     ROOT version
 * @param fgeant    GEANT3 version
 * @param turl      UrL we're at 
 * @param guid      Global universal identifier 
 * @param fperiod   Period
 * @param fpass     Pass
 * @param fname     Full name 
 * 
 * @return 
 */
Bool_t UpdateTag(TString faliroot, 
		 TString froot, 
		 TString fgeant, 
		 TString turl, 
		 TString guid,
		 TString fperiod, 
		 TString fpass, 
		 TString fname) 
{
  Info("", "Updating tags (%s,%s,%s,%s,%s,%s,%s,%s",
       faliroot.Data(), froot.Data(), fgeant.Data(), 
       turl.Data(), guid.Data(), fperiod.Data(), 
       fpass.Data(),fname.Data());

  const TString tagPattern = "tag.root";

  // --- Open the working directory ----------------------------------
  TSystemDirectory dir(".", gSystem->pwd());
  TIter            next(dir.GetListOfFiles());
  TSystemFile*     file = 0;

  // --- Add all files matching *pattern* to the chain ---------------
  while ((file = static_cast<TSystemFile*>(next()))) {
    TString name(file->GetName());
    if (!name.Contains(tagPattern)) continue;
    
    // --- Open file matching pattern --------------------------------
    TFile*      f     = TFile::Open(name,"read") ;
    if (!f) { 
      continue;
    }
    Info("", "Updating tags in %s", name.Data());

    // --- Find the tree ---------------------------------------------
    AliRunTag*  tag   = 0x0;
    AliFileTag* flTag = 0x0;
    TTree*      fTree = (TTree *)f->Get("T");
    if (!fTree) { 
      f->Close(); 
      continue; 
    }
    fTree->SetBranchAddress("AliTAG",&tag);
   
    // --- Defining new tag objects ----------------------------------
    AliRunTag* newTag = 0x0;
    TTree      ttag("T","A Tree with event tags");
    TBranch*   btag = ttag.Branch("AliTAG", &newTag);
    btag->SetCompressionLevel(9);
    // --- disassociate the tree with underlying directory -----------
    ttag.SetDirectory(0);
      
    Printf(">>>>> Found %d entries....",fTree->GetEntries());

    for (Int_t iTagFiles = 0; iTagFiles < fTree->GetEntries(); iTagFiles++) {
      fTree->GetEntry(0);
      newTag = new AliRunTag(*tag);
      newTag->SetAlirootVersion(faliroot);
      newTag->SetRootVersion(froot);
      newTag->SetGeant3Version(fgeant);
      newTag->SetLHCPeriod(fperiod);
      newTag->SetReconstructionPass(fpass);
      newTag->SetProductionName(fname);
      Printf("Found %d file tags",newTag->GetNFiles());
      for(Int_t j = 0; j < newTag->GetNFiles(); j++) {
	flTag = (AliFileTag *) newTag->GetFileTag(j);
	flTag->SetTURL(turl);
	flTag->SetGUID(guid);
      }
      ttag.Fill();

      delete tag;
      delete newTag;
    }//tag file loop 
    
    // --- Close the input file --------------------------------------
    f->Close();

    // --- Overwrite the file ----------------------------------------
    TFile* ftag = TFile::Open(name, "recreate");
    ftag->cd();
    ttag.Write();
    ftag->Close();
    Info("", "Overwrote %s with new tags", name.Data());
  }//directory loop
  return kTRUE;
}

// 
// EOF
// 
