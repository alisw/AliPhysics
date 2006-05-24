TChain* CreateESDChainFromDir(const char* aDataDir, Int_t aRuns = 20, Bool_t aAddHeader = kTRUE)
{
  if (!aDataDir)
    return 0;

  TChain* chain = new TChain("esdTree");
  TChain* chaingAlice = 0;

  if (aAddHeader != kFALSE)
    chaingAlice = new TChain("TE");

  TString execDir(gSystem->pwd());
  TSystemDirectory* baseDir = new TSystemDirectory(".", aDataDir);
  TList* dirList            = baseDir->GetListOfFiles();
  Int_t nDirs               = dirList->GetEntries();
  gSystem->cd(execDir);

  Int_t count = 0;

  for (Int_t iDir=0; iDir<nDirs; ++iDir)
  {
    TSystemFile* presentDir = (TSystemFile*) dirList->At(iDir);
    if (!presentDir || !presentDir->IsDirectory() || strcmp(presentDir->GetName(), ".") == 0 || strcmp(presentDir->GetName(), "..") == 0)
      continue;

    if (count++ == aRuns)
      break;

    TString presentDirName(aDataDir);
    presentDirName += "/";
    presentDirName += presentDir->GetName();

    chain->Add(presentDirName + "/AliESDs.root/esdTree");

    if (aAddHeader != kFALSE)
      chaingAlice->Add(presentDirName + "/galice.root/TE");
  }

  if (aAddHeader != kFALSE)
    chain->AddFriend(chaingAlice);

  return chain;
}

TChain* CreateESDChainFromList(const char* listFile, Int_t aRuns = 20, Bool_t aAddHeader = kTRUE)
{
  if (!listFile)
    return 0;

  TChain* chain = new TChain("esdTree");
  TChain* chaingAlice = 0;

  if (aAddHeader != kFALSE)
    chaingAlice = new TChain("TE");

  // Open the input stream
  ifstream in;
  in.open(listFile);

  Int_t count = 0;

  // Read the input list of files and add them to the chain
  TString esdfile;
  while(in.good()) {
    in >> esdfile;
    if (!esdfile.Contains("root")) continue; // protection

    if (count++ == aRuns)
      break;

      // add esd file
    chain->Add(esdfile);

      // add header
    if (aAddHeader != kFALSE)
    {
      esdfile.ReplaceAll("AliESDs", "galice");
      chaingAlice->Add(esdfile + "/TE");
    }
  }

  in.close();

  chain->Lookup();

  if (aAddHeader != kFALSE)
  {
    chaingAlice->Lookup();
    chain->AddFriend(chaingAlice);
  }

  return chain;
}
