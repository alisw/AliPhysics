TChain* CreateESDChain(const char* aDataDir, Int_t aRuns = 20, Bool_t aAddHeader = kTRUE)
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
