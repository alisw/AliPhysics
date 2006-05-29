/* $Id$ */

// Helper macros for creating chains

TChain* CreateESDChainFromDir(const char* aDataDir, Int_t aRuns = 20, Int_t offset = 0)
{
  // creates chain of files in a given directory. The structure is expected as:
  // <aDataDir>/<dir0>/AliESDs.root
  // <aDataDir>/<dir1>/AliESDs.root
  // ...

  if (!aDataDir)
    return 0;

  TChain* chain = new TChain("esdTree");
  TChain* chaingAlice = 0;

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

    if (offset > 0)
    {
      --offset;
      continue;
    }

    if (count++ == aRuns)
      break;

    TString presentDirName(aDataDir);
    presentDirName += "/";
    presentDirName += presentDir->GetName();

    chain->Add(presentDirName + "/AliESDs.root/esdTree");
  }

  return chain;
}

TChain* CreateESDChainFromList(const char* listFile, Int_t aRuns = 20)
{
  // Creates a chain from a file which contains a list of ESD files

  if (!listFile)
    return 0;

  TChain* chain = new TChain("esdTree");
  TChain* chaingAlice = 0;

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
  }

  in.close();

  chain->Lookup();

  return chain;
}
