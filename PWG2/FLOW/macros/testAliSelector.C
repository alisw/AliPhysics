/* $Id$ */

//
// This is an example how to run an AliSelector ...
// this script runs the AliSelectorFlow on a chain of ESD files, some debug information is printed
//

// parameters are
//   dataDir: the directory containing subdirectories that contain the ESD files
//   nRuns: the number of files that should be processed
//   offset: the directory to start with

// from CreateESDChain.C - instead of  #include "CreateESDChain.C"
TChain* CreateESDChain(const char* aDataDir = "ESDfiles.txt", Int_t aRuns = 20, Int_t offset = 0) ;
void LookupWrite(TChain* chain, const char* target) ;
//

void testAliSelector(const Char_t* dataDir=".", Int_t nRuns = 10, Int_t offset = 0)
{
 // load needed libraries
  gSystem->Load("libESD");
  
 // PWG2 libraries 
  gSystem->Load("libPWG2.so");
  gSystem->Load("libPWG2flow.so");
  //char *loadFrom = gSystem->ExpandPathName("$ALICE_ROOT/PWG2/FLOW/AliFlow_Pure.so") ;
  //gSystem->Load(loadFrom);

 // Selector
  char *loadSelector = gSystem->ExpandPathName("$ALICE_ROOT/PWG2/FLOW/AliSelectorFlow.cxx+") ;
  //gInterpreter->AddIncludePath(inc) ;

 // create the TChain. CreateESDChain() is defined in CreateESDChain.C
  TChain* chain = CreateESDChain(dataDir, nRuns, offset);

  cout << " * " << chain->GetEntries() << " * " << endl ;

 // enable debugging
  AliLog::SetClassDebugLevel("AliSelectorFlow", AliLog::kInfo);

 // run selector on chain
  Long64_t result = chain->Process(loadSelector);

  if (result != 0)
  {
    printf("ERROR: Executing process failed with %d.\n", result);
    return;
  }

  printf("Execution complete.\n");
}

// Helper macros for creating chains
// from: CreateESDChain.C,v 1.10 jgrosseo Exp

TChain* CreateESDChain(const char* aDataDir, Int_t aRuns, Int_t offset)
{
  // creates chain of files in a given directory or file containing a list.
  // In case of directory the structure is expected as:
  // <aDataDir>/<dir0>/AliESDs.root
  // <aDataDir>/<dir1>/AliESDs.root
  // ...

  if (!aDataDir)
    return 0;

  Long_t id, size, flags, modtime;
  if (gSystem->GetPathInfo(aDataDir, &id, &size, &flags, &modtime))
  {
    printf("%s not found.\n", aDataDir);
    return 0;
  }

  TChain* chain = new TChain("esdTree");
  TChain* chaingAlice = 0;

  if (flags & 2)
  {
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
  }
  else
  {
    // Open the input stream
    ifstream in;
    in.open(aDataDir);

    Int_t count = 0;

    // Read the input list of files and add them to the chain
    TString esdfile;
    while(in.good()) {
      in >> esdfile;
      if (!esdfile.Contains("root")) continue; // protection

      if (offset > 0)
      {
        --offset;
        continue;
      }

      if (count++ == aRuns)
        break;

        // add esd file
      chain->Add(esdfile);
    }

    in.close();
  }

  return chain;
}

void LookupWrite(TChain* chain, const char* target)
{
  // looks up the chain and writes the remaining files to the text file target

  chain->Lookup();

  TObjArray* list = chain->GetListOfFiles();
  TIterator* iter = list->MakeIterator();
  TObject* obj = 0;

  ofstream outfile;
  outfile.open(target);

  while ((obj = iter->Next()))
    outfile << obj->GetTitle() << "#AliESDs.root" << endl;

  outfile.close();

  delete iter;
}
