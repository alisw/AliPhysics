////////////////////////////////////////////////////////////////////////////////
//macro to run AliSelectorLYZ for LYZ analysis on ESDs                        //
//run first with "Bool_t firstrun = kTRUE" by running .x makeSelectorLYZ()    //
//for the first run over the data (integrated flow)                           //
//run second with "Bool_t firstrun = kFALSE" by running .x makeSelectorLYZ(0) //
//for the second run over the data (differential flow)                        //
////////////////////////////////////////////////////////////////////////////////


// from CreateESDChain.C - instead of  #include "CreateESDChain.C"
TChain* CreateESDChain(const char* aDataDir = "ESDfiles.txt", Int_t aRuns = 20, Int_t offset = 0) ;
void LookupWrite(TChain* chain, const char* target) ;



void makeSelectorLYZ(Bool_t firstrun = kTRUE, Bool_t usesum = kTRUE, const Char_t* dataDir="/data/alice1/simili/lcgHijing/4", Int_t nRuns = -1, Int_t offset = 0)
//void makeSelectorLYZ(Bool_t firstrun = kTRUE, Bool_t usesum = kTRUE, const Char_t* dataDir="/data/alice/kolk/TestLYZData/data/cent5", Int_t nRuns = 1, Int_t offset = 0)
//void makeSelectorLYZ(Bool_t firstrun = kTRUE, Bool_t usesum = kTRUE, const Char_t* dataDir="/data/alice/simili/lcgGevSim3x/LDL/0", Int_t nRuns = -1, Int_t offset = 0)

//void makeSelectorLYZ(Bool_t firstrun = kTRUE, Bool_t usesum = kFALSE, const Char_t* dataDir="/data/alice/simili/lcgHijing/5",Int_t nRuns = -1, Int_t offset = 0)
//void makeSelectorLYZ(Bool_t firstrun = kTRUE, Bool_t usesum = kFALSE, const Char_t* dataDir="/data/alice/kolk/TestLYZData/data/cent5", Int_t nRuns =1, Int_t offset = 0)



{
 // include path (to find the .h files when compiling)
 gSystem->AddIncludePath("-I$ALICE_ROOT/include") ;
 gSystem->AddIncludePath("-I$ROOTSYS/include") ;
 gSystem->AddIncludePath("-I$ALICE_ROOT/PWG2/FLOW") ;

 // load needed libraries
 gSystem->Load("libESD");
  
 // Flow libraries 
 gSystem->Load("libPWG2flow.so");
 gROOT->LoadMacro("AliFlowLYZConstants.cxx+");
 gROOT->LoadMacro("AliFlowLYZHist1.cxx+");
 gROOT->LoadMacro("AliFlowLYZHist2.cxx+");
 gROOT->LoadMacro("AliFlowLeeYangZerosMaker.cxx+");
 gROOT->LoadMacro("AliSelectorLYZ.cxx+");
 
 // create the TChain. CreateESDChain() is defined in CreateESDChain.C
 TChain* chain = CreateESDChain(dataDir, nRuns, offset);
 cout << " * " << chain->GetEntries() << " * " << endl ;

 // enable debugging
 AliLog::SetClassDebugLevel("AliSelectorLYZ", AliLog::kInfo);

 // run selector on chain
 AliSelectorLYZ* selector = new AliSelectorLYZ; //create new selector object
 selector->SetFirstRunLYZ(firstrun);            //set to first or second run
 selector->SetUseSumLYZ(usesum);                //set to sum gen.function or product gen.function
 Long64_t result = chain->Process(selector);    //proces the selector on the chain of data
  
 if (result != 0)  {
   cout<<"ERROR: Executing process failed with "<<result<<endl;
   return; }

 cout<<"Execution complete."<<endl<<endl;
 delete selector;
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
