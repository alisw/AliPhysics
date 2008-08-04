// from CreateESDChain.C - instead of  #include "CreateESDChain.C"
TChain* CreateESDChain(const char* aDataDir = "ESDfiles.txt", Int_t aRuns = 20, Int_t offset = 0) ;
void LookupWrite(TChain* chain, const char* target) ;

void runAliAnalysisTaskCumulants(Int_t nRuns =-1, TString type = "ESD", const Char_t* dataDir="/data/alice2/ante/ab2", Int_t offset = 0) 
//void runAliAnalysisTaskCumulants(Int_t nRuns = -1, TString type = "ESD", const Char_t* dataDir="/data/alice2/ab2", Int_t offset = 0) 
//void runAliAnalysisTaskCumulants(Int_t nRuns = -1, TString type = "ESD", const Char_t* dataDir="/data/alice2/LHyquid3_rot", Int_t offset = 0)
//void runAliAnalysisTaskCumulants(Int_t nRuns = 2, TString type = "MC", const Char_t* dataDir="/Users/snelling/alice_data/TherminatorFix", Int_t offset = 0) 
{
  TStopwatch timer;
  timer.Start();

  // include path (to find the .h files when compiling)
  gSystem->AddIncludePath("-I$ALICE_ROOT/include") ;
  gSystem->AddIncludePath("-I$ROOTSYS/include") ;

  // load needed libraries
  gSystem->Load("libTree.so");
  gSystem->Load("libESD.so");
  cerr<<"libESD loaded..."<<endl;
  gSystem->Load("libANALYSIS.so");
  cerr<<"libANALYSIS.so loaded..."<<endl;
  gSystem->Load("libPWG2flow.so");
  cerr<<"libPWG2flow.so loaded..."<<endl;
  
  // create the TChain. CreateESDChain() is defined in CreateESDChain.C
  TChain* chain = CreateESDChain(dataDir, nRuns, offset);
  cout<<"chain ("<<chain<<")"<<endl;

  //____________________________________________//
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");
   
  //AliVEventHandler* esdH = new AliESDInputHandler;
  //mgr->SetInputEventHandler(esdH);  
  //AliMCEventHandler *mc = new AliMCEventHandler();
  //mgr->SetMCtruthEventHandler(mc);
  
  if (type == "ESD"){
  AliVEventHandler* esdH = new AliESDInputHandler;
  mgr->SetInputEventHandler(esdH); }
   
  if (type == "AOD"){
  AliVEventHandler* aodH = new AliAODInputHandler;
  mgr->SetInputEventHandler(aodH); }
  
  if (type == "MC"){
  AliMCEventHandler *mc = new AliMCEventHandler();
  mgr->SetMCtruthEventHandler(mc); }
  
  //____________________________________________//
  // 1st cumulant task
  AliAnalysisTaskCumulants *task1 = new AliAnalysisTaskCumulants("TaskCumulants");
  task1->SetAnalysisType(type);
  mgr->AddTask(task1);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput1 = 
    mgr->CreateContainer("cchain1",TChain::Class(),AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *coutput1 = 
    mgr->CreateContainer("cobj1", TList::Class(),AliAnalysisManager::kOutputContainer,"outputFromCumulantAnalysisESD.root");

  //____________________________________________//
  mgr->ConnectInput(task1,0,cinput1);
  mgr->ConnectOutput(task1,0,coutput1);

  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("local",chain);

  timer.Stop();
  timer.Print();
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
	  cerr<<presentDirName<<endl;
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
1
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

