//=========================================================================//
//                                                                         //
//                   A template runPmdTask for PMD analysis                //
//            You can copy it and add features according to your need      //
//                                                                         //
//                               Satyajit Jena                             //
//                               sjena@cern.ch                             //
//                                13/04/2012                               //
//                                                                         //
//=========================================================================//

const char * incollection = "file.txt";
//______________________________________________________________________________
void runPmdTask(Bool_t         isGrid = 0,
		Bool_t           isMC = 0,
		const char *gridmode  = "test") {

 
  if (isGrid) {
    Printf("Strating the Grid Job ");
    Printf("Grid Mode %s",gridmode);
    if (isMC)  Printf("It is Data Type Run");
    if (!isMC) Printf("It is MC Type Run");
  }
    
  // Load the needed libraries most of 
  // them already loaded by aliroot
  //--------------------------------------
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libXMLIO");
  gSystem->Load("libPhysics");

  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");

  // Use AliRoot includes to compile our task                                   
  gROOT->ProcessLine(".include $ALICE_INSTALL/include");
  
  AliAnalysisManager *mgr = new AliAnalysisManager("PMDAnalysis");
  
  if(isGrid) {
    gROOT->LoadMacro("CreateAlienHandler.C");
    AliAnalysisGrid *alienHandler = CreateAlienHandler(gridmode);  
    if (!alienHandler) return;
    mgr->SetGridHandler(alienHandler);
  }
    
  
 
  //  TChain *  chain;
  //  if(!isGrid) {
  //  gROOT->LoadMacro("CreateESDChain.C"); // use it if you know it
  //  chain = CreateESDChain("file.txt", 10);
  //  }
 
  TChain *  chain = new TChain("esdTree"); 
  if(!isGrid) {
    ifstream file_collect(incollection);
    TString line;
    while (line.ReadLine(file_collect) ) {
      chain->Add(line.Data());
    }
  }
    
  if(isMC) {
    AliVEventHandler* esdH = new AliESDInputHandler();
    mgr->SetInputEventHandler(esdH);
    AliMCEventHandler *mc = new AliMCEventHandler();
    //mc->SetReadTR(kTRUE);
    mgr->SetMCtruthEventHandler(mc);
    gROOT->LoadMacro("$ALICE_ROOT/OADB/macros/AddTaskPhysicsSelection.C");
    AliPhysicsSelectionTask* physicsSelTask = AddTaskPhysicsSelection(isMC);
    physicsSelTask->GetPhysicsSelection()->SetAnalyzeMC();

    gROOT->LoadMacro("$ALICE_ROOT/OADB/macros/AddTaskCentrality.C");
    AliCentralitySelectionTask *centralityTask = AddTaskCentrality();
  } else {
    AliVEventHandler* esdH = new AliESDInputHandler();
    mgr->SetInputEventHandler(esdH);
    gROOT->LoadMacro("$ALICE_ROOT/OADB/macros/AddTaskPhysicsSelection.C");
    AliPhysicsSelectionTask* physicsSelTask = AddTaskPhysicsSelection(isMC);
    
    gROOT->LoadMacro("$ALICE_ROOT/OADB/macros/AddTaskCentrality.C");
    AliCentralitySelectionTask *centralityTask = AddTaskCentrality();
  }

  gROOT->LoadMacro("AliPMDAnalysisTaskPbPb.cxx++g");   
  gROOT->LoadMacro("AddAliPMDAnalysisTaskPbPb.C");
  AddAliPMDAnalysisTaskPbPb("MyTask",isMC);
  
  // Enable debug printouts
  mgr->SetDebugLevel(0);
  if (!mgr->InitAnalysis())
    return;
  
 mgr->PrintStatus();
 
 if(isGrid)
   mgr->StartAnalysis("grid");
 else  
   mgr->StartAnalysis("local",chain);
 
};


//__________________________________________________________________________________________
// This helper macros creates a chain of ESD files for you. Source can be either a text
// file with the file paths or a directory. In the latter case all ESD files in all subdirectories
// are considered.
//
// Inspired by: Jan.Fiete.Grosse-Oetringhaus@cern.ch

TChain* CreateESDChain(const char* aDataDir = "ESDfiles.txt", Int_t aRuns = 20, Int_t offset = 0, Bool_t addFileName = kFALSE, Bool_t addFriend = kFALSE, const char* check = 0)
{
  // creates chain of files in a given directory or file containing a list.
  // In case of directory the structure is expected as:
  // <aDataDir>/<dir0>/AliESDs.root
  // <aDataDir>/<dir1>/AliESDs.root
  // ...
  //
  // if addFileName is true the list only needs to contain the directories that contain the AliESDs.root files
  // if addFriend is true a file AliESDfriends.root is expected in the same directory and added to the chain as friend
  // if check is != 0 the files that work are written back into the textfile with the name check

  if (!aDataDir)
    return 0;

  Long_t id, size, flags, modtime;
  if (gSystem->GetPathInfo(aDataDir, &id, &size, &flags, &modtime))
  {
    printf("%s not found.\n", aDataDir);
    return 0;
  }

  TChain* chain = new TChain("esdTree");
  TChain* chainFriend = 0;
  
  if (addFriend)
    chainFriend = new TChain("esdFriendTree");

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

    ofstream outfile;
    if (check)
      outfile.open(check);

    Int_t count = 0;

    // Read the input list of files and add them to the chain
    TString line;
    while (in.good())
    {
      in >> line;

      if (line.Length() == 0)
        continue;

      if (offset > 0)
      {
        offset--;
        continue;
      }

      if (count++ == aRuns)
        break;

      TString esdFile(line);

      if (addFileName)
        esdFile += "/AliESDs.root";
        
      TString esdFileFriend(esdFile);
      esdFileFriend.ReplaceAll("AliESDs.root", "AliESDfriends.root");
        
      if (check)
      {
        TFile* file = TFile::Open(esdFile);
        if (!file)
          continue;
        file->Close();
        
        if (chainFriend)
        {
          TFile* file = TFile::Open(esdFileFriend);
          if (!file)
            continue;
          file->Close();
        }
        
        outfile << line.Data() << endl;
        printf("%s\n", line.Data());
      }        
        
        // add esd file
      chain->Add(esdFile);

        // add file
      if (chainFriend)
        chainFriend->Add(esdFileFriend);
    }

    in.close();
    
    if (check)
      outfile.close();
  }
  
  if (chainFriend)
    chain->AddFriend(chainFriend);

  return chain;
}

void ChainToTextFile(TChain* chain, const char* target)
{
  // write a text list of the files in the chain
  
  TObjArray* list = chain->GetListOfFiles();
  TIterator* iter = list->MakeIterator();
  TObject* obj = 0;

  ofstream outfile;
  outfile.open(target);

  while ((obj = iter->Next())) {
    TString fileName(obj->GetTitle());
    
    outfile << fileName.Data() << endl;
  }

  outfile.close();

  delete iter;
} 

TObjArray* Chain2List(TChain* chain)
{
  // returns a TObjArray of TObjStrings of the file names in the chain

  TObjArray* result = new TObjArray;

  for (Int_t i=0; i<chain->GetListOfFiles()->GetEntries(); i++)
    result->Add(new TObjString(chain->GetListOfFiles()->At(i)->GetTitle()));

  return result;
}

void LookupWrite(TChain* chain, const char* target)
{
  // looks up the chain and writes the remaining files to the text file target

  chain->Lookup();

  ChainToTextFile(chain, target);
}

TChain* CreateChain(const char* treeName, const char* aDataDir, Int_t aRuns, Int_t offset = 0)
{
  // creates chain of files in a given directory or file containing a list.

  if (!treeName || !aDataDir)
    return 0;

  TChain* chain = new TChain(treeName);
  
  // Open the input stream
  ifstream in;
  in.open(aDataDir);

  Int_t count = 0;

  // Read the input list of files and add them to the chain
  TString line;
  while(in.good()) 
  {
    in >> line;
      
    if (line.Length() == 0)
      continue;      
    
    if (offset > 0)
    {
      --offset;
      continue;
    }

    if (count++ == aRuns)
      break;

    chain->Add(line);
  }

  in.close();

  return chain;
}
