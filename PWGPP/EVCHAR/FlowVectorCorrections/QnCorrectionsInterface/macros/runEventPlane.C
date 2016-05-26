//TChain* CreateESDChain(const char* aDataDir, Int_t aRuns, Int_t offset);
//TChain* CreateAODChain(const char* aDataDir, Int_t aRuns, Int_t offset);

void runEventPlane(Bool_t isESD, Bool_t onGrid)
{
  // Load common libraries
  gSystem->Load("libCore.so");  
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  // gSystem->Load("libNetx.so");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  // Use AliRoot includes to compile our task
  gROOT->ProcessLine(".include $ROOTSYS/include");
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  gROOT->ProcessLine(".include $ALICE_PHYSICS/include");

  gSystem->Load("libPWGPPevcharQn.so");
  gSystem->Load("libPWGPPevcharQnInterface.so");

  //Bool_t onGrid=kFALSE;
  //Bool_t onGrid=kTRUE;

  // Select input format
  //Bool_t isESD=kTRUE;
  Bool_t isAOD=!isESD;

  //Create and configure the alien handler plugin




  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    mgr = new AliAnalysisManager("EventPlaneCorrections");

    //Error("AliAnalysisManager", "No analysis manager found.");
    //return 0;
  }


  if(isESD) {
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddESDHandler.C");
    AliVEventHandler* esdH = AddESDHandler();
  }
  if(isAOD) {
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddAODHandler.C");
    AliVEventHandler* aodH = AddAODHandler();
  }

  TChain* chain=0X0;
  if(onGrid){
    gROOT->LoadMacro("CreateAlienHandlerEventPlane.C");
    AliAnalysisGrid *alienHandler = CreateAlienHandler();
    if (!alienHandler) return;
    mgr->SetGridHandler(alienHandler);
  }
  else {

    if(isESD) {
      gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/CreateESDChain.C");
      //chain = CreateESDChain("/u/jonderw/file137231.txt", 1, 0);
      //chain = CreateESDChain("/u/jonderw/files.txt", 1, 0);
      chain = CreateESDChain("/u/jonderw/files.txt", 1, 0);
    }
    if(isAOD) {
      gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/CreateAODChain.C");
      //chain = CreateAODChain("/u/jonderw/aodfiles.txt", 1, 0);
      chain = CreateAODChain("/mnt/home/jonderw/alicesoftware/alice/alien/data/2010/LHC10h/000138275/ESDs/pass2/AOD160/files.txt", 1, 0);
    }
  }

  if(isESD){

    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
    AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection();

    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
    AliMultSelectionTask * task = AddTaskMultSelection(kFALSE); // user mode:

    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");
    AliCentralitySelectionTask *taskCentrality = AddTaskCentrality();

  }

  // uncomment if par file is used (may have to comment out includes at beginning of file)
  //if (!AliAnalysisAlien::SetupPar("PWGPPevcharQnInterface.par")) return;

  //gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/EVCHAR/EventPlaneFramework/macros/AddTask_ep.C");
  gROOT->LoadMacro("AddTask_ep.C");
  AddTask_ep();

  // Enable debug printouts
  mgr->SetDebugLevel(5);

  if (!mgr->InitAnalysis()) return 0x0;
  mgr->PrintStatus();

  if(onGrid) mgr->StartAnalysis("grid");
  else mgr->StartAnalysis("local",chain);



}




// Helper macros for creating chains
// from: CreateESDChain.C,v 1.10 jgrosseo Exp

//TChain* CreateESDChain(const char* aDataDir, Int_t aRuns, Int_t offset)
//{
//  // creates chain of files in a given directory or file containing a list.
//  // In case of directory the structure is expected as:
//  // <aDataDir>/<dir0>/AliESDs.root
//  // <aDataDir>/<dir1>/AliESDs.root
//  // ...
//
//  if (!aDataDir)
//    return 0;
//
//  Long_t id, size, flags, modtime;
//  if (gSystem->GetPathInfo(aDataDir, &id, &size, &flags, &modtime))
//  {
//    printf("%s not found.\n", aDataDir);
//    return 0;
//  }
//
//  TChain* chain = new TChain("esdTree");
//  TChain* chaingAlice = 0;
//
//  if (flags & 2)
//  {
//    TString execDir(gSystem->pwd());
//    TSystemDirectory* baseDir = new TSystemDirectory(".", aDataDir);
//    TList* dirList            = baseDir->GetListOfFiles();
//    Int_t nDirs               = dirList->GetEntries();
//    gSystem->cd(execDir);
//
//    Int_t count = 0;
//
//    for (Int_t iDir=0; iDir<nDirs; ++iDir)
//    {
//      TSystemFile* presentDir = (TSystemFile*) dirList->At(iDir);
//      if (!presentDir || !presentDir->IsDirectory() || strcmp(presentDir->GetName(), ".") == 0 || strcmp(presentDir->GetName(), "..") == 0)
//        continue;
//
//      if (offset > 0)
//      {
//        --offset;
//        continue;
//      }
//
//      if (count++ == aRuns)
//        break;
//
//      TString presentDirName(aDataDir);
//      presentDirName += "/";
//      presentDirName += presentDir->GetName();	  
//      chain->Add(presentDirName + "/AliESDs.root/esdTree");
//      //  cerr<<presentDirName<<endl;
//    }
//
//  }
//  else
//  {
//    // Open the input stream
//    ifstream in;
//    in.open(aDataDir);
//
//    Int_t count = 0;
//
//    // Read the input list of files and add them to the chain
//    TString esdfile;
//    while(in.good()) {
//      in >> esdfile;
//      if (!esdfile.Contains("root")) continue; // protection
//
//      if (offset > 0)
//      {
//        --offset;
//        continue;
//      }
//
//      if (count++ == aRuns)
//        break;
//
//      // add esd file
//      chain->Add(esdfile);
//    }
//
//    in.close();
//  }
//
//  return chain;
//}
//
//
//TChain* CreateAODChain(const char* aDataDir, Int_t aRuns, Int_t offset)
//{
//  // creates chain of files in a given directory or file containing a list.
//  // In case of directory the structure is expected as:
//  // <aDataDir>/<dir0>/AliAODs.root
//  // <aDataDir>/<dir1>/AliAODs.root
//  // ...
//
//  if (!aDataDir)
//    return 0;
//
//  Long_t id, size, flags, modtime;
//  if (gSystem->GetPathInfo(aDataDir, &id, &size, &flags, &modtime))
//  {
//    printf("%s not found.\n", aDataDir);
//    return 0;
//  }
//
//  TChain* chain = new TChain("aodTree");
//  TChain* chaingAlice = 0;
//
//  if (flags & 2)
//  {
//    TString execDir(gSystem->pwd());
//    TSystemDirectory* baseDir = new TSystemDirectory(".", aDataDir);
//    TList* dirList            = baseDir->GetListOfFiles();
//    Int_t nDirs               = dirList->GetEntries();
//    gSystem->cd(execDir);
//
//    Int_t count = 0;
//
//    for (Int_t iDir=0; iDir<nDirs; ++iDir)
//    {
//      TSystemFile* presentDir = (TSystemFile*) dirList->At(iDir);
//      if (!presentDir || !presentDir->IsDirectory() || strcmp(presentDir->GetName(), ".") == 0 || strcmp(presentDir->GetName(), "..") == 0)
//        continue;
//
//      if (offset > 0)
//      {
//        --offset;
//        continue;
//      }
//
//      if (count++ == aRuns)
//        break;
//
//      TString presentDirName(aDataDir);
//      presentDirName += "/";
//      presentDirName += presentDir->GetName();	  
//      chain->Add(presentDirName + "/AliAODs.root/aodTree");
//      //  cerr<<presentDirName<<endl;
//    }
//
//  }
//  else
//  {
//    // Open the input stream
//    ifstream in;
//    in.open(aDataDir);
//
//    Int_t count = 0;
//
//    // Read the input list of files and add them to the chain
//    TString aodfile;
//    while(in.good()) {
//      in >> aodfile;
//      if (!aodfile.Contains("root")) continue; // protection
//
//      if (offset > 0)
//      {
//        --offset;
//        continue;
//      }
//
//      if (count++ == aRuns)
//        break;
//
//      // add aod file
//      chain->Add(aodfile);
//    }
//
//    in.close();
//  }
//
//  return chain;
//}
//
