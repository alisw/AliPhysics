
//--------------------------------------------------------------------------
// Example macro for running locally AliAnalysisTaskNorm task
//
//
//--------------------------------------------------------------------------

enum {kLocal, kInteractif_xml, kInteractif_ESDList};
Int_t GetMode(TString inputFileName);
TChain* CreateChainFromCollection(const char *xmlfile);
TChain* CreateChainFromAODFile(const char *rootfile,Bool_t isESD);
TChain* CreateChainFromESDList(const char *esdList);
TChain* CreateChain(TString inputFileName,Bool_t isESD);

void runNorm(TString inputFileName = "AliAOD.root", Int_t nEvents = 1e6, TString beamConf = "p-Pb", Bool_t isESD = kFALSE, Bool_t isMC = kFALSE, Int_t debugLevel = 0)
{

  TStopwatch timer;
  timer.Start();
  
  // Check runing mode
  Int_t mode = GetMode(inputFileName);
  if(mode < 0){
    Error("runAnalysis","Please provide either an ESD/AOD root file or a collection of ESDs/AODs.");
    return;
  }
  
  // Load common libraries
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libPhysics");
  gSystem->Load("libSTEERBase");
  //  gSystem->Load("libSTEER");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libEventMixing");
  gSystem->Load("libCORRFW");
  gSystem->Load("libPWGmuon"); 
  gSystem->Load("libPWGmuondep"); 
  
  gSystem->AddIncludePath(Form("-I\"%s/include\" -I\"%s/include\"", gSystem->ExpandPathName("$ALICE_ROOT"),gSystem->ExpandPathName("$ALICE_PHYSICS")));
  gROOT->ProcessLine(Form(".include %s/include %s/include", gSystem->ExpandPathName("$ALICE_ROOT"), gSystem->ExpandPathName("$ALICE_PHYSICS")));

  // Create input chain
  TChain* chain = CreateChain(inputFileName,isESD);
  if (!chain) return;
  
  // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("MuonTask");
  
  if ( isESD) {
    // ESD input handler
    AliESDInputHandler* esdH = new AliESDInputHandler();
    esdH->SetReadFriends(kFALSE);
    mgr->SetInputEventHandler(esdH);
  }
  else {
    // AOD input handler
    AliAODInputHandler* aodH = new AliAODInputHandler();
    mgr->SetInputEventHandler(aodH);
  }

  TString dataType = mgr->GetInputEventHandler()->GetDataType();
  Info("runLocal",Form("Manager with %s",dataType.Data()));

 // Enable MC event handler for ESDs
  if ( isMC && isESD ){
    AliVEventHandler* handler = new AliMCEventHandler;
    mgr->SetMCtruthEventHandler(handler);
  }
  
  // event selection and centrality framework
  if ( isESD ){
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
    AliPhysicsSelectionTask* physicsSelection = AddTaskPhysicsSelection(isMC);
    if ( !physicsSelection ) {
      Error("runLocal","AliPhysicsSelectionTask not created!");
      return;
    }
  }
  if ( isESD ){
    // event centrality
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");
    AliCentralitySelectionTask *taskCentrality = AddTaskCentrality();
    if ( !taskCentrality ) {
      Error("runLocal on ESD","AliCentralitySelectionTask not created!");
      return;
    }
    if ( isMC ) taskCentrality->SetMCInput();     
    //taskCentrality->SetPass(1); // remember to set the pass you are processing!!!
  }
  /*else {
    //Only on full AOD, it is possible to reprocess the centrality framework
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");
    AliCentralitySelectionTask *taskCentrality = AddTaskCentrality(kTRUE,kTRUE);
    if ( !taskCentrality ) {
      Error("runLocal on AOD","AliCentralitySelectionTask not created!");
      return;
    }
    if ( isMC ) taskCentrality->SetMCInput();    
    }*/
  
  // Example analysis
  // Create task
  //Analysis task is on working directory
  //  gROOT->LoadMacro("AliAnalysisTaskNorm.cxx+g");
  //  gROOT->LoadMacro("AddTaskNorm.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/../src/PWG/muon/AliAnalysisTaskNorm.cxx+g");
  gROOT->LoadMacro("$ALICE_PHYSICS/../src/PWG/muon/AddTaskNorm.C");
  AliAnalysisTaskNorm* task = AddTaskNorm(isESD,isMC,beamConf);
  if (!task) {
    Error("runAnalysis","AliAnalysisTaskNorm not created!");
    return;
  }
  task->SetDebugLevel(debugLevel);
  
  // Enable debug printouts
  mgr->SetDebugLevel(debugLevel);

  // start local analysis
  if (mgr->InitAnalysis()) {
    mgr->PrintStatus();
    if(debugLevel>=2) mgr->SetNSysInfo(100);
    mgr->StartAnalysis("local", chain, nEvents);
  }

  if(debugLevel>=2){
    mgr->ProfileTask("Task");
  }

  timer.Stop();
  timer.Print();
}

//______________________________________________________________________________
Int_t GetMode(TString inputFileName)
{
  if ( inputFileName.EndsWith(".xml") ) return kInteractif_xml;
  else if ( inputFileName.EndsWith(".txt") ) return kInteractif_ESDList;
  else if ( inputFileName.EndsWith(".root") ) return kLocal;
  return -1;
}

//______________________________________________________________________________
TChain* CreateChainFromCollection(const char *xmlfile)
{
  // Create a chain from the collection of tags.
  TGridCollection* coll = gGrid->OpenCollection(xmlfile);
  if (!coll) {
    ::Error("CreateChainFromTags", "Cannot create an AliEn collection from %s", xmlfile);
    return NULL;
  }
  
  TGridResult* tagResult = coll->GetGridResult("",kFALSE,kFALSE);
  AliTagAnalysis *tagAna = new AliTagAnalysis("ESD");
  tagAna->ChainGridTags(tagResult);
  
  AliRunTagCuts      *runCuts = new AliRunTagCuts();
  AliLHCTagCuts      *lhcCuts = new AliLHCTagCuts();
  AliDetectorTagCuts *detCuts = new AliDetectorTagCuts();
  AliEventTagCuts    *evCuts  = new AliEventTagCuts();
  
  // Check if the cuts configuration file was provided
  if (!gSystem->AccessPathName("ConfigureCuts.C")) {
    gROOT->LoadMacro("ConfigureCuts.C");
    ConfigureCuts(runCuts, lhcCuts, detCuts, evCuts);
  }
  
  TChain *chain = tagAna->QueryTags(runCuts, lhcCuts, detCuts, evCuts);
  if (!chain || !chain->GetNtrees()) return NULL;
  chain->ls();
  return chain;
}

//______________________________________________________________________________
TChain* CreateChainFromAODFile(const char *rootfile, Bool_t isESD)
{
  // Create a chain using the root file.
  TChain* chain = 0;
  if ( !isESD) chain = new TChain("aodTree");
  else chain = new TChain("esdTree");
  chain->Add(rootfile);
  if (!chain->GetNtrees()) return NULL;
  chain->ls();
  return chain;
}

//______________________________________________________________________________
TChain* CreateChainFromESDList(const char *esdList)
{
  // Create a chain using tags from the run list.
  TChain* chain = new TChain("esdTree");
  ifstream inFile(esdList);
  TString inFileName;
  if (inFile.is_open()) {
    while (! inFile.eof() ) {
      inFileName.ReadLine(inFile,kFALSE);
      if(!inFileName.EndsWith(".root")) continue;
      chain->Add(inFileName.Data());
    }
  }
  inFile.close();
  if (!chain->GetNtrees()) return NULL;
  chain->ls();
  return chain;
}

//______________________________________________________________________________
TChain* CreateChain(TString inputFileName, Bool_t isESD)
{
  printf("*******************************\n");
  printf("*** Getting the Chain       ***\n");
  printf("*******************************\n");
  Int_t mode = GetMode(inputFileName);
  if(mode == kInteractif_xml) return CreateChainFromCollection(inputFileName.Data());
  else if (mode == kInteractif_ESDList) return CreateChainFromESDList(inputFileName.Data());
  else if (mode == kLocal) return CreateChainFromAODFile(inputFileName.Data(),isESD);
  else return NULL;
}

