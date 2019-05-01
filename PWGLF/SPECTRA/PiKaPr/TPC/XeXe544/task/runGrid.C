void runGrid() {
  // Load common libraries

  Bool_t local = kTRUE; // kTRUE when running locally. kFALSE when running on grid 
  
  gSystem->Load("libCore.so");  
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  
  // Use AliRoot includes to compile our task
  // gROOT->ProcessLine(".include $ALICE_ROOT/include");
  
  // gROOT->ProcessLine(".include $ALICE_PHYSICS/include");
  gSystem->AddIncludePath("-I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");

  // Create and configure the alien handler plugin
  gROOT->LoadMacro("CreateAlienHandler.C");
  AliAnalysisGrid *alienHandler = CreateAlienHandler();  
  if (!alienHandler) return;

  // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("testAnalysis");
  AliESDInputHandler* esdH = new AliESDInputHandler();
  esdH->SetNeedField(kTRUE);
  mgr->SetInputEventHandler(esdH);

  // Connect plug-in to the analysis manager
  if (!local) mgr->SetGridHandler(alienHandler);

  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask * physSelTask = AddTaskPhysicsSelection(kTRUE); //kTRUE while running on MC,kFALSE for data

  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
  AliMultSelectionTask * multSelTask = AddTaskMultSelection();

  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  AddTaskPIDResponse();

  AliMCEventHandler* handler = new AliMCEventHandler;
  handler->SetReadTR(kFALSE);
  mgr->SetMCtruthEventHandler(handler); // comment while running on Data

  gROOT->LoadMacro("AliAnalysisTaskTpcSpectra.cxx++g");   
  AliAnalysisTaskTpcSpectra* task = new AliAnalysisTaskTpcSpectra("SpectraTask");
  task->SelectCollisionCandidates(AliVEvent::kINT7);

  mgr->AddTask(task);  

  // Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput = mgr->CreateContainer("chist", TList::Class(),AliAnalysisManager::kOutputContainer, "PiKpOutput.root");

  // Connect input/output
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput);

  // Enable debug printouts
  mgr->SetDebugLevel(20);

  if (!mgr->InitAnalysis())
    return;

  mgr->PrintStatus();
  // Start analysis in grid.

  if (local) {
    TChain* chain = new TChain("esdTree");
    // add a few files to the chain (change this so that your local files are added)
    chain->Add("root_archive.zip#AliESDs.root");
    // start the analysis locally, reading the events from the tchain
    mgr->StartAnalysis("local", chain);
  } else {
    mgr->StartAnalysis("grid");
  }
};
