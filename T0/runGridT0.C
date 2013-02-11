{
  Int_t debugLevel = 2;
  // Load common libraries
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
  gSystem->Load("libOADB.so");
  // Use AliRoot includes to compile our task
  gROOT->ProcessLine(".include $ALICE_ROOT/include");

  // Create and configure the alien handler plugin
  gROOT->LoadMacro("CreateAlienHandlerT0.C");
  AliAnalysisGrid *alienHandler = CreateAlienHandler();  
  if (!alienHandler) return;

  // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("testAnalysis");

  // Connect plug-in to the analysis manager
  mgr->SetGridHandler(alienHandler);


  AliESDInputHandler* esdH = new AliESDInputHandler();
  mgr->SetInputEventHandler(esdH);
  
  gROOT->ProcessLine(".L $ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection();
  
  Bool_t isMC=false;
  
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");
    AliCentralitySelectionTask *taskCentrality = AddTaskCentrality();
  
  
  gROOT->LoadMacro("AliT0HIanalysisTask.cxx+g");   
  AliT0HIanalysisTask *task = new AliT0HIanalysisTask("TaskT0");
  task->SetDebugLevel(debugLevel);
  // if you use the following line, your task only gets the selected events
  task->SelectCollisionCandidates(AliVEvent::kINT7 );
   mgr->AddTask(task);


  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += ":Alla_histograms";
  // Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput = mgr->CreateContainer("chist", TList::Class(),    AliAnalysisManager::kOutputContainer, outputFileName);

  // Connect input/output
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput);
  mgr->SetSkipTerminate(kFALSE);
  // Enable debug printouts
  mgr->SetDebugLevel(debugLevel);

  if (!mgr->InitAnalysis())
    return;

  mgr->PrintStatus();
  // Start analysis in grid.
   mgr->StartAnalysis("grid");

};
