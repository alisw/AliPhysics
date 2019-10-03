void analyzeJpsiME(TString tag="./"){

  TStopwatch timer;
  timer.Start();
  
  //_____________Setting up libraries_______________
  gSystem->Load("libTree");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libVMC");
  gSystem->Load("libESD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libPWGDQdielectron");
  gSystem->AddIncludePath("-I$ALICE_ROOT/PWGDQ/dielectron/ -I$ALICE_ROOT/ANALYSIS -I$ALICE_ROOT/CORRFW");

  //_____________Load Macro_____________________________
  gROOT->LoadMacro("ConfigJpsi2eeData.C");    // user Config file
  
  // Make the analysis manager__________________________
  AliAnalysisManager *mgr       = new AliAnalysisManager("AnalysisManager");

  // Event handler setting
  Int_t nEventsToLoad = 5;    // pool depth
  AliMultiEventInputHandler *inputHandler  = new AliMultiEventInputHandler(nEventsToLoad,0);  // 0 for ESD, 1 for AOD
  mgr->SetInputEventHandler(inputHandler);
  
  // event pool setting
  AliEventPoolOTF *myPool = new AliEventPoolOTF("event pool","ESD");
  myPool->SetTagDirectory(tag.Data());
  myPool->SetMultiplicityBin(0,100,100);  // (min, max, bin width)
  myPool->Init();

  mgr->SetEventPool(myPool);            // link event pool with manager
  inputHandler->SetEventPool(myPool);   // link event pool with input handler

  // Analysis Task Jpsi->e+e-___________________________
  AliAnalysisTaskDielectronME *task = new AliAnalysisTaskDielectronME("DielectronTaskME");
  task->RequireFreshBuffer();         // refresh mode in the pool
  task->SetPoolDepth(nEventsToLoad);  // set pool depth in the task
  mgr->AddTask(task);

  // add dielectron analysis with different cuts to the task
  for (Int_t i=0; i<nDie; ++i){       //nDie defined in config file
    AliDielectron *jpsi=ConfigJpsi2ee(i);
    task->AddDielectron(jpsi);
  }

  // Create containers for input/output_________________
  AliAnalysisDataContainer *coutput1     = mgr->CreateContainer("jpsi_treeEff",TTree::Class(),AliAnalysisManager::kExchangeContainer,"jpsi_Effdefault");
  AliAnalysisDataContainer *cOutputHist1 = mgr->CreateContainer("QA", TList::Class(), AliAnalysisManager::kOutputContainer, "jpsi_QA.root");
  AliAnalysisDataContainer *cOutputHist2 = mgr->CreateContainer("CF", TList::Class(), AliAnalysisManager::kOutputContainer, "jpsi_CF.root");
  AliAnalysisDataContainer *cOutputHist3 = mgr->CreateContainer("jpsi_Eff_EventStat",TH1D::Class(),AliAnalysisManager::kOutputContainer,"jpsi_Eff.root");
  
  mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 0, coutput1 );
  mgr->ConnectOutput(task, 1, cOutputHist1);
  mgr->ConnectOutput(task, 2, cOutputHist2); 
  mgr->ConnectOutput(task, 3, cOutputHist3);

  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  TChain *chain = NULL;     // null pointer to avoid entering grid mode. temporary solution
  mgr->StartAnalysis("mix",chain);
  //mgr->Terminate();

  timer.Stop();
  timer.Print();
}
