AliAnalysisTask *AddTask_shin_pPb(){
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_shin_pPb", "No analysis manager found.");
    return 0;
  }


  //Do we have an MC handler?
  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
  
  //Get the current train configuration
  //  TString trainConfig=gSystem->Getenv("CONFIG_FILE");
  
  //set config file name
  TString configBasePath("$ALICE_ROOT/PWGDQ/dielectron/macrosLMEE/");
  TString configFile("Config_shin_pPb.C");
  TString configFilePath(configBasePath+configFile);

  // TString list=gSystem->Getenv("LIST");
  //create task and add it to the manager
  AliAnalysisTaskMultiDielectron *task=new AliAnalysisTaskMultiDielectron("MultiDieData");
  if (!hasMC ) task->UsePhysicsSelection();
  task->SetTriggerMask(AliVEvent::kINT7); 
  mgr->AddTask(task);

  
  //load dielectron configuration file
  gROOT->LoadMacro(configFilePath.Data());
  
  //add dielectron analysis with different cuts to the task
   for (Int_t i=0; i<nDie; ++i){ //nDie defined in config file
    AliDielectron *dile=Config_shin_pPb(i);
    if (!dile) continue;
    task->AddDielectron(dile);
  }
  
  //Add event filter
  AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
  eventCuts->SetRequireVertex();
  eventCuts->SetMinVtxContributors(1);
  eventCuts->SetVertexZ(-10.,10.);
  task->SetEventFilter(eventCuts);
  
  //create output container
  TString containerName = "hayashi_lowmass.root";
  AliAnalysisDataContainer *coutput1 =
    mgr->CreateContainer("tree_lowmass",
                         TTree::Class(),
                         AliAnalysisManager::kExchangeContainer,
                         containerName.Data());
  
  AliAnalysisDataContainer *cOutputHist1 =
    mgr->CreateContainer("Histos_diel_lowmass",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         containerName.Data());

  AliAnalysisDataContainer *cOutputHist2 =
    mgr->CreateContainer("CF_diel_lowmass",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         containerName.Data());
  
  AliAnalysisDataContainer *cOutputHist3 =
    mgr->CreateContainer("EventStat",
                         TH1D::Class(),
                         AliAnalysisManager::kOutputContainer,
                         containerName.Data());
  
  mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 0, coutput1 );
  mgr->ConnectOutput(task, 1, cOutputHist1);
  mgr->ConnectOutput(task, 2, cOutputHist2);
  mgr->ConnectOutput(task, 3, cOutputHist3);
  
  return task;
}
