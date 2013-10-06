AliAnalysisTask* AddTask_jpsi_Default(TString prod="")
{
  //get the current analysis manager

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_jpsi_JPsi", "No analysis manager found.");
    return 0;
  }

  //Do we have an MC handler?
  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
  //AOD input?
  Bool_t isAOD=mgr->GetInputEventHandler()->IsA()==AliAODInputHandler::Class();
  
  //Get the current train configuration
  TString trainConfig=gSystem->Getenv("CONFIG_FILE");
  
  //set config file name
  TString configFile("$ALICE_ROOT/PWGDQ/dielectron/macrosJPSI/ConfigDefault.C");
  
  TString list=gSystem->Getenv("LIST");
  if( list.IsNull()) list=prod;
  
  //create task and add it to the manager
  AliAnalysisTaskMultiDielectron *task=new AliAnalysisTaskMultiDielectron("JpsiDefault");
  if (!hasMC) task->UsePhysicsSelection();
  if (list.Contains("LHC11d")) task->SetTriggerMask(AliVEvent::kEMCEJE+AliVEvent::kEMC7+AliVEvent::kEMCEGA);
  if (list.Contains("LHC11h")) {
    task->SetTriggerMask(AliVEvent::kMB+AliVEvent::kCentral+AliVEvent::kSemiCentral);
    printf("AddTask_jpsi_Default: LHC11h data detected. Using trigger classes: %d\n",AliVEvent::kMB+AliVEvent::kCentral+AliVEvent::kSemiCentral);
  }
  if (list.Contains("LHC12h")) task->SetTriggerMask(AliVEvent::kAnyINT);
  mgr->AddTask(task);

  //load dielectron configuration file
  gROOT->LoadMacro(configFile.Data());
  
  //add dielectron analysis with different cuts to the task
  for (Int_t i=0; i<nDie; ++i){ //nDie defined in config file
    AliDielectron *jpsi=ConfigDefault(i);
    if (!jpsi) continue;
    task->AddDielectron(jpsi);
  }
  
  //Add event filter
  AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
  if (isAOD) eventCuts->SetVertexType(AliDielectronEventCuts::kVtxAny);
  eventCuts->SetRequireVertex();
  eventCuts->SetMinVtxContributors(1);
  eventCuts->SetVertexZ(-10.,10.);
  task->SetEventFilter(eventCuts);

  //   task->SetTriggerOnV0AND();
//   if ( trainConfig=="pp" ) task->SetRejectPileup();
  
  //create output container
  AliAnalysisDataContainer *coutput1 =
    mgr->CreateContainer("jpsi_Default_tree",
                         TTree::Class(),
                         AliAnalysisManager::kExchangeContainer,
                         "jpsi_Default_default");
  
  AliAnalysisDataContainer *cOutputHist1 =
    mgr->CreateContainer("jpsi_Default_QA",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         "jpsi_Default.root");

  AliAnalysisDataContainer *cOutputHist2 =
    mgr->CreateContainer("jpsi_Default_CF",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         "jpsi_Default.root");
  
  AliAnalysisDataContainer *cOutputHist3 =
    mgr->CreateContainer("jpsi_Default_EventStat",
                         TH1D::Class(),
                         AliAnalysisManager::kOutputContainer,
                         "jpsi_Default.root");
  
  mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 0, coutput1 );
  mgr->ConnectOutput(task, 1, cOutputHist1);
  mgr->ConnectOutput(task, 2, cOutputHist2);
  mgr->ConnectOutput(task, 3, cOutputHist3);
  
  return task;
}
