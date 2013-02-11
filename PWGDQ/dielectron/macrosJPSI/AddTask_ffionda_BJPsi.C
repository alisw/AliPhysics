AliAnalysisTask *AddTask_ffionda_BJPsi(Bool_t isMC=kFALSE){
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_ffionda_BJPsi", "No analysis manager found.");
    return 0;
  }

  //Do we have an MC handler?
  Bool_t hasMC=isMC;
    //(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
  //Do we have an AOD handler?
  Bool_t isAOD=(mgr->GetInputEventHandler()->IsA()==AliAODInputHandler::Class() ? kTRUE : kFALSE);
  
  //set config file name
  TString configFile("$TRAIN_ROOT/jbook_jpsi/ConfigJpsi_jb_PbPb.C");
  TString trainRoot=gSystem->Getenv("TRAIN_ROOT");                                                                                   
  if (trainRoot.IsNull()) configFile="$ALICE_ROOT/PWGDQ/dielectron/macrosJPSI/ConfigBJpsi_ff_PbPb.C";

  
  //create task and add it to the manager
  AliAnalysisTaskMultiDielectron *task=new AliAnalysisTaskMultiDielectron("MultiDieData");
  task->SetTriggerMask(AliVEvent::kMB+AliVEvent::kCentral+AliVEvent::kSemiCentral);
  if (!hasMC) task->UsePhysicsSelection();
  mgr->AddTask(task);
  

  //load dielectron configuration file
  TString checkconfig="ConfigJpsi_ff_PbPb";
  if (!gROOT->GetListOfGlobalFunctions()->FindObject(checkconfig.Data()))
    gROOT->LoadMacro(configFile.Data());

  //add dielectron analysis with different cuts to the task
  for (Int_t i=0; i<2; ++i){ //nDie defined in config file
    AliDielectron *jpsi=ConfigBJpsi_ff_PbPb(i,hasMC);
    if (jpsi ) task->AddDielectron(jpsi);
    if (jpsi ) printf("add: %s\n",jpsi->GetName());
  }
  
  //Add event filter
  AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
  if(!hasMC) eventCuts->SetRequireVertex();
  if (isAOD) eventCuts->SetVertexType(AliDielectronEventCuts::kVtxAny);
  eventCuts->SetMinVtxContributors(1);
  eventCuts->SetVertexZ(-10.,10.);
  eventCuts->SetCentralityRange(0.0,80.0);
  task->SetEventFilter(eventCuts);

  
  //create output container
  TString containerName = "JPSI.root";
  AliAnalysisDataContainer *coutput1 =
    mgr->CreateContainer("ffionda_tree",
                         TTree::Class(),
                         AliAnalysisManager::kExchangeContainer,
			 containerName.Data());
  
  AliAnalysisDataContainer *cOutputHist1 =
    mgr->CreateContainer("ffionda_QA",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
			 containerName.Data());

  AliAnalysisDataContainer *cOutputHist2 =
    mgr->CreateContainer("ffionda_CF",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
			 containerName.Data());
  
  AliAnalysisDataContainer *cOutputHist3 =
    mgr->CreateContainer("ffionda_EventStat",
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
