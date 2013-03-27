AliAnalysisTask *AddTask_mkoehler_CCbar(Bool_t hasMC_aod = kFALSE){
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTask_mkoehler_CCbar", "No analysis manager found.");
    return NULL;
  }
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTask_mkoehler_CCbar", "This task requires an input event handler");
    return NULL;
  }

  //Do we have an MC handler?
  Bool_t hasMC=(mgr->GetMCtruthEventHandler()!=0x0);
  
  TString configFile("ConfigCCbar_mk_pp.C");
  if (hasMC) configFile="$ALICE_ROOT/PWGDQ/dielectron/macros/ConfigJpsi2eeEff.C";//check later
  Bool_t isAOD=mgr->GetInputEventHandler()->IsA()==AliAODInputHandler::Class();

  //create task and add it to the manager
  AliAnalysisTaskMultiDielectron *task=new AliAnalysisTaskMultiDielectron("MultiDie");
  mgr->AddTask(task);
  
  //load dielectron configuration file
  TString checkconfig="ConfigCCbar_mk_pp";
  if (!gROOT->GetListOfGlobalFunctions()->FindObject(checkconfig.Data()))
    gROOT->LoadMacro(configFile.Data());
  
  //add dielectron analysis with different cuts to the task
  for (Int_t i=0; i<nDie; ++i){ //nDie defined in config file
    AliDielectron *jpsi=ConfigCCbar_mk_pp(i,isAOD);
    if (isAOD) jpsi->SetHasMC(hasMC_aod);
    if (jpsi) task->AddDielectron(jpsi);
  }

  //Add event filter
  AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
  if(isAOD) eventCuts->SetVertexType(AliDielectronEventCuts::kVtxAny);
  eventCuts->SetRequireVertex();
  eventCuts->SetMinVtxContributors(1);
  eventCuts->SetVertexZ(-10.,10.);
  task->SetEventFilter(eventCuts);

  // pileup rejection
//    task->SetRejectPileup();
   task->SetTriggerMask(AliVEvent::kEMC7);
   task->UsePhysicsSelection();
  //----------------------
  //create data containers
  //----------------------


  //create output container
  
  AliAnalysisDataContainer *cOutputHist1 =
    mgr->CreateContainer("mkoehler_pp_QA", TList::Class(), AliAnalysisManager::kOutputContainer,
                         "LMEEoutput.root");
  
  AliAnalysisDataContainer *cOutputHist2 =
    mgr->CreateContainer("mkoehler_pp_CF", TList::Class(), AliAnalysisManager::kOutputContainer,
                         "LMEEoutput.root");

  AliAnalysisDataContainer *cOutputHist3 =
    mgr->CreateContainer("mkoehler_pp_EventStat", TH1D::Class(), AliAnalysisManager::kOutputContainer,
                        "LMEEoutput.root");
  
  mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, cOutputHist1);
  mgr->ConnectOutput(task, 2, cOutputHist2);
  mgr->ConnectOutput(task, 3, cOutputHist3);
  
  return task;
}
