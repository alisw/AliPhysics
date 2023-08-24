AliAnalysisCuts* CreateEventFilter(Bool_t isAOD);

//__________________________________________________________________________________________
AliAnalysisTask* AddTask_iarsene_ComputeLumi(){
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_iarsene_ComputeLumi", "No analysis manager found.");
    return 0;
  }

  Bool_t isAOD=mgr->GetInputEventHandler()->IsA()==AliAODInputHandler::Class();
  
  //create task and add it to the manager
  AliAnalysisTaskComputeLumi *task=new AliAnalysisTaskComputeLumi("ComputeLumi");
  
  task->SetEventFilter(CreateEventFilter(isAOD));
  
  mgr->AddTask(task);
    
  AliAnalysisDataContainer *cEventStats = mgr->CreateContainer("histos", TList::Class(),
                                           AliAnalysisManager::kOutputContainer, "ComputeLumi.root");
    
  mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, cEventStats);
  
  return task;
}

//_______________________________________________________________________________________________
AliAnalysisCuts* CreateEventFilter(Bool_t isAOD) {
  //
  // Event wise cuts
  //
  AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
  if(isAOD) eventCuts->SetVertexType(AliDielectronEventCuts::kVtxAny);
  eventCuts->SetRequireVertex();
  eventCuts->SetMinVtxContributors(1);
  eventCuts->SetVertexZ(-10.,10.);
  return eventCuts;
}
