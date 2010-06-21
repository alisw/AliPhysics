AliAnalysisTaskVertexESD *AddTaskVertexESD(Bool_t readMC=kFALSE) 
{
  //
  // Task for validation of the primary vertices (SPD,TPC,ITS+TPC)
  //
  // andrea.dainese@pd.infn.it
  //


  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTask", "No analysis manager to connect to.");
    return NULL;
  }   

  // Create the task
  AliAnalysisTaskVertexESD *taskVtxESD = new AliAnalysisTaskVertexESD("VertexESD");
  taskVtxESD->SetReadMC(readMC);
  taskVtxESD->SetFillNtuple(kFALSE);
  taskVtxESD->SetRerecoVertexITSTPCHalfEvent(kFALSE);
  taskVtxESD->SetFillNtupleBeamSpot(kTRUE);
  taskVtxESD->SetRerecoVertexITSTPC(kTRUE);
  AliLog::SetClassDebugLevel("AliAnalysisTaskVertexESD",10);
  // Add to the manager
  mgr->AddTask(taskVtxESD);

  //
  // Create containers for input/output
  AliAnalysisDataContainer *cInputVtxESD = mgr->CreateContainer("cInputVtxESD",TChain::Class(),AliAnalysisManager::kInputContainer);

  AliAnalysisDataContainer *cOutputVtxESD = mgr->CreateContainer("cOutputVtxESD",TList::Class(),AliAnalysisManager::kOutputContainer, Form("%s:Vertex_Performance", mgr->GetCommonFileName()));


  // Attach input
  mgr->ConnectInput(taskVtxESD,0,mgr->GetCommonInputContainer());
  // Attach output
  mgr->ConnectOutput(taskVtxESD,1,cOutputVtxESD);
  
  return taskVtxESD;
}
