AliAlignmentDataFilterITS *AddTaskAlignmentDataFilterITS() 
{
  //
  // Task for the extraction of the ITS alignment data (AliTrackPoints)
  //
  // andrea.dainese@pd.infn.it
  //


  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskCompareHF", "No analysis manager to connect to.");
    return NULL;
  }   

  // Create the task
  AliAlignmentDataFilterITS *taskFilter = new AliAlignmentDataFilterITS("filterITS");
  AliLog::SetClassDebugLevel("AliAlignmentDataFilterITS",10);
  // configuration via AliITSRecoParam (should be taken from OCDB)
  AliITSRecoParam *itsRecoParam = AliITSRecoParam::GetCosmicTestParam();
  itsRecoParam->SetAlignFilterUseLayer(0,kTRUE);
  itsRecoParam->SetAlignFilterUseLayer(1,kTRUE);
  itsRecoParam->SetAlignFilterUseLayer(2,kTRUE);
  itsRecoParam->SetAlignFilterUseLayer(3,kTRUE);
  itsRecoParam->SetAlignFilterUseLayer(4,kTRUE);
  itsRecoParam->SetAlignFilterUseLayer(5,kTRUE);
  taskFilter->SetITSRecoParam(itsRecoParam);
  //taskFilter->SetOnlySPDFO();
  taskFilter->SetGeometryFileName("alien:///alice/cern.ch/user/d/dainesea/geometry.root");
  // Add to the manager
  mgr->AddTask(taskFilter);

  //
  // Create containers for input/output
  AliAnalysisDataContainer *cInput = mgr->CreateContainer("cInput",TChain::Class(),AliAnalysisManager::kInputContainer);

  AliAnalysisDataContainer *cOutput0= mgr->CreateContainer("cOutput0",TTree::Class(),AliAnalysisManager::kOutputContainer,"AliTrackPoints.root");
  AliAnalysisDataContainer *cOutput1= mgr->CreateContainer("cOutput1",TList::Class(),AliAnalysisManager::kOutputContainer,"AliTrackPoints.root");


  // Attach input
  mgr->ConnectInput(taskFilter,0,mgr->GetCommonInputContainer());
  // Attach output
  mgr->ConnectOutput(taskFilter,0,cOutput0);
  mgr->ConnectOutput(taskFilter,1,cOutput1);
  
  return taskFilter;
}
