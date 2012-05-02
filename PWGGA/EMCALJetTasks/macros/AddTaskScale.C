// $Id$

AliAnalysisTaskScale* AddTaskScale(
  const char *nTracks        = "Tracks",
  const char *nClusters      = "CaloClustersCorr",
  const char *outfilename    = "AnalysisResults.root"
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskScale", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskScale", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  TString name(Form("Scale_%s", nClusters));
  AliAnalysisTaskScale *scaletask = new AliAnalysisTaskScale(name);
  scaletask->SetTracksName(nTracks);
  scaletask->SetClustersName(nClusters);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(scaletask);

  // Create containers for input/output
  mgr->ConnectInput (scaletask, 0, mgr->GetCommonInputContainer() );
  AliAnalysisDataContainer *coscale = mgr->CreateContainer(name,
                                                           TList::Class(),
                                                           AliAnalysisManager::kOutputContainer,
                                                           outfilename);
  mgr->ConnectOutput(scaletask,1,coscale);

  return scaletask;
}
