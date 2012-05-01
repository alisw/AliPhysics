// $Id$
AliAnalysisTaskScale* AddTaskScale(
  const char *ntracks        = "Tracks",
  const char *nclusters      = "CaloClustersOut",
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

  AliAnalysisTaskScale *scale = new AliAnalysisTaskScale("scale");
  scale->SetTracksName(ntracks);
  scale->SetClustersName(nclusters);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(scale);

  // Create containers for input/output
  mgr->ConnectInput (scale, 0, mgr->GetCommonInputContainer() );
  AliAnalysisDataContainer *coscale = mgr->CreateContainer(Form("scale_histos"),TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s", outfilename));
  mgr->ConnectOutput(scale,1,coscale);

  return scale;
}
