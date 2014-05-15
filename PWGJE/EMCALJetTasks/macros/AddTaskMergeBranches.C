// $Id$

AliJetModelMergeBranches* AddTaskMergeBranches(
  const char     *tracksName1   = "Tracks",
  const char     *tracksName2   = "Tracks2",
  const char     *suffix        = "Emb",
  const char     *taskName      = "JetModelMergeBranches"
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskMergeBranches", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskMergeBranches", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  AliJetModelMergeBranches *jetMerge = new AliJetModelMergeBranches(taskName);
  jetMerge->SetTracksName(tracksName1);
  jetMerge->SetTracksMergeName(tracksName2);
  jetMerge->SetClusName("");
  jetMerge->SetSuffix(suffix);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(jetMerge);
    
  // Create containers for input/output
  mgr->ConnectInput (jetMerge, 0, mgr->GetCommonInputContainer() );

  return jetMerge;
}
