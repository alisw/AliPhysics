// $Id$

AliJetModelCopyTracks* AddTaskModelCopyTracks(
					     const char     *tracksName1   = "Tracks",
					     const char     *tracksName2   = "Tracks2",
					     Int_t           massType      = AliJetModelCopyTracks::kMassless,
					     const char     *taskName      = "JetModelCopyTracks"
					    )
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskModelCopyTracks", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskModelCopyTracks", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  AliJetModelCopyTracks *copyTask = new AliJetModelCopyTracks(taskName);
  copyTask->SetTracksName(tracksName1);
  copyTask->SetTracksOutName(tracksName2);
  copyTask->SetParticleMassType(massType);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------
  mgr->AddTask(copyTask);
    
  // Create containers for input/output
  mgr->ConnectInput (copyTask, 0, mgr->GetCommonInputContainer() );

  TString contName = taskName;
  contName += "_histos";
  TString outputfile = Form("%s",AliAnalysisManager::GetCommonFileName());
  AliAnalysisDataContainer *outc = mgr->CreateContainer(contName.Data(),
							TList::Class(),
							AliAnalysisManager::kOutputContainer,
							outputfile);
  mgr->ConnectOutput(copyTask, 1, outc);

  return copyTask;
}
