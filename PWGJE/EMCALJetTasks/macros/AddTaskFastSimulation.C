// $Id$

AliJetFastSimulation* AddTaskFastSimulation(
					    const char     *tracksName1   = "Tracks",
					    const char     *tracksName2   = "Tracks2",
					    const char     *taskName      = "JetFastSimulation"
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

  AliJetFastSimulation *jetFastSim = new AliJetFastSimulation(taskName);
  jetFastSim->AddParticleContainer(tracksName1);
  jetFastSim->SetTracksOutName(tracksName2);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------
  mgr->AddTask(jetFastSim);
    
  // Create containers for input/output
  mgr->ConnectInput (jetFastSim, 0, mgr->GetCommonInputContainer() );

  TString contName = taskName;
  contName += "_histos";
  TString outputfile = Form("%s",AliAnalysisManager::GetCommonFileName());
  AliAnalysisDataContainer *outc = mgr->CreateContainer(contName.Data(),
							TList::Class(),
							AliAnalysisManager::kOutputContainer,
							outputfile);
  mgr->ConnectOutput(jetFastSim, 1, outc);

  return jetFastSim;
}
