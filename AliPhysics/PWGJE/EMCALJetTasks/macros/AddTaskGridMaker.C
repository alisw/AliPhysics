// $Id$

AliEmcalPicoTrackInGridMaker* AddTaskGridMaker(
					       const char *inname       = "PicoTracks",
					       const char *taskName     = "AliEmcalPicoTrackInGridMaker"
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalPicoTrackFromJetMaker", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskEmcalPicoTrackFromJetMaker", "This task requires an input event handler");
    return NULL;
  }

  TString wagonName = Form("%s",taskName);
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  AliEmcalPicoTrackInGridMaker *eTask = new AliEmcalPicoTrackInGridMaker(taskName);
  eTask->AddParticleContainer(inname);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------
  mgr->AddTask(eTask);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer();
  mgr->ConnectInput(eTask, 0, cinput1 );

  //Connect output
  TString contName(wagonName);
  TString outputfile = Form("%s",AliAnalysisManager::GetCommonFileName());
  Printf("outputfile: %s",outputfile.Data());
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contName.Data(), TList::Class(),AliAnalysisManager::kOutputContainer,outputfile);
  mgr->ConnectOutput(eTask,1,coutput1);
  
  return eTask;
}
