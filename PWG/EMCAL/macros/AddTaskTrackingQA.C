AliEmcalTrackingQATask* AddTaskTrackingQA(const char *nDetLev      = "PicoTracks",
                                          const char *nGenLev      = "mcparticles")
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskTrackingQA", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskTrackingQA", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  TString name(Form("AliEmcalTrackingQATask_%s_%s", nGenLev, nDetLev));
  AliEmcalTrackingQATask *qaTask = new AliEmcalTrackingQATask(name);
  qaTask->SetDetectorLevelName(nDetLev);
  if (nGenLev && strcmp(nGenLev,"")!=0) qaTask->SetGeneratorLevelName(nGenLev);
  qaTask->SetVzRange(-10,10);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(qaTask);
    
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;

  TString contName(name);
  contName += "_histos";
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contName.Data(), 
							    TList::Class(),AliAnalysisManager::kOutputContainer,
							    Form("%s", AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectInput  (qaTask, 0,  cinput1 );
  mgr->ConnectOutput (qaTask, 1, coutput1 );
    
  return qaTask;
}
