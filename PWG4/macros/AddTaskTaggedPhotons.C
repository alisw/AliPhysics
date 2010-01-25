AliAnalysisTaskTaggedPhotons* AddTaskTaggedPhotons(Bool_t bPHOS = kFALSE)
{
  // Creates an tagged photons task, 
  // configures it and adds it to the analysis manager.

  // Get the pointer to the existing analysis manager via the static access method.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskTaggedPhotons", "No analysis manager to connect to.");
    return NULL;
  }  

  // Check the analysis type using the event handlers connected to the analysis manager.
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskTaggedPhotons", "This task requires an input event handler");
    return NULL;
  }
  
  // Add task
  
  TString det;
  if(bPHOS)det = "PHOS";
  else det = "EMCAL";
  
  AliAnalysisTaskTaggedPhotons * task = new AliAnalysisTaskTaggedPhotons(Form("Tagged%s",det.Data()));
  task->SetPHOS(bPHOS);
  mgr->AddTask(task);



  TString outputfile = AliAnalysisManager::GetCommonFileName();                                                              
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(Form("tagged%s",det.Data()), TList::Class(),
							    AliAnalysisManager::kOutputContainer, Form("%s:PWG4_tagged%s",outputfile.Data(),det.Data()));
  mgr->ConnectInput  (task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (task, 0, mgr->GetCommonOutputContainer());

  mgr->ConnectOutput (task,    1, coutput2 );
  
  return task;
}


