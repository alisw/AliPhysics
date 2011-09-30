AliAnalysisTaskSpectraAOD * AddTaskSpectraAOD(const char * outfilename)
{
  // TODO: add some parameters to set the centrality for this task, and maybe the name of the task
  // TODO: shall I use the same file and different dirs for the different centralities?

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskSpectraAOD", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskSpectraAOD", "This task requires an input event handler");
    return NULL;
  }
  TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  
  if (inputDataType != "AOD") {
    Printf("ERROR! This task can only run on AODs!");
  }

  // Configure analysis
  //===========================================================================
    
   
  // Set I/O
  AliAnalysisDataContainer *cinput0 = mgr->GetCommonInputContainer();
  AliAnalysisTaskSpectraAOD *task = new AliAnalysisTaskSpectraAOD("TaskAODSpectra");
  mgr->AddTask(task);
 
  // Set the cuts
  AliSpectraAODVertexCuts * vcuts = new AliSpectraAODVertexCuts("VertexCuts");
  AliSpectraAODTrackCuts  * tcuts = new AliSpectraAODTrackCuts ("TracksCuts");
  tcuts->SetTrackType(6);
  tcuts->SetEta(1.);
  task->SetVertexCuts(vcuts);
  task->SetTrackCuts (tcuts);


  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutputpt1 = mgr->CreateContainer("chistpt", AliSpectraAODHistoManager::Class(),  AliAnalysisManager::kOutputContainer, outfilename);
  AliAnalysisDataContainer *coutputpt2 = mgr->CreateContainer("cvcutpt", AliSpectraAODVertexCuts::Class(),    AliAnalysisManager::kOutputContainer, outfilename);
  AliAnalysisDataContainer *coutputpt3 = mgr->CreateContainer("ctcutpt", AliSpectraAODTrackCuts::Class(),     AliAnalysisManager::kOutputContainer, outfilename);

  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutputpt1);
  mgr->ConnectOutput(task, 2, coutputpt2);
  mgr->ConnectOutput(task, 3, coutputpt3);

  return task;
}   
