AliMeanVertexCalibTask *AddTaskMeanVertexCalib(){
  //
  //AddTask for MeanVertex Task to run with pass0
  //Author: D.Caffarri davide.caffarri@pd.infn.it
  //

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskMeanVertex", "No analysis manager to connect to.");
    return NULL;
  }


  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskMeanVertex", "This task requires an input event handler");
    return NULL;
  }  

  
  AliMeanVertexCalibTask *meanVertexTask = new AliMeanVertexCalibTask("AliMeanVertexCalibTask");
  
  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
  if (!cinput1) cinput1 = mgr->CreateContainer("cchain",TChain::Class(), AliAnalysisManager::kInputContainer);
  
  
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("MeanVertex", TList::Class(),AliAnalysisManager::kOutputContainer, "AliESDfriends_v1.root"); 
  
  mgr->ConnectInput(meanVertexTask,0,cinput1);
  mgr->ConnectOutput(meanVertexTask,1,coutput1);
  
  return meanVertexTask;
  
}
