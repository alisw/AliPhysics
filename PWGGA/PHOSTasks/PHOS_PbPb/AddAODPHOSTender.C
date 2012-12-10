AliPHOSTenderTask* AddAODPHOSTender()
{
  //Add a task with PHOS tender which works with AOD to the analysis train
  //Author: D.Peressounko

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPHOSPi0Flow", "No analysis manager to connect to");
    return NULL;
  }
  
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskPHOSPi0Flow", "This task requires an input event handler");
    return NULL;
  }

  AliPHOSTenderTask * tenderTask = new AliPHOSTenderTask("AODPHOSTender") ;
  mgr->AddTask(tenderTask);

  AliPHOSTenderSupply *PHOSSupply=new AliPHOSTenderSupply("PHOStender");
  PHOSSupply->SetReconstructionPass(1) ;
  tenderTask->SetPHOSTenderSupply(PHOSSupply);

  AliAnalysisDataContainer *cinput   = mgr->GetCommonInputContainer();

  // Connect input/output
  mgr->ConnectInput(tenderTask , 0, cinput);

  return tenderTask;
}
