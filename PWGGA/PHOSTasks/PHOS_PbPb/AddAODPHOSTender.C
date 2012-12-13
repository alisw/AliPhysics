AliPHOSTenderTask* AddAODPHOSTender()
{
  //Add a task with PHOS tender which works with AOD to the analysis train
  //Author: D.Peressounko

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddAODPHOSTender", "No analysis manager to connect to");
    return NULL;
  }
  
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddAODPHOSTender", "This task requires an input event handler");
    return NULL;
  }

  // input must be AOD
  TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if( "AOD" != inputDataType )
    ::Error("AddAODPHOSTender", Form("AOD input data required, input data is of type: %s", inputDataType.Data()));

  // create and add task
  AliPHOSTenderTask * tenderTask = new AliPHOSTenderTask("AODPHOSTender") ;
  AliPHOSTenderSupply *PHOSSupply=new AliPHOSTenderSupply("PHOStender");
  PHOSSupply->SetReconstructionPass(1) ;
  tenderTask->SetPHOSTenderSupply(PHOSSupply);

  mgr->AddTask(tenderTask);

  // Connect input/output
  mgr->ConnectInput(tenderTask , 0, mgr->GetCommonInputContainer());

  return tenderTask;
}
