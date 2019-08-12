AliAnalysisTaskSECleanupVertexingHF *AddTaskCleanupVertexingHF()
{
  //
  // AddTask for the AliAnalysisTaskSECleanupVertexingHF to delete secondary vertex of the candidates
  // Get the pointer to the existing analysis manager via the static access method.
  //=================================================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskD0Distr", "No analysis manager to connect to.");
    return NULL;
  }   
  TString filename="";
  filename = AliAnalysisManager::GetCommonFileName();
  filename += ":PWG3_D2H_Delete";

  TString name = "DeleteTask";
  AliAnalysisTaskSECleanupVertexingHF *secVertDelete = new AliAnalysisTaskSECleanupVertexingHF(name.Data());

  mgr->AddTask(secVertDelete);


//Create containers for input/output
 name = "cinputDeleteMY";
 AliAnalysisDataContainer *cinputDelete = mgr->CreateContainer(name,TChain::Class(),
                                                            AliAnalysisManager::kInputContainer);

//TString out1name="nEntriesD0_delete";
//AliAnalysisDataContainer *coutputmassD01 = mgr->CreateContainer(out1name,TH1F::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //nev




 //mgr->ConnectInput(secVertDelete,0,cinputDelete);
 mgr->ConnectInput(secVertDelete,0,mgr->GetCommonInputContainer());
// mgr->ConnectOutput(secVertDelete,1,coutputmassD01);
  return secVertDelete;
}
