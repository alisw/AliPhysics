AliAnalysisTaskSEDplus *AddTaskDplus()
{
  //                                                                                                                                    
  // Test macro for the AliAnalysisTaskSE for heavy-flavour candidates                                                                  
  // association with MC truth (using MC info in AOD)                                                                                   
                                                                                               
  //                                                                                                                                    


  // Get the pointer to the existing analysis manager via the static access method.                                                     
  //==============================================================================                                                      
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskDplus", "No analysis manager to connect to.");
    return NULL;
  }


  // Aanalysis task                                                                                                                     
  AliAnalysisTaskSEDplus *dplusTask = new AliAnalysisTaskSEDplus("DplusAnalysis");
  dplusTask->SetDebugLevel(2);
  mgr->AddTask(dplusTask);

  //                                                                                                                                    
  // Create containers for input/output                                                                                                 
  AliAnalysisDataContainer *cinputDplus = mgr->CreateContainer("cinputDplus",TChain::Class(),
                                                          AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *coutputDplus = mgr->CreateContainer("coutputDplus",TList::Class(),
                                                           AliAnalysisManager::kOutputContainer,
                                                           "InvMassDplus.root");
  mgr->ConnectInput(dplusTask,0,mgr->GetCommonInputContainer());

  mgr->ConnectOutput(dplusTask,1,coutputDplus);

  return dplusTask;
}
