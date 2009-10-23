AliAnalysisTaskSEDplus *AddTaskDplus(Bool_t storeNtuple=kFALSE)
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
  AliAnalysisTaskSEDplus *dplusTask = new AliAnalysisTaskSEDplus("DplusAnalysis",storeNtuple);
  dplusTask->SetDebugLevel(0);
  mgr->AddTask(dplusTask);

  //                                                                                                                                    
  // Create containers for input/output                                                                                                 
  AliAnalysisDataContainer *cinputDplus = mgr->CreateContainer("cinputDplus",TChain::Class(),
                                                          AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *coutputDplus = mgr->CreateContainer("coutputDplus",TList::Class(),
                                                           AliAnalysisManager::kOutputContainer,
                                                           "InvMassDplus.root");
  if(storeNtuple){
    AliAnalysisDataContainer *coutputDplus2 = mgr->CreateContainer("coutputDplus2",TNtuple::Class(),
                                                           AliAnalysisManager::kOutputContainer,
								 "InvMassDplus_nt1.root");
    AliAnalysisDataContainer *coutputDplus3 = mgr->CreateContainer("coutputDplus3",TNtuple::Class(),
                                                           AliAnalysisManager::kOutputContainer,
								 "InvMassDplus_nt2.root");
    coutputDplus2->SetSpecialOutput();
    coutputDplus3->SetSpecialOutput();
  }
  mgr->ConnectInput(dplusTask,0,mgr->GetCommonInputContainer());

  mgr->ConnectOutput(dplusTask,1,coutputDplus);
  
  if(storeNtuple){
    mgr->ConnectOutput(dplusTask,2,coutputDplus2);
    mgr->ConnectOutput(dplusTask,3,coutputDplus3);
  }
  return dplusTask;
}
