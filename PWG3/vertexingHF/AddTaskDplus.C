AliAnalysisTaskSEDplus *AddTaskDplus(Bool_t storeNtuple=kFALSE)
{
  //                                                                                                                                    
  // Test macro for the AliAnalysisTaskSE for D+ candidates 

  //Invariant mass histogram and                                                 
  // association with MC truth (using MC info in AOD)                                                                                   
  //  R. Bala, bala@to.infn.it                                                                                                                                  
  // Get the pointer to the existing analysis manager via the static access method.                                                     
  //==============================================================================                                                      
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskDplus", "No analysis manager to connect to.");
    return NULL;
  }


  // Aanalysis task                                                                                                                     
  AliAnalysisTaskSEDplus *dplusTask = new AliAnalysisTaskSEDplus("DplusAnalysis",storeNtuple);
  dplusTask->SetReadMC(kTRUE);
  dplusTask->SetDoLikeSign(kTRUE);
  dplusTask->SetDebugLevel(0);
  mgr->AddTask(dplusTask);

  //                                                                                                                                    
  // Create containers for input/output                                                                                                 
  AliAnalysisDataContainer *cinputDplus = mgr->CreateContainer("cinputDplus",TChain::Class(),
                                                          AliAnalysisManager::kInputContainer);
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += ":PWG3_D2H_InvMassDplus";
  AliAnalysisDataContainer *coutputDplus = mgr->CreateContainer("coutputDplus",TList::Class(),
                                                           AliAnalysisManager::kOutputContainer,
								outputfile.Data());
  if(storeNtuple){
    AliAnalysisDataContainer *coutputDplus2 = mgr->CreateContainer("coutputDplus2",TNtuple::Class(),
                                                           AliAnalysisManager::kOutputContainer,
								 "InvMassDplus_nt1.root");

    coutputDplus2->SetSpecialOutput();
  }
  mgr->ConnectInput(dplusTask,0,mgr->GetCommonInputContainer());

  mgr->ConnectOutput(dplusTask,1,coutputDplus);
  
  if(storeNtuple){
    mgr->ConnectOutput(dplusTask,2,coutputDplus2);
  }
  return dplusTask;
}
