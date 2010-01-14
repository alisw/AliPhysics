AliAnalysisTaskSEDs *AddTaskDs()
{
  //                                                                           
  // Test macro for the AliAnalysisTaskSE for Ds candidates 


  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskDs", "No analysis manager to connect to.");
    return NULL;
  }


  // Aanalysis task                                                                                                                     
  AliAnalysisTaskSEDs *dsTask = new AliAnalysisTaskSEDs("DsAnalysis");
  dsTask->SetReadMC(kTRUE);
  dsTask->SetDebugLevel(0);
  mgr->AddTask(dsTask);

  //                                                                                                                                    
  // Create containers for input/output                                                                                                 
  AliAnalysisDataContainer *cinputDs = mgr->CreateContainer("cinputDs",TChain::Class(),
							    AliAnalysisManager::kInputContainer);

  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += ":PWG3_D2H_InvMassDs";
  AliAnalysisDataContainer *coutputDs = mgr->CreateContainer("coutputDs",TList::Class(),
							     AliAnalysisManager::kOutputContainer,
							     outputfile.Data());

  mgr->ConnectInput(dsTask,0,mgr->GetCommonInputContainer());

  mgr->ConnectOutput(dsTask,1,coutputDs);
  
  return dsTask;
}
