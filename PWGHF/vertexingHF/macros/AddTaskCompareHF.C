AliAnalysisTaskSECompareHF *AddTaskCompareHF()
{
  //
  // Test macro for the AliAnalysisTaskSE for heavy-flavour candidates
  // association with MC truth (using MC info in AOD)
  // A.Dainese, andrea.dainese@lnl.infn.it
  //


  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskCompareHF", "No analysis manager to connect to.");
    return NULL;
  }   

  
  // Aanalysis task    
  AliAnalysisTaskSECompareHF *hfTask = new AliAnalysisTaskSECompareHF("CompareHFAnalysis");
  hfTask->SetDebugLevel(0);
  mgr->AddTask(hfTask);
  
  //
  // Create containers for input/output
  AliAnalysisDataContainer *cinputCmp = mgr->CreateContainer("cinput",TChain::Class(), 
							  AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *coutputCmp1 = mgr->CreateContainer("coutputCmp1",TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   "CmpHF.root");
  AliAnalysisDataContainer *coutputCmp2 = mgr->CreateContainer("coutputCmp2",TNtuple::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   "CmpHFnt.root");
  coutputCmp2->SetSpecialOutput();

  mgr->ConnectInput(hfTask,0,mgr->GetCommonInputContainer());

  mgr->ConnectOutput(hfTask,1,coutputCmp1);
  mgr->ConnectOutput(hfTask,2,coutputCmp2);

  return hfTask;
}
