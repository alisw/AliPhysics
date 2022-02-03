AliAnalysisTaskUpc2Pi2E *AddTaskUpc2Pi2E( 
    const char* outputFileName = 0,
    const char* folderName = "2Pi2ETree1"
){

  //--- get the current analysis manager ---//
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
      Error("AddTask_Upc2Pi2E", "No analysis manager found.");
      return 0;
   }
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    Error("AddTask_Upc2Pi2E", "This task requires an input event handler");
    return 0;
  }
	
  TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  Bool_t isMC = kFALSE;
  if(mgr->GetMCtruthEventHandler()) isMC = kTRUE;
  
  // Create tasks
  AliAnalysisTaskUpc2Pi2E *task = new AliAnalysisTaskUpc2Pi2E(inputDataType.Data(), isMC);
  mgr->AddTask(task);
  task->SetIsMC(isMC);

   // Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  if (!outputFileName) outputFileName = AliAnalysisManager::GetCommonFileName();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("2Pi2ETree1", TTree::Class(),AliAnalysisManager::kOutputContainer,Form("%s:%s", outputFileName, folderName));
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("Histograms", TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:%s", outputFileName, folderName));

//  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("4PiTree", TTree::Class(), AliAnalysisManager::kOutputContainer,Form("%s:4PiCentral", AliAnalysisManager::GetCommonFileName()));
//  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("Histograms", TList::Class(), AliAnalysisManager::kOutputContainer,Form("%s:4PiCentral", AliAnalysisManager::GetCommonFileName()));
//  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("4PiTree1", TTree::Class(),AliAnalysisManager::kOutputContainer,Form("%s:%s", outputFileName, folderName));
  AliAnalysisDataContainer *coutput4 = NULL;
  if (isMC) coutput4 = mgr->CreateContainer("MCTree", TTree::Class(), AliAnalysisManager::kOutputContainer,Form("%s:Rho0Central", AliAnalysisManager::GetCommonFileName()));

  // Connect input/output
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput1);
  mgr->ConnectOutput(task, 2, coutput2);
//  mgr->ConnectOutput(task, 3, coutput3);
  if (isMC) mgr->ConnectOutput(task, 4, coutput4);
  

return task;
}
