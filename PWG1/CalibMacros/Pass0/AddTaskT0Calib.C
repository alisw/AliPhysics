//_____________________________________________________________________________
AliAnalysisTask  *AddTaskT0Calib(Int_t runNumber)
{
  //
  // add calibration task
  //
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskT0Calib", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // check the input handler
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskT0Calib", "This task requires an input event handler");
    return NULL;
  }  

  // set TPC OCDB parameters
  //ConfigOCDB(runNumber);

  // setup task
 AliT0CalibOffsetChannelsTask  *task1 = new AliT0CalibOffsetChannelsTask("CalibObjectsTrain1");
   mgr->AddTask(task1);

  //  AliT0AnalysisTaskQA * task2 = new AliT0AnalysisTaskQA("QA task");
  //    mgr->AddTask(task2);

  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
  if (!cinput1) cinput1 = mgr->CreateContainer("cchain",TChain::Class(), 
                                      AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("T0Calib",TObjArray::Class(), AliAnalysisManager::kOutputContainer, "AliESDfriends_v1.root");  

  mgr->ConnectInput(task1,0,cinput1);
  mgr->ConnectOutput(task1,1,coutput1);
  return task1;
}
