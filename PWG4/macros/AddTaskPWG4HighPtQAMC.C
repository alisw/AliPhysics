//DEFINITION OF A FEW CONSTANTS

AliPWG4HighPtQAMC* AddTaskPWG4HighPtQAMC()
{
  // Creates HighPtQAMC analysis task and adds it to the analysis manager.
  
  // A. Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskPWG4HighPtQMC", "No analysis manager to connect to.");
    return NULL;
  }  

  // B. Check the analysis type using the event handlers connected to the analysis
  //    manager. The availability of MC handler can also be checked here.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddPWG4TaskHighPtQAMC", "This task requires an input event handler");
    return NULL;
  }  
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  const char *analysisType = "ESD";//"TPC"

  // C. Create the task, add it to manager.
  //===========================================================================
 
  //CREATE THE  CUTS -----------------------------------------------
  //Use AliESDtrackCuts
  AliESDtrackCuts *trackCuts = new AliESDtrackCuts("AliESDtrackCuts","Standard Cuts");
  //Standard Cuts
  trackCuts=trackCuts->GetStandardITSTPCTrackCuts2009(kTRUE);//Primary Track Selection
  trackCuts->SetEtaRange(-0.9,0.9);
  trackCuts->SetPtRange(0.15, 1e10);
  trackCuts->SetRequireITSRefit(kFALSE);
  
  AliESDtrackCuts *trackCutsITS = new AliESDtrackCuts("AliESDtrackCuts","Standard Cuts with ITSrefit");
  //Standard Cuts
  trackCuts=trackCuts->GetStandardITSTPCTrackCuts2009(kTRUE);//Primary Track Selection
  trackCuts->SetEtaRange(-0.9,0.9);
  trackCuts->SetPtRange(0.15, 1e10);
  trackCuts->SetRequireITSRefit(kTRUE);

  //Create the task
  AliPWG4HighPtQAMC *taskPWG4QAMC = new AliPWG4HighPtQAMC("AliPWG4HighPtQAMC");
  taskPWG4QAMC->SetCuts(trackCuts);
  taskPWG4QAMC->SetCutsITS(trackCutsITS);
  
 
  // E. Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================

  //------ input data ------
  //  AliAnalysisDataContainer *cinput0  = mgr->GetCommonInputContainer();
  printf("Create output containers \n");
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += ":PWG4_HighPtQAMC"; 
  //char *outputfile = "outputAliPWG4HighPtQAMCTestTrain.root";
  AliAnalysisDataContainer *cout_hist0 = mgr->CreateContainer("qa_histsMC", TList::Class(), AliAnalysisManager::kOutputContainer,outputfile);
  AliAnalysisDataContainer *cout_hist2 = mgr->CreateContainer("qa_histsMCITS", TList::Class(), AliAnalysisManager::kOutputContainer,outputfile);  

  mgr->AddTask(taskPWG4QAMC);
  mgr->ConnectInput(taskPWG4QAMC,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskPWG4QAMC,0,cout_hist0);
  mgr->ConnectOutput(taskPWG4QAMC,1,cout_hist2);

  // Return task pointer at the end
  return taskPWG4QAMC;
}
