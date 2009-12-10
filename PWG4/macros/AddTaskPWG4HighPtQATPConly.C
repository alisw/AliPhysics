//DEFINITION OF A FEW CONSTANTS

AliPWG4HighPtQATPConly* AddTaskPWG4HighPtQATPConly()//<some_parameters>)
{
  // Creates a proton analysis task and adds it to the analysis manager.
  
  // A. Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskPWG4HighPtSpectra", "No analysis manager to connect to.");
    return NULL;
  }  

  // B. Check the analysis type using the event handlers connected to the analysis
  //    manager. The availability of MC handler can also be checked here.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskPWG4HighPtTPConly", "This task requires an input event handler");
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
  trackCuts->SetAcceptKinkDaughters(kFALSE);//
  trackCuts->SetRequireTPCRefit(kTRUE);
  trackCuts->SetEtaRange(-0.9,0.9);//-0.5,0.5);//
  trackCuts->SetMaxCovDiagonalElements(2,2,0.5,0.5,2);//
  trackCuts->SetPtRange(0.15, 1e10);//
  trackCuts->SetMinNClustersTPC(50);//
  trackCuts->SetMaxChi2PerClusterTPC(3.5);//
  //trackCuts->SetRequireITSRefit(kTRUE);
  trackCuts->SetMaxDCAToVertexXY(2.4);
  trackCuts->SetMaxDCAToVertexZ(3.2);
  trackCuts->SetDCAToVertex2D(kTRUE);
  
  AliESDtrackCuts *trackCutsITS = new AliESDtrackCuts("AliESDtrackCuts","Standard Cuts with ITSrefit");
  trackCutsITS->SetAcceptKinkDaughters(kFALSE);//
  trackCutsITS->SetRequireTPCRefit(kTRUE);//
  trackCutsITS->SetEtaRange(-0.9,0.9);//-0.5,0.5);//
  trackCutsITS->SetMaxCovDiagonalElements(2,2,0.5,0.5,2);//
  trackCutsITS->SetPtRange(0.15, 1e10);//
  trackCutsITS->SetMinNClustersTPC(50);//
  trackCutsITS->SetMaxChi2PerClusterTPC(3.5);//
  trackCutsITS->SetRequireITSRefit(kTRUE);
  trackCutsITS->SetMaxDCAToVertexXY(2.4);
  trackCutsITS->SetMaxDCAToVertexZ(3.2);
  trackCutsITS->SetDCAToVertex2D(kTRUE); 

  //Create the task
  AliPWG4HighPtQATPConly *taskPWG4QA = new AliPWG4HighPtQATPConly("AliPWG4HighPtQATPConly");
  taskPWG4QA->SetCuts(trackCuts);
  taskPWG4QA->SetCutsITS(trackCutsITS);
  taskPWG4QA->SelectTrigger(AliAnalysisHelperJetTasks::kMB1); 
  
 
  // E. Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================

  //------ input data ------
  //  AliAnalysisDataContainer *cinput0  = mgr->GetCommonInputContainer();
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += ":PWG4_HighPtQATPConly"; 
  //char *outputfile = "outputAliPWG4HighPtQATPConlyTestTrain.root";
  AliAnalysisDataContainer *cout_hist0 = mgr->CreateContainer("qa_hists", TList::Class(), AliAnalysisManager::kOutputContainer,outputfile);
  AliAnalysisDataContainer *cout_hist1 = mgr->CreateContainer("qa_histsTPC", TList::Class(), AliAnalysisManager::kOutputContainer,outputfile);
  AliAnalysisDataContainer *cout_hist2 = mgr->CreateContainer("qa_histsITS", TList::Class(), AliAnalysisManager::kOutputContainer,outputfile);  

 
  mgr->AddTask(taskPWG4QA);

  mgr->ConnectInput(taskPWG4QA,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskPWG4QA,0,cout_hist0);
  mgr->ConnectOutput(taskPWG4QA,1,cout_hist1);
  mgr->ConnectOutput(taskPWG4QA,2,cout_hist2);

  // Return task pointer at the end
  return taskPWG4QA;
}
