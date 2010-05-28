//DEFINITION OF A FEW CONSTANTS

AliPWG4HighPtQATPConly* AddTaskPWG4HighPtQATPConly(int cuts=1)//1: Standard Cuts 2009 2: GetStandardITSTPCTrackCuts2009
{
  // Creates HighPtQATPConly analysis task and adds it to the analysis manager.
  
  // A. Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskPWG4HighPtQATPConly", "No analysis manager to connect to.");
    return NULL;
  }  

  // B. Check the analysis type using the event handlers connected to the analysis
  //    manager. The availability of MC handler can also be checked here.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskPWG4HighPtQATPConly", "This task requires an input event handler");
    return NULL;
  }  
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  const char *analysisType = "ESD";//"TPC"

  // C. Create the task, add it to manager.
  //===========================================================================
 
  //CREATE THE  CUTS -----------------------------------------------
  //Use AliESDtrackCuts
  AliESDtrackCuts *trackCuts = new AliESDtrackCuts("AliESDtrackCuts","Standard Cuts");
  if(cuts==1) {
    //Standard Cuts
    trackCuts->SetAcceptKinkDaughters(kFALSE);
    trackCuts->SetRequireTPCStandAlone(kTRUE); 
    trackCuts->SetRequireTPCRefit(kTRUE);
    trackCuts->SetMinNClustersTPC(70);
    trackCuts->SetEtaRange(-0.9,0.9);
    trackCuts->SetMaxCovDiagonalElements(2,2,0.5,0.5,2);
    trackCuts->SetPtRange(0.15, 1e10);
    trackCuts->SetMaxChi2PerClusterTPC(3.5);
    trackCuts->SetMaxDCAToVertexXY(2.4);
    trackCuts->SetMaxDCAToVertexZ(3.2);
    trackCuts->SetDCAToVertex2D(kTRUE);
  }
  else if(cuts==2) {
    trackCuts=trackCuts->GetStandardITSTPCTrackCuts2009(kTRUE);//Primary Track Selection
    trackCuts->SetEtaRange(-0.9,0.9);
    trackCuts->SetPtRange(0.15, 1e10);
    trackCuts->SetRequireITSRefit(kFALSE);
  }

  AliESDtrackCuts *trackCutsITS = new AliESDtrackCuts("AliESDtrackCuts","Standard Cuts with ITSrefit");
  if(cuts==1) {
    trackCutsITS->SetAcceptKinkDaughters(kFALSE);
    trackCutsITS->SetRequireTPCRefit(kTRUE);
    trackCutsITS->SetEtaRange(-0.9,0.9);
    trackCutsITS->SetMaxCovDiagonalElements(2,2,0.5,0.5,2);
    trackCutsITS->SetPtRange(0.15, 1e10);
    trackCutsITS->SetMinNClustersTPC(70);
    trackCutsITS->SetMaxChi2PerClusterTPC(3.5);
    trackCutsITS->SetRequireITSRefit(kTRUE);
    trackCutsITS->SetMaxDCAToVertexXY(2.4);
    trackCutsITS->SetMaxDCAToVertexZ(3.2);
    trackCutsITS->SetDCAToVertex2D(kTRUE); 
  }
 else if(cuts==2) {
   trackCutsITS=trackCutsITS->GetStandardITSTPCTrackCuts2009(kTRUE);//Primary Track Selection
   trackCutsITS->SetEtaRange(-0.9,0.9);
   trackCutsITS->SetPtRange(0.15, 1e10);
 }

//Create the task
  AliPWG4HighPtQATPConly *taskPWG4QA = new AliPWG4HighPtQATPConly(Form("AliPWG4HighPtQATPConly%d",cuts));
  taskPWG4QA->SetCuts(trackCuts);
  taskPWG4QA->SetCutsITS(trackCutsITS);
  taskPWG4QA->SetMaxCosmicAngle(0.008);
  taskPWG4QA->SetCutType(cuts);

 
  // E. Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================

  //------ input data ------
  TString outputfile = "";
  outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += Form(":PWG4_HighPtQATPConly%d",cuts);
  
  AliAnalysisDataContainer *cout_hist0;
  AliAnalysisDataContainer *cout_hist1;
  AliAnalysisDataContainer *cout_hist2;
  AliAnalysisDataContainer *cout_hist3;
  AliAnalysisDataContainer *cout_cuts0;
  AliAnalysisDataContainer *cout_cuts1;

  cout_hist0 = mgr->CreateContainer(Form("qa_histsCuts%d",cuts), TList::Class(), AliAnalysisManager::kOutputContainer,outputfile);
  cout_hist1 = mgr->CreateContainer(Form("qa_histsTPCCuts%d",cuts), TList::Class(), AliAnalysisManager::kOutputContainer,outputfile);
  cout_hist2 = mgr->CreateContainer(Form("qa_histsITSCuts%d",cuts), TList::Class(), AliAnalysisManager::kOutputContainer,outputfile);  
  cout_hist3 = mgr->CreateContainer(Form("qa_histsCosmicsCuts%d",cuts), TList::Class(), AliAnalysisManager::kOutputContainer,outputfile);  
  cout_cuts0 = mgr->CreateContainer(Form("qa_trackCuts%d",cuts), AliESDtrackCuts::Class(), AliAnalysisManager::kParamContainer,outputfile);
  cout_cuts1 = mgr->CreateContainer(Form("qa_trackCutsITS%d",cuts), AliESDtrackCuts::Class(), AliAnalysisManager::kParamContainer,outputfile);

  //Add task to manager
  mgr->AddTask(taskPWG4QA);

  //Connect input containter to manager
  mgr->ConnectInput(taskPWG4QA,0,mgr->GetCommonInputContainer());

  //Connect output containers to manager
  mgr->ConnectOutput(taskPWG4QA,0,cout_hist0);
  mgr->ConnectOutput(taskPWG4QA,1,cout_hist1);
  mgr->ConnectOutput(taskPWG4QA,2,cout_hist2);
  mgr->ConnectOutput(taskPWG4QA,3,cout_hist3);
  mgr->ConnectOutput(taskPWG4QA,4,cout_cuts0);
  mgr->ConnectOutput(taskPWG4QA,5,cout_cuts1);

  // Return task pointer at the end
  return taskPWG4QA;
}
