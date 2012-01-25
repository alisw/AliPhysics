//DEFINITION OF A FEW CONSTANTS
const Double_t maxDeltaTheta = 0.01;
const Double_t ptMin = 5.;

AliPWG4CosmicCandidates* AddTaskPWG4CosmicCandidates(int cuts=1)
{
  // Creates HighPtQATPConly analysis task and adds it to the analysis manager.
  
  // A. Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskPWG4CosmicCandidates", "No analysis manager to connect to.");
    return NULL;
  }  

  // B. Check the analysis type using the event handlers connected to the analysis
  //    manager. The availability of MC handler can also be checked here.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskPWG4CosmicCandidates", "This task requires an input event handler");
    return NULL;
  }  

  TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); 

  // C. Create the task, add it to manager.
  //===========================================================================
 
  //CREATE THE  CUTS -----------------------------------------------
  //Use AliESDtrackCuts
  AliESDtrackCuts *trackCuts = new AliESDtrackCuts("AliESDtrackCuts","Standard Cuts");
  if(cuts==1) trackCuts=trackCuts->GetStandardITSTPCTrackCuts2009(kTRUE);//Primary Track Selection
  trackCuts->SetEtaRange(-0.9,0.9);
  trackCuts->SetPtRange(0.15, 1e10);


  //Create the task
  cout << "Create the task AliPWG4CosmicCandidates" << endl;
  AliPWG4CosmicCandidates *taskPWG4CC = new AliPWG4CosmicCandidates(Form("AliPWG4CosmicCandidates%d",cuts));
  taskPWG4CC->SetCuts(trackCuts);
  taskPWG4CC->SetMaxCosmicAngle(maxDeltaTheta);//0.008);
  taskPWG4CC->SetPtMin(ptMin);
  taskPWG4CC->SelectCollisionCandidates();

  //Add task to manager
  mgr->AddTask(taskPWG4CC);

  // E. Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  //  TString commonFileName = AliAnalysisManager::GetCommonFileName();
  //  commonFileName += ":PWG4_CosmicCandidates"; 
  AliAnalysisDataContainer *cout_hist1 = mgr->CreateContainer(Form("cosmic_hists%d",cuts), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:PWG4_CosmicCandidates%d", AliAnalysisManager::GetCommonFileName(),cuts));
  AliAnalysisDataContainer *cout_cuts0 = mgr->CreateContainer(Form("cosmic_cuts%d",cuts), AliESDtrackCuts::Class(), AliAnalysisManager::kParamContainer, Form("%s:PWG4_CosmicCandidates%d",AliAnalysisManager::GetCommonFileName(),cuts));
  
  
  //Connect input containter to manager
  mgr->ConnectInput(taskPWG4CC,0,mgr->GetCommonInputContainer());

  //Connect output containers to manager
  mgr->ConnectOutput(taskPWG4CC,1,cout_hist1);
  mgr->ConnectOutput(taskPWG4CC,2,cout_cuts0);
  
  // Return task pointer at the end
  return taskPWG4CC;
}
