AliAnalysisTaskCheckAODTracks *AddTaskCheckAODTracks(TString suffix="",
						     Bool_t readMC=kFALSE,
						     Bool_t useMCtruthForPID=kFALSE,
						     Double_t minCent=-1.,
						     Double_t maxCent=999.,
						     TString estim="V0M")
{

  // Creates, configures and attaches to the train the task for tracking checks
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskCheckAODTracks", "No analysis manager to connect to.");
    return NULL;
  }   
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskCheckAODTracks", "This task requires an input event handler");
    return NULL;
  }   
  
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if(type.Contains("ESD")){
    ::Error("AddTaskCheckAODTracks", "This task requires to run on AOD");
    return NULL;
  }
  
  //Bool_t isMC=kFALSE;
  //if (mgr->GetMCtruthEventHandler()) isMC=kTRUE;
  
  // Add MC handler (for kinematics)
  if(readMC){
    AliMCEventHandler* handler = (AliMCEventHandler*)mgr->GetMCtruthEventHandler();
    if (!handler) {
      ::Error("AddTaskCheckAODTracks","Macro called with readMC=true but MC handler not present");
      return 0;
    }
  }
  // Create and configure the task
  AliAnalysisTaskCheckAODTracks *tasktr = new AliAnalysisTaskCheckAODTracks();
  tasktr->SetReadMC(readMC);
  tasktr->SetUseMCtruthForPID(useMCtruthForPID);
  if(minCent>=0 && maxCent<=100) tasktr->SetCentralityInterval(minCent,maxCent,estim);
  mgr->AddTask(tasktr);
  
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += ":CheckAODTracks";

  TString listname="clistCheckAODTracks";
  listname+=suffix.Data();
  
  TString treeFileName="TreeOfAODTracks.root";
  treeFileName += ":CheckAODTracksTree";
  treeFileName += suffix.Data();

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(listname,
							    TList::Class(),
							    AliAnalysisManager::kOutputContainer,
							    outputFileName.Data() );
  
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(Form("ctreeAODTracks%s",suffix.Data()),
							    TTree::Class(),
							    AliAnalysisManager::kOutputContainer,
							    treeFileName.Data() );
  coutput2->SetSpecialOutput();

  mgr->ConnectInput(tasktr, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(tasktr, 1, coutput1);
  mgr->ConnectOutput(tasktr, 2, coutput2);
  return tasktr;
}   

