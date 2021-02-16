AliAnalysisTaskFluctMCTEPOS * AddTaskMC(const char * outfilename, bool  isMonteCarlo) {

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPhysicsSelection", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
    if (!mgr->GetMCtruthEventHandler()) {
    ::Error("AddTaskPhysicsSelection", "This task requires an input event handler");
    return NULL;
  }
	
    //TString inputDataType = mgr->GetMCtruthEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  
  // if (inputDataType != "AOD") {
  //   Printf("ERROR! This task can only run on AODs!");
  // }

  // Configure analysis
  //===========================================================================
  
    // TString lAnalysisCut      = "no";    
    //Int_t iMCAnalysis = 0;
   
   
  Bool_t isMC = 0;


   
  char taskName[15];
  sprintf(taskName,"example_task");
  
  AliAnalysisTaskFluctMCTEPOS *task = new AliAnalysisTaskFluctMCTEPOS(taskName);
  //task->SelectCollisionCandidates(AliVEvent::kINT7); // if physics selection performed in UserExec(), this line should be commented
  
   task-> SetIsMonteCarlo		(isMonteCarlo);
   // task-> SetGenerator("EPOS-LHC");  
    
  mgr->AddTask(task);
  AliAnalysisDataContainer *cinput0 = mgr->GetCommonInputContainer();
  
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(outfilename, TList::Class(),AliAnalysisManager::kOutputContainer, Form("%s:lambdak0Luke", AliAnalysisManager::GetCommonFileName()));
  
  mgr->ConnectInput (task, 0, cinput0);
  mgr->ConnectOutput(task,1,coutput1);
  
  return task;
}   

