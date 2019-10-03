AliAnalysisTaskLukeAOD * AddTaskLukeAOD(const char * outfilename, bool	isMonteCarlo, double cutCosPa, double cutcTauMin, double cutNcTauMax, double cutBetheBloch, double	cutMinNClustersTPC, double cutRatio, double	cutEta, double cutRapidity, double cutArmenteros) {


  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPhysicsSelection", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskPhysicsSelection", "This task requires an input event handler");
    return NULL;
  }
	
  TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  
  if (inputDataType != "AOD") {
    Printf("ERROR! This task can only run on AODs!");
  }

  // Configure analysis
  //===========================================================================
  // Int_t nbMinTPCclusters = 80;
  // Int_t lCollidingSystems = 1; 
  // TString fAnalysisType = "ESD";
    //    TString lAnalysisPidMode  = "withPID";
    // TString lAnalysisCut      = "no";    
    //Int_t iMCAnalysis = 0;
     
 	char taskName[15];
	sprintf(taskName,"example_task");
	
	AliAnalysisTaskLukeAOD * task = new AliAnalysisTaskLukeAOD(taskName);
	task->SelectCollisionCandidates(AliVEvent::kMB); // if physics selection performed in UserExec(), this line should be commented
	
	task->SetIsMonteCarlo		(isMonteCarlo);
	task->SetCutCosPa			(cutCosPa);
	task->SetCutcTauMin			(cutcTauMin);
	task->SetCutNcTauMax		(cutNcTauMax);
	task->SetCutBetheBloch		(cutBetheBloch);
	task->SetCutMinNClustersTPC	(cutMinNClustersTPC);
	task->SetCutRatio			(cutRatio)	;
	task->SetCutEta				(cutEta);
	task->SetCutRapidity		(cutRapidity);
	task->SetCutArmenteros		(cutArmenteros);
	
	
	//task->SetOption("sample");
   
    mgr->AddTask(task);
	AliAnalysisDataContainer *cinput0 = mgr->GetCommonInputContainer();
   
	AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(outfilename, TList::Class(),AliAnalysisManager::kOutputContainer, Form("%s:lambdak0Luke", AliAnalysisManager::GetCommonFileName()));
    
	mgr->ConnectInput (task, 0, cinput0);
    mgr->ConnectOutput(task,1,coutput1);
    
  return task;
}   


