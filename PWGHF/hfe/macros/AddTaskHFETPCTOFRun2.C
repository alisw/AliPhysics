//AliAnalysisTask AddTaskHFETPCTOFRun2(
AliAnalysisHFETPCTOF AddTaskHFETPCTOFRun2( ///Zaida recommendation

			Bool_t 	isMC 			= kFALSE, 
			Bool_t 	isAOD 			= kFALSE,
			Bool_t 	isPP 			= kFALSE,
			char * period			= "c"
)
{
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	
	if (!mgr) {
	::Error("AddTaskHFETPCTOFRun2", "No analysis manager to connect to.");
	return NULL;
	}
	
	if (!mgr->GetInputEventHandler()) {
	::Error("AddTaskHFETPCTOFRun2", "This task requires an input event handler");
	return NULL;
	}
	
	//_______________________
	//Config Task
	gROOT->LoadMacro("ConfigHFETPCTOF.C");
	AliAnalysisHFETPCTOF *task = ConfigHFETPCTOF(isMC,isAOD,isPP);
	
	//_______________________
	//Trigger
	if(!isMC)
	{
		task->SelectCollisionCandidates(AliVEvent::kINT7); //Selecting Minumum Bias events (selected randomlly)
		//task->SelectCollisionCandidates(AliVEvent::kAny);
    }
	
	mgr->AddTask(task);
	
	//added to run on grid
    //==================================================================================
    TString containerName = mgr->GetCommonFileName();
    containerName += ":HFE_Camila";
    //containerName += Form("_%d_%d_%d_%d",triggerIndex,configIndex,centralityIndex,EMCalThreshould);
    
    
    //Create containers for input/output
    AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput = mgr->CreateContainer("hist", TList::Class(),    AliAnalysisManager::kOutputContainer, containerName.Data());
    // end of modifications to run on grid
    //==================================================================================
	
	
	//Create containers for input/output
	//AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
	//AliAnalysisDataContainer *coutput = mgr->CreateContainer("hist",  TList::Class(),    AliAnalysisManager::kOutputContainer, "OutPutDataMC.root");

	//Connect input/output
	mgr->ConnectInput(task, 0, cinput);
	mgr->ConnectOutput(task, 1, coutput);
	
	return task;
}
