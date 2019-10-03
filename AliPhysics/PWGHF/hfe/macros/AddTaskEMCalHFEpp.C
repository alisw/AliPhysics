AliAnalysisTaskEMCalHFEpA *AddTaskEMCalHFEpp(

			Bool_t 	isMC 			= kFALSE, 
			Int_t 	triggerIndex 	= 0, 
			Int_t 	configIndex 	= 0, 
			Int_t 	centralityIndex = 0, 
			Bool_t 	isAOD 			= kFALSE,
			Bool_t isEMCal 			= kFALSE,
			Bool_t isTrigger 		= kFALSE,
			char * period			= "b",
			Int_t EMCalThreshould 	= 0,
			Bool_t isTender = kFALSE
		)
{
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	
	if (!mgr) {
	::Error("AddTaskEMCalHFEpp", "No analysis manager to connect to.");
	return NULL;
	}
	
	if (!mgr->GetInputEventHandler()) {
	::Error("AddTaskEMCalHFEpp", "This task requires an input event handler");
	return NULL;
	}
	
	//_______________________
	//Config Task
	//gROOT->LoadMacro("ConfigEMCalHFEpp.C");
	gROOT->LoadMacro("$ALICE_PHYSICS/PWGHF/hfe/macros/configs/pp/ConfigEMCalHFEpp.C");
	AliAnalysisTaskEMCalHFEpA *task = ConfigEMCalHFEpp(isMC,triggerIndex,configIndex,centralityIndex,isAOD,isEMCal,isTrigger, EMCalThreshould, isTender);
	
	//_______________________
	//Trigger
	if(!isMC)
	{
		if(triggerIndex == 0) task->SelectCollisionCandidates(AliVEvent::kINT7);
		if(triggerIndex == 1) task->SelectCollisionCandidates(AliVEvent::kEMC7);
		if(triggerIndex == 2) task->SelectCollisionCandidates(AliVEvent::kEMCEGA);
		if(triggerIndex == 3) task->SelectCollisionCandidates(AliVEvent::kEMCEJE); //Jet Trigger
		//if(triggerIndex == 4) task->SelectCollisionCandidates(AliVEvent::kEMC8);
		
	}
	
	mgr->AddTask(task);

	TString containerName = mgr->GetCommonFileName();
        containerName += ":HFE_EMCal_pp";
        containerName += Form("_%d_%d_%d_%d",triggerIndex,configIndex,centralityIndex,EMCalThreshould);
	
	//Create containers for input/output
	AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
	AliAnalysisDataContainer *coutput = mgr->CreateContainer(Form("chist_%d_%d_%d_%d",triggerIndex,configIndex,centralityIndex, EMCalThreshould), TList::Class(),    AliAnalysisManager::kOutputContainer, containerName.Data());

	//Connect input/output
	mgr->ConnectInput(task, 0, cinput);
	mgr->ConnectOutput(task, 1, coutput);
	
	return task;
}
