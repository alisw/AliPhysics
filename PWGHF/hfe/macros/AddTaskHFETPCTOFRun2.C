AliAnalysisTask AddTaskHFETPCTOFRun2(

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
	//gROOT->LoadMacro("ConfigHFETPCTOF.C");
	gROOT->LoadMacro("$ALICE_PHYSICS/PWGHF/hfe/macros/configs/pp/ConfigHFETPCTOF.C");
	AliAnalysisHFETPCTOF *task = ConfigHFETPCTOF(isMC,isAOD,isPP);
	
	//_______________________
	//Trigger
	if(!isMC)
	{
		task->SelectCollisionCandidates(AliVEvent::kINT7); //Selecting Minumum Bias events (selected randomlly)
		//task->SelectCollisionCandidates(AliVEvent::kAny);
    }
	
	mgr->AddTask(task);
	
	//Create containers for input/output
	AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
	AliAnalysisDataContainer *coutput = mgr->CreateContainer("hist",  TList::Class(),    AliAnalysisManager::kOutputContainer, "OutPutDataMC.root");

	//Connect input/output
	mgr->ConnectInput(task, 0, cinput);
	mgr->ConnectOutput(task, 1, coutput);
	
	return task;
}
