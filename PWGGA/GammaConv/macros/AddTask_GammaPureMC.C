void AddTask_GammaPureMC(Int_t isK0 = 0 ) {

	// ================== GetAnalysisManager ===============================
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr) {
		Error("AddTask_GammaPureMC", "No analysis manager found.");
		return ;
	}

	// ================== GetInputEventHandler =============================
	AliVEventHandler *inputHandler=mgr->GetInputEventHandler();

	AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

	
	//================================================
	//========= Add task to the ANALYSIS manager =====
	//================================================
	//            find input container
	AliAnalysisTaskGammaPureMC *task=NULL;
	task= new AliAnalysisTaskGammaPureMC("GammaPureMC");
    // if no k0 desired, set to isK0
    task->SetIsK0(isK0);

	
	//connect containers
	AliAnalysisDataContainer *coutput =
		mgr->CreateContainer("GammaPureMC", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:GammaPureMC",AliAnalysisManager::GetCommonFileName()));
		
	mgr->AddTask(task);
	mgr->ConnectInput(task,0,cinput);
	mgr->ConnectOutput(task,1,coutput);
	
	return;
	
}
