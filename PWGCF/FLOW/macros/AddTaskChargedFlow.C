void AddTaskChargedFlow(
	Double_t	fEtaCut 			= 0.8,
	Double_t	fVtxCut				= 10.0,
	Double_t	fMinPt				= 0.2,
	Double_t	fMaxPt				= 3.0,
	Int_t			IsSample			= 10
 ) {

	// ================== GetAnalysisManager ===============================
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr) {
		Error("AddTask_ChargedFlow", "No analysis manager found.");
		return ;
	}

	// ================== GetInputEventHandler =============================
	AliVEventHandler *inputHandler=mgr->GetInputEventHandler();

	AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

	
	//================================================
	//========= Add task to the ANALYSIS manager =====
	//================================================
	//            find input container
	AliAnalysisTaskChargedFlow *task=NULL;
	task= new AliAnalysisTaskChargedFlow("ChargedFlow");
  task->SetEtaCut(fEtaCut);
  task->SetVtxCut(fVtxCut);
  task->SetMinPt(fMinPt);
  task->SetMaxPt(fMaxPt);
  task->SetIsSample(IsSample);

	
	//connect containers
	AliAnalysisDataContainer *coutput =
		mgr->CreateContainer("ChargedFlow", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:ChargedFlow",AliAnalysisManager::GetCommonFileName()));
		
	mgr->AddTask(task);
	mgr->ConnectInput(task,0,cinput);
	mgr->ConnectOutput(task,1,coutput);
	
	return;
	
}
