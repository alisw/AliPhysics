AliAnalysisTaskInclusivef0f2* AddTaskInclusivef0f2(
	const char* taskname = "test",
	const char* option = "LHC16lAODMCCSysZSysTrk",
	bool ismc = kFALSE){


	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr) return 0x0;
	if (!mgr->GetInputEventHandler())  return 0x0;

	AliAnalysisTaskInclusivef0f2* taskInclusivef0f2 =
		new AliAnalysisTaskInclusivef0f2(taskname, Form("%s_%s",taskname,option) );
//	taskInclusivef0f2 -> SetIsMC(1);
	if( !taskInclusivef0f2 ) return 0x0;

	mgr->AddTask(taskInclusivef0f2);

	AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
	AliAnalysisDataContainer *coutputInclusivef0f2 = mgr->CreateContainer(Form("%s_%s",taskname,option),
		TList::Class(), AliAnalysisManager::kOutputContainer,"AnalysisResults.root");

	mgr->ConnectInput(taskInclusivef0f2, 0, cinput);
//	mgr->ConnectOutput(taskInclusivef0f2, 1, coutputInclusivef0f2);
	mgr->ConnectOutput(taskInclusivef0f2,1,mgr->CreateContainer("output", TList::Class(), AliAnalysisManager::kOutputContainer,"AnalysisResults.root"));

	return taskInclusivef0f2;
}
