AliAnalysisTaskRidge* AddTaskRidge(
	const char* taskname = "test",
        const char* option = "LHC16lAODSysZSysTrk",
        bool ismc = kFALSE,
	const char* suffix = "" ){

	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr) return 0x0;
	if (!mgr->GetInputEventHandler())  return 0x0;

	AliAnalysisTaskRidge* taskRidge =
		new AliAnalysisTaskRidge(taskname, Form("%s_%s",taskname,option) );
	if( !taskRidge ) return 0x0;

	mgr->AddTask(taskRidge);

	AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
	AliAnalysisDataContainer *coutputRidge = mgr->CreateContainer(Form("%s_%s",taskname,option),
		AliDirList::Class(), AliAnalysisManager::kOutputContainer,"AnalysisResults.root");

	mgr->ConnectInput(taskRidge, 0, cinput);
	mgr->ConnectOutput(taskRidge,1,mgr->CreateContainer(Form("output%s%s",suffix,option), AliDirList::Class(), AliAnalysisManager::kOutputContainer,"AnalysisResults.root"));

	return taskRidge;
}
