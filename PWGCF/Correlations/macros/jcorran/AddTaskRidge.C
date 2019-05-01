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
		TList::Class(), AliAnalysisManager::kOutputContainer,"AnalysisResults.root");


//	TGrid::Connect("alien://");
//	taskRidge->SetEfficiencyFile("/alien/alice/cern.ch/user/j/junlee/Efficiency_RIDGE/EffOut.root?ALICE::CERN::se01");
//	taskRidge->SetEfficiencyFile("/alien/alice/cern.ch/user/j/junlee/Efficiency_RIDGE/EffOut.root");
	TGrid::Connect("alien://");
	taskRidge->SetEfficiencyFile("alien:///alice/cern.ch/user/j/junlee/Efficiency_RIDGE/EffOut.root?ALICE::CERN::se01");

	mgr->ConnectInput(taskRidge, 0, cinput);
	mgr->ConnectOutput(taskRidge,1,mgr->CreateContainer(Form("output%s",suffix), TList::Class(), AliAnalysisManager::kOutputContainer,"AnalysisResults.root"));

	return taskRidge;


}
