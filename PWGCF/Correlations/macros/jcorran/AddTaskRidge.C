AliAnalysisTaskRidge* AddTaskRidge(
	const char* taskname = "test",
        const char* option = "LHC16lAODSysZSysTrk",
        bool ismc = kFALSE,
	const char* suffix = "" ){

	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr) return 0x0;
	if (!mgr->GetInputEventHandler())  return 0x0;

	TGrid::Connect("alien://");

	AliAnalysisTaskRidge* taskRidge =
		new AliAnalysisTaskRidge(taskname, Form("%s_%s",taskname,option) );
	if( !taskRidge ) return 0x0;

	mgr->AddTask(taskRidge);

	AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
	AliAnalysisDataContainer *coutputRidge = mgr->CreateContainer(Form("%s_%s",taskname,option),
//		AliDirList::Class(), AliAnalysisManager::kOutputContainer,"AnalysisResults.root");
		AliDirList::Class(), AliAnalysisManager::kOutputContainer,"AnalysisResults.root");

//	TGrid::Connect("alien://");
//	taskRidge->SetEfficiencyFile("/alien/alice/cern.ch/user/j/junlee/Efficiency_RIDGE/EffOut.root?ALICE::CERN::se01");
//	taskRidge->SetEfficiencyFile("/alien/alice/cern.ch/user/j/junlee/Efficiency_RIDGE/EffOut.root");
//	TGrid::Connect("alien:");
//	taskRidge->SetEfficiencyFile("alien:///alice/cern.ch/user/j/junlee/Efficiency_RIDGE/EffOut.root?ALICE::CERN::se01");
//	TGrid::Connect("alien://");
//	taskRidge->fefficiencyFile = (TFile*)TFile::Open("alien:///alice/cern.ch/user/j/junlee/Efficiency_RIDGE/EffOut.root","READ");
	taskRidge->SetEfficiencyFile("alien:///alice/cern.ch/user/j/junlee/Efficiency_RIDGE/EffOut.root");
	taskRidge->SetEfficiency3DFile("alien:///alice/cern.ch/user/j/junlee/Efficiency_RIDGE/Eff3DOut.root");

//	taskRidge->SetEfficiencyFile("./EffOut.root");

	mgr->ConnectInput(taskRidge, 0, cinput);
//	mgr->ConnectOutput(taskRidge,1,mgr->CreateContainer(Form("output%s",suffix), AliDirList::Class(), AliAnalysisManager::kOutputContainer,"AnalysisResults.root"));
	mgr->ConnectOutput(taskRidge,1,mgr->CreateContainer(Form("output%s",suffix), AliDirList::Class(), AliAnalysisManager::kOutputContainer,"AnalysisResults.root"));

	return taskRidge;
}
