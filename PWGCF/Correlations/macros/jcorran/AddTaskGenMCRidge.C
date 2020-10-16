//AliAnalysisTaskGenMCRidge* AddTaskGenMCRidge(
AliAnalysisTask* AddTaskGenMCRidge(
	const char* taskname = "test",
        const char* option = "LHC",
        bool ismc = kTRUE,
        const char* suffix = "" ){

	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if( !mgr ) return NULL;
	if( !mgr->GetMCtruthEventHandler()) return NULL;
//	if (!mgr->GetInputEventHandler())  return 0x0;

	AliAnalysisTaskGenMCRidge* taskRidgeMC = new AliAnalysisTaskGenMCRidge(taskname, Form("%s_%s",taskname,option));
	if( !taskRidgeMC ) return NULL;

	mgr->AddTask( taskRidgeMC );

	AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
        AliAnalysisDataContainer *coutputRidge = mgr->CreateContainer(Form("%s_%s",taskname,option),
                AliDirList::Class(), AliAnalysisManager::kOutputContainer,"AnalysisResults.root");

        mgr->ConnectInput(taskRidgeMC, 0, cinput);
        mgr->ConnectOutput(taskRidgeMC,1,mgr->CreateContainer(Form("output%s%s",suffix,option), AliDirList::Class(), AliAnalysisManager::kOutputContainer,"AnalysisResults.root"));

        return taskRidgeMC;
}
