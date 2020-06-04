AliAnalysisTaskGenMCRidge* AddTaskGenMCRidge(
	const char* taskname = "test",
        const char* option = "LHC16lAODSysZSysTrk",
        bool ismc = kFALSE,
        const char* suffix = "" ){


	AliAnalysisTaskGenMCRidge* taskRidgeMC = new AliAnalysisTaskGenMCRidge("AliAnalysisTaskGenMCRidge");

	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if( !mgr ) return NULL;

	if(!mgr->GetMCtruthEventHandler()) return NULL;

	mgr->AddTask( taskRidgeMC );

        mgr->ConnectInput(taskRidgeMC, 0, cinput);
        mgr->ConnectOutput(taskRidgeMC,1,mgr->CreateContainer(Form("output%s%s",suffix,option), AliDirList::Class(), AliAnalysisManager::kOutputContainer,"AnalysisResults.root"));

        return taskRidgeMC;
}
