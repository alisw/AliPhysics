	AliAnalysisPseudoRapidityDensity* AddTaskPseudoRapidityDensity(const char* taskname, const char* option){
		
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        return 0x0;
    }
    if (!mgr->GetInputEventHandler()) {
        return 0x0;
    }

		AliAnalysisPseudoRapidityDensity *task = new AliAnalysisPseudoRapidityDensity(taskname, option);
		AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
		AliAnalysisDataContainer *coutput = mgr->CreateContainer("output", TList::Class(), AliAnalysisManager::kOutputContainer,"AnalysisResults.root");
		mgr->AddTask(task);
		mgr->ConnectInput(task, 0, cinput);
		mgr->ConnectOutput(task, 1, coutput);
		return task;
	}

