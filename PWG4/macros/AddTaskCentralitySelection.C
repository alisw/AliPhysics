AliCentralitySelectionTask *AddTaskCentralitySelection(const char* percentilefile1, const char* percentilefile2){


    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if(!mgr){
	::Error("AddTaskCentralitySelection", "No analysis manager to connect ot.");
	return NULL;
    }
    if(!mgr->GetInputEventHandler()){
        ::Error("AddTaskCentralitySelection", "This task requires an input event handler.");
	return NULL;
    }


    AliCentralitySelectionTask *task = new AliCentralitySelectionTask("CentralitySelection");

    if(percentilefile1) task->SetPercentileFile(percentilefile1);
    if(percentilefile2) task->SetPercentileFile2(percentilefile2);

    mgr->AddTask(task);

    mgr->ConnectInput (task, 0, mgr->GetCommonInputContainer());

    return task;
}
