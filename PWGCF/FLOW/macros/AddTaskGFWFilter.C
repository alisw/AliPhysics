AliGFWFilterTask *AddTaskGFWFilter(const char *name) {
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        return 0x0;
    }
    if (!mgr->GetInputEventHandler()) {
        return 0x0;
    }
    AliGFWFilterTask* task = new AliGFWFilterTask("AliGFWFilter");
    if(!task) return 0x0;
    mgr->AddTask(task);
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task,1,mgr->CreateContainer("GFWFilterQA", TList::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root:GFWTrackFilter"));
    return task;
}
