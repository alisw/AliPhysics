AliAnalysisTaskEmeanCalib* AddTaskEmeanCalib(Bool_t usePhysicsSelection) {
    
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddTaskEmeanCalib", "No analysis manager to connect to.");
        return NULL;
    }
    
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(
                                                              "clist",
                                                              TList::Class(),
                                                              AliAnalysisManager::kOutputContainer,
                                                              mgr->GetCommonFileName() );
    
    AliAnalysisTaskEmeanCalib *task = new AliAnalysisTaskEmeanCalib("AliAnalysisTaskEmeanCalib");
    mgr->AddTask(task);
    
    mgr->ConnectInput (task, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task, 1, coutput1);
    
    return task;
}
