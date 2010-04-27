
AliAnalysisTaskSE* AddSDDPoints() {
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
    Error("AddTaskESDFilter", "No analysis manager to connect to.");
    return NULL;
    }   

    AliAnalysisTaskSDDRP *task= new AliAnalysisTaskSDDRP();
    mgr->AddTask(task);
    mgr->SetDebugLevel(2);
    
    
    AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("coutputRP",TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:SDD_Performance", mgr->GetCommonFileName()));
    
    
    
    mgr->ConnectInput(task,  0, cinput1);
    mgr->ConnectOutput(task, 1, coutput1);
    
    return task;
    
}   



