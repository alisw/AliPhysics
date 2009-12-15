
AliAnalysisTaskSE* AddSDDPoints(Int_t run) {
    gROOT->LoadMacro("AliAnalysisTaskSDDRP.cxx++g");
    
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
    Error("AddTaskESDFilter", "No analysis manager to connect to.");
    return NULL;
    }   

    AliAnalysisTaskSDDRP *task= new AliAnalysisTaskSDDRP();
    task->SetRunNumber(run);
    mgr->AddTask(task);
    mgr->SetDebugLevel(2);
    
    
    AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("coutputRP",TList::Class(),AliAnalysisManager::kOutputContainer,"SDD.Performance.root");
    
    
    
    mgr->ConnectInput(task,  0, cinput1);
    mgr->ConnectOutput(task, 1, coutput1);
    
    return task;
    
}   



