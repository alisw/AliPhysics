AliAnalysisTaskDensity* AddTaskDensity(TString name = "name")
{
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        return 0x0;
    }
    
    if (!mgr->GetInputEventHandler()) {
        return 0x0;
    }
    TString fileName = AliAnalysisManager::GetCommonFileName();
    fileName += ":DensityTask";     
    
    // now we create an instance of your task
    AliAnalysisTaskDensity* task = new AliAnalysisTaskDensity(name.Data());   
    if(!task) return 0x0;
    task->SelectCollisionCandidates(AliVEvent::kINT7);    
    mgr->AddTask(task);
    
    // task needs input
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    
    // output
    mgr->ConnectOutput(task,1,mgr->CreateContainer("PtSubEventSamples", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));

    return task;
}
