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
    
    AliAnalysisTaskDensity* task = new AliAnalysisTaskDensity(name.Data());   
    if(!task) return 0x0;
    task->SelectCollisionCandidates(AliVEvent::kINT7+AliVEvent::kCentral+AliVEvent::kSemiCentral);    
    mgr->AddTask(task);
    
    // input
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    
    // output
    mgr->ConnectOutput(task,1,mgr->CreateContainer("PtSubEventSamples", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task,2,mgr->CreateContainer("QAAliEventCuts", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task,3,mgr->CreateContainer("QATrackCuts", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));

    return task;
}
