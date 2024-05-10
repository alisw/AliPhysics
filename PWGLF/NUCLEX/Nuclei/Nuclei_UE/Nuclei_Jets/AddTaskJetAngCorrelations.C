AliAnalysisTaskJetAngCorrelations* AddTaskJetAngCorrelations(TString name = "JetAngCorrTask")
{
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        return 0x0;
    }
    if (!mgr->GetInputEventHandler()) {
        return 0x0;
    }
    
    TString fileName = AliAnalysisManager::GetCommonFileName();
    fileName += ":Jets";

    AliAnalysisTaskJetAngCorrelations* task = new AliAnalysisTaskJetAngCorrelations(name.Data());
    if(!task) return 0x0;
    task->SelectCollisionCandidates(AliVEvent::kINT7);
    
    mgr->AddTask(task);
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task,1,mgr->CreateContainer("Jets", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task,2,mgr->CreateContainer("QAPlots", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));

    return task;
}