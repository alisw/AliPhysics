AliAnalysisTaskCorrForFlow* AddCorrForFlowTask(TString name = "name", const char* suffix = "")
{
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        return 0x0;
    }
    if (!mgr->GetInputEventHandler()) {
        return 0x0;
    }
    TString fileName = AliAnalysisManager::GetCommonFileName();
    fileName += ":CorrForFlow";
    TString taskName = Form("%S_%s",name.Data(),suffix);
    AliAnalysisTaskCorrForFlow* task = new AliAnalysisTaskCorrForFlow(taskName.Data());
    if(!task) return 0x0;
    task->SelectCollisionCandidates(AliVEvent::kAnyINT);
    mgr->AddTask(task);
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task,1,mgr->CreateContainer(Form("Charged_%s",taskname.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    return task;
}
