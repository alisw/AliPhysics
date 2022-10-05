AliAnalysisTaskJetQ* AddTaskJetQ(TString name = "FlowTask", TString weightsFile = "", const char* suffix = "")
{
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        return 0x0;
    }
    if (!mgr->GetInputEventHandler()) {
        return 0x0;
    }

    AliAnalysisTaskJetQ* task = new AliAnalysisTaskJetQ(name);
    if(!task) return 0x0;

    mgr->AddTask(task);
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task,1,mgr->CreateContainer("OutputJW", TList::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root"));
    mgr->ConnectOutput(task,2,mgr->CreateContainer("MaxPt", TH1D::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root"));

    return task;
}
