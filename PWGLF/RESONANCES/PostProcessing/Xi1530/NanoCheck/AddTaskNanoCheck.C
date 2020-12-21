AliAnalysisTaskNanoCheck* AddTaskNanoCheck(const char *taskname = "NanoCheck"
                                     , const char *option = "MB"
                                     , int nmix=10
                                     , const char* suffix = "MB")
{
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        return 0x0;
    }
    if (!mgr->GetInputEventHandler()) {
        return 0x0;
    }
    TString foption = option;
    Bool_t IsMC = kFALSE;
    if(foption.Contains("MC"))
        IsMC = kTRUE;
    AliAnalysisTaskNanoCheck *taskNanoCheck = new AliAnalysisTaskNanoCheck(Form("%s%s", taskname,suffix), IsMC);
    taskNanoCheck->fEventCuts.fCentralityFramework = 1;
    taskNanoCheck->fEventCuts.SetMaxVertexZposition(10);
    
    if(!taskNanoCheck) return 0x0;
    mgr->AddTask(taskNanoCheck);
    
    
    AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
    mgr->ConnectInput(taskNanoCheck, 0, cinput);
    
    AliAnalysisDataContainer *coutputNanoCheck = mgr->CreateContainer(Form("%s%s", taskname, suffix), TList::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
    mgr->ConnectOutput(taskNanoCheck, 1, coutputNanoCheck);

    return taskNanoCheck;
}

