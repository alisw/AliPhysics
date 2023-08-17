AliAnalysisTaskMultspec_MCpred *AddTaskMultspec_MCpred()
{

    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AliAnalysisTaskMultspec_MCpred", "No analysis manager to connect to.");
        return NULL;
    }

    if (!mgr->GetInputEventHandler()) printf("Caution: No input handler detected!");

    AliAnalysisTaskMultspec_MCpred *taskMCPred = new AliAnalysisTaskMultspec_MCpred("taskMCPred");

    mgr->AddTask(taskMCPred);
    TString outputFileName = AliAnalysisManager::GetCommonFileName();
    outputFileName += ":PWGLF_Multspec_MCpred";
    Printf("Set OutputFileName : \n %s\n", outputFileName.Data() );

    AliAnalysisDataContainer *coutputList = mgr->CreateContainer("fHistos_misc", TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName );

    mgr->ConnectInput (taskMCPred, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(taskMCPred, 1, coutputList);

    return taskMCPred;
}
