AliAnalysisTaskCheckAODdAODMatching* AddTaskCheckAODdAODMatching()
{
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr)
        ::Error("AddTaskDs", "No analysis manager to connect to.");

    AliAnalysisTaskCheckAODdAODMatching* taskMathing = new AliAnalysisTaskCheckAODdAODMatching("taskCheckAODmatching");
    mgr->AddTask(taskMathing);

    // Create containers for input/output
    AliAnalysisDataContainer *cinput = mgr->CreateContainer("cinputAODmatching", TChain::Class(), AliAnalysisManager::kInputContainer);

    TString outputfile = AliAnalysisManager::GetCommonFileName();
    outputfile += ":AOD_dAOD_Matching";
    AliAnalysisDataContainer *coutput = mgr->CreateContainer("coutputAODmatching", TList::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
    AliAnalysisDataContainer *coutputTree = mgr->CreateContainer("coutputAODtreeMismatch", TTree::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());

    mgr->ConnectInput(taskMathing, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(taskMathing, 1, coutput);
    mgr->ConnectOutput(taskMathing, 2, coutputTree);

    return taskMathing;
}
