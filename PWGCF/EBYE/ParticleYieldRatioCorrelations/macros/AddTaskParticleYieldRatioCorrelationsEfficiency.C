AliAnalysisTaskParticleYieldRatioCorrelationsEfficiency *AddTaskParticleYieldRatioCorrelationsEfficiency(
    TString name = "name", int iTask = 0, const char* suffix = "")
{
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr)
    {
        return 0x0;
    }
    if (!mgr->GetInputEventHandler())
    {
        return 0x0;
    }
    TString fileName = AliAnalysisManager::GetCommonFileName();
    fileName += ":output_";
    TString SuffixAdd;
    SuffixAdd.Form("%s", suffix);
    name += "_";
    name += SuffixAdd;
    fileName += name;
    AliAnalysisTaskParticleYieldRatioCorrelationsEfficiency *task = new AliAnalysisTaskParticleYieldRatioCorrelationsEfficiency(name.Data());
    if (!task)
        return 0x0;
    task->SelectCollisionCandidates(AliVEvent::kINT7);

    mgr->AddTask(task);
    mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
    cout << "fileName: " << fileName << endl;
    AliAnalysisDataContainer *coutput_list = mgr->CreateContainer(name, TList::Class(), AliAnalysisManager::kOutputContainer, fileName);
    mgr->ConnectOutput(task, 1, coutput_list);
    return task;
}
