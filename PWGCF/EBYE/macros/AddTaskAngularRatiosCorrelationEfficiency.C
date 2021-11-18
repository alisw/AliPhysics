AliAnalysisTaskAngularRatiosCorrelationEfficiency *AddTaskAngularRatiosCorrelationEfficiency(
    TString name = "name",
    int iTask = 0,
    Bool_t MCGen = kFALSE,
    Bool_t pbpb = kTRUE,
    int nPhiBins = 16,
    int nVertexBins = 1,
    int nPBins = 5,
    int nTPCcrossedRows = 70,
    int nSigma = 2,
    int nCentrClasses = 4,
    int nEtaClasses = 16,
    int nSorts = 8,
    int minCent = 0,
    int maxCent = 80,
    UInt_t filterbit = 96,
    Float_t minP = 0.2,
    Float_t maxP = 2.0,
    Float_t Vertexmin = -8,
    Float_t Vertexmax = 8)
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
    fileName += name;
    AliAnalysisTaskAngularRatiosCorrelationEfficiency *task = new AliAnalysisTaskAngularRatiosCorrelationEfficiency(name.Data());
    if (!task)
        return 0x0;
    task->SelectCollisionCandidates(AliVEvent::kINT7);
    task->SetFilterBit(filterbit);
    task->SetParams(iTask, nPhiBins, nVertexBins, nPBins, minCent, maxCent,
                    minP, maxP, Vertexmin, Vertexmax, MCGen, pbpb,
                    nSigma, nTPCcrossedRows);

    mgr->AddTask(task);
    mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
    cout << "fileName: " << fileName << endl;
    AliAnalysisDataContainer *coutput_list = mgr->CreateContainer(name, TList::Class(), AliAnalysisManager::kOutputContainer, fileName);
    mgr->ConnectOutput(task, 1, coutput_list);
    return task;
}
