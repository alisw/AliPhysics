AliAnalysisTaskAngularRatiosCorrelation *AddTaskAngularRatiosCorrelation(
    TString name = "name", int iTask = 0, TString InputFileName = "")
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
    AliAnalysisTaskAngularRatiosCorrelation *task = new AliAnalysisTaskAngularRatiosCorrelation(name.Data());
    if (!task)
        return 0x0;

    task->SelectCollisionCandidates(AliVEvent::kINT7);

    mgr->AddTask(task);
    mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

    // in case of non-local run, establish connection to ALiEn for loading the file
    if (InputFileName.Contains("alien://"))
    {
        gGrid->Connect("alien://");
    }
    TFile *InputFile = TFile::Open(InputFileName.Data());
    if (!InputFile)
        printf("Could not open input file\n");
    TList *fList = (TList *)InputFile->Get(Form("List%d", iTask));
    AliAnalysisDataContainer *cinput_list = mgr->CreateContainer("Efficiency", TList::Class(), AliAnalysisManager::kInputContainer);
    cinput_list->SetData(fList);
    mgr->ConnectInput(task, 1, cinput_list);

    AliAnalysisDataContainer *coutput_list = mgr->CreateContainer(name, TList::Class(), AliAnalysisManager::kOutputContainer, fileName);
    mgr->ConnectOutput(task, 1, coutput_list);

    return task;
}
