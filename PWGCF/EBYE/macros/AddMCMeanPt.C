AliAnalysisMCMeanPt *AddMCMeanPt(
    int trackBit = 768,
    float zvtxcut1 = -10,
    float zvtxcut2 = 10,
    float pt_min = 0.15,
    float pt_max = 1.0,
    int MaxTPCclus = 70.,
    TString period = "def",
    TString suffixName = "a")
{
    // an instance of the class to call the function
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr)
    {
        return 0x0;
    }

    TString TaskMCMeanPt;
    TaskMCMeanPt.Form("Taskfor_%d_%s", trackBit, suffixName.Data());
    AliAnalysisMCMeanPt *task = new AliAnalysisMCMeanPt(TaskMCMeanPt);
    if (!task)
        return 0x0;

    mgr->AddTask(task);

    task->Setzvtxcut(zvtxcut1, zvtxcut2);
    task->SettrackBit(trackBit);
    task->SetMaxTPCCluster(MaxTPCclus);
    task->Setptrange(pt_min, pt_max);


    TString finDirname = suffixName.Data();
    TString fileName = mgr->GetCommonFileName();

    printf("container name: %s\n", fileName.Data());
    printf("file name: %s\n", finDirname.Data());

    mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

    // same for the output
    mgr->ConnectOutput(task, 1, mgr->CreateContainer(Form("fList_%s_MC_FB%d_%s", period.Data(), trackBit, suffixName.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));

    return task;
}
