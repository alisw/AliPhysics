AliAnalysisTaskK1 *AddTaskK1(const char *taskname = "K1", const char *option = "MB_Mix", int nmix = 10, const char *suffix = "MB")
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
    TString foption = option;
    Bool_t IsMC = kFALSE;
    if (foption.Contains("MC"))
        IsMC = kTRUE;
    AliAnalysisTaskK1 *taskK1 = new AliAnalysisTaskK1(Form("%s%s", taskname, suffix), IsMC);
    taskK1->fEventCut.fCentralityFramework = 1;
    taskK1->fEventCut.SetMaxVertexZposition(10);
    taskK1->fEventCut.SelectOnlyInelGt0(false);
    std::cout << Form("AddTaskK1:: Option: %s", option) << std::endl;
    if (foption.Contains("MC"))
    {
        std::cout << "AliAnalysisTaskK1:: MC mode " << std::endl;
        if (foption.Contains("Gen"))
        {
            taskK1->SetIsPrimaryMC(kFALSE); // default: kTRUE
            std::cout << "<GENERAL PURPOSE MC>" << std::endl;
        }
    }
    if (foption.Contains("pA"))
    {
        taskK1->SetK1RapidityCutLow(0); // default: -0.5
        std::cout << "AliAnalysisTaskK1:: pA mode" << std::endl;
    }
    if (foption.Contains("Ap"))
    {
        taskK1->SetK1RapidityCutHigh(0); // default: 0.5
        std::cout << "AliAnalysisTaskK1:: Ap mode" << std::endl;
    }
    if (foption.Contains("HM"))
    {
        taskK1->SetHighMult(kTRUE);                              // default: kFALSE
        taskK1->fEventCut.fTriggerMask = AliVEvent::kHighMultV0; // default: kINT7
        std::cout << "AliAnalysisTaskK1:: HighMultV0 mode " << std::endl;
    }
    // Mixing
    if (foption.Contains("Mix"))
    {
        taskK1->SetnMix(nmix);
        taskK1->SetMixing(kTRUE);
        std::cout << Form("Event Mix mode: %d", nmix) << std::endl;
    }
    // Skip Filling Histograms
    if (foption.Contains("NoHist"))
    {
        taskK1->SetSkipFillHistos(true);
        std::cout << "Skip filling histograms (for Tree study)" << std::endl;
    }
    // Fill Tree
    if (foption.Contains("Tree"))
    {
        taskK1->SetFillTree(true);
        std::cout << "Fill Tree" << std::endl;
    }
    if (!taskK1)
        return 0x0;

    mgr->AddTask(taskK1);

    TString output = "AnalysisResults.root";
    AliAnalysisDataContainer *coutputK1 = mgr->CreateContainer(Form("%s%s", taskname, suffix), TList::Class(), AliAnalysisManager::kOutputContainer, output.Data());
    AliAnalysisDataContainer *coutputK1RTree = mgr->CreateContainer(Form("%s%s_NanoTree", taskname, suffix), TTree::Class(), AliAnalysisManager::kOutputContainer, output.Data());

    mgr->ConnectInput(taskK1, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(taskK1, 1, coutputK1);
    mgr->ConnectOutput(taskK1, 2, coutputK1RTree);

    return taskK1;
}
