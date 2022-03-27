AliAnalysisTaskXi1530PbPb *AddTaskXi1530PbPb(
    const char *taskname = "Xi1530",
    const char *option = "MB_Mix",
    int nmix = 10,
    const char *suffix = "MB") {
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr)
        return 0x0;
    if (!mgr->GetInputEventHandler())
        return 0x0;

    TString foption = option;
    Bool_t IsMC = kFALSE;
    if (foption.Contains("MC"))
        IsMC = kTRUE;
    AliAnalysisTaskXi1530PbPb *taskXi1530PbPb =
        new AliAnalysisTaskXi1530PbPb(Form("%s%s", taskname, suffix), IsMC);
    taskXi1530PbPb->fEventCuts.fCentralityFramework = 1;
    taskXi1530PbPb->fEventCuts.SetMaxVertexZposition(10);
    taskXi1530PbPb->fEventCuts.SelectOnlyInelGt0(false);
    std::cout << "AddTaskXi1530PbPb:: Option: " << option << std::endl;
    if (foption.Contains("MC")) {
        std::cout << "AliAnalysisTaskXi1530PbPb:: MC mode " << std::endl;
        if (foption.Contains("Gen")) {
            taskXi1530PbPb->SetIsPrimaryMC(kFALSE);  // default: kTRUE
            std::cout << "<GENERAL PURPOSE MC>" << std::endl;
        }
    }
    if (foption.Contains("pA")) {
        taskXi1530PbPb->SetXi1530RapidityCutLow(0);  // default: -0.5
        std::cout << "AliAnalysisTaskXi1530PbPb:: pA mode " << std::endl;
    }
    if (foption.Contains("Ap")) {
        taskXi1530PbPb->SetXi1530RapidityCutHigh(0);  // default: 0.5
        std::cout << "AliAnalysisTaskXi1530PbPb:: Ap mode " << std::endl;
    }
    if (foption.Contains("HM")) {
        taskXi1530PbPb->SetHighMult(kTRUE);  // default: kFALSE
        taskXi1530PbPb->fEventCuts.fTriggerMask =
            AliVEvent::kHighMultV0;  // default: kINT7
        std::cout << "AliAnalysisTaskXi1530PbPb:: HighMultV0 mode " << std::endl;
    }
    if (foption.Contains("Tree")) {
        taskXi1530PbPb->SetFillnTuple();
        std::cout << "AliAnalysisTaskXi1530PbPb:: Fill nTuple" << std::endl;
    }
    if (foption.Contains("Open")) {
        taskXi1530PbPb->SetCutOpen();
        std::cout << "AliAnalysisTaskXi1530PbPb:: Cut open for Tree study" << std::endl;
    }
    // Mixing
    if (foption.Contains("Mix")) {
        taskXi1530PbPb->SetMixing(kTRUE);
        auto tasks = (TObjArray *)mgr->GetTasks();
        auto IsMixer = kFALSE;
        for (int ntask = 0; ntask < tasks->GetEntries(); ntask++) {
            // If we have a TrackMixer task, use it.
            if ((tasks->At(ntask))->InheritsFrom(AliAnalysisTaskTrackMixer::Class())) {
                taskXi1530PbPb->SetMixerTask((AliAnalysisTaskTrackMixer *)tasks->At(ntask));
                ((AliAnalysisTaskTrackMixer *)tasks->At(ntask))->SetnMix(nmix);
                IsMixer = kTRUE;
                break;
            }
        }
        if (!IsMixer) {
            // If we don't have a TrackMixer task, use Built-in mixer.
            taskXi1530PbPb->SetUseBuiltinMixer(kTRUE);
            std::cout << "Built-in Mixer mode" << std::endl;
        }
        taskXi1530PbPb->SetnMix(nmix);
        std::cout << "Event Mix mode: " << nmix << "times" << std::endl;
    }

    if (!taskXi1530PbPb)
        return 0x0;
    mgr->AddTask(taskXi1530PbPb);

    AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
    mgr->ConnectInput(taskXi1530PbPb, 0, cinput);

    AliAnalysisDataContainer *coutputXi1530PbPb = mgr->CreateContainer(
        Form("%s%s", taskname, suffix), TList::Class(),
        AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
    mgr->ConnectOutput(taskXi1530PbPb, 1, coutputXi1530PbPb);

    AliAnalysisDataContainer *coutputXi1530PbPbTree = mgr->CreateContainer(
        Form("%s%s_tree", taskname, suffix), TTree::Class(), AliAnalysisManager::kOutputContainer,
        "AnalysisResults.root");
    coutputXi1530PbPbTree->SetSpecialOutput();
    mgr->ConnectOutput(taskXi1530PbPb, 2, coutputXi1530PbPbTree);

    return taskXi1530PbPb;
}
