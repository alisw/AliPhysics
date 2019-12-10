AliAnalysisTaskSigma1385PM* AddTaskSigma1385(
    const char* taskname = "Sigma1385",
    const char* option = "MB_Mix",
    int nmix = 10,
    const char* suffix = "MB") {
    AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        return 0x0;
    }
    if (!mgr->GetInputEventHandler()) {
        return 0x0;
    }
    TString foption = option;
    Bool_t IsMC = kFALSE;
    if (foption.Contains("MC"))
        IsMC = kTRUE;
    AliAnalysisTaskSigma1385PM* taskSigma1385 =
        new AliAnalysisTaskSigma1385PM(Form("%s%s", taskname, suffix), IsMC);
    taskSigma1385->fEventCuts.fCentralityFramework = 1;
    taskSigma1385->fEventCuts.SetMaxVertexZposition(10);
    std::cout << "AddTaskSigma1385:: Option: " << option << std::endl;
    if (foption.Contains("MC")) {
        std::cout << "AliAnalysisTaskSigma1385PM:: MC mode " << std::endl;
        if (foption.Contains("Gen")) {
            taskSigma1385->SetIsPrimaryMC(kFALSE);  // default: kTRUE
            std::cout << "<GENERAL PURPOSE MC>" << std::endl;
        }
    }
    if (foption.Contains("pA")) {
        taskSigma1385->SetSigmaStarRapidityCutLow(0);  // default: -0.5
        std::cout << "AliAnalysisTaskSigma1385PM:: pA mode " << std::endl;
    }
    if (foption.Contains("Ap")) {
        taskSigma1385->SetSigmaStarRapidityCutHigh(0);  // default: 0.5
        std::cout << "AliAnalysisTaskSigma1385PM:: Ap mode " << std::endl;
    }
    if (foption.Contains("Mix")) {
        taskSigma1385->SetMixing(kTRUE);  // default: kFALSE
        std::cout << "AliAnalysisTaskSigma1385PM:: Event Mix(" << nmix
                  << ") mode " << std::endl;
    }
    if (foption.Contains("HM")) {
        taskSigma1385->SetHighMult(kTRUE); // default: kFALSE
        taskSigma1385->fEventCuts.fTriggerMask =
            AliVEvent::kHighMultV0;  // default: kINT7
        std::cout << "AliAnalysisTaskSigma1385PM:: HighMultV0 mode "
                  << std::endl;
    }
    if (foption.Contains("INEL")) {
        taskSigma1385->fEventCuts.fCentralityFramework = 0;
        taskSigma1385->fEventCuts.SelectOnlyInelGt0(false);
        taskSigma1385->SetINEL(kTRUE);  // default: kFALSE
        std::cout << "AliAnalysisTaskSigma1385PM:: Inelastic mode " << std::endl;
    }
    if (foption.Contains("Study")) {
        taskSigma1385->SetCutOpen();
        taskSigma1385->SetFillnTuple(true);
        std::cout << "AliAnalysisTaskSigma1385PM:: Cut Study mode "
                  << std::endl;
    }
    taskSigma1385->SetnMix(nmix);

    if (!taskSigma1385)
        return 0x0;
    mgr->AddTask(taskSigma1385);

    AliAnalysisDataContainer* cinput = mgr->GetCommonInputContainer();
    mgr->ConnectInput(taskSigma1385, 0, cinput);

    AliAnalysisDataContainer* coutputSigma1385 = mgr->CreateContainer(
        Form("%s%s", taskname, suffix), TList::Class(),
        AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
    mgr->ConnectOutput(taskSigma1385, 1, coutputSigma1385);

    AliAnalysisDataContainer* coutputSigma1385Tuple = mgr->CreateContainer(
        Form("%s%s_tree",taskname, suffix), TNtupleD::Class(), AliAnalysisManager::kOutputContainer,
        "AnalysisResults.root");
    coutputSigma1385Tuple->SetSpecialOutput();
    mgr->ConnectOutput(taskSigma1385, 2, coutputSigma1385Tuple);

    return taskSigma1385;
}
