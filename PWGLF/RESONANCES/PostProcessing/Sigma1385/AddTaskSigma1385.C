AliAnalysisSigma1385* AddTaskSigma1385(
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
    AliAnalysisSigma1385* taskSigma1385 =
        new AliAnalysisSigma1385(Form("%s%s", taskname, suffix), IsMC);
    taskSigma1385->fEventCuts.fCentralityFramework = 1;
    taskSigma1385->fEventCuts.SetMaxVertexZposition(10);
    std::cout << "AddTaskSigma1385:: Option: " << option << std::endl;
    if (foption.Contains("MC")) {
        std::cout << "AliAnalysisSigma1385:: MC mode " << std::endl;
        if (foption.Contains("Gen")) {
            taskSigma1385->SetIsPrimaryMC(kFALSE);  // default: kTRUE
            std::cout << "<GENERAL PURPOSE MC>" << std::endl;
        }
    }
    if (foption.Contains("pA")) {
        taskSigma1385->SetSigmaStarRapidityCutLow(0);  // default: -0.5
        std::cout << "AliAnalysisSigma1385:: pA mode " << std::endl;
    }
    if (foption.Contains("Ap")) {
        taskSigma1385->SetSigmaStarRapidityCutHigh(0);  // default: 0.5
        std::cout << "AliAnalysisSigma1385:: Ap mode " << std::endl;
    }
    if (foption.Contains("Mix")) {
        taskSigma1385->SetMixing(kTRUE);  // default: kFALSE
        std::cout << "AliAnalysisSigma1385:: Event Mix(" << nmix
                  << ") mode " << std::endl;
    }
    if (foption.Contains("HM")) {
        taskSigma1385->fEventCuts.fTriggerMask =
            AliVEvent::kHighMultV0;  // default: kINT7
        std::cout << "AliAnalysisSigma1385:: HighMultV0 mode "
                  << std::endl;
    }
    if (foption.Contains("Study")) {
        taskSigma1385->SetCutOpen();
        taskSigma1385->SetFillnTuple(true);
        std::cout << "AliAnalysisSigma1385:: Cut Study mode "
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
        "Sigma1385", TNtupleD::Class(), AliAnalysisManager::kOutputContainer,
        "AnalysisResults.root");
    coutputSigma1385Tuple->SetSpecialOutput();
    mgr->ConnectOutput(taskSigma1385, 2, coutputSigma1385Tuple);

    return taskSigma1385;
}
