AliAnalysisTaskSENonPromptLc *AddTaskNonPromptLc(int decCh = AliAnalysisTaskSENonPromptLc::kLctopKpi,
                                                 bool readMC = false,
                                                 TString fileName = "",
                                                 TString suffix = "",
                                                 bool createMLtree = true, // for the moment only possibility true
                                                 TString cutObjName = "AnalysisCuts",
                                                 int AODProtection = 0)
{
    // \brief: AddTask for AliAnalysisTaskSENonPromptLc
    // \authors:
    // F. Grosa, fabrizio.grosa@cern.ch
    //==============================================================================

    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr)
        ::Error("AddTaskNonPromptLc", "No analysis manager to connect to.");

    TFile *fileCuts = TFile::Open(fileName.Data());
    if (!fileCuts || (fileCuts && !fileCuts->IsOpen()))
        ::Fatal("AddTaskNonPromptLc", "Cut file not found on Grid: analysis will not start!\n");
    else
        printf("Cut file correctly found\n");

    AliRDHFCuts *analysisCuts = NULL;
    if(decCh == AliAnalysisTaskSENonPromptLc::kLctopKpi)
    {
        analysisCuts = (AliRDHFCutsLctopKpi*)fileCuts->Get(cutObjName.Data());
        suffix.Prepend("LctopKpi");
    }
    else if(decCh == AliAnalysisTaskSENonPromptLc::kLctopK0s)
    {
        analysisCuts = (AliRDHFCutsLctoV0*)fileCuts->Get(cutObjName.Data());
        suffix.Prepend("LctopK0s");
    }
    else if(decCh == AliAnalysisTaskSENonPromptLc::kLctopiL)
    {
        analysisCuts = (AliRDHFCutsLctoV0*)fileCuts->Get(cutObjName.Data());
        suffix.Prepend("LctopiL");
    }
    else
    {
        ::Fatal("AddTaskNonPromptLc", "Specific AliRDHFCuts not found!\n");
        return NULL;
    }

    // Analysis Task
    AliAnalysisTaskSENonPromptLc *npLcTask = new AliAnalysisTaskSENonPromptLc("NonPromptLcAnalysis", decCh, analysisCuts, createMLtree);
    npLcTask->SetAODMismatchProtection(AODProtection);
    npLcTask->SetReadMC(readMC);
    if (createMLtree && readMC)
        npLcTask->SetFillOnlySignalInMLtree();
    mgr->AddTask(npLcTask);

    // Create containers for input/output
    TString name = Form("cinputNonPrompt%s", suffix.Data());
    AliAnalysisDataContainer *cinputNonPromptLc = mgr->CreateContainer(name, TChain::Class(), AliAnalysisManager::kInputContainer);
    TString outputfile = AliAnalysisManager::GetCommonFileName();
    outputfile += ":PWGHF_D2H_NonPrompt";
    outputfile += suffix.Data();

    name = Form("coutputCutsNonPrompt%s", suffix.Data());
    AliAnalysisDataContainer *coutputNonPromptLcCuts = mgr->CreateContainer(name, TList::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
    name = Form("coutputNonPrompt%s", suffix.Data());
    AliAnalysisDataContainer *coutputNonPromptLc = mgr->CreateContainer(name, TList::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
    name = Form("coutputNormNonPrompt%s", suffix.Data());
    AliAnalysisDataContainer *coutputNonPromptLcNorm = mgr->CreateContainer(name, AliNormalizationCounter::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
    AliAnalysisDataContainer *coutputNonPromptLcML = nullptr;
    if (createMLtree)
    {
        name = Form("coutputMLNonPrompt%s", suffix.Data());
        coutputNonPromptLcML = mgr->CreateContainer(name, TTree::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
    }

    mgr->ConnectInput(npLcTask, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(npLcTask, 1, coutputNonPromptLc);
    mgr->ConnectOutput(npLcTask, 2, coutputNonPromptLcCuts);
    mgr->ConnectOutput(npLcTask, 3, coutputNonPromptLcNorm);
    if (createMLtree)
        mgr->ConnectOutput(npLcTask, 4, coutputNonPromptLcML);

    return npLcTask;
}
