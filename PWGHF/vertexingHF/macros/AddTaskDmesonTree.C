AliAnalysisTaskSEDmesonTree *AddTaskDmesonTree(int decCh = AliAnalysisTaskSEDmesonTree::kD0toKpi,
                                               bool readMC = false,
                                               TString fileName = "",
                                               TString suffix = "",
                                               bool createMLtree = true, // for the moment only possibility true
                                               TString cutObjName = "AnalysisCuts",
                                               int AODProtection = 0)
{
    // \brief: AddTask for AliAnalysisTaskSEDmesonTree
    // \authors:
    // F. Grosa, fabrizio.grosa@cern.ch
    //==============================================================================

    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr)
        ::Error("AddTaskTreeLc", "No analysis manager to connect to.");

    TFile *fileCuts = TFile::Open(fileName.Data());
    if (!fileCuts || (fileCuts && !fileCuts->IsOpen()))
        ::Fatal("AddTaskTreeLc", "Cut file not found on Grid: analysis will not start!\n");
    else
        printf("Cut file correctly found\n");

    AliRDHFCuts *analysisCuts = NULL;
    if(decCh == AliAnalysisTaskSEDmesonTree::kD0toKpi)
    {
        analysisCuts = (AliRDHFCutsD0toKpi*)fileCuts->Get(cutObjName.Data());
        suffix.Prepend("D0toKpi");
    }
    else if(decCh == AliAnalysisTaskSEDmesonTree::kDplustoKpipi)
    {
        analysisCuts = (AliRDHFCutsDplustoKpipi*)fileCuts->Get(cutObjName.Data());
        suffix.Prepend("DplustoKpipi");
    }
    else if(decCh == AliAnalysisTaskSEDmesonTree::kDstartoD0pi)
    {
        analysisCuts = (AliRDHFCutsDStartoKpipi*)fileCuts->Get(cutObjName.Data());
        suffix.Prepend("DstartoD0pi");
    }
    else
    {
        ::Fatal("AddTaskDmesonTree", "Specific AliRDHFCuts not found!\n");
        return NULL;
    }

    // Analysis Task
    AliAnalysisTaskSEDmesonTree *dMesonTask = new AliAnalysisTaskSEDmesonTree("DmesonTreeAnalysis", decCh, analysisCuts, createMLtree);
    dMesonTask->SetAODMismatchProtection(AODProtection);
    dMesonTask->SetReadMC(readMC);
    if (createMLtree && readMC)
        dMesonTask->SetFillOnlySignalInMLtree();
    mgr->AddTask(dMesonTask);

    // Create containers for input/output
    TString name = Form("cinputTree%s", suffix.Data());
    AliAnalysisDataContainer *cinputTreeD = mgr->CreateContainer(name, TChain::Class(), AliAnalysisManager::kInputContainer);
    TString outputfile = AliAnalysisManager::GetCommonFileName();
    outputfile += ":PWGHF_D2H_Tree";
    outputfile += suffix.Data();

    name = Form("coutputCutsTree%s", suffix.Data());
    AliAnalysisDataContainer *coutputTreeDCuts = mgr->CreateContainer(name, TList::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
    name = Form("coutputTree%s", suffix.Data());
    AliAnalysisDataContainer *coutputTreeD = mgr->CreateContainer(name, TList::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
    name = Form("coutputNormTree%s", suffix.Data());
    AliAnalysisDataContainer *coutputTreeDNorm = mgr->CreateContainer(name, AliNormalizationCounter::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
    AliAnalysisDataContainer *coutputTreeDML = nullptr;
    if (createMLtree)
    {
        name = Form("coutputMLTree%s", suffix.Data());
        coutputTreeDML = mgr->CreateContainer(name, TTree::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
    }

    mgr->ConnectInput(dMesonTask, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(dMesonTask, 1, coutputTreeD);
    mgr->ConnectOutput(dMesonTask, 2, coutputTreeDCuts);
    mgr->ConnectOutput(dMesonTask, 3, coutputTreeDNorm);
    if (createMLtree)
        mgr->ConnectOutput(dMesonTask, 4, coutputTreeDML);

    return dMesonTask;
}
