AliAnalysisTaskSEDstarPolarization *AddTaskDstarPolarization(bool readMC = false,
                                                             TString fileName = "",
                                                             TString suffix = "",
                                                             TString cutObjName = "AnalysisCuts",
                                                             int AODProtection = 0)
{
    // \brief: AddTask for AliAnalysisTaskSEDstarPolarization
    // \authors:
    // F. Grosa, fabrizio.grosa@cern.ch
    // S. Kundu, sourav.kundu@cern.ch
    //==============================================================================

    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr)
        ::Error("AddTaskDstarPolarization", "No analysis manager to connect to.");

    TFile *fileCuts = TFile::Open(fileName.Data());
    if (!fileCuts || (fileCuts && !fileCuts->IsOpen()))
        ::Fatal("AddTaskDstarPolarization", "Cut file not found on Grid: analysis will not start!\n");
    analysisCuts = (AliRDHFCutsDStartoKpipi*)fileCuts->Get(cutObjName.Data());

    // Analysis Task
    AliAnalysisTaskSEDstarPolarization *dMesonTask = new AliAnalysisTaskSEDstarPolarization("DstarPolarizationAnalysis", analysisCuts);
    dMesonTask->SetAODMismatchProtection(AODProtection);
    dMesonTask->SetReadMC(readMC);
    mgr->AddTask(dMesonTask);

    // Create containers for input/output
    TString name = Form("cinputDstarPolarization%s", suffix.Data());
    AliAnalysisDataContainer *cinput = mgr->CreateContainer(name, TChain::Class(), AliAnalysisManager::kInputContainer);
    TString outputfile = AliAnalysisManager::GetCommonFileName();
    outputfile += ":PWGHF_D2H_DstarPolarization";
    outputfile += suffix.Data();

    name = Form("coutputCutsDstarPolarization%s", suffix.Data());
    AliAnalysisDataContainer *coutputCuts = mgr->CreateContainer(name, TList::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
    name = Form("coutputDstarPolarization%s", suffix.Data());
    AliAnalysisDataContainer *coutput = mgr->CreateContainer(name, TList::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
    name = Form("coutputNormDstarPolarization%s", suffix.Data());

    mgr->ConnectInput(dMesonTask, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(dMesonTask, 1, coutput);
    mgr->ConnectOutput(dMesonTask, 2, coutputCuts);

    return dMesonTask;
}
