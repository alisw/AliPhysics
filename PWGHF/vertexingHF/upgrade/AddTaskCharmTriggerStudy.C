AliAnalysisTaskSECharmTriggerStudy *AddTaskCharmTriggerStudy(int system = AliAnalysisTaskSECharmTriggerStudy::kpp,
                                                             bool enable2prongs = true,
                                                             bool enable3prongsDplus = true,
                                                             bool enable3prongsDs = true,
                                                             bool enable3prongsLc = true,
                                                             bool enableDstars = true,
                                                             bool enableCascades = false,
                                                             bool enableBplus = false,
                                                             bool enableB0 = false,
                                                             bool enableBs = false,
                                                             bool enableLb = false,
                                                             bool fillOnlySignal = false,
                                                             TString cutfilename = "",
                                                             bool applyCuts = false,
                                                             TString suffix = "")
{
    //
    // Test macro for the AliAnalysisTaskSE for charm-trigger studies

    // Get the pointer to the existing analysis manager via the static access method.
    //============================================================================
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr)
    {
        ::Error("AliAnalysisTaskSECharmTriggerStudy", "No analysis manager to connect to.");
        return NULL;
    }
    if (!mgr->GetInputEventHandler())
    {
        ::Error("AliAnalysisTaskSECharmTriggerStudy", "This task requires an input event handler");
        return NULL;
    }
    TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
    if (type.Contains("ESD"))
    {
        ::Error("AliAnalysisTaskSECharmTriggerStudy", "This task requires to run on AOD");
        return NULL;
    }

    //cuts files
    AliRDHFCutsD0toKpi* cutsD0 = nullptr;
    AliRDHFCutsDplustoKpipi* cutsDplus = nullptr;
    AliRDHFCutsDStartoKpipi* cutsDstar = nullptr;
    AliRDHFCutsDstoKKpi* cutsDs = nullptr;
    AliRDHFCutsLctopKpi* cutsLc = nullptr;
    AliRDHFCutsLctoV0* cutsLctoV0 = nullptr;
    if(cutfilename != "")
    {
        TFile* infilecuts = TFile::Open(cutfilename);
        if(infilecuts)
        {
            cutsD0     = (AliRDHFCutsD0toKpi*)infilecuts->Get("D0toKpiCuts");
            cutsDplus  = (AliRDHFCutsDplustoKpipi*)infilecuts->Get("DplustoKpipiCuts");
            cutsDstar  = (AliRDHFCutsDStartoKpipi*)infilecuts->Get("DstoKKpiCuts");
            cutsDs     = (AliRDHFCutsDstoKKpi*)infilecuts->Get("LctopKpiCuts");
            cutsLc     = (AliRDHFCutsLctopKpi*)infilecuts->Get("DstartoKpipiCuts");
            cutsLctoV0 = (AliRDHFCutsLctoV0*)infilecuts->Get("LctoV0bachCuts");
        }
    }

    TList* listofcuts = new TList();
    listofcuts->SetName("fListCuts");
    listofcuts->SetOwner(true);
    if(cutsD0)
        listofcuts->Add(cutsD0);
    if(cutsDplus)
        listofcuts->Add(cutsDplus);
    if(cutsDstar)
        listofcuts->Add(cutsDstar);
    if(cutsDs)
        listofcuts->Add(cutsDs);
    if(cutsLc)
        listofcuts->Add(cutsLc);
    if(cutsLctoV0)
        listofcuts->Add(cutsLctoV0);

    // Analysis task
    AliAnalysisTaskSECharmTriggerStudy *chTask = new AliAnalysisTaskSECharmTriggerStudy("CharmTriggerStudyTask", listofcuts);
    chTask->Enable2Prongs(enable2prongs);
    chTask->Enable3Prongs(enable3prongsDplus, enable3prongsDs, enable3prongsLc);
    chTask->EnableDstars(enableDstars);
    chTask->EnableCascades(enableCascades);
    chTask->EnableBeauty3Prongs(enableBplus);
    chTask->EnableBeauty4Prongs(enableB0, enableBs, enableLb);
    chTask->SetFillOnlySignal(fillOnlySignal);
    chTask->SetSystem(system);
    chTask->ApplyCuts(applyCuts);
    mgr->AddTask(chTask);

    // Create containers for input/output
    TString contname = Form("cinputChTrigger%s", suffix.Data());
    AliAnalysisDataContainer *cinputcont = mgr->CreateContainer(contname.Data(), TChain::Class(), AliAnalysisManager::kInputContainer);

    TString outputfile = AliAnalysisManager::GetCommonFileName();
    TString outputdirname = Form("%s:PWGHF_D2H_CharmTrigger_%s", outputfile.Data(), suffix.Data());

    contname = Form("coutputChTrigger%s", suffix.Data());
    AliAnalysisDataContainer *coutputlist = mgr->CreateContainer(contname.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, outputdirname.Data());

    contname = Form("coutputChTriggerRecoTree%s", suffix.Data());
    AliAnalysisDataContainer *coutputrecotree = mgr->CreateContainer(contname.Data(), TTree::Class(), AliAnalysisManager::kOutputContainer, outputdirname.Data());
    coutputrecotree->SetSpecialOutput();

    contname = Form("coutputChTriggerGenTree%s", suffix.Data());
    AliAnalysisDataContainer *coutputgentree = mgr->CreateContainer(contname.Data(), TTree::Class(), AliAnalysisManager::kOutputContainer, outputdirname.Data());
    coutputgentree->SetSpecialOutput();

    contname = Form("coutputChTriggerCuts%s", suffix.Data());
    AliAnalysisDataContainer *coutputcuts    = mgr->CreateContainer(contname.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, outputdirname.Data());

    mgr->ConnectInput(chTask, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(chTask, 1, coutputlist);
    mgr->ConnectOutput(chTask, 2, coutputrecotree);
    mgr->ConnectOutput(chTask, 3, coutputgentree);
    mgr->ConnectOutput(chTask, 4, coutputcuts);

    return chTask;
}
