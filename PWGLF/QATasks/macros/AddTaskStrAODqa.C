AliAnalysisTaskStrAODqa *AddTaskStrAODqa(bool isMC=kTRUE, bool IsOOBPileUpRem=kTRUE, TString suffix="" )
{

    // analysis manager
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) { ::Error("AddTaskStrAODqa", "No analysis manager to connect to."); return NULL; }
    if (!mgr->GetInputEventHandler()) { ::Error("AddTaskStrAODqa", "This task requires an input event handler"); return NULL; }

    // Create the task and add it to the manager
    TString combinedName=Form("StrAODqa_Task_%s", suffix.Data());
    AliAnalysisTaskStrAODqa *mytask = new AliAnalysisTaskStrAODqa(combinedName);
    mytask->SetMC(isMC);
    mytask->SetOOBPU(IsOOBPileUpRem);
    mgr->AddTask(mytask);

    // output file name
    TString outputFileName = AliAnalysisManager::GetCommonFileName();
    outputFileName += ":PWGLF_StrAODqa";
    if (IsOOBPileUpRem)    outputFileName += "_OOBPileUpRemoved";
    printf("Set OutputFileName : \n %s\n", outputFileName.Data() );

    //output containers

    AliAnalysisDataContainer *coutput_0, *coutput_1, *coutput_2, *coutput_3, *coutput_4, *coutput_5, *coutput_6, *coutput_7;
    coutput_0 = mgr->CreateContainer(Form("chists_eve_%s",combinedName), TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName );
    coutput_1 = mgr->CreateContainer(Form("chists_V0_%s",combinedName), TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName );
    coutput_2 = mgr->CreateContainer(Form("chists_Casc_%s",combinedName), TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName );
    coutput_3 = mgr->CreateContainer(Form("AliEventCuts_%s",combinedName), TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);

    //connecting input and output
    mgr->ConnectInput (mytask, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(mytask, 1, coutput_0);
    mgr->ConnectOutput(mytask, 2, coutput_1);
    mgr->ConnectOutput(mytask, 3, coutput_2);
    mgr->ConnectOutput(mytask, 4, coutput_3);
    return mytask;

}
