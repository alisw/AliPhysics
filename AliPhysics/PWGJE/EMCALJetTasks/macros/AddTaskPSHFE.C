    
AliAnalysisTaskPSHFE* AddTaskPSHFE(const char* taskname, Bool_t trkCutsStrong=kFALSE, Bool_t SSCuts=kFALSE, Bool_t UseNonSignalEvents=kFALSE)
{
    //==============================================================================
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr)
    {
        ::Error("AddTaskEmcalJetSample", "No analysis manager to connect to.");
        return NULL;
    }

    // Check the analysis type using the event handlers connected to the analysis manager.
    //==============================================================================
    if (!mgr->GetInputEventHandler())
    {
        ::Error("AddTaskEmcalJetSample", "This task requires an input event handler");
        return NULL;
    }

    AliAnalysisTaskPSHFE* PSHFEtask = new AliAnalysisTaskPSHFE(taskname);

    PSHFEtask->SetElectronTrackCuts(trkCutsStrong);
    PSHFEtask->SetSSCutBool(SSCuts);
    PSHFEtask->SetUseNonSignalEvents(UseNonSignalEvents);
    mgr->AddTask(PSHFEtask);  

    TString contname = "";

    if(trkCutsStrong){contname+=TString("_Strong");}else{contname+=TString("_Weak");}
    if(SSCuts){contname+=TString("_SS");}else{contname+=TString("_NoSS");}
    if(UseNonSignalEvents){contname+=TString("_NoSig");}else{contname+=TString("_Sig");}
    // create containers for input/output
    AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(Form("Min-Bias%s",contname.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s", AliAnalysisManager::GetCommonFileName()));

    AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(Form("EMCal7%s",contname.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s", AliAnalysisManager::GetCommonFileName()));

    AliAnalysisDataContainer *coutput3 = mgr->CreateContainer(Form("EMCalEGA%s",contname.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s", AliAnalysisManager::GetCommonFileName()));
    
    AliAnalysisDataContainer *coutput4 = mgr->CreateContainer(Form("EMCalJet%s",contname.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s", AliAnalysisManager::GetCommonFileName()));



    // connect input/output
    mgr->ConnectInput(PSHFEtask, 0, cinput);
    mgr->ConnectOutput(PSHFEtask, 1, coutput1);
    mgr->ConnectOutput(PSHFEtask, 2, coutput2);
    mgr->ConnectOutput(PSHFEtask, 3, coutput3);
    mgr->ConnectOutput(PSHFEtask, 4, coutput4);

    return PSHFEtask;
}
