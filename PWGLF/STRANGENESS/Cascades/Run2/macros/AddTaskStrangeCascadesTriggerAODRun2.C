AliAnalysisTaskStrangeCascadesTriggerAODRun2* AddTaskStrangeCascadesTriggerAODRun2(TString name = "name", TString lExtraOutputName = "", Bool_t lSaveV0s = kFALSE, Bool_t lSaveRsn = kFALSE, Bool_t lSavePrimaries = kFALSE, Bool_t lUseEventMixing = kFALSE)
{
    // Creates, configures and attaches to the train a cascades check task.
    // Get the pointer to the existing analysis manager via the static access method.
    //==============================================================================
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddTaskStrangeCascadesTriggerAODRun2", "No analysis manager to connect to.");
        return 0x0;
    }
    
    // Check the analysis type using the event handlers connected to the analysis manager.
    //==============================================================================
    if (!mgr->GetInputEventHandler()) {
        ::Error("AddTaskStrangeCascadesTriggerAODRun2", "This task requires an input event handler");
        return 0x0;
    }
    TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
    
    // by default, a file is open for writing. here, we get the filename
    TString outputFileName = AliAnalysisManager::GetCommonFileName();
    
    outputFileName += ":PWGLF_StrVsMult";
    if (mgr->GetMCtruthEventHandler()) outputFileName += "_MC";
    outputFileName += lExtraOutputName.Data();
    Printf("Set OutputFileName : \n %s\n", outputFileName.Data() );
    
    // Create and configure the task
    AliAnalysisTaskStrangeCascadesTriggerAODRun2* task = new AliAnalysisTaskStrangeCascadesTriggerAODRun2(name.Data(), lSaveV0s, lSaveRsn, lSavePrimaries, lUseEventMixing);   
    if(!task) return 0x0;
    
    mgr->AddTask(task);
    
    AliAnalysisDataContainer *coutputList           = 0x0;
    AliAnalysisDataContainer *coutputTreeCascade    = 0x0;
    AliAnalysisDataContainer *coutputTreeV0         = 0x0;
    AliAnalysisDataContainer *coutputTreeRsn        = 0x0;
    AliAnalysisDataContainer *coutputTreeRsnMixed   = 0x0;
    AliAnalysisDataContainer *coutputTreePrimTrack  = 0x0;
    
    coutputList = mgr->CreateContainer("coutputList", TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName.Data());
    
    coutputTreeCascade = mgr->CreateContainer("fTreeCascade", TTree::Class(), AliAnalysisManager::kOutputContainer, outputFileName.Data());
    coutputTreeCascade->SetSpecialOutput();
    
    if( lSaveV0s )
    {
        coutputTreeV0 = mgr->CreateContainer("fTreeV0", TTree::Class(), AliAnalysisManager::kOutputContainer, outputFileName.Data());
        coutputTreeV0->SetSpecialOutput();
    }
    
    if( lSaveRsn )
    {
        coutputTreeRsn = mgr->CreateContainer("fTreeRsn", TTree::Class(), AliAnalysisManager::kOutputContainer, outputFileName.Data());
        coutputTreeRsn->SetSpecialOutput();
        
        if( lUseEventMixing )
        {
            coutputTreeRsnMixed = mgr->CreateContainer("fTreeRsnMixed", TTree::Class(), AliAnalysisManager::kOutputContainer, outputFileName.Data());
            coutputTreeRsnMixed->SetSpecialOutput();
        }
    }
    
    if( lSavePrimaries )
    {
        coutputTreePrimTrack = mgr->CreateContainer("fTreePrimTracks", TTree::Class(), AliAnalysisManager::kOutputContainer, outputFileName.Data());
        coutputTreePrimTrack->SetSpecialOutput();
    }
    
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    
    mgr->ConnectOutput(task, 1, coutputList);
    mgr->ConnectOutput(task, 2, coutputTreeCascade);
    if( lSaveV0s                    ) mgr->ConnectOutput(task, 3, coutputTreeV0);
    if( lSaveRsn                    ) mgr->ConnectOutput(task, 4, coutputTreeRsn);
    if( lSaveRsn && lUseEventMixing ) mgr->ConnectOutput(task, 5, coutputTreeRsnMixed);
    if( lSavePrimaries              ) mgr->ConnectOutput(task, 6, coutputTreePrimTrack);
    
    return task;
}
