AliAnalysisTaskStrangeCascadesTriggerAODRun2MC* AddTaskStrangeCascadesTriggerAODRun2MC(TString name = "name", TString lExtraOutputName = "", Bool_t lSaveCascades = kTRUE, Bool_t lSaveV0s = kFALSE, Bool_t lSaveTracks = kFALSE, Bool_t lSaveMCGenPart = kFALSE)
{
    // Creates, configures and attaches to the train a cascades check task.
    // Get the pointer to the existing analysis manager via the static access method.
    //==============================================================================
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddTaskStrangeCascadesTriggerAODRun2MC", "No analysis manager to connect to.");
        return 0x0;
    }
    
    // Check the analysis type using the event handlers connected to the analysis manager.
    //==============================================================================
    if (!mgr->GetInputEventHandler()) {
        ::Error("AddTaskStrangeCascadesTriggerAODRun2MC", "This task requires an input event handler");
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
    AliAnalysisTaskStrangeCascadesTriggerAODRun2MC* task = new AliAnalysisTaskStrangeCascadesTriggerAODRun2MC(name.Data(), lSaveCascades, lSaveV0s, lSaveTracks, lSaveMCGenPart);   
    if(!task) return 0x0;
    
    mgr->AddTask(task);
    
    AliAnalysisDataContainer *coutputList           = 0x0;
    AliAnalysisDataContainer *coutputTree           = 0x0;
    
    coutputList = mgr->CreateContainer("coutputList", TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName.Data());
    
    coutputTree = mgr->CreateContainer("fOutputTree", TTree::Class(), AliAnalysisManager::kOutputContainer, outputFileName.Data());
    coutputTree->SetSpecialOutput();
    
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    
    mgr->ConnectOutput(task, 1, coutputList);
    mgr->ConnectOutput(task, 2, coutputTree);
    
    return task;
}
