AliAnalysisTaskhStrCorr *AddTaskhStrCorr( TString lCustomName = "", const TString lMasterJobSessionFlag = "")
{
    // Creates, configures and attaches to the train a cascades check task.
    // Get the pointer to the existing analysis manager via the static access method.
    //==============================================================================
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddTaskhStrCorr", "No analysis manager to connect to.");
        return NULL;
    }
    
    // Check the analysis type using the event handlers connected to the analysis manager.
    //==============================================================================
    if (!mgr->GetInputEventHandler()) {
        ::Error("AddTaskhStrCorr", "This task requires an input event handler");
        return NULL;
    }
    TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
    
    // Create and configure the task
    AliAnalysisTaskhStrCorr *taskv0extract = new AliAnalysisTaskhStrCorr("taskhStrCorr");
    
    mgr->AddTask(taskv0extract);
    
    TString outputFileName = AliAnalysisManager::GetCommonFileName();
    
    outputFileName += ":PWGLF";
    outputFileName += "_";
    outputFileName += lCustomName.Data();
    
    Printf("Set OutputFileName : \n %s\n", outputFileName.Data() );
    
    AliAnalysisDataContainer *coutputList = mgr->CreateContainer("clist",
                                                                 TList::Class(),
                                                                 AliAnalysisManager::kOutputContainer,
                                                                 outputFileName );
    
    mgr->ConnectInput( taskv0extract, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(taskv0extract, 1, coutputList);
    
    return taskv0extract;
}   
