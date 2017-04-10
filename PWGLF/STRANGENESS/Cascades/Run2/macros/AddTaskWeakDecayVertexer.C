AliAnalysisTaskWeakDecayVertexer *AddTaskWeakDecayVertexer( TString lExtraOptions = "", const TString lMasterJobSessionFlag = "")
{
    // Creates, configures and attaches to the train a cascades check task.
    // Get the pointer to the existing analysis manager via the static access method.
    //==============================================================================
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddTaskWeakDecayVertexer", "No analysis manager to connect to.");
        return NULL;
    }
    
    // Check the analysis type using the event handlers connected to the analysis manager.
    //==============================================================================
    if (!mgr->GetInputEventHandler()) {
        ::Error("AddTaskWeakDecayVertexer", "This task requires an input event handler");
        return NULL;
    }
    TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
    
    // Create and configure the task
    AliAnalysisTaskWeakDecayVertexer *taskWDvertexer = new AliAnalysisTaskWeakDecayVertexer("taskWDvertexer", lExtraOptions);
    mgr->AddTask(taskWDvertexer);
    TString outputFileName = AliAnalysisManager::GetCommonFileName();
    
    outputFileName += ":PWGLF_WDVertexer";
    //if (mgr->GetMCtruthEventHandler()) outputFileName += "_MC";
    
    Printf("Set OutputFileName : \n %s\n", outputFileName.Data() );
    
    AliAnalysisDataContainer *coutputList = mgr->CreateContainer("cListVertexer",
                                                                 TList::Class(),
                                                                 AliAnalysisManager::kOutputContainer,
                                                                 outputFileName );
    
    //Recommendation: Tree as a single output slot
    mgr->ConnectInput (taskWDvertexer, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(taskWDvertexer, 1, coutputList);
    
    return taskWDvertexer;
}   
