AliAnalysisTaskMCPredictions *AddTaskMCPredictions( Int_t lNSmallBinning = 1000, Int_t lNLargeBinning = 4000, Int_t lRebinFactor = 1, const TString lMasterJobSessionFlag = "")
{
    // Creates, configures and attaches to the train a cascades check task.
    // Get the pointer to the existing analysis manager via the static access method.
    //==============================================================================
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddTaskMCPredictions", "No analysis manager to connect to.");
        return NULL;
    }
    
    // Check the analysis type using the event handlers connected to the analysis manager.
    //==============================================================================
    if (!mgr->GetInputEventHandler()) {
        printf("Caution: No input handler detected!");
    }
    //TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
    
    // Create and configure the task
    AliAnalysisTaskMCPredictions *taskMCPred = new AliAnalysisTaskMCPredictions("taskMCPred");

    mgr->AddTask(taskMCPred);
    TString outputFileName = AliAnalysisManager::GetCommonFileName();
    
    outputFileName += ":PWGLF_MCPredictions";
    
    Printf("Set OutputFileName : \n %s\n", outputFileName.Data() );
    
    AliAnalysisDataContainer *coutputList = mgr->CreateContainer("cList",
                                                                 TList::Class(),
                                                                 AliAnalysisManager::kOutputContainer,
                                                                 outputFileName );
    
    //Recommendation: Tree as a single output slot
    mgr->ConnectInput (taskMCPred, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(taskMCPred, 1, coutputList);
    
    return taskMCPred;
}   
