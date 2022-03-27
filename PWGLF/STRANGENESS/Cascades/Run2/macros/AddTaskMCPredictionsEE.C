AliAnalysisTaskMCPredictionsEE *AddTaskMCPredictionsEE(TString lOutputDir = "", const TString lMasterJobSessionFlag = "")
{
    // Creates, configures and attaches to the train a cascades check task.
    // Get the pointer to the existing analysis manager via the static access method.
    //==============================================================================
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddTaskMCPredictionsEE", "No analysis manager to connect to.");
        return NULL;
    }
    
    // Check the analysis type using the event handlers connected to the analysis manager.
    //==============================================================================
    if (!mgr->GetInputEventHandler()) {
        printf("Caution: No input handler detected!");
    }
    //TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
    
    // Create and configure the task
    AliAnalysisTaskMCPredictionsEE *taskMCPredEE = new AliAnalysisTaskMCPredictionsEE("taskMCPredEE");

    mgr->AddTask(taskMCPredEE);
    TString outputFileName = AliAnalysisManager::GetCommonFileName();
    
    outputFileName += ":PWGLF_MCPredictionsEE";
    
    Printf("Set OutputFileName : \n %s\n", outputFileName.Data() );
    
    AliAnalysisDataContainer *cOutputT= mgr->CreateContainer("treeZDCmc",TTree::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
    
    //Recommendation: Tree as a single output slot
    mgr->ConnectInput (taskMCPredEE, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(taskMCPredEE, 1, cOutputT);
    
    return taskMCPredEE;
}   
