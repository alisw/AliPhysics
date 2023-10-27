AliAnalysisTaskMCPredictions2pc *AddTaskMCPredictions2pc( Int_t lNSmallBinning = 1000, Int_t lNLargeBinning = 2000, Int_t lRebinFactor = 1, Int_t lNEtaBins = 80, TString lOutputDir = "", const TString lMasterJobSessionFlag = "")
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
    AliAnalysisTaskMCPredictions2pc *taskMCPred2pc = new AliAnalysisTaskMCPredictions2pc("taskMCPred2pc", lNSmallBinning, lNLargeBinning, lRebinFactor, lNEtaBins);

    mgr->AddTask(taskMCPred2pc);
    TString outputFileName = AliAnalysisManager::GetCommonFileName();
    
    outputFileName += ":PWGLF_MCPredictions2pc";
    
    Printf("Set OutputFileName : \n %s\n", outputFileName.Data() );
    
    AliAnalysisDataContainer *coutputList = mgr->CreateContainer(Form("cList%s", lOutputDir.Data()),
                                                                 TList::Class(),
                                                                 AliAnalysisManager::kOutputContainer,
                                                                 outputFileName );
  
    //Recommendation: Tree as a single output slot
    mgr->ConnectInput (taskMCPred2pc, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(taskMCPred2pc, 1, coutputList);
    
    return taskMCPred2pc;
}   
