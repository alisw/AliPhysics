AliAnalysisTaskMCPredictions *AddTaskMCPredictions( Int_t lNSmallBinning = 1000, Int_t lNLargeBinning = 2000, Int_t lRebinFactor = 1, Int_t lNBBins = 1, Int_t lNNpartBins = 1, Int_t lNEtaBins = 800, TString lOutputDir = "", const TString lMasterJobSessionFlag = "")
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
    AliAnalysisTaskMCPredictions *taskMCPred = new AliAnalysisTaskMCPredictions("taskMCPred", lNSmallBinning, lNLargeBinning, lRebinFactor, lNBBins, lNNpartBins, lNEtaBins);

    mgr->AddTask(taskMCPred);
    TString outputFileName = AliAnalysisManager::GetCommonFileName();
    
    outputFileName += ":PWGLF_MCPredictions";
    
    Printf("Set OutputFileName : \n %s\n", outputFileName.Data() );
    
    AliAnalysisDataContainer *coutputList = mgr->CreateContainer(Form("cList%s", lOutputDir.Data()),
                                                                 TList::Class(),
                                                                 AliAnalysisManager::kOutputContainer,
                                                                 outputFileName );
  AliAnalysisDataContainer *coutputTree = mgr->CreateContainer(Form("cTree%s", lOutputDir.Data()),
                                                               TTree::Class(),
                                                               AliAnalysisManager::kOutputContainer,
                                                               outputFileName );
  
    //Recommendation: Tree as a single output slot
    mgr->ConnectInput (taskMCPred, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(taskMCPred, 1, coutputList);
    mgr->ConnectOutput(taskMCPred, 2, coutputTree);
    
    return taskMCPred;
}   
