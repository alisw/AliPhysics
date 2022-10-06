AliAnalysisTaskMCPredictionsStrgVsMultVsZDC *AddTaskMCPredictionsStrgVsMultVsZDC(Float_t lCenterOfMassEnergy = 13000, Bool_t kDoPythia = kTRUE, Bool_t kDoEPOS = kFALSE, TString lOutputDir = "", const TString lMasterJobSessionFlag = "")
{
    // Creates, configures and attaches to the train a cascades check task.
    // Get the pointer to the existing analysis manager via the static access method.
    //==============================================================================
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddTaskMCPredictionsStrgVsMultVsZDC", "No analysis manager to connect to.");
        return NULL;
    }
    
    // Check the analysis type using the event handlers connected to the analysis manager.
    //==============================================================================
    if (!mgr->GetInputEventHandler()) {
        printf("Caution: No input handler detected!");
    }
    //TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
    
    // Create and configure the task
    AliAnalysisTaskMCPredictionsStrgVsMultVsZDC *taskMCPred = new AliAnalysisTaskMCPredictionsStrgVsMultVsZDC("taskMCPred", lCenterOfMassEnergy, kDoPythia, kDoEPOS);

    mgr->AddTask(taskMCPred);
    TString outputFileName = AliAnalysisManager::GetCommonFileName();
    
    outputFileName += ":PWGLF_MCPredictionsStrgVsMultVsZDC";
    
    Printf("Set OutputFileName : \n %s\n", outputFileName.Data() );
    
    AliAnalysisDataContainer *coutputList = mgr->CreateContainer(Form("cList%s", lOutputDir.Data()),
                                                                 TList::Class(),
                                                                 AliAnalysisManager::kOutputContainer,
                                                                 outputFileName );
    
    //Recommendation: Tree as a single output slot
    mgr->ConnectInput (taskMCPred, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(taskMCPred, 1, coutputList);
    
    return taskMCPred;
}   
