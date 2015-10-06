AliMultSelectionTask *AddTaskMultSelection( Bool_t lCalibration = kTRUE, const TString lMasterJobSessionFlag = "")
{
    // Creates, configures and attaches to the train a Multiplicity Selection Task
    // Get the pointer to the existing analysis manager via the static access method.
    //==============================================================================
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddTaskMultSelection", "No analysis manager to connect to.");
        return NULL;
    }
    
    // Check the analysis type using the event handlers connected to the analysis manager.
    //==============================================================================
    if (!mgr->GetInputEventHandler()) {
        ::Error("AddTaskMultSelection", "This task requires an input event handler");
        return NULL;
    }
    TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
    
    // Create and configure the task
    AliMultSelectionTask *taskMultSelection = new AliMultSelectionTask("taskMultSelection");
    mgr->AddTask(taskMultSelection);
    TString outputFileName = AliAnalysisManager::GetCommonFileName();
    
    outputFileName += ":MultSelection";
    if (mgr->GetMCtruthEventHandler()) outputFileName += "_MC";
    
    Printf("Set OutputFileName : \n %s\n", outputFileName.Data() );
    
    AliAnalysisDataContainer *coutputList = mgr->CreateContainer("cListMultSelection",
                                                                 TList::Class(),
                                                                 AliAnalysisManager::kOutputContainer,
                                                                 outputFileName );
    
    if ( lCalibration ){
        AliAnalysisDataContainer *coutputTree = mgr->CreateContainer("cCalibrationTree",
                                                                     TTree::Class(),
                                                                     AliAnalysisManager::kOutputContainer,
                                                                     outputFileName );
    }
    
    //This one you should merge in file-resident ways...
    coutputTree->SetSpecialOutput();
    
    //Recommendation: Tree as a single output slot
    mgr->ConnectInput (taskMultSelection, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(taskMultSelection, 1, coutputList);
    
    if ( lCalibration ) mgr->ConnectOutput(taskMultSelection, 2, coutputTree);
    
    return taskMultSelection;
}   
