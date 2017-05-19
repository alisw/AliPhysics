AliAnalysisTaskPPVsMultCrossCheckMC *AddTaskPPVsMultCrossCheckMC( const TString lMasterJobSessionFlag = "")
{
    // Creates, configures and attaches to the train a cascades check task.
    // Get the pointer to the existing analysis manager via the static access method.
    //==============================================================================
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AliAnalysisTaskPPVsMultCrossCheckMC", "No analysis manager to connect to.");
        return NULL;
    }
    
    // Check the analysis type using the event handlers connected to the analysis manager.
    //==============================================================================
    if (!mgr->GetInputEventHandler()) {
        printf("Caution: No input handler detected!");
    }
    //TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
    
    // Create and configure the task
    AliAnalysisTaskPPVsMultCrossCheckMC *taskAuxiliary = new AliAnalysisTaskPPVsMultCrossCheckMC("taskAuxiliary");
      taskAuxiliary->SetSelectedTriggerClass(AliVEvent::kINT7);
      taskAuxiliary->SetUseMultSelection(kTRUE);
      taskAuxiliary->SetSkipPS(kTRUE);
    mgr->AddTask(taskAuxiliary);
    TString outputFileName = AliAnalysisManager::GetCommonFileName();
    
    outputFileName += ":PWGLF_PPVsMultXCheck";
    if (mgr->GetMCtruthEventHandler()) outputFileName += "_MC";
    
    Printf("Set OutputFileName : \n %s\n", outputFileName.Data() );
    
    AliAnalysisDataContainer *coutputList = mgr->CreateContainer("cList",
                                                                 TList::Class(),
                                                                 AliAnalysisManager::kOutputContainer,
                                                                 outputFileName );
    
    //Recommendation: Tree as a single output slot
    mgr->ConnectInput (taskAuxiliary, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(taskAuxiliary, 1, coutputList);
    
    return taskAuxiliary;
}   
