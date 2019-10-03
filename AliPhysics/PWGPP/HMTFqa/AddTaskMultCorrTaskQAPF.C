const char *anatype = "ESD";
AliAnalysisMultCorrTaskQAPF *AddTaskMultCorrTaskQAPF( const TString lMasterJobSessionFlag = "")
//void AddTaskMultCorrTaskQA()
{
    
    // Creates, configures and attaches to the train a cascades check task.
    // Get the pointer to the existing analysis manager via the static access method.
    //==============================================================================
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddTaskMultCorrTaskQAPF", "No analysis manager to connect to.");
        return NULL;
    }
    // Check the analysis type using the event handlers connected to the analysis manager.
    //==============================================================================
    if (!mgr->GetInputEventHandler()) {
        ::Error("AddTaskMultCorrTaskQAPF", "This task requires an input event handler");
        return NULL;
    }
    TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
    
    // Create and configure the task
    //gROOT->LoadMacro("AliAnalysisMultCorrTaskQA.cxx++g"); //-- me
    AliAnalysisMultCorrTaskQAPF *taskAuxiliary = new AliAnalysisMultCorrTaskQAPF("taskAuxiliary");
    mgr->AddTask(taskAuxiliary); // comment this line in order to work
    
    //-------------------------------------------------------------------------------
    TString outputFileName = AliAnalysisManager::GetCommonFileName();
    
    outputFileName += ":PWGLF_PPVsMult";
    if (mgr->GetMCtruthEventHandler()) outputFileName += "_MC";
    
    Printf("Set OutputFileName : \n %s\n", outputFileName.Data() );
    
    AliAnalysisDataContainer *coutputList = mgr->CreateContainer("cList",
                                                                 TList::Class(),
                                                                 AliAnalysisManager::kOutputContainer,
                                                                 outputFileName );
    
    mgr->ConnectInput (taskAuxiliary, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(taskAuxiliary, 1, coutputList);
    
    return taskAuxiliary;

}
