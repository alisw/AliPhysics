AliAnalysisTaskCentralJpsi_DG *AddTaskCentralJpsi_DG(Bool_t neutral, const char* suffix = "")
{
    // get the current analysis manager
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        Error("AddTaskCentralJpsi_DG", "No analysis manager found.");
        return 0;
    }
    // check the analysis type using the event handlers connected to the analysis manager
    if (!mgr->GetInputEventHandler()) {
        Error("AddTaskCentralJpsi_DG", "This task requires an input event handler");
        return 0;
    }
    // check if this is an analysis of MC data 
    Bool_t isMC = kFALSE;
    if(mgr->GetMCtruthEventHandler()) isMC = kTRUE;

    // create the task and configure it
    TString combinedName;
    combinedName.Form("fOutputList%s", suffix);
    AliAnalysisTaskCentralJpsi_DG* task = new AliAnalysisTaskCentralJpsi_DG(combinedName.Data());
    task->SetIsMC(isMC);
    task->SetIsNeutral(neutral);
    // add the task to the analysis manager
    mgr->AddTask(task);

    // get the name of the file that is created and create subfolder
    TString fileName = AliAnalysisManager::GetCommonFileName();
    fileName += ":AnalysisOutput";      

    // create container for output
    AliAnalysisDataContainer *coutput = mgr->CreateContainer(combinedName, TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data());
    // connect input and output containers with the manager
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task,1,coutput);

    return task;
}
