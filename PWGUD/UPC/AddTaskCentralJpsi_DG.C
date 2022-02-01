AliAnalysisTaskCentralJpsi_DG *AddTaskCentralJpsi_DG()
{
    TString name = "CentralJpsi2018";
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
    Bool_t isMC = kFALSE;
    if(mgr->GetMCtruthEventHandler()) isMC = kTRUE;

    // create task
    AliAnalysisTaskCentralJpsi_DG* task = new AliAnalysisTaskCentralJpsi_DG(name.Data());
    task->SetIsMC(isMC);
    // add task to the manager
    mgr->AddTask(task);

    // get the name of the file that is created and create subfolder
    TString fileName = AliAnalysisManager::GetCommonFileName();
    fileName += ":AnalysisOutput";      

    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    // same for the output
    mgr->ConnectOutput(task,1,mgr->CreateContainer("fTreeJpsi", TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task,2,mgr->CreateContainer("fTreeJpsiMCGen", TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task,3,mgr->CreateContainer("fOutputList", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));

    return task;
}