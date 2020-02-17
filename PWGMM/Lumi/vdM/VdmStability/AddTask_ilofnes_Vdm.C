AliAnalysisTaskVdmStability* AddTask_ilofnes_Vdm(TString name = "name") {
    
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  
    // resolve the name of the output file
    TString fileName = AliAnalysisManager::GetCommonFileName();
    fileName += ":MyTask";      // create a subfolder in the file
    
    // now we create an instance of your task
    AliAnalysisTaskVdmStability* task = new AliAnalysisTaskVdmStability(name.Data());

    // add your task to the manager
    mgr->AddTask(task);
    
    
    //Add event filter
    //define evnt cuts ...
    //task->SetEventFilter(eventCuts);
    
    //Pileup rejection
    //task->SetRejectPileup(kFALSE);

    // connect the manager to your task
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());

    // same for the output
    mgr->ConnectOutput(task,1,mgr->CreateContainer("EventStatHist", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task,2,mgr->CreateContainer("EventTree", TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    
    // important: return a pointer to your task
    return task;
    
}
