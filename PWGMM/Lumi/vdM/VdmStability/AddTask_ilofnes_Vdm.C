AliAnalysisTaskVdmStability* AddTask_ilofnes_Vdm(TString name = "name", char *year = "16", Bool_t fillTTree = false, Int_t nRuns = 0) {
    
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  
    // resolve the name of the output file
    TString fileName = AliAnalysisManager::GetCommonFileName();
    fileName += ":MyTask";      // create a subfolder in the file
    
    // now we create an instance of your task
    AliAnalysisTaskVdmStability* task = new AliAnalysisTaskVdmStability(name.Data());
    if (year == "16") task->SetNRuns(650);//max
    if (year == "17") task->SetNRuns(864);//max
    if (year == "18") task->SetNRuns(1000);//790 CB
    if (nRuns > 0) task->SetNRuns(nRuns);
    task->SetNCases(21);
    task->SetFillTTree(fillTTree);

    // add your task to the manager
    mgr->AddTask(task);

    // connect the manager to your task
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());

    // same for the output
    mgr->ConnectOutput(task,1,mgr->CreateContainer("EventStatHist", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task,2,mgr->CreateContainer("EventTree", TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    
    // important: return a pointer to your task
    return task;
    
}
