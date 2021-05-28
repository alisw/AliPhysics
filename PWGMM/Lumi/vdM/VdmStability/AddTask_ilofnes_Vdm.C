AliAnalysisTaskVdmStability* AddTask_ilofnes_Vdm2(TString name = "name", char *year = "16", Bool_t fillTTree = false, Float_t diffMin = 5.5, Float_t diffMax = 11.5, Float_t sumMin = 11.5, Float_t sumMax = 17.5) {
    
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  
    // resolve the name of the output file
    TString fileName = AliAnalysisManager::GetCommonFileName();
    fileName += ":MyTask";      // create a subfolder in the file
    
    // now we create an instance of your task
    AliAnalysisTaskVdmStability* task = new AliAnalysisTaskVdmStability(name.Data());
    if (year == "16") task->SetNRuns(650);//max
    if (year == "17") task->SetNRuns(864);//max
    if (year == "18") task->SetNRuns(1000);//790 CB
    task->SetNCases(25);
    task->SetFillTTree(fillTTree);
    
    task->SetDiffTimingCut(diffMin,diffMax);
    task->SetSumTimingCut(sumMin,sumMax);

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
