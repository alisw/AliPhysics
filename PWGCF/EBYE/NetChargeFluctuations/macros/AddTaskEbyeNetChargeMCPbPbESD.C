
AliAnalysisTaskEbyeNetChargeMCPbPbESD* AddTaskEbyeNetChargeMCPbPbESD()
{
    // get the manager via the static access member. since it's static, you don't need
    // an instance of the class to call the function
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        return 0x0;
    }
    
    // now we create an instance for task
    AliAnalysisTaskEbyeNetChargeMCPbPbESD* task = new AliAnalysisTaskEbyeNetChargeMCPbPbESD("TaskEbyE");
    if(!task) return 0x0;
    
    // add your task to the manager
    mgr->AddTask(task);
    
    // your task needs input: here we connect the manager to your task
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    
    // same for the output
    mgr->ConnectOutput(task,1,mgr->CreateContainer("fOutputList", TList::Class(), AliAnalysisManager::kOutputContainer, AliAnalysisManager::GetCommonFileName()));
    mgr->ConnectOutput(task,2,mgr->CreateContainer("fTreeMCgen", TTree::Class(), AliAnalysisManager::kOutputContainer, AliAnalysisManager::GetCommonFileName()));
    
    
    return task;
    
}
