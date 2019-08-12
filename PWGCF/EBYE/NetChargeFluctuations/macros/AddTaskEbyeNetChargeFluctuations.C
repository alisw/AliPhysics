
AliAnalysisTaskEbyeNetChargeFluctuations* AddTaskEbyeNetChargeFluctuations(Int_t MCthere)
{
    // get the manager via the static access member. since it's static, you don't need
    // an instance of the class to call the function
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        return 0x0;
    }
    
    // now we create an instance for task
    AliAnalysisTaskEbyeNetChargeFluctuations* task = new AliAnalysisTaskEbyeNetChargeFluctuations("TaskEbyE");
    if(!task) return 0x0;
  //   task->SelectCollisionCandidates(AliVEvent::kINT7);
    
    if(MCthere){
        task->SetAnalysisType(MCthere);
    }
    
    // add your task to the manager
    mgr->AddTask(task);
    
    // your task needs input: here we connect the manager to your task
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    
    // same for the output
    mgr->ConnectOutput(task,1,mgr->CreateContainer("fOutputList", TList::Class(), AliAnalysisManager::kOutputContainer, AliAnalysisManager::GetCommonFileName()));
    mgr->ConnectOutput(task,2,mgr->CreateContainer("fTree", TTree::Class(), AliAnalysisManager::kOutputContainer, AliAnalysisManager::GetCommonFileName()));
    mgr->ConnectOutput(task,3,mgr->CreateContainer("fTreeMCrec", TTree::Class(), AliAnalysisManager::kOutputContainer, AliAnalysisManager::GetCommonFileName()));
    mgr->ConnectOutput(task,4,mgr->CreateContainer("fTreeMCallrec", TTree::Class(), AliAnalysisManager::kOutputContainer, AliAnalysisManager::GetCommonFileName()));
    mgr->ConnectOutput(task,5,mgr->CreateContainer("fTreeMCgen", TTree::Class(), AliAnalysisManager::kOutputContainer, AliAnalysisManager::GetCommonFileName()));
    mgr->ConnectOutput(task,6,mgr->CreateContainer("fTreeMCallgen", TTree::Class(), AliAnalysisManager::kOutputContainer, AliAnalysisManager::GetCommonFileName()));
    
    return task;
    
}
