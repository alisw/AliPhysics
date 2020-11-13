// TODO LIST
// TODO: You're all set!

AliAnalysisTaskPhiCount* AddMyTask( Bool_t MCFlag, Bool_t PhiFlag, Bool_t KaonFlag, TString name = "name" )
{
    // Analysis Manager
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    
    // Filename
    TString fileName = AliAnalysisManager::GetCommonFileName();
    
    // Task initialisation
    AliAnalysisTaskPhiCount* task = new AliAnalysisTaskPhiCount(name.Data());
    
    // Checks
    if (!mgr)                           return 0x0;
    if (!mgr->GetInputEventHandler())   return 0x0;
    if (!task)                          return 0x0;
    
    // task Selection
    task->SelectCollisionCandidates(AliVEvent::kAnyINT);
    
    // MC option
    task-> fSetMCFlag(MCFlag);
    task-> fSetPhiFlag(PhiFlag);
    task-> fSetKaonFlag(KaonFlag);
    
    // Input
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    
    // Output
    // - // TLists
    
    mgr->ConnectOutput(task,1,mgr->CreateContainer("fAnalysisOutputList",   TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task,2,mgr->CreateContainer("fQCOutputList",         TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    
    // - // TTrees
    
    if ( PhiFlag )  mgr->ConnectOutput(task,3,mgr->CreateContainer("OutContainer3", TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    if ( KaonFlag ) mgr->ConnectOutput(task,4,mgr->CreateContainer("OutContainer4", TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    if ( MCFlag )
    {
        if ( PhiFlag )  mgr->ConnectOutput(task,5,mgr->CreateContainer("OutContainer5", TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
        if ( KaonFlag ) mgr->ConnectOutput(task,6,mgr->CreateContainer("OutContainer6", TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    }
    
    // Add task
    mgr->AddTask(task);
    
    return task;
}
