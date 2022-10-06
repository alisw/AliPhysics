AliAnalysisTaskTRDtriggerTracks* AddTask_TRDtriggerTracks(TString name = "name")
{
    // get the manager via the static access member. since it's static, you don't need
    // to create an instance of the class here to call the function
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        return 0x0;
    }
    // get the input event handler, again via a static method. 
    // this handler is part of the managing system and feeds events
    // to your task
    if (!mgr->GetInputEventHandler()) {
        return 0x0;
    }
    // by default, a file is open for writing. here, we get the filename
    TString fileName = AliAnalysisManager::GetCommonFileName();
    fileName += ":MyTask";      // create a subfolder in the file
    // now we create an instance of your task
    AliAnalysisTaskTRDtriggerTracks* task = new AliAnalysisTaskTRDtriggerTracks(name.Data());   
    if(!task) return 0x0;
    task->SelectCollisionCandidates(AliVEvent::kAny);

    mgr->AddTask(task);

    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task,1,mgr->CreateContainer("MyOutputContainer", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task, 2,	mgr->CreateContainer("tree", TTree::Class(),AliAnalysisManager::kOutputContainer,mgr->GetCommonFileName()));
    return task;
}
