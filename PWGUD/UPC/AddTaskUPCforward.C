AliAnalysisTaskUPCforward* AddTaskUPCforward( TString name = "name", Int_t _fSetSingleMuonPt = 0, const char* suffix = "" )
{
    // TString name = "PolarisationSys";
    TString combinedName;
    combinedName.Form("PolarisationSys_%s", suffix);


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
    // fileName += ":MyTask";      // create a subfolder in the file
    // now we create an instance of your task
    // AliAnalysisTaskUPCforward* task = new AliAnalysisTaskUPCforward(combinedName.Data(), _fSetSingleMuonPt);
    AliAnalysisTaskUPCforward* task = new AliAnalysisTaskUPCforward(Form("MyOutputContainer%s", suffix), _fSetSingleMuonPt);
    if(!task) return 0x0;
    // task->SelectCollisionCandidates(AliVEvent::kAnyINT);     // Physics Selection used by "everybody" but NOT in UPC
    // add your task to the manager
    mgr->AddTask(task);
    // your task needs input: here we connect the manager to your task
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    // same for the output
    // mgr->ConnectOutput(task,1,mgr->CreateContainer(combinedName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task,1,mgr->CreateContainer(Form("MyOutputContainer%s", suffix), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    // in the end, this macro returns a pointer to your task. this will be convenient later on
    // when you will run your analysis in an analysis train on grid
    return task;
}
