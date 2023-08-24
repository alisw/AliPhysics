AliAnalysisTaskJPsiMC_DG *AddTaskJPsiMC_DG(Bool_t Neutral = kFALSE)
{
    TString name = "CentralJPsi2018";
    // get the manager via the static access member. since it's static, you don't need
    // to create an instance of the class here to call the function
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        return 0x0; // returns zero if mgr is a null pointer
    }
    // get the input event handler, again via a static method. 
    // this handler is part of the managing system and feeds events
    // to your task
    if (!mgr->GetInputEventHandler()) {
        return 0x0; //again return zero if the input is a null pointer
    }
    //if (!mgr->GetMCtruthEventHandler()) {
    //    return 0x0;
    //}
    // by default, a file is open for writing. here, we get the filename
    TString fileName = AliAnalysisManager::GetCommonFileName();
    fileName += ":AnalysisOutput";      // create a subfolder in the file
    // now we create an instance of your task
    AliAnalysisTaskJPsiMC_DG* task = new AliAnalysisTaskJPsiMC_DG(name.Data());   
    if(!task) return 0x0;
    // set if neutral pions or not
    task->SetNeutralPions(Neutral);

    //task->SelectCollisionCandidates(AliVEvent::kAnyINT);
    // add your task to the manager
    mgr->AddTask(task);
    // your task needs input: here we connect the manager to your task
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    // same for the output
    mgr->ConnectOutput(task,1,mgr->CreateContainer("fTreeJPsiMCRec", TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task,2,mgr->CreateContainer("fTreeJPsiMCGen", TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task,3,mgr->CreateContainer("fOutputList", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    // in the end, this macro returns a pointer to your task. this will be convenient later on
    // when you will run your analysis in an analysis train on grid
    return task;
}
