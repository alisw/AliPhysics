//class AliAnalysisDataContainer;

AliAnalysisTaskStronglyIntensiveCorrTree* AddMyTask(TString name               = "AliAnalysisTaskStronglyIntensiveCorrTree",
					  UInt_t collisionCandidates = AliVEvent::kINT7,
					  Bool_t isMC = kFALSE,
					  Int_t filterBit = 128,
					  Double_t ptMin=0.2,
					  Double_t ptMax=5
					  )
{
    // get the manager via the static access member. since it's static, you don't need
    // an instance of the class to call the function
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
    fileName +=":fb_";
    fileName += filterBit;      // create a subfolder in the file
    
    // now we create an instance of your task 
    AliAnalysisTaskStronglyIntensiveCorrTree* task = new AliAnalysisTaskStronglyIntensiveCorrTree(name.Data());   
    //task->SetParameterName(param,...)
    if(!task) return 0x0;
    
    task->SelectCollisionCandidates(collisionCandidates); //tigger
    // add your task to the manager
    task->SetTrackBit(filterBit);
    task->SetMCStatus(isMC);
    task->SetPtRange(ptMin,ptMax);
    mgr->AddTask(task);
    // your task needs input: here we connect the manager to your task
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    // same for the output
    mgr->ConnectOutput(task,1,mgr->CreateContainer("fTree", TTree::Class(),
						   AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task,2,mgr->CreateContainer("OutputList", TList::Class(),
						   AliAnalysisManager::kOutputContainer, fileName.Data()));
    if(isMC) mgr->ConnectOutput(task,3,mgr->CreateContainer("fTreeMCPrim", TTree::Class(),
						   AliAnalysisManager::kOutputContainer, fileName.Data()));
    if(isMC) mgr->ConnectOutput(task,4,mgr->CreateContainer("fTreeMC", TTree::Class(),
						   AliAnalysisManager::kOutputContainer, fileName.Data()));
    // in the end, this macro returns a pointer to your task. this will be convenient later on
    // when you will run your analysis in an analysis train on grid
    return task;
}
