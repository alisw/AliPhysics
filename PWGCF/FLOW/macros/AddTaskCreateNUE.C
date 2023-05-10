AliAnalysisTaskCreateNUE* AddTaskCreateNUE(
    TString fPeriod         ="LHC17",
    Int_t   fSystFlag       =0,
    Bool_t  fUseHM          =true,
    TString	uniqueID        = "Default"
)
{
        TString name            = "MyNUEtask";
		    Double_t	fMinPt			= 0.2;
		    Double_t	fMaxPt			= 3.0;
        Bool_t fUseCuts          =true;
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
    //fileName += ":NUE";      // create a subfolder in the file
    // now we create an instance of your task
    AliAnalysisTaskCreateNUE* task = new AliAnalysisTaskCreateNUE(name.Data());   
    if(!task) return 0x0;
    task->SetPeriod(fPeriod);
    task->SetSystFlag(fSystFlag);
    task->SetMinPt(fMinPt);
	  task->SetMaxPt(fMaxPt);
    task->SetUseHM(fUseHM);
    task->SetUseCuts(fUseCuts);
    //Minimal bias Trigger
    //task->SelectCollisionCandidates(AliVEvent::kINT7);
    //Hign Multiplicity Trigger
    if(fUseHM){
        task->SelectCollisionCandidates(AliVEvent::kHighMultV0);
    }
    else{
        task->SelectCollisionCandidates(AliVEvent::kINT7);
    }
    // add your task to the manager
    mgr->AddTask(task);
    // your task needs input: here we connect the manager to your task
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    // same for the output
    mgr->ConnectOutput(task,1,mgr->CreateContainer(Form("Efficiency_%s", uniqueID.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    // in the end, this macro returns a pointer to your task. this will be convenient later on
    // when you will run your analysis in an analysis train on grid
    return task;
}
