///////////////////////////////////////////////////////////////////
//                                                               //            
// AddTaskTrans                                                  //
//                                                               //
///////////////////////////////////////////////////////////////////
class AliAnalysisDataContainer;

AliAnalysisTaskTransTask* AddTaskTrans(TString name = "name")
{
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        return 0x0;
    }
    if (!mgr->GetInputEventHandler()) {
        return 0x0;
    }
    Bool_t isMC = kFALSE;
    if(mgr->GetMCtruthEventHandler())isMC = kTRUE;
    // by default, a file is open for writing. here, we get the filename
    TString fileName = AliAnalysisManager::GetCommonFileName();
    fileName += ":CtrueTask";      // create a subfolder in the file
    // now we create an instance of your task
    AliAnalysisTaskTransTask* task = new AliAnalysisTaskTransTask(name.Data());   
    if(!task) return 0x0;
    task->SetIsMC(isMC);
    // add your task to the manager
    mgr->AddTask(task);
    // your task needs input: here we connect the manager to your task
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    // same for the output
    mgr->ConnectOutput(task,1,
		       mgr->CreateContainer("fAnaTree", TTree::Class(),
					    AliAnalysisManager::kOutputContainer,
					    fileName.Data()));
    mgr->ConnectOutput(task,2,
		       mgr->CreateContainer("fOutputList", TList::Class(),
					    AliAnalysisManager::kOutputContainer,
					    fileName.Data()));
    // in the end, this macro returns a pointer to your task. this will be convenient later on
    // when you will run your analysis in an analysis train on grid
    return task;
}
