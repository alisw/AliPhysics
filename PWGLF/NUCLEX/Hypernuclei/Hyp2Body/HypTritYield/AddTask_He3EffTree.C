AliAnalysisTaskHe3EffTree* AddTask_He3EffTree(TString name = "name")
{
    // get the manager via the static access member. since it's static, you don't need
    // to create an instance of the class here to call the function
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        return 0x0;
    }
    if (!mgr->GetInputEventHandler()) {
        return 0x0;
    }

    TString fileName = AliAnalysisManager::GetCommonFileName();

    AliAnalysisTaskHe3EffTree* task = new AliAnalysisTaskHe3EffTree(name.Data());   
    if(!task) return 0x0;
    task->SelectCollisionCandidates(AliVEvent::kINT7 | AliVEvent::kTRD | AliVEvent::kHighMultV0 | AliVEvent::kHighMultSPD);

    mgr->AddTask(task);

    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());

	AliAnalysisDataContainer *coutput1 =
	mgr->CreateContainer("histograms", TList::Class(),AliAnalysisManager::kOutputContainer,mgr->GetCommonFileName());
	AliAnalysisDataContainer *coutput2 =
	mgr->CreateContainer("tree", TTree::Class(),AliAnalysisManager::kOutputContainer,mgr->GetCommonFileName());
	AliAnalysisDataContainer *coutput3 =
	mgr->CreateContainer("treeGen", TTree::Class(),AliAnalysisManager::kOutputContainer,mgr->GetCommonFileName());
	mgr->ConnectOutput(task, 1, coutput1);
	mgr->ConnectOutput(task, 2, coutput2);
	mgr->ConnectOutput(task, 3, coutput3);
	return task;
}
