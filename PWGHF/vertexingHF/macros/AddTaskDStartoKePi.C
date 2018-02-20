AliAnalysisTask *AddTaskDStartoKePi(Int_t trigger = 0)
{
    cout << "----------------------------------------------------------------------------------------------------------------------------------" << endl;
    cout << " in add macro " << endl;
    //get the current analysis manager
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr)
    {
	Error("AddTask_TRDMatching", "No analysis manager found.");
	return 0;
    }

    //========= Add task to the ANALYSIS manager =====
    AliAnalysisTaskDStartoKePi *task = new AliAnalysisTaskDStartoKePi("DStarTask");
    mgr->AddTask(task); // <-

    //----------------------------------
    TString containerName = "D0_output";

    TString results = "AnalysisResults.root";

    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(containerName, TList::Class(),AliAnalysisManager::kOutputContainer,results.Data());
    AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("events", TTree::Class(),AliAnalysisManager::kOutputContainer,results.Data());

    mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task, 1, coutput1); // TList
    mgr->ConnectOutput(task, 2, coutput2); // TTree
    //----------------------------------

    return task;
}
