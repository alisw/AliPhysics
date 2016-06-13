class AliAnalysisDataContainer;
AliAnalysisBGMonitorQATrain* AddTaskBGMonitorQATrain()
{
    //
    //This macro configures the task for the Beam Gas Monitoring QA
    //==============================================================================
    std::cout << "marker 1 for debuging" << std::endl;
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddTaskBGMonitorQATrain", "No analysis manager to connect to.");
        return NULL;
    }
    std::cout << "marker 2 for debuging" << std::endl;
    if (!mgr->GetInputEventHandler()) {
        ::Error("AddTaskBGMonitorQATrain", "This task requires an input event handler");
        return NULL;
    }
    std::cout << "marker 3 for debuging" << std::endl;
    // Create and configure the task
    AliAnalysisBGMonitorQATrain *taskBGMonitorQATrain = new AliAnalysisBGMonitorQATrain("taskBGMonitorQATrain");
    std::cout << "marker 4 for debuging" << std::endl;
    taskBGMonitorQATrain->SetDebugLevel(0);
    std::cout << "marker 5 for debuging" << std::endl;
    mgr->AddTask(taskBGMonitorQATrain);
    std::cout << "marker 6 for debuging" << std::endl;
    mgr->ConnectInput(taskBGMonitorQATrain,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(taskBGMonitorQATrain,1,mgr->CreateContainer("cOutputH", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:BeamGasMon", mgr->GetCommonFileName())));
    mgr->ConnectOutput(taskBGMonitorQATrain,0,mgr->CreateContainer("cOutputT", TTree::Class(), AliAnalysisManager::kOutputContainer, Form("%s:BeamGasMon", mgr->GetCommonFileName())));
    std::cout << "marker 7 for debuging" << std::endl;
    return taskBGMonitorQATrain;
}