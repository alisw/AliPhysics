class AliAnalysisDataContainer;
AliAnalysisBGMonitorQATrain* AddTaskBGMonitorQATrain(Bool_t UseTree = kFALSE)
{
    //
    //This macro configures the task for the Beam Gas Monitoring QA
    //==============================================================================
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddTaskBGMonitorQATrain", "No analysis manager to connect to.");
        return NULL;
    }
    if (!mgr->GetInputEventHandler()) {
        return 0x0;
    }
    if (!mgr->GetInputEventHandler()) {
        ::Error("AddTaskBGMonitorQATrain", "This task requires an input event handler");
        return NULL;
    }
    // Create and configure the task
    AliAnalysisBGMonitorQATrain* taskBGMonitorQA = new AliAnalysisBGMonitorQATrain("taskBGMonitorQA");
    if(!taskBGMonitorQA) return 0x0;
    mgr->AddTask(taskBGMonitorQA);
    mgr->ConnectInput(taskBGMonitorQA,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(taskBGMonitorQA,1,mgr->CreateContainer("cOutputH", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:BeamGasMon", mgr->GetCommonFileName())));
    mgr->ConnectOutput(taskBGMonitorQA,0,mgr->CreateContainer("cOutputT", TTree::Class(), AliAnalysisManager::kOutputContainer, Form("%s:BeamGasMon", mgr->GetCommonFileName())));
    return taskBGMonitorQA;
}