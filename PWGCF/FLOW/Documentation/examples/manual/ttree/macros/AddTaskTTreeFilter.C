// addtask macro for AliAnalysisTaskTTreeFilter
// author: redmer alexander bertens (rbertens@cern.ch)
//
// note that this addtask macro is not necessary to run
// the example macro in this folder, but is added for 
// clarity

AliAnalysisTaskTTreeFilter* AddTaskTTreeFilter(
        TString outfile = "myFilteredTree.root",
        UInt_t trigger = AliVEvent::kAnyINT
        ) 
{
    // check validity of framework setup 
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddTaskTTreeFilter", "No analysis manager to connect to.");
        return NULL;
    }  
    if (!mgr->GetInputEventHandler())  {
        ::Error("AddTaskTTreeFilter", "This task requires an input event handler");
        return NULL;
    }
  

    // initialize task, connect it to the manager and return it
    AliAnalysisTaskTTreeFilter* filter = new AliAnalysisTaskTTreeFilter("filter");
    // add the task to the manager
    mgr->AddTask(filter);
    // set the trigger selection
    filter->SelectCollisionCandidates(trigger);

    // get the common input container from the analysis manager
    AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  
    AliAnalysisDataContainer *coutput = mgr->CreateContainer(
        "Tree",
        TTree::Class(),
        AliAnalysisManager::kOutputContainer,
        outfile.Data());
    // connect the input data to the flow event task
    mgr->ConnectInput(filter, 0, cinput);
    // and connect the output to the flow event task
    mgr->ConnectOutput(filter, 1, coutput);

    return filter;
}
