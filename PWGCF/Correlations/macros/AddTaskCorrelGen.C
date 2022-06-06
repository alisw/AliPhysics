
/// Create an instance of this class and add it to the analysis manager. Often, this function is called
/// by an AddTask C macro. However, by compiling the code, it ensures that we do not
/// have to deal with difficulties caused by CINT.
///
/// \param suffix additional suffix that can be added at the end of the task name
/// \return pointer to the new AddTaskCorrelGen task
AliAnalysisTaskCorrelGen *AddTaskCorrelGen(TString suffix="")
{
    // Get the pointer to the existing analysis manager via the static access method.
    // Since it's static, you don't need an instance of the class to call the function
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr)
    {
        ::Error("AliAnalysisTaskCorrelGen::AddTaskCorrelGen", "No analysis manager to connect to.");
        return 0;
    }

    // Get the input event handler, again via a static method.
    // This handler is part of the managing system and feeds events to your task
    AliVEventHandler *handler = mgr->GetInputEventHandler();
    if (!handler)
    {
        ::Error("AliAnalysisTaskCorrelGen::AddTaskCorrelGen", "This task requires an input event handler");
        return 0;
    }

    // Initialize the task and do some settings
    TString name("AliAnalysisTaskCorrelGen");

    if (!suffix.IsNull())
    {
        name += "_";
        name += suffix;
    }

    AliAnalysisTaskCorrelGen *correlTask = new AliAnalysisTaskCorrelGen(name);
    if (!correlTask)
        return 0;

    
    // Final settings, pass to manager and set the containers
    mgr->AddTask(correlTask);

    // Create containers for input/output
    AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
    TString contname(name);
    contname += "_histos";
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname.Data(),
                                                              TList::Class(), AliAnalysisManager::kOutputContainer,
                                                              Form("%s", AliAnalysisManager::GetCommonFileName()));
    mgr->ConnectInput(correlTask, 0, cinput1);
    mgr->ConnectOutput(correlTask, 1, coutput1);

    return correlTask;
}
