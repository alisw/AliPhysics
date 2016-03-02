AliAnalysisTask* AddTaskAnalysisTaskExampleV()
{
  //get the analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTask", "No analysis manager to connect to.");
    return NULL;
  }  

  // check if there is some input handler
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTask", "This task requires an input event handler");
    return NULL;
  }

  //attach the actual task
  AliAnalysisTask* task = new AliAnalysisTaskExampleV("example");
  mgr->AddTask(task);

  //Create the in/out containers
  AliAnalysisDataContainer *input = mgr->GetCommonInputContainer(); //this is standard
  AliAnalysisDataContainer *output1 = mgr->CreateContainer("outputList",TList::Class(),AliAnalysisManager::kOutputContainer);

  // Connect to the input and output containers
  mgr->ConnectInput(task,0,input); 
  mgr->ConnectOutput(task,0,output1);

  return task;
}
