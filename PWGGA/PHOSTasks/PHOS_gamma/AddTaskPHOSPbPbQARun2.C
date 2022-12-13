AliAnalysisTaskPHOSPbPbQARun2* AddTaskPHOSPbPbQARun2( const char* filename="AnalysisResults_PHOSPbPbQA.root",
					              const char* containername=NULL)
{
  //Add PHOS PbPb QA task to the PWGPP QA train.
  //See PHOSPbPb.C how to run it locally or standalone.

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPHOSPbPbQA", "No analysis manager to connect to");
    return NULL;
  }
  
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskPHOSPbPbQA", "This task requires an input event handler");
    return NULL;
  }

  AliAnalysisTaskPHOSPbPbQARun2* task = new AliAnalysisTaskPHOSPbPbQARun2("PbPbQA");
  mgr->AddTask(task);

  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

  if(!containername) containername = "PHOSPbPbQAResults";
  if(!filename)     filename =  mgr->GetCommonFileName();

  mgr->ConnectOutput(task, 1, mgr->CreateContainer(containername,TList::Class(), AliAnalysisManager::kOutputContainer, filename));
  return task;
}
