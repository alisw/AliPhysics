AliAnalysisTaskLMeePureMC* AddTask_LMeePureMC() {

  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_LMeePureMC", "No analysis manager found.");
    return ;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();

  //================================================
  //========= Add task to the ANALYSIS manager =====
  //================================================
  //            find input container
  AliAnalysisTaskLMeePureMC *task=NULL;
  task= new AliAnalysisTaskLMeePureMC("LMeePureMC");

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task,1,mgr->CreateContainer("LMeePureMC", TList::Class(), AliAnalysisManager::kOutputContainer, AliAnalysisManager::GetCommonFileName()));

  return task;

}
