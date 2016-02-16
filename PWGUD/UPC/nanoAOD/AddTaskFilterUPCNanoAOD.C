///
/// Macro to keep, from an input AOD, only the events that do
/// satisfy couple of conditions
///
///
/// \author Michal Broz
///

AliAnalysisTask* AddTaskFilterUPCNanoAOD(Bool_t withSPDtracklets)
{
  
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskFilterUPCNanoAOD", "No analysis manager to connect to.");
    return 0x0;
  }
  
  AliInputEventHandler* input = mgr->GetInputEventHandler();
  
  if (!input)
  {
    ::Error("AddTaskFilterUPCNanoAOD", "This task requires an input event handler");
    return 0x0;
  }
  
  TString inputDataType = input->GetDataType(); // can be "ESD" or "AOD"

  if (inputDataType != "AOD")
  {
    ::Error("AddTaskFilterUPCNanoAOD", "This task requires an AOD input event handler");
    return 0x0;
  }
  
  AliAODHandler* aodHandler = dynamic_cast<AliAODHandler*>(mgr->GetOutputEventHandler());
  if (!aodHandler)
  {
    ::Error("AddTaskFilterUPCNanoAOD", "This task requires an AOD output event handler");
    return 0x0;
  }
  
  aodHandler->SetCreateNonStandardAOD();
  
  AliAnalysisTask* task = new AliAnalysisTaskFilterUPCNanoAOD(withSPDtracklets);
    
  mgr->AddTask(task);
    
  // Connect input/output
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 0, mgr->GetCommonOutputContainer());

  return task;
}
