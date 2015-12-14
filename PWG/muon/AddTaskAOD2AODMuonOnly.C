///
/// Macro to keep, from an input AOD (full or already muon one), only the events that do
/// satisfy one of these two conditions :
///
/// 1 - at least one muon track
/// 2 - at least one of the muon L0 inputs (0MSL or 0MUL or 0MLL or MSH) is present
///
/// Warning : the triggerinputs bits must be set according to the trigger setup for the required
/// period. The actual numbers for those bits, following a given runlist, can be found using the
/// TriggerInputsForMuonEventCuts macro
///
/// The default values here are for PbPb2015
///
/// \author Laurent Aphecetche
///

AliAnalysisTask* AddTaskAOD2AODMuonOnly(const char* triggerinputs,Bool_t withSPDtracklets, Int_t mcMode)
{
  
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskAOD2AODMuonOnly", "No analysis manager to connect to.");
    return 0x0;
  }
  
  AliInputEventHandler* input = mgr->GetInputEventHandler();
  
  if (!input)
  {
    ::Error("AddTaskAOD2AODMuonOnly", "This task requires an input event handler");
    return 0x0;
  }
  
  TString inputDataType = input->GetDataType(); // can be "ESD" or "AOD"

  if (inputDataType != "AOD")
  {
    ::Error("AddTaskAOD2AODMuonOnly", "This task requires an AOD input event handler");
    return 0x0;
  }
  
  AliAODHandler* aodHandler = dynamic_cast<AliAODHandler*>(mgr->GetOutputEventHandler());
  if (!aodHandler)
  {
    ::Error("AddTaskAOD2AODMuonOnly", "This task requires an AOD output event handler");
    return 0x0;
  }
  
  aodHandler->SetCreateNonStandardAOD();

  AliMuonEventCuts* eventCuts = new AliMuonEventCuts("L0cutter","");
    
  eventCuts->SetTrigClassPatterns("0MSL|0MUL|0MSH|0MLL",triggerinputs);
  
  AliAnalysisTask* task = new AliAnalysisTaskAOD2MuonAOD(mcMode,withSPDtracklets,eventCuts);
    
  mgr->AddTask(task);
    
  // Connect input/output
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 0, mgr->GetCommonOutputContainer());

  return task;
}
