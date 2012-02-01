AliAnalysisTaskUnicor *AddTaskUnicor()
{
  // Creates a unicor analysis task and adds it to the analysis manager 

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {Error("AddTaskUnicor", "no analysis manager"); return NULL;}  
  if (!mgr->GetInputEventHandler()) {Error("AddTaskUnicor", "no input event handler"); return NULL;}
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if (type!="ESD") {Error("AddTaskUnicor","handler not of ESD type"); return NULL;}
  AliAnalysisTaskUnicor *mytask = new AliAnalysisTaskUnicor();
  mgr->AddTask(mytask);
  mgr->ConnectInput (mytask,0,mgr->GetCommonInputContainer());
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += ":PWG2UNICOR"; 
  AliAnalysisDataContainer *coutpt = mgr->CreateContainer("unilis", TList::Class(),
  							  AliAnalysisManager::kOutputContainer,
  							  outputfile);
  mgr->ConnectOutput(mytask,1,coutpt);
  return mytask;
}
