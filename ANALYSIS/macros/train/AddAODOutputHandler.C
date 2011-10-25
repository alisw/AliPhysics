AliVEventHandler* AddAODOutputHandler()
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) 
  {
    ::Error("AddAODOutputHandler", "No analysis manager to connect to.");
    return NULL;
  }

  AliAODHandler* handler = new AliAODHandler();
  handler->SetOutputFileName("AliAOD.root");
  mgr->SetOutputEventHandler(handler);
  AliAnalysisDataContainer* cout_aod = mgr->GetCommonOutputContainer();
  cout_aod->SetSpecialOutput();
  
  return handler;
}
