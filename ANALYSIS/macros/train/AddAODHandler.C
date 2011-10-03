AliVEventHandler* AddAODHandler()
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) 
  {
    ::Error("AddAODHandler", "No analysis manager to connect to.");
    return NULL;
  }

  AliAODInputHandler* handler = new AliAODInputHandler();
  mgr->SetInputEventHandler(handler);
  
  return handler;
}
