AliMultiInputEventHandler* AddMultiHandler()
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) 
  {
    ::Error("AddMultiHandler", "No analysis manager to connect to.");
    return NULL;
  }
  
  AliMultiInputEventHandler *handler = new AliMultiInputEventHandler();
  mgr->SetInputEventHandler(handler);
  
  return handler;
}
