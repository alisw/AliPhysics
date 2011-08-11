AliVEventHandler* AddESDHandler()
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) 
  {
    ::Error("AddESDHandler", "No analysis manager to connect to.");
    return NULL;
  }

  AliESDInputHandler *esdHandler = new AliESDInputHandler();
  mgr->SetInputEventHandler(esdHandler);
  
  return esdHandler;
}
