AliVEventHandler* AddMCHandler(Bool_t readTrackRef = kFALSE)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) 
  {
    ::Error("AddMCHandler", "No analysis manager to connect to.");
    return NULL;
  }

  AliMCEventHandler* handler = new AliMCEventHandler();
  mgr->SetMCtruthEventHandler(handler);
  handler->SetReadTR(readTrackRef);

  return handler;
}
