AliMCEventHandler* AddMCHandler(Bool_t readTrackRef = kFALSE)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) 
  {
    ::Error("AddMCHandler", "No analysis manager to connect to.");
    return NULL;
  }

  AliMCEventHandler* handler = new AliMCEventHandler();
  handler->SetReadTR(readTrackRef);
  
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();
  if (inputHandler && (inputHandler->IsA() == AliMultiInputEventHandler::Class())) {
    AliMultiInputEventHandler *multiInputHandler=(AliMultiInputEventHandler*)inputHandler;
    multiInputHandler->AddInputEventHandler(handler);
  } else {
    mgr->SetMCtruthEventHandler(handler);
  }
  
  return handler;
}
