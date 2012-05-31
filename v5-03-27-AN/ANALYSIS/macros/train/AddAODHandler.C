AliAODInputHandler* AddAODHandler()
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) 
  {
    ::Error("AddAODHandler", "No analysis manager to connect to.");
    return NULL;
  }
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();

  AliAODInputHandler* handler = new AliAODInputHandler();
  
  if (inputHandler && (inputHandler->IsA() == AliMultiInputEventHandler::Class())) {
    AliMultiInputEventHandler *multiInputHandler=(AliMultiInputEventHandler*)inputHandler;
    multiInputHandler->AddInputEventHandler(handler);
  } else {
    if (!inputHandler) {
      mgr->SetInputEventHandler(handler);
    } else {
      ::Error("AddAODHandler", "inputHandler is NOT null. AOD handler was NOT added !!!");
    }
  }
  
  return handler;
}
