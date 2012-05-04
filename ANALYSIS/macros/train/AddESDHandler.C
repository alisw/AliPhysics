AliESDInputHandler* AddESDHandler()
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) 
  {
    ::Error("AddESDHandler", "No analysis manager to connect to.");
    return NULL;
  }

  AliESDInputHandler *esdHandler = new AliESDInputHandler();
  
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();
  if (inputHandler && (inputHandler->IsA() == AliMultiInputEventHandler::Class())) {
    AliMultiInputEventHandler *multiInputHandler=(AliMultiInputEventHandler*)inputHandler;
    multiInputHandler->AddInputEventHandler(esdHandler);
  } else {
    if (!inputHandler) {
      mgr->SetInputEventHandler(esdHandler);
    } else {
      ::Error("AddESDHandler", "inputHandler is NOT null. ESD handler was NOT added !!!");
      return NULL;
    }
  }
  
  return esdHandler;
}
