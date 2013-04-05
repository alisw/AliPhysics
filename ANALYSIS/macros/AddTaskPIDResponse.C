AliAnalysisTask *AddTaskPIDResponse(Bool_t isMC=kFALSE, Bool_t autoMCesd=kTRUE,
                                    Bool_t tuneOnData=kFALSE, Int_t recoPass=2,
                                    Bool_t cachePID=kFALSE, TString detResponse="",
                                    Bool_t useTPCEtaCorrection = kTRUE)
{
// Macro to connect a centrality selection task to an existing analysis manager.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPIDResponse", "No analysis manager to connect to.");
    return 0x0;
  }

  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();

  //case of multi input event handler (needed for mixing)
  if (inputHandler->IsA() == AliMultiInputEventHandler::Class()) {
    printf("========================================================================================\n");
    printf("PIDResponse: AliMultiInputEventHandler detected, initialising AliPIDResponseInputHandler\n");
    printf("========================================================================================\n");
    AliMultiInputEventHandler *multiInputHandler=(AliMultiInputEventHandler*)inputHandler;
    
    AliPIDResponseInputHandler *pidResponseIH = new AliPIDResponseInputHandler();
    multiInputHandler->AddInputEventHandler(pidResponseIH);

    if (autoMCesd &&
        multiInputHandler->GetFirstInputEventHandler()->IsA()==AliESDInputHandler::Class() &&
        multiInputHandler->GetFirstMCEventHandler()
       ) isMC=kTRUE;
    pidResponseIH->SetIsMC(isMC);

    return 0x0;
  }

  // standard with task
  printf("========================================================================================\n");
  printf("PIDResponse: Initialising AliAnalysisTaskPIDResponse\n");
  printf("========================================================================================\n");
  
  if ( autoMCesd && (inputHandler->IsA() == AliESDInputHandler::Class()) ) {
    isMC=mgr->GetMCtruthEventHandler()!=0x0;
  }

  AliAnalysisTaskPIDResponse *pidTask = new AliAnalysisTaskPIDResponse("PIDResponseTask");
//   pidTask->SelectCollisionCandidates(AliVEvent::kMB);
  pidTask->SetIsMC(isMC);
  if(isMC&&tuneOnData) {
    pidTask->SetTuneOnData(kTRUE,recoPass);
    // tuning on MC is by default active on TPC and TOF, to enable it only on one of them use:
    // pidTask->SetTuneOnDataMask(AliPIDResponse::kDetTPC);   
    // pidTask->SetTuneOnDataMask(AliPIDResponse::kDetTOF);   
  }
  pidTask->SetCachePID(cachePID);
  pidTask->SetSpecialDetectorResponse(detResponse);
  pidTask->SetUseTPCEtaCorrection(useTPCEtaCorrection);
  mgr->AddTask(pidTask);
  
//   AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("PIDResponseQA",
//     TList::Class(), AliAnalysisManager::kOutputContainer,
//     "PIDResponseQA.root");
  
  mgr->ConnectInput(pidTask, 0, mgr->GetCommonInputContainer());
//   mgr->ConnectOutput(pidTask,1,coutput1);
  
  return pidTask;
}
