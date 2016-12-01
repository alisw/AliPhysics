AliAnalysisTask *AddTaskPIDResponse(Bool_t isMC=kFALSE, Bool_t autoMCesd=kTRUE,
                                    Bool_t tuneOnData=kTRUE, TString recoPass="2",
                                    Bool_t cachePID=kFALSE, TString detResponse="",
                                    Bool_t useTPCEtaCorrection = kTRUE,/*Please use default value! Otherwise splines can be off*/
                                    Bool_t useTPCMultiplicityCorrection = kTRUE,/*Please use default value! Otherwise splines can be off*/
                                    Int_t  recoDataPass = -1)
{
// Macro to connect a centrality selection task to an existing analysis manager.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPIDResponse", "No analysis manager to connect to.");
    return 0x0;
  }

  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();

  //case of multi input event handler (needed for mixing)
  //WARNING: most probably this is not fully supported!
  if (inputHandler->IsA() == AliMultiInputEventHandler::Class()) {
    printf("========================================================================================\n");
    printf("PIDResponse: AliMultiInputEventHandler detected, initialising AliPIDResponseInputHandler\n");
    printf("========================================================================================\n");
    AliMultiInputEventHandler *multiInputHandler=(AliMultiInputEventHandler*)inputHandler;

    AliPIDResponseInputHandler *pidResponseIH = new AliPIDResponseInputHandler();
    multiInputHandler->AddInputEventHandler(pidResponseIH);

    pidResponseIH->SetIsMC(isMC);

    return 0x0;
  }

  // standard with task
  printf("========================================================================================\n");
  printf("PIDResponse: Initialising AliAnalysisTaskPIDResponse\n");

  AliAnalysisTaskPIDResponse *pidTask = new AliAnalysisTaskPIDResponse("PIDResponseTask");
//   pidTask->SelectCollisionCandidates(AliVEvent::kMB);
  pidTask->SetIsMC(isMC);
  if(isMC){
    if (tuneOnData) {
      // get numerical reco pass number
      Int_t recoPassNumber=0;
      TString recoPassName=recoPass;
      if (recoPass.IsDigit()) {
        recoPassNumber=recoPass.Atoi();
        recoPassName="";
      }
      else {
        if (recoPass.Contains("pass1") ) {
          recoPassNumber=1;
        } else if (recoPass.Contains("pass2") ) {
          recoPassNumber=2;
        } else if (recoPass.Contains("pass3") ) {
          recoPassNumber=3;
        } else if (recoPass.Contains("pass4") ) {
          recoPassNumber=4;
        } else if (recoPass.Contains("pass5") ) {
          recoPassNumber=5;
        }

      }

      printf("             Using MC with tune on data.\n");
      printf("             !!! ATTENTION ATTENTION ATTENTION !!!\n");
      printf("             You MUST make sure the reco pass set (%s: %d) corresponds to the one this MC was produced for!\n",recoPass.Data(), recoPassNumber);
      pidTask->SetTuneOnData(kTRUE,recoPassNumber, recoPassName);
      // tuning on MC is by default active on TPC and TOF, to enable it only on one of them use:
      // pidTask->SetTuneOnDataMask(AliPIDResponse::kDetTPC);
      // pidTask->SetTuneOnDataMask(AliPIDResponse::kDetTOF);
    } else {
      printf("             !!! ATTENTION ATTENTION ATTENTION !!!\n");
      printf("             You are using MC without the tune on data option.\n");
      printf("             NOTE that this is not supported any longer!.\n");
      printf("             !!! ATTENTION ATTENTION ATTENTION !!!\n");
    }
  }
  pidTask->SetCachePID(cachePID);
  pidTask->SetSpecialDetectorResponse(detResponse);
  pidTask->SetUseTPCEtaCorrection(useTPCEtaCorrection);
  pidTask->SetUseTPCMultiplicityCorrection(useTPCMultiplicityCorrection);
  pidTask->SetUserDataRecoPass(recoDataPass);
  mgr->AddTask(pidTask);

//   AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("PIDResponseQA",
//     TList::Class(), AliAnalysisManager::kOutputContainer,
//     "PIDResponseQA.root");

  mgr->ConnectInput(pidTask, 0, mgr->GetCommonInputContainer());
//   mgr->ConnectOutput(pidTask,1,coutput1);
  printf("========================================================================================\n");

  return pidTask;
}
