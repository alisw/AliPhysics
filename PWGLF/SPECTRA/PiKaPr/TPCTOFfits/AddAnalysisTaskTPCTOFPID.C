AliAnalysisTaskTPCTOFPID *
AddAnalysisTaskTPCTOFPID(Bool_t mcFlag = kFALSE, Bool_t mcTuneFlag = kFALSE, Bool_t pbpbFlag = kFALSE)
{

  /* check analysis manager */
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddAnalysisTaskTPCTOFPID", "cannot get analysis manager");
    return NULL;
  }

  /* check input event handler */
  if (!mgr->GetInputEventHandler()) {
    Error("AddAnalysisTaskTPCTOFPID", "cannot get input event handler");
    return NULL;
  }
  
  /* check input data type */
  TString str = mgr->GetInputEventHandler()->GetDataType();
  if (str.CompareTo("ESD")) {
    Error("AddAnalysisTaskTPCTOFPID", "input data type is not \"ESD\"");
    return NULL;
  }

  /* check MC truth event handler */
  if (mcFlag) {
    if (!mgr->GetMCtruthEventHandler()) {
      Error("AddAnalysisTaskTPCTOFPID", "cannot get MC truth event handler");
      return NULL;
    }
    AliMCEventHandler * handler = (AliMCEventHandler *) mgr->GetMCtruthEventHandler();
    handler->SetReadTR(kTRUE);
  }
  
  /* get common input data container */
  AliAnalysisDataContainer *inputc = mgr->GetCommonInputContainer();
  if (!inputc) {
    Error("AddAnalysisTaskTPCTOFPID", "cannot get common input container");
    return NULL;
  }
  

  /*  create task and connect input/output */
  AliAnalysisTaskTPCTOFPID *task = new AliAnalysisTaskTPCTOFPID(mcFlag);
  mgr->ConnectInput(task, 0, inputc);
  AliAnalysisDataContainer *outcont = mgr->CreateContainer("PIDTree",TTree::Class(), AliAnalysisManager::kOutputContainer,AliAnalysisManager::GetCommonFileName());
  mgr->ConnectOutput(task,1,outcont);
  AliAnalysisDataContainer *outcont2 = mgr->CreateContainer("StatHist",TH1D::Class(), AliAnalysisManager::kOutputContainer,AliAnalysisManager::GetCommonFileName());
  mgr->ConnectOutput(task,2,outcont2);
  /* setup task */
  task->SetMCFlag(mcFlag);
  task->SetMCTuneFlag(mcTuneFlag);
  task->SetPbPbFlag(pbpbFlag);
  task->SelectCollisionCandidates(AliVEvent::kAny);
  task->SetVertexSelectionFlag(kTRUE);
  task->SetVertexCut(15.0);
  task->SetRapidityCut(1.0);
  /* setup TOF calib */
  task->GetTOFcalib()->SetRemoveMeanT0(!mcFlag);
  task->GetTOFcalib()->SetCalibrateTOFsignal(!mcFlag);
  task->GetTOFcalib()->SetCorrectTExp(kFALSE);
  /* setup resolution */
  Double_t timeReso = 85.;
  if (mcFlag && !mcTuneFlag) timeReso = 80.;
  task->SetTimeResolution(timeReso);
  task->GetESDpid()->GetTOFResponse().SetTimeResolution(timeReso);
  task->GetTOFT0maker()->SetTimeResolution(timeReso);
  
  /* return task */
  mgr->AddTask(task);
  return task;
  
}
