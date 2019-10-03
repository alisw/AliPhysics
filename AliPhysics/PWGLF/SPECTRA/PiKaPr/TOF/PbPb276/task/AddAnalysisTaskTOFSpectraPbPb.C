AliAnalysisTaskTOFSpectraPbPb *
AddAnalysisTaskTOFSpectraPbPb(Bool_t mcFlag = kFALSE, Bool_t mcTuneFlag = kFALSE, Bool_t pbpbFlag = kFALSE)
{

  /* check analysis manager */
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddAnalysisTaskTOFSpectraPbPb", "cannot get analysis manager");
    return NULL;
  }

  /* check input event handler */
  if (!mgr->GetInputEventHandler()) {
    Error("AddAnalysisTaskTOFSpectraPbPb", "cannot get input event handler");
    return NULL;
  }
  
  /* check input data type */
  TString str = mgr->GetInputEventHandler()->GetDataType();
  if (str.CompareTo("ESD")) {
    Error("AddAnalysisTaskTOFSpectraPbPb", "input data type is not \"ESD\"");
    return NULL;
  }

  /* check MC truth event handler */
  if (mcFlag) {
    if (!mgr->GetMCtruthEventHandler()) {
      Error("AddAnalysisTaskTOFSpectraPbPb", "cannot get MC truth event handler");
      return NULL;
    }
  }
  
  /* get common input data container */
  AliAnalysisDataContainer *inputc = mgr->GetCommonInputContainer();
  if (!inputc) {
    Error("AddAnalysisTaskTOFSpectraPbPb", "cannot get common input container");
    return NULL;
  }
  
  /* setup output event handler */
  AliAODHandler *outputh = (AliAODHandler *)mgr->GetOutputEventHandler();
  outputh->SetCreateNonStandardAOD();
  outputh->SetOutputFileName("TOFSpectraPbPb.root");

  /*  create task and connect input/output */
  AliAnalysisTaskTOFSpectraPbPb *task = new AliAnalysisTaskTOFSpectraPbPb();
  mgr->ConnectInput(task, 0, inputc);

  /* setup task */
  task->SetMCFlag(mcFlag);
  task->SetMCTuneFlag(mcTuneFlag);
  task->SetPbPbFlag(pbpbFlag);
  task->SelectCollisionCandidates(AliVEvent::kMB);
  task->SetVertexSelectionFlag(kTRUE);
  task->SetVertexCut(10.0);
  task->SetRapidityCut(0.5);
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
  /* setup track cuts */
  AliESDtrackCuts *trackCuts = new AliESDtrackCuts;
  trackCuts->SetMaxChi2PerClusterTPC(4.5);
  trackCuts->SetAcceptKinkDaughters(kFALSE);
  trackCuts->SetRequireTPCRefit(kTRUE);
  trackCuts->SetRequireITSRefit(kTRUE);
  trackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					 AliESDtrackCuts::kAny);
  trackCuts->SetMaxDCAToVertexZ(2);
  trackCuts->SetDCAToVertex2D(kFALSE);
  trackCuts->SetRequireSigmaToVertex(kFALSE);
  trackCuts->SetPtRange(0.15, 10.);
  trackCuts->SetEtaRange(-0.9, 0.9);
  task->SetTrackCuts(trackCuts);
  
  /* return task */
  return task;
  
}
