AliTOFAnalysisTaskCalibPass0 *
AddTOFAnalysisTaskCalibPass0()
{

  /* check analysis manager */
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddAnalysisTaskEventTime", "cannot get analysis manager");
    return NULL;
  }

  /* check input event handler */
  if (!mgr->GetInputEventHandler()) {
    Error("AddAnalysisTaskEventTime", "cannot get input event handler");
    return NULL;
  }
  
  /* check input data type */
  TString str = mgr->GetInputEventHandler()->GetDataType();
  if (str.CompareTo("ESD")) {
    Error("AddAnalysisTaskEventTime", "input data type is not \"ESD\"");
    return NULL;
  }

  /* get common input data container */
  AliAnalysisDataContainer *inputc = mgr->GetCommonInputContainer();
  if (!inputc) {
    Error("AddAnalysisTaskEventTime", "cannot get common input container");
    return NULL;
  }
  
  /* create output data container */
  AliAnalysisDataContainer *outputc1 = mgr->CreateContainer("TOFHistos", TList::Class(), AliAnalysisManager::kOutputContainer, "AliESDfriends_v1.root");
  if (!outputc1) {
    Error("AddAnalysisTaskEventTime", "cannot create output container \"Histos\"");
    return NULL;
  }

  /*  create task and connect input/output */
  AliTOFAnalysisTaskCalibPass0 *task = new AliTOFAnalysisTaskCalibPass0();

  // adding the task
  mgr->AddTask(task);

  mgr->ConnectInput(task, 0, inputc);
  mgr->ConnectOutput(task, 1, outputc1);

  /* setup task */
  task->SetEventSelectionFlag(kTRUE);
  task->SetVertexSelectionFlag(kTRUE);
  task->SetVertexCut(25.0);
  /* setup TOF calib */
  task->GetTOFcalib()->SetRemoveMeanT0(kFALSE);
  task->GetTOFcalib()->SetCalibrateTOFsignal(kTRUE);
  task->GetTOFcalib()->SetCorrectTExp(kFALSE);
  /* setup track cuts */
  AliESDtrackCuts *trackCuts = task->GetTrackCuts();
  trackCuts->SetPtRange(0.5, 10.);
  trackCuts->SetEtaRange(-1.0, 1.0);
  trackCuts->SetRequireITSRefit(kTRUE);
  trackCuts->SetMinNClustersITS(1);
  trackCuts->SetRequireTPCRefit(kTRUE);
  trackCuts->SetMinNClustersTPC(70);
  trackCuts->SetMaxChi2PerClusterTPC(4.);
  trackCuts->SetAcceptKinkDaughters(kFALSE);
  trackCuts->SetMaxDCAToVertexZ(3.2);
  trackCuts->SetMaxDCAToVertexXY(2.4);
  trackCuts->SetDCAToVertex2D(kTRUE);

  /* return task */
  return task;

}
