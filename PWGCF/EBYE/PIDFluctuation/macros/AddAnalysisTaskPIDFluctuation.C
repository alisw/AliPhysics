//__________________________________________________________________

AliAnalysisTaskPIDFluctuation *
AddAnalysisTaskPIDFluctuation(Int_t aodFilterBit, Float_t ptMin, Float_t ptMax, Float_t etaMin, Float_t etaMax)
{

  /* init analysis name */
  TString analysisName = "PIDFluctuation";
  analysisName += "_";
  analysisName += GetAODFilterBitName(aodFilterBit);
  analysisName += "_";
  analysisName += Form("pt_%.1f_%.1f", ptMin, ptMax);
  analysisName += "_";
  analysisName += Form("eta_%.1f_%.1f", etaMin, etaMax);

  /* check analysis manager */
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("", "cannot get analysis manager");
    return NULL;
  }

  /* check input event handler */
  if (!mgr->GetInputEventHandler()) {
    Error("", "cannot get input event handler");
    return NULL;
  }
  
  /* get common input data container */
  AliAnalysisDataContainer *inputc = mgr->GetCommonInputContainer();
  if (!inputc) {
    Error("", "cannot get common input container");
    return NULL;
  }
   
  /* create output data container */
  AliAnalysisDataContainer *outputc1 = mgr->CreateContainer(analysisName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, "PIDFluctuation.root");
  if (!outputc1) {
    Error("", "cannot create output container \"Histos\"");
    return NULL;
  }

  /*  create task and connect input/output */
  AliAnalysisTaskPIDFluctuation *task = new AliAnalysisTaskPIDFluctuation(analysisName.Data());
  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, inputc);
  mgr->ConnectOutput(task, 1, outputc1);

  /* setup task */
  task->SetESDtrackCuts(GetESDtrackCuts(aodFilterBit));
  task->SetAODfilterBit(aodFilterBit);
  task->SetEtaRange(etaMin, etaMax);
  task->SetPtRange(ptMin, ptMax);

  task->Dump();
  return task;
}

//__________________________________________________________________

AliESDtrackCuts *
GetESDtrackCuts(Int_t type)
{
  AliESDtrackCuts *trackCuts;
  switch (type) {
  case AliAODTrack::kTrkGlobal:
    trackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010();
    break;
  case AliAODTrack::kTrkTPCOnlyConstrained:
  case AliAODTrack::kTrkTPCOnly:
    trackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
    trackCuts->SetMinNClustersTPC(70);
    break;
  }
  return trackCuts;
}

//__________________________________________________________________

Char_t *
GetAODFilterBitName(Int_t type)
{
  switch (type) {
  case AliAODTrack::kTrkGlobal:
    return "TrkGlobal";
    break;
  case AliAODTrack::kTrkTPCOnlyConstrained:
    return "TrkTPCOnlyConstrained";
    break;
  case AliAODTrack::kTrkTPCOnly:
    return "TrkTPCOnly";
    break;
  }
}
