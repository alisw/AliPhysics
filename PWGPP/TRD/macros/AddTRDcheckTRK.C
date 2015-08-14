void AddTRDcheckTRK(AliAnalysisDataContainer **ci)
{
  Info("AddTRDcheckTRK",  "[0]=\"%s\" [1]=\"%s\" [2]=\"%s\" [3]=\"%s\"", ci[0]->GetName(), ci[1]->GetName(), ci[2]->GetName(), ci[3]->GetName());

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr) return;

  //AliLog::SetClassDebugLevel("AliTRDcheckTRK", 3);
  // global settings for tracking
  AliTRDcheckTRK::SetKalmanUpdate(kTRUE);
  AliTRDcheckTRK::SetKalmanStep(0.2);
  AliTRDcheckTRK::SetClRecalibrate(kTRUE);
  AliTRDtrackerV1::SetBetheBloch(AliTRDtrackerV1::kGeant);
/*  AliTRDtransform::SetVd(.2);
  AliTRDtransform::SetT0(.2);
  AliTRDtransform::SetExB(.2);*/
  char bb[10];
  switch(AliTRDtrackerV1::GetBetheBloch()){
  case AliTRDtrackerV1::kGeant:
    snprintf(bb, 10, "Geant"); break;
  case AliTRDtrackerV1::kSolid:
    snprintf(bb, 10, "Solid"); break;
  case AliTRDtrackerV1::kGas:
    snprintf(bb, 10, "Gas"); break;
  }
  Info("AddTRDcheckTRK", Form("Tracking settings:\n"
    "  BetheBloch    [%s]\n"
    "  KalmanUpdate  [%c]\n"
    "  KalmanStep    [%f]\n"
    "  ClRecalibrate [%c]\n"
    , bb
    , AliTRDcheckTRK::HasKalmanUpdate()?'y':'n'
    , AliTRDcheckTRK::GetKalmanStep()
    , AliTRDcheckTRK::HasClRecalibrate()?'y':'n'
  ));

  AliTRDcheckTRK *trk(NULL);;
  mgr->AddTask(trk = new AliTRDcheckTRK((char*)"TRDcheckTRK"));
  trk->SetDebugLevel(0);
  trk->SetMCdata((Bool_t)mgr->GetMCtruthEventHandler());
  trk->SetFriends(kTRUE);
  mgr->ConnectInput( trk, 0, mgr->GetCommonInputContainer()); // connect main (ESD) container
  mgr->ConnectInput( trk, 1, ci[0]);                          // connect barrel tracks container
  mgr->ConnectInput( trk, 2, ci[1]);                          // connect event info container
  mgr->ConnectInput( trk, 3, ci[2]);                          // connect onl.tracklets container
  mgr->ConnectInput( trk, 4, ci[3]);                          // connect clusters container

  mgr->ConnectOutput(trk, 1, mgr->CreateContainer(trk->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("%s:TRD_Performance",mgr->GetCommonFileName())));
}

