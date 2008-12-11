AliAnalysisVertexingHF* ConfigVertexingHF() {

  printf("Call to AliAnalysisVertexingHF parameters setting :\n");
  vHF = new AliAnalysisVertexingHF();
 
  //--- switch-off candidates finding (default: all on)
  //vHF->SetD0toKpiOff();
  //vHF->SetJPSItoEleOff();
  //vHF->Set3ProngOff();
  vHF->Set4ProngOff();
  //--- secondary vertex with KF?
  //vHF->SetSecVtxWithKF();
  //--- set cuts for single-track selection
  
  AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts","default");
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  esdTrackCuts->SetMinNClustersITS(5);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
   					AliESDtrackCuts::kBoth);
  esdTrackCuts->SetMinDCAToVertexXY(0.);
  esdTrackCuts->SetPtRange(0.3,1.e10);
  AliAnalysisFilter *trkFilter = new AliAnalysisFilter("trackFilter");
  trkFilter->AddCuts(esdTrackCuts);
  vHF->SetTrackFilter(trkFilter);
  //--- set cuts for candidates selection
  vHF->SetD0toKpiCuts(0.7,999999.,1.1,0.,0.,999999.,999999.,999999.,0.);
  vHF->SetBtoJPSICuts(0.350);
  vHF->SetDplusCuts(0.2,0.,0.,0.,0.,0.01,0.06,0.,0.,0.8);
  vHF->SetDsCuts(0.2,0.,0.,0.,0.,0.005,0.06,0.,0.,0.8,0.,0.1,0.1);
  vHF->SetLcCuts(0.2,0.,0.,0.,0.,0.01,0.06,0.,0.,0.8);
  //--- set this if you want to reconstruct primary vertex candidate by
  //    candidate using other tracks in the event (for pp, broad 
  //    interaction region)
  vHF->SetRecoPrimVtxSkippingTrks();
  //--- OR set this if you want to remove the candidate daughters from 
  //    the primary vertex, without recostructing it from scratch
  //vHF->SetRmTrksFromPrimVtx();

  //--- check the settings
  vHF->PrintStatus();
  //--- verbose
  //AliLog::SetClassDebugLevel("AliAnalysisVertexingHF",1);

 
  return vHF;
}



