AliAnalysisVertexingHF* ConfigVertexingHF() {

  printf("Call to AliAnalysisVertexingHF parameters setting :\n");
  vHF = new AliAnalysisVertexingHF();
 
  //--- switch-off candidates finding (default: all on)
  //vHF->SetD0toKpiOff();
  //vHF->SetJPSItoEleOff();
  //vHF->Set3ProngOff();
  vHF->SetLikeSignOn(); // like-sign pairs and triplets
  //vHF->Set4ProngOff();
  //vHF->SetDstarOff();
  vHF->SetFindVertexForDstar(kFALSE);
  //--- secondary vertex with KF?
  //vHF->SetSecVtxWithKF();
  //  vHF->SetCascadesOff();
  vHF->SetFindVertexForCascades(kFALSE);

  //--- set cuts for single-track selection  
  //     displaced tracks
  AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts","default");
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  esdTrackCuts->SetMinNClustersITS(4);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					 AliESDtrackCuts::kAny);
  esdTrackCuts->SetMinDCAToVertexXY(0.);
  esdTrackCuts->SetPtRange(0.3,1.e10);
  AliAnalysisFilter *trkFilter = new AliAnalysisFilter("trackFilter");
  trkFilter->AddCuts(esdTrackCuts);
  vHF->SetTrackFilter(trkFilter);
  //     D* soft pion tracks
  AliESDtrackCuts *esdTrackCutsSoftPi = new AliESDtrackCuts("AliESDtrackCuts","default");
  esdTrackCutsSoftPi->SetMinNClustersITS(4);
  AliAnalysisFilter *trkFilterSoftPi = new AliAnalysisFilter("trackFilterSoftPi");
  trkFilterSoftPi->AddCuts(esdTrackCutsSoftPi);
  vHF->SetTrackFilterSoftPi(trkFilterSoftPi);
  //--- set cuts for candidates selection
  AliRDHFCutsD0toKpi *cutsD0toKpi = new AliRDHFCutsD0toKpi("CutsD0toKpi");
  Float_t cutsArrayD0toKpi[9]={0.3,999999.,1.1,0.,0.,999999.,999999.,999999.,0.};
  cutsD0toKpi->SetCuts(9,cutsArrayD0toKpi);
  cutsD0toKpi->AddTrackCuts(esdTrackCuts);
  vHF->SetCutsD0toKpi(cutsD0toKpi);
  AliRDHFCutsJpsitoee *cutsJpsitoee = new AliRDHFCutsJpsitoee("CutsJpsitoee");
  Float_t cutsArrayJpsitoee[9]={0.350,100000.,1.1,0.,0.,100000.,100000.,100000000.,-1.1};
  cutsJpsitoee->SetCuts(9,cutsArrayJpsitoee);
  cutsJpsitoee->AddTrackCuts(esdTrackCuts);
  vHF->SetCutsJpsitoee(cutsJpsitoee);
  AliRDHFCutsDplustoKpipi *cutsDplustoKpipi = new AliRDHFCutsDplustoKpipi("CutsDplustoKpipi");
  Float_t cutsArrayDplustoKpipi[12]={0.2,0.4,0.4,0.,0.,0.01,0.06,0.02,0.,0.85,0.,10000000000.};
  cutsDplustoKpipi->SetCuts(12,cutsArrayDplustoKpipi);
  cutsDplustoKpipi->AddTrackCuts(esdTrackCuts);
  vHF->SetCutsDplustoKpipi(cutsDplustoKpipi);
  AliRDHFCutsDstoKKpi *cutsDstoKKpi = new AliRDHFCutsDstoKKpi("CutsDstoKKpi");
  Float_t cutsArrayDstoKKpi[14]={0.2,0.4,0.4,0.,0.,0.005,0.06,0.,0.,0.85,0.,1000.,0.1,0.1};
  cutsDstoKKpi->SetCuts(14,cutsArrayDstoKKpi);
  cutsDstoKKpi->AddTrackCuts(esdTrackCuts);
  vHF->SetCutsDstoKKpi(cutsDstoKKpi);
  AliRDHFCutsLctopKpi *cutsLctopKpi = new AliRDHFCutsLctopKpi("CutsLctopKpi");
  Float_t cutsArrayLctopKpi[12]={0.2,0.4,0.4,0.,0.,0.01,0.06,0.02,0.,0.85,0.,10000000000.};
  cutsLctopKpi->SetCuts(12,cutsArrayLctopKpi);
  cutsLctopKpi->AddTrackCuts(esdTrackCuts);
  vHF->SetCutsLctopKpi(cutsLctopKpi);
  AliRDHFCutsD0toKpipipi *cutsD0toKpipipi = new AliRDHFCutsD0toKpipipi("CutsD0toKpipipi");
  Float_t cutsArrayD0toKpipipi[9]={0.2,0.04,0.00,0.01,0.02,0.8,0.,0.1,0.};
  cutsD0toKpipipi->SetCuts(9,cutsArrayD0toKpipipi);
  cutsD0toKpipipi->AddTrackCuts(esdTrackCuts);
  vHF->SetCutsD0toKpipipi(cutsD0toKpipipi);
  AliRDHFCutsD0toKpi *cutsD0fromDstar = new AliRDHFCutsD0toKpi("CutsD0fromDstar");
  Float_t cutsArrayD0fromDstar[9]={0.3,999999.,1.1,0.,0.,999999.,999999.,999999.,0.};
  cutsD0fromDstar->SetCuts(9,cutsArrayD0fromDstar);
  cutsD0fromDstar->AddTrackCuts(esdTrackCuts);
  vHF->SetCutsD0fromDstar(cutsD0fromDstar);
  AliRDHFCutsLctoV0 *cutsLctoV0 = new AliRDHFCutsLctoV0("CutsLctoV0");
  Float_t cutsArrayLctoV0[9]={4.0,4.0,2.0,2.0,0.0,0.0,0.0,1000.,1000.};
  cutsLctoV0->SetCuts(9,cutsArrayLctoV0);
  cutsLctoV0->AddTrackCuts(esdTrackCuts);
  vHF->SetCutsLctoV0(cutsLctoV0);
  // 
  // to be removed:
  vHF->SetD0toKpiCuts(0.3,999999.,1.1,0.,0.,999999.,999999.,999999.,0.);
  vHF->SetBtoJPSICuts(0.350,999999.,1.1,0.,0.,999999.,999999.,999999.,0.);
  vHF->SetDplusCuts(0.2,0.4,0.4,0.,0.,0.01,0.06,0.02,0.,0.85,0,1e6);
  vHF->SetDsCuts(0.2,0.4,0.4,0.,0.,0.005,0.06,0.,0.,0.85,0.,0.1,0.1,0.1);
  vHF->SetLcCuts(0.2,0.4,0.4,0.,0.,0.01,0.06,0.02,0.,0.85,0,1e6);
  vHF->SetD0to4ProngsCuts(0.2,0.04,0.00,0.01,0.02,0.8,0.,0.1,0.);
  vHF->SetDstarCuts(0.3, 0.1, 0.05, 100000000000.0, 0.5);
  vHF->SetD0fromDstarCuts(0.3,999999.,1.1,0.,0.,999999.,999999.,999999.,0.);
  vHF->SetLctoV0Cuts(4.0,4.0,2.0,2.0,0.0,0.0,0.0,1000,1000);
  //--- set this if you want to reconstruct primary vertex candidate by
  //    candidate using other tracks in the event (for pp, broad 
  //    interaction region)
  //vHF->SetRecoPrimVtxSkippingTrks();
  //--- OR set this if you want to remove the candidate daughters from 
  //    the primary vertex, without recostructing it from scratch
  //vHF->SetRmTrksFromPrimVtx();

  //--- check the settings
  vHF->PrintStatus();
  //--- verbose
  //  AliLog::SetClassDebugLevel("AliAnalysisVertexingHF",2);

 
  return vHF;
}


