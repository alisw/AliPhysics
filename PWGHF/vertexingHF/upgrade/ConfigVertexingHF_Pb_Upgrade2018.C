AliAnalysisVertexingHF* ConfigVertexingHF() {

  printf("MYCONFIGPID Call to AliAnalysisVertexingHF parameters setting :\n");
  vHF = new AliAnalysisVertexingHF();
    //Set Reduce Size dAOD
  vHF->SetMakeReducedRHF(kTRUE);
 
  //--- switch-off candidates finding (default: all on)
  //vHF->SetD0toKpiOff();
  // vHF->SetD0toKpiOn();
  vHF->SetJPSItoEleOff();
  //vHF->Set3ProngOff();
  //vHF->SetLikeSignOn(); // like-sign pairs and triplets
  //  vHF->SetLikeSign3prongOff();
  //vHF->Set4ProngOff();
  // vHF->SetDstarOn();
  vHF->SetFindVertexForDstar(kFALSE);
  //--- secondary vertex with KF?
  //vHF->SetSecVtxWithKF();
  //Cascade
  vHF->SetCascadesOn();                 
  vHF->SetFindVertexForCascades(kFALSE);
   vHF->SetMassCutBeforeVertexing(kTRUE); // PbPb
  vHF->SetV0TypeForCascadeVertex(AliRDHFCuts::kAllV0s); //All V0s 0, Offline 1, OnTheFly 2 //as in Config_Pb_AllCent
  vHF->SetUseProtonPIDforLambdaC2V0();

  vHF->SetMassCutBeforeVertexing(kTRUE); // PbPb

 
  
//--- set cuts for single-track selection  
  //     displaced tracks
  AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts","default");
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetMinNClustersTPC(70);
 esdTrackCuts->SetRequireITSRefit(kTRUE);
  //esdTrackCuts->SetMinNClustersITS(4);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,	 AliESDtrackCuts::kAny);
  // |d0|>30 micron for pt<2GeV, no cut above 2
  esdTrackCuts->SetMinDCAToVertexXYPtDep("0.003*TMath::Max(0.,(1-TMath::Floor(TMath::Abs(pt)/2.)))");
  esdTrackCuts->SetMaxDCAToVertexXY(1.);  
  esdTrackCuts->SetMaxDCAToVertexZ(1.);
  esdTrackCuts->SetPtRange(0.4,1.e10);
  esdTrackCuts->SetEtaRange(-0.8,+0.8);
  AliAnalysisFilter *trkFilter = new AliAnalysisFilter("trackFilter");
  trkFilter->AddCuts(esdTrackCuts);
  vHF->SetTrackFilter(trkFilter);
  //     D* soft pion tracks
  AliESDtrackCuts *esdTrackCutsSoftPi = new AliESDtrackCuts("AliESDtrackCuts","default");
  esdTrackCutsSoftPi->SetMinNClustersITS(2);
  esdTrackCutsSoftPi->SetMaxDCAToVertexXY(1.);  
  esdTrackCutsSoftPi->SetMaxDCAToVertexZ(1.);
  esdTrackCutsSoftPi->SetPtRange(0.1,1.e10);
  esdTrackCutsSoftPi->SetEtaRange(-0.8,+0.8);
  AliAnalysisFilter *trkFilterSoftPi = new AliAnalysisFilter("trackFilterSoftPi");
  trkFilterSoftPi->AddCuts(esdTrackCutsSoftPi);
  vHF->SetTrackFilterSoftPi(trkFilterSoftPi);
  //--- set cuts for candidates selection
  Int_t nptbins=2; Float_t ptlimits[2]={0.,1000000.};
  AliRDHFCutsD0toKpi *cutsD0toKpi = new AliRDHFCutsD0toKpi("CutsD0toKpi");
  cutsD0toKpi->SetStandardCutsPbPb2010();
   cutsD0toKpi->SetUseCentrality(kFALSE);
  // cutsD0toKpi->SetMinCentrality(-10);
  //cutsD0toKpi->SetMaxCentrality(110);
  cutsD0toKpi->SetUseSpecialCuts(kFALSE);
  cutsD0toKpi->SetMinPtCandidate(0.);
  cutsD0toKpi->SetUsePID(kFALSE);
  cutsD0toKpi->SetUsePhysicsSelection(kFALSE);
  cutsD0toKpi->SetMaxVtxZ(1.e6);
  cutsD0toKpi->SetTriggerClass("");
  Float_t cutsArrayD0toKpi[11]={0.4,999999.,1.1,0.,0.,999999.,999999.,0.,0.5,-1,0.};
  cutsD0toKpi->SetPtBins(nptbins,ptlimits);
  cutsD0toKpi->SetCuts(11,cutsArrayD0toKpi);
  cutsD0toKpi->AddTrackCuts(esdTrackCuts);
  vHF->SetCutsD0toKpi(cutsD0toKpi);
  AliRDHFCutsJpsitoee *cutsJpsitoee = new AliRDHFCutsJpsitoee("CutsJpsitoee");
  Float_t cutsArrayJpsitoee[9]={0.350,100000.,1.1,0.,0.,100000.,100000.,100000000.,-1.1};
  cutsJpsitoee->SetCuts(9,cutsArrayJpsitoee);
  cutsJpsitoee->AddTrackCuts(esdTrackCuts);
  vHF->SetCutsJpsitoee(cutsJpsitoee);
  AliRDHFCutsDplustoKpipi *cutsDplustoKpipi = new AliRDHFCutsDplustoKpipi("CutsDplustoKpipi");
  cutsDplustoKpipi->SetStandardCutsPbPb2010();
  cutsDplustoKpipi->SetUseCentrality(kFALSE);
  cutsDplustoKpipi->SetUsePID(kFALSE);
  Float_t cutsArrayDplustoKpipi[14]={0.25,0.3,0.3,0.,0.,0.01,0.05,0.05,0.,0.7,0.,10000000000.,0.,-1.};
  cutsDplustoKpipi->SetPtBins(nptbins,ptlimits);
  cutsDplustoKpipi->SetCuts(14,cutsArrayDplustoKpipi);
  cutsDplustoKpipi->AddTrackCuts(esdTrackCuts);
  cutsDplustoKpipi->SetMinPtCandidate(2.);
  vHF->SetCutsDplustoKpipi(cutsDplustoKpipi);
  AliRDHFCutsDstoKKpi *cutsDstoKKpi = new AliRDHFCutsDstoKKpi("CutsDstoKKpi");
  cutsDstoKKpi->SetStandardCutsPbPb2010();
   cutsDstoKKpi->SetUseCentrality(kFALSE);
  cutsDstoKKpi->SetUsePID(kFALSE);
  Float_t cutsArrayDstoKKpi[20]={0.35,0.3,0.3,0.,0.,0.005,0.06,0.,0.,0.7,0.,100000.,0.035,0.0001,-1.,1.,0.,0.,0.,-1.};
  cutsDstoKKpi->SetPtBins(nptbins,ptlimits);
  cutsDstoKKpi->SetCuts(20,cutsArrayDstoKKpi);
  cutsDstoKKpi->AddTrackCuts(esdTrackCuts);
  cutsDstoKKpi->SetMinPtCandidate(2.);
  vHF->SetCutsDstoKKpi(cutsDstoKKpi);
  AliRDHFCutsLctopKpi *cutsLctopKpi = new AliRDHFCutsLctopKpi("CutsLctopKpi");
  cutsLctopKpi->SetStandardCutsPbPb2010();
  cutsLctopKpi->SetUseCentrality(kFALSE);
  cutsLctopKpi->SetUsePID(kFALSE);
  Float_t cutsArrayLctopKpi[13]={0.13,0.4,0.4,0.,0.,0.,0.06,0.,0.,0.,0.,0.05,0.4};
  cutsLctopKpi->SetPtBins(nptbins,ptlimits);
  cutsLctopKpi->SetCuts(13,cutsArrayLctopKpi);
  cutsLctopKpi->AddTrackCuts(esdTrackCuts);
  cutsLctopKpi->SetMinPtCandidate(2.);
  vHF->SetCutsLctopKpi(cutsLctopKpi);
  AliRDHFCutsD0toKpipipi *cutsD0toKpipipi = new AliRDHFCutsD0toKpipipi("CutsD0toKpipipi");
  Float_t cutsArrayD0toKpipipi[9]={0.2,0.04,0.00,0.01,0.02,0.8,0.,0.1,0.};
  cutsD0toKpipipi->SetCuts(9,cutsArrayD0toKpipipi);
  cutsD0toKpipipi->AddTrackCuts(esdTrackCuts);
  vHF->SetCutsD0toKpipipi(cutsD0toKpipipi);


 
   AliRDHFCutsDStartoKpipi *cutsDStartoKpipi = new AliRDHFCutsDStartoKpipi("CutsDStartoKpipi");
  Float_t cutsArrayDStartoKpipi[16]={0.1,999999.,1.1,0.,0.,999999.,999999.,999999.,0.7,0.03, 0.1, 0.05, 100000000000.0, 0.5,-1.,1.}; // first 9 for D0 from D*, next 5 for D*, last 2 for D0 again
  cutsDStartoKpipi->SetUsePID(kFALSE);
  cutsDStartoKpipi->SetCuts(16,cutsArrayDStartoKpipi);
  cutsDStartoKpipi->AddTrackCuts(esdTrackCuts);
  cutsDStartoKpipi->AddTrackCutsSoftPi(esdTrackCutsSoftPi);
  cutsDStartoKpipi->SetMinPtCandidate(2.);
  vHF->SetCutsDStartoKpipi(cutsDStartoKpipi);
 
  //--------------------------------------------------------

  AliRDHFCutsLctoV0 *cutsLctoV0 = new AliRDHFCutsLctoV0("CutsLctoV0");
  /*
  Float_t cutsArrayLctoV0[9]={4.0,4.0,2.0,2.0,0.0,0.0,0.0,1000.,1000.};
  cutsLctoV0->SetCuts(9,cutsArrayLctoV0);
  */
  //last dummy
  Float_t cutsArrayLctoV0[21]={1.0,1.0,0.05,0.05,0.0,0.0,0.0,1000.,1000.,0.99,3.,1000.,0.,0.,0.,0.,9999.,-9999.,-9999.,-9999.,0.0};//from the Config_Pb_AllCent_NOLS
 
  cutsLctoV0->SetCuts(21,cutsArrayLctoV0);
  cutsLctoV0->AddTrackCuts(esdTrackCuts);
    AliESDtrackCuts *esdV0daughterTrackCuts = new AliESDtrackCuts("AliESDtrackCutsForV0D","default cuts for V0 daughters");
  esdV0daughterTrackCuts->SetRequireTPCRefit(kTRUE);
  esdV0daughterTrackCuts->SetMinNClustersTPC(30);
  esdV0daughterTrackCuts->SetRequireITSRefit(kFALSE);
  esdV0daughterTrackCuts->SetMinDCAToVertexXY(0.);
  esdV0daughterTrackCuts->SetPtRange(0.05,1.e10);
  esdV0daughterTrackCuts->SetEtaRange(-1.1,+1.1);
  esdV0daughterTrackCuts->SetAcceptKinkDaughters(kTRUE);
  esdV0daughterTrackCuts->SetRequireSigmaToVertex(kFALSE);
  cutsLctoV0->AddTrackCutsV0daughters(esdV0daughterTrackCuts);
  vHF->SetCutsLctoV0(cutsLctoV0);
  // 
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
