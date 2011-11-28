AliAnalysisVertexingHF* ConfigVertexingHF() {

  printf("Call to AliAnalysisVertexingHF parameters setting :\n");
  vHF = new AliAnalysisVertexingHF();
 
  //--- switch-off candidates finding (default: all on)
  //vHF->SetD0toKpiOff();
  vHF->SetJPSItoEleOff();
  //vHF->Set3ProngOff();
  vHF->SetLikeSignOn(); // like-sign pairs and triplets
  vHF->Set4ProngOff();
  //vHF->SetDstarOff();
  vHF->SetFindVertexForDstar(kFALSE);
  //--- secondary vertex with KF?
  //vHF->SetSecVtxWithKF();
  vHF->SetCascadesOff();
  vHF->SetFindVertexForCascades(kFALSE);
  vHF->SetMassCutBeforeVertexing(kTRUE); // PbPb

  //--- set cuts for single-track selection  
  //     displaced tracks
  AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts","default");
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetMinNClustersTPC(70);
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  //esdTrackCuts->SetMinNClustersITS(4);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					 AliESDtrackCuts::kAny);
  // |d0|>75 micron for pt<2GeV, no cut above 2
  esdTrackCuts->SetMinDCAToVertexXYPtDep("0.0075*TMath::Max(0.,(1-TMath::Floor(TMath::Abs(pt)/2.)))");
  esdTrackCuts->SetMaxDCAToVertexXY(1.);  
  esdTrackCuts->SetMaxDCAToVertexZ(1.);
  esdTrackCuts->SetPtRange(0.7,1.e10);
  esdTrackCuts->SetEtaRange(-0.8,+0.8);
  AliAnalysisFilter *trkFilter = new AliAnalysisFilter("trackFilter");
  trkFilter->AddCuts(esdTrackCuts);
  vHF->SetTrackFilter(trkFilter);
  //     D* soft pion tracks
  AliESDtrackCuts *esdTrackCutsSoftPi = new AliESDtrackCuts("AliESDtrackCuts","default");
  esdTrackCutsSoftPi->SetMinNClustersITS(4);
  esdTrackCutsSoftPi->SetMaxDCAToVertexXY(1.);  
  esdTrackCutsSoftPi->SetMaxDCAToVertexZ(1.);
  esdTrackCutsSoftPi->SetPtRange(0.2,1.e10);
  esdTrackCutsSoftPi->SetEtaRange(-0.8,+0.8);
  AliAnalysisFilter *trkFilterSoftPi = new AliAnalysisFilter("trackFilterSoftPi");
  trkFilterSoftPi->AddCuts(esdTrackCutsSoftPi);
  vHF->SetTrackFilterSoftPi(trkFilterSoftPi);
  //--- set cuts for candidates selection
  Int_t nptbins=2; Float_t ptlimits[2]={0.,1000000.};
  AliRDHFCutsD0toKpi *cutsD0toKpi = new AliRDHFCutsD0toKpi("CutsD0toKpi");
  cutsD0toKpi->SetStandardCutsPbPb2010();
  cutsD0toKpi->SetMinPtCandidate(0.);
  cutsD0toKpi->SetUsePID(kFALSE);
  cutsD0toKpi->SetUsePhysicsSelection(kFALSE);
  cutsD0toKpi->SetMaxVtxZ(1.e6);
  cutsD0toKpi->SetTriggerClass("");
  Float_t cutsArrayD0toKpi[11]={0.2,999999.,1.1,0.,0.,999999.,999999.,0.,0.5,-1,0.};
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
  cutsDplustoKpipi->SetUsePID(kFALSE);
  Float_t cutsArrayDplustoKpipi[14]={0.2,0.3,0.3,0.,0.,0.01,0.06,0.05,1.,0.85,0.,10000000000.,0.,0.};
  cutsDplustoKpipi->SetPtBins(nptbins,ptlimits);
  cutsDplustoKpipi->SetCuts(14,cutsArrayDplustoKpipi);
  cutsDplustoKpipi->AddTrackCuts(esdTrackCuts);
  cutsDplustoKpipi->SetMinPtCandidate(3.);
  vHF->SetCutsDplustoKpipi(cutsDplustoKpipi);
  AliRDHFCutsDstoKKpi *cutsDstoKKpi = new AliRDHFCutsDstoKKpi("CutsDstoKKpi");
  cutsDstoKKpi->SetStandardCutsPbPb2010();
  cutsDstoKKpi->SetUsePID(kFALSE);
  Float_t cutsArrayDstoKKpi[16]={0.2,0.4,0.4,0.,0.,0.005,0.045,0.,0.,0.9,0.,100000.,0.035,0.05,-1.,1.};
  cutsDstoKKpi->SetPtBins(nptbins,ptlimits);
  cutsDstoKKpi->SetCuts(16,cutsArrayDstoKKpi);
  cutsDstoKKpi->AddTrackCuts(esdTrackCuts);
  cutsDstoKKpi->SetMinPtCandidate(4.);
  vHF->SetCutsDstoKKpi(cutsDstoKKpi);
  AliRDHFCutsLctopKpi *cutsLctopKpi = new AliRDHFCutsLctopKpi("CutsLctopKpi");
  cutsLctopKpi->SetStandardCutsPbPb2010();
  cutsLctopKpi->SetUsePID(kFALSE);
  Float_t cutsArrayLctopKpi[13]={0.13,0.9,1.,0.,0.,0.01,0.04,0.006,1.,0.5,0.,0.05,0.4};
  cutsLctopKpi->SetPtBins(nptbins,ptlimits);
  cutsLctopKpi->SetCuts(13,cutsArrayLctopKpi);
  cutsLctopKpi->AddTrackCuts(esdTrackCuts);
  cutsLctopKpi->SetMinPtCandidate(4.);
  vHF->SetCutsLctopKpi(cutsLctopKpi);
  AliRDHFCutsD0toKpipipi *cutsD0toKpipipi = new AliRDHFCutsD0toKpipipi("CutsD0toKpipipi");
  Float_t cutsArrayD0toKpipipi[9]={0.2,0.04,0.00,0.01,0.02,0.8,0.,0.1,0.};
  cutsD0toKpipipi->SetCuts(9,cutsArrayD0toKpipipi);
  cutsD0toKpipipi->AddTrackCuts(esdTrackCuts);
  vHF->SetCutsD0toKpipipi(cutsD0toKpipipi);
  AliRDHFCutsDStartoKpipi *cutsDStartoKpipi = new AliRDHFCutsDStartoKpipi("CutsDStartoKpipi");
  cutsDStartoKpipi->SetStandardCutsPbPb2010();
  cutsDStartoKpipi->SetUsePID(kFALSE);


 // D* pt dependent cuts ------------------------------------------

  AliRDHFCutsDStartoKpipi *cutsDStartoKpipi = new AliRDHFCutsDStartoKpipi("CutsDStartoKpipi");
  
  const Int_t nvars=16;
  const Int_t nptbins=2;
  
  Float_t* ptbins;
  ptbins=new Float_t[nptbins+1];
  ptbins[0]=0.;
  ptbins[1]=4.;
  ptbins[2]=999.;
  
  cutsDStartoKpipi->SetPtBins(nptbins+1,ptbins);
  
  Float_t** rdcutsvalmine;
  rdcutsvalmine=new Float_t*[nvars];
  for(Int_t iv=0;iv<nvars;iv++){
    rdcutsvalmine[iv]=new Float_t[nptbins];
  }
  //0-4
  rdcutsvalmine[0][0]=0.15;
  rdcutsvalmine[1][0]=0.07;
  rdcutsvalmine[2][0]=0.8;
  rdcutsvalmine[3][0]=0.8;
  rdcutsvalmine[4][0]=0.8;
  rdcutsvalmine[5][0]=0.08;
  rdcutsvalmine[6][0]=0.08;
  rdcutsvalmine[7][0]=-0.00002;
  rdcutsvalmine[8][0]=0.72;
  rdcutsvalmine[9][0]=0.15;
  rdcutsvalmine[10][0]=0.03;
  rdcutsvalmine[11][0]=0.2;
  rdcutsvalmine[12][0]=100.;
  rdcutsvalmine[13][0]=0.5;
  rdcutsvalmine[14][0]=-1.;
  rdcutsvalmine[15][0]=0.;
  //4-999
  rdcutsvalmine[0][1]=0.24;
  rdcutsvalmine[1][1]=0.07;
  rdcutsvalmine[2][1]=0.8;
  rdcutsvalmine[3][1]=0.8;
  rdcutsvalmine[4][1]=0.8;
  rdcutsvalmine[5][1]=0.1;
  rdcutsvalmine[6][1]=0.1;
  rdcutsvalmine[7][1]=0.00001;
  rdcutsvalmine[8][1]=0.72;
  rdcutsvalmine[9][1]=0.15;
  rdcutsvalmine[10][1]=0.03;
  rdcutsvalmine[11][1]=0.2;
  rdcutsvalmine[12][1]=100.;
  rdcutsvalmine[13][1]=0.5;
  rdcutsvalmine[14][1]=-1.;
  rdcutsvalmine[15][1]=0.;

  cutsDStartoKpipi->SetCuts(nvars,nptbins,rdcutsvalmine);
 
  cutsDStartoKpipi->AddTrackCuts(esdTrackCuts);
  cutsDStartoKpipi->AddTrackCutsSoftPi(esdTrackCutsSoftPi);
  cutsDStartoKpipi->SetMinPtCandidate(3.);
  vHF->SetCutsDStartoKpipi(cutsDStartoKpipi);

  //--------------------------------------------------------

  AliRDHFCutsLctoV0 *cutsLctoV0 = new AliRDHFCutsLctoV0("CutsLctoV0");
  Float_t cutsArrayLctoV0[9]={4.0,4.0,2.0,2.0,0.0,0.0,0.0,1000.,1000.};
  cutsLctoV0->SetCuts(9,cutsArrayLctoV0);
  cutsLctoV0->AddTrackCuts(esdTrackCuts);
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


