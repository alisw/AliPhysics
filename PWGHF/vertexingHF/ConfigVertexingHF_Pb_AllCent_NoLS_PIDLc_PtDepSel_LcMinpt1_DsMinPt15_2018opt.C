AliAnalysisVertexingHF* ConfigVertexingHF() {

  printf("Call to AliAnalysisVertexingHF parameters setting (2018 PbPb file):\n");
  AliAnalysisVertexingHF *vHF = new AliAnalysisVertexingHF();
  //Set Reduce Size dAOD
  vHF->SetMakeReducedRHF(kTRUE);

  //--- switch-off candidates finding (default: all on)
  //vHF->SetD0toKpiOff();
  vHF->SetJPSItoEleOff();
  //  vHF->Set3ProngOff();
  vHF->SetLikeSignOff(); // like-sign pairs and triplets
  vHF->SetLikeSign3prongOff();
  vHF->Set4ProngOff();
  //  vHF->SetDstarOff();
  vHF->SetFindVertexForDstar(kFALSE);
  //--- secondary vertex with KF?
  //vHF->SetSecVtxWithKF();
  //vHF->SetCascadesOff();
  vHF->SetFindVertexForCascades(kFALSE);
  vHF->SetMassCutBeforeVertexing(kTRUE); // PbPb
  vHF->SetV0TypeForCascadeVertex(AliRDHFCuts::kOnlyOfflineV0s);

   //set PID
   vHF->SetUseKaonPIDfor3Prong(kFALSE);
   vHF->SetUseProtonAndPionPIDforLambdaC();
   vHF->SetUseProtonPIDforLambdaC2V0();
   //  vHF->SetnSigmaTOFforKaonSel(5., 5.);
   // vHF->SetnSigmaTPCforKaonSel(5., 5.);
   vHF->SetnSigmaTOFforProtonSel(40.,40.);
   vHF->SetnSigmaTPCforProtonSel(5., 5.);
   vHF->SetnSigmaTPCforPionSel(40., 40.);
   vHF->SetnSigmaTOFforPionSel(40.,40.);
   vHF->SetMaxMomForTPCPid(9999999999.);
   vHF->SetUseTPCPID(kTRUE);
   vHF->SetUseTOFPID(kFALSE);
   // vHF->SetUseTPCPIDOnlyIfNoTOF(kTRUE);


  //--- set cuts for single-track selection  
  //     displaced tracks 
  AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts","default");
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetMinNClustersTPC(50);
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  //esdTrackCuts->SetMinNClustersITS(4);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					 AliESDtrackCuts::kAny);
  // |d0|>25 micron for pt<2GeV, no cut above 2
  esdTrackCuts->SetMinDCAToVertexXYPtDep("0.0025*TMath::Max(0.,(1-TMath::Floor(TMath::Abs(pt)/2.)))");
  esdTrackCuts->SetMaxDCAToVertexXY(1.);
  esdTrackCuts->SetMaxDCAToVertexZ(1.);
  esdTrackCuts->SetPtRange(0.4,1.e10);
  esdTrackCuts->SetEtaRange(-0.8,+0.8);
  AliAnalysisFilter *trkFilter = new AliAnalysisFilter("trackFilter");
  trkFilter->AddCuts(esdTrackCuts);
  vHF->SetTrackFilter(trkFilter);

  //   displaced tracks for 20% most central events 2 prongs
  AliESDtrackCuts *esdTrackCuts1 = new AliESDtrackCuts("AliESDtrackCuts","default");
  esdTrackCuts1->SetRequireTPCRefit(kTRUE);
  esdTrackCuts1->SetMinNClustersTPC(50);
  esdTrackCuts1->SetRequireITSRefit(kTRUE);
  //esdTrackCuts->SetMinNClustersITS(4);
  esdTrackCuts1->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					 AliESDtrackCuts::kAny);
  // |d0|>25 micron for pt<2GeV, no cut above 2
  esdTrackCuts1->SetMinDCAToVertexXYPtDep("0.0025*TMath::Max(0.,(1-TMath::Floor(TMath::Abs(pt)/2.)))");
  esdTrackCuts1->SetMaxDCAToVertexXY(1.);
  esdTrackCuts1->SetMaxDCAToVertexZ(1.);
  esdTrackCuts1->SetPtRange(0.5,1.e10);
  esdTrackCuts1->SetEtaRange(-0.8,+0.8);
  AliAnalysisFilter *trkFilter1 = new AliAnalysisFilter("trackFilter");
  trkFilter1->AddCuts(esdTrackCuts1);
  vHF->SetTrackFilter2prongPbCentral(20.,trkFilter1); // for centrality 0-20%

  //     displaced tracks for 20% most central events 3 prongs
  AliESDtrackCuts *esdTrackCuts2 = new AliESDtrackCuts("AliESDtrackCuts","default");
  esdTrackCuts2->SetRequireTPCRefit(kTRUE);
  esdTrackCuts2->SetMinNClustersTPC(50);
  esdTrackCuts2->SetRequireITSRefit(kTRUE);
  //esdTrackCuts->SetMinNClustersITS(4);
  esdTrackCuts2->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					 AliESDtrackCuts::kAny);
  // |d0|>60 micron for pt<2GeV, no cut above 2
  esdTrackCuts2->SetMinDCAToVertexXYPtDep("0.0060*TMath::Max(0.,(1-TMath::Floor(TMath::Abs(pt)/2.)))");
  esdTrackCuts2->SetMaxDCAToVertexXY(1.);
  esdTrackCuts2->SetMaxDCAToVertexZ(1.);
  esdTrackCuts2->SetPtRange(0.6,1.e10);
  esdTrackCuts2->SetEtaRange(-0.8,+0.8);
  AliAnalysisFilter *trkFilter2 = new AliAnalysisFilter("trackFilter");
  trkFilter2->AddCuts(esdTrackCuts2);
  vHF->SetTrackFilter3prongPbCentral(20.,trkFilter2); // for centrality 0-20%

  //     D* soft pion tracks
  AliESDtrackCuts *esdTrackCutsSoftPi = new AliESDtrackCuts("AliESDtrackCuts","default");
  esdTrackCutsSoftPi->SetRequireITSRefit(kTRUE);
  esdTrackCutsSoftPi->SetMinNClustersITS(2);
  esdTrackCutsSoftPi->SetMaxDCAToVertexXY(1.);  
  esdTrackCutsSoftPi->SetMaxDCAToVertexZ(1.);
  esdTrackCutsSoftPi->SetPtRange(0.1,1.e10);
  esdTrackCutsSoftPi->SetEtaRange(-0.8,+0.8);
  AliAnalysisFilter *trkFilterSoftPi = new AliAnalysisFilter("trackFilterSoftPi");
  trkFilterSoftPi->AddCuts(esdTrackCutsSoftPi);
  vHF->SetTrackFilterSoftPi(trkFilterSoftPi);

  //     bachelor track cuts
  AliESDtrackCuts *esdTrackCutsBach = new AliESDtrackCuts("AliESDtrackCuts","default");
  esdTrackCutsBach->SetRequireTPCRefit(kTRUE);
  esdTrackCutsBach->SetMinNClustersTPC(50);
  esdTrackCutsBach->SetRequireITSRefit(kTRUE);
  esdTrackCutsBach->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					 AliESDtrackCuts::kAny);
  esdTrackCutsBach->SetMaxDCAToVertexXY(1.);  
  esdTrackCutsBach->SetMaxDCAToVertexZ(1.);
  esdTrackCutsBach->SetPtRange(0.5,1.e10);
  esdTrackCutsBach->SetEtaRange(-0.8,+0.8);
  AliAnalysisFilter *trkFilterBach = new AliAnalysisFilter("trackFilterBachelor");
  trkFilterBach->AddCuts(esdTrackCutsBach);
  vHF->SetTrackFilterBachelor(trkFilterBach);


  //--- set cuts for candidates selection
  Int_t nptbins=2; Float_t ptlimits[2]={0.,1000000.};

  const Int_t nptbinsD0=2; 
  Float_t ptlimitsD0[nptbinsD0+1]={0.,5.,1000000.};
  Float_t** cutsArrayD0toKpi;
  cutsArrayD0toKpi=new Float_t*[11];
  for(Int_t iv=0;iv<11;iv++){
    cutsArrayD0toKpi[iv]=new Float_t[nptbinsD0];
  }
  //0-5
  cutsArrayD0toKpi[0][0]=0.25;  //D0 inv mass window
  cutsArrayD0toKpi[1][0]=999999.;
  cutsArrayD0toKpi[2][0]=1.1;
  cutsArrayD0toKpi[3][0]=0.;
  cutsArrayD0toKpi[4][0]=0.;
  cutsArrayD0toKpi[5][0]=999999.;
  cutsArrayD0toKpi[6][0]=999999.;
  cutsArrayD0toKpi[7][0]=0.;  // d0xd0
  cutsArrayD0toKpi[8][0]=0.5;
  cutsArrayD0toKpi[9][0]=-1.;
  cutsArrayD0toKpi[10][0]=0.;
  //5-inf
  cutsArrayD0toKpi[0][1]=0.4;  //D0 inv mass window
  cutsArrayD0toKpi[1][1]=999999.;
  cutsArrayD0toKpi[2][1]=1.1;
  cutsArrayD0toKpi[3][1]=0.;
  cutsArrayD0toKpi[4][1]=0.;
  cutsArrayD0toKpi[5][1]=999999.;
  cutsArrayD0toKpi[6][1]=999999.;
  cutsArrayD0toKpi[7][1]=0.0001; // d0xd0
  cutsArrayD0toKpi[8][1]=0.5;
  cutsArrayD0toKpi[9][1]=-1.;
  cutsArrayD0toKpi[10][1]=0.;

 
  AliRDHFCutsD0toKpi *cutsD0toKpi = new AliRDHFCutsD0toKpi("CutsD0toKpi");
  cutsD0toKpi->SetStandardCutsPbPb2010();
  cutsD0toKpi->SetMinPtCandidate(0.);
  cutsD0toKpi->SetUsePID(kFALSE);
  cutsD0toKpi->SetUseTrackSelectionWithFilterBits(kFALSE);
  cutsD0toKpi->SetUseSpecialCuts(kFALSE);
  cutsD0toKpi->SetUsePhysicsSelection(kFALSE);
  cutsD0toKpi->SetMaxCentrality(90.);
  cutsD0toKpi->SetMaxVtxZ(1.e6);
  cutsD0toKpi->SetTriggerClass("");
  cutsD0toKpi->SetPtBins(nptbinsD0+1,ptlimitsD0);
  cutsD0toKpi->SetCuts(11,nptbinsD0,cutsArrayD0toKpi);
  cutsD0toKpi->AddTrackCuts(esdTrackCuts);
  vHF->SetCutsD0toKpi(cutsD0toKpi);
  AliRDHFCutsJpsitoee *cutsJpsitoee = new AliRDHFCutsJpsitoee("CutsJpsitoee");
  Float_t cutsArrayJpsitoee[9]={0.350,100000.,1.1,0.,0.,100000.,100000.,100000000.,-1.1};
  cutsJpsitoee->SetCuts(9,cutsArrayJpsitoee);
  cutsJpsitoee->AddTrackCuts(esdTrackCuts);
  vHF->SetCutsJpsitoee(cutsJpsitoee);


  const Int_t nptbinsDp=2; 
  Float_t ptlimitsDp[nptbinsDp+1]={0.,4.,1000000.};
  Float_t** cutsArrayDplustoKpipi;
  cutsArrayDplustoKpipi=new Float_t*[14];
  for(Int_t iv=0;iv<14;iv++){
    cutsArrayDplustoKpipi[iv]=new Float_t[nptbinsDp];
  }
  //0-4
  cutsArrayDplustoKpipi[0][0]=0.2;  
  cutsArrayDplustoKpipi[1][0]=0.3;
  cutsArrayDplustoKpipi[2][0]=0.3;
  cutsArrayDplustoKpipi[3][0]=0.;
  cutsArrayDplustoKpipi[4][0]=0.;
  cutsArrayDplustoKpipi[5][0]=0.01;
  cutsArrayDplustoKpipi[6][0]=0.05;
  cutsArrayDplustoKpipi[7][0]=0.05;
  cutsArrayDplustoKpipi[8][0]=0.;
  cutsArrayDplustoKpipi[9][0]=0.95;
  cutsArrayDplustoKpipi[10][0]=0.;
  cutsArrayDplustoKpipi[11][0]=10000000000.;
  cutsArrayDplustoKpipi[12][0]=3.;
  cutsArrayDplustoKpipi[13][0]=0.;
  //4-inf
  cutsArrayDplustoKpipi[0][1]=0.25;  
  cutsArrayDplustoKpipi[1][1]=0.3;
  cutsArrayDplustoKpipi[2][1]=0.3;
  cutsArrayDplustoKpipi[3][1]=0.;
  cutsArrayDplustoKpipi[4][1]=0.;
  cutsArrayDplustoKpipi[5][1]=0.01;
  cutsArrayDplustoKpipi[6][1]=0.05;
  cutsArrayDplustoKpipi[7][1]=0.05;
  cutsArrayDplustoKpipi[8][1]=0.;
  cutsArrayDplustoKpipi[9][1]=0.95;
  cutsArrayDplustoKpipi[10][1]=0.;
  cutsArrayDplustoKpipi[11][1]=10000000000.;
  cutsArrayDplustoKpipi[12][1]=3.;
  cutsArrayDplustoKpipi[13][1]=0.;

  AliRDHFCutsDplustoKpipi *cutsDplustoKpipi = new AliRDHFCutsDplustoKpipi("CutsDplustoKpipi");
  cutsDplustoKpipi->SetStandardCutsPbPb2010();
  cutsDplustoKpipi->SetUseTrackSelectionWithFilterBits(kFALSE);
  cutsDplustoKpipi->SetUsePID(kFALSE);
  //  Float_t cutsArrayDplustoKpipi[14]={0.2,0.3,0.3,0.,0.,0.01,0.05,0.05,0.,0.95,0.,10000000000.,3.,0.};
  cutsDplustoKpipi->SetPtBins(nptbinsDp+1,ptlimitsDp);
  cutsDplustoKpipi->SetCuts(14,nptbinsDp,cutsArrayDplustoKpipi);
  cutsDplustoKpipi->AddTrackCuts(esdTrackCuts);
  cutsDplustoKpipi->SetMinPtCandidate(2.);
  vHF->SetCutsDplustoKpipi(cutsDplustoKpipi);

  const Int_t nptbinsDs=2; 
  Float_t ptlimitsDs[nptbinsDs+1]={0.,4.,1000000.};
  Float_t** cutsArrayDstoKKpi;
  cutsArrayDstoKKpi=new Float_t*[20];
  for(Int_t iv=0;iv<20;iv++){
    cutsArrayDstoKKpi[iv]=new Float_t[nptbinsDp];
  }
  //0-4
  cutsArrayDstoKKpi[0][0]=0.22;
  cutsArrayDstoKKpi[1][0]=0.3;
  cutsArrayDstoKKpi[2][0]=0.3;
  cutsArrayDstoKKpi[3][0]=0.;
  cutsArrayDstoKKpi[4][0]=0.;
  cutsArrayDstoKKpi[5][0]=0.  ;
  cutsArrayDstoKKpi[6][0]=0.06;
  cutsArrayDstoKKpi[7][0]=0.02;
  cutsArrayDstoKKpi[8][0]=0.;
  cutsArrayDstoKKpi[9][0]=0.94;
  cutsArrayDstoKKpi[10][0]=0.;
  cutsArrayDstoKKpi[11][0]=100000.;
  cutsArrayDstoKKpi[12][0]=0.02;
  cutsArrayDstoKKpi[13][0]=0.0001;
  cutsArrayDstoKKpi[14][0]=-1.;
  cutsArrayDstoKKpi[15][0]=1.;
  cutsArrayDstoKKpi[16][0]=0.;
  cutsArrayDstoKKpi[17][0]=0.;
  cutsArrayDstoKKpi[18][0]=0.;
  cutsArrayDstoKKpi[19][0]=-1.;
  //4-inf
  cutsArrayDstoKKpi[0][1]=0.3;
  cutsArrayDstoKKpi[1][1]=0.3;
  cutsArrayDstoKKpi[2][1]=0.3;
  cutsArrayDstoKKpi[3][1]=0.;
  cutsArrayDstoKKpi[4][1]=0.;
  cutsArrayDstoKKpi[5][1]=0.  ;
  cutsArrayDstoKKpi[6][1]=0.06;
  cutsArrayDstoKKpi[7][1]=0.02;
  cutsArrayDstoKKpi[8][1]=0.;
  cutsArrayDstoKKpi[9][1]=0.92;
  cutsArrayDstoKKpi[10][1]=0.;
  cutsArrayDstoKKpi[11][1]=100000.;
  cutsArrayDstoKKpi[12][1]=0.02;
  cutsArrayDstoKKpi[13][1]=0.0001;
  cutsArrayDstoKKpi[14][1]=-1.;
  cutsArrayDstoKKpi[15][1]=1.;
  cutsArrayDstoKKpi[16][1]=0.;
  cutsArrayDstoKKpi[17][1]=0.;
  cutsArrayDstoKKpi[18][1]=0.;
  cutsArrayDstoKKpi[19][1]=-1.;

  AliRDHFCutsDstoKKpi *cutsDstoKKpi = new AliRDHFCutsDstoKKpi("CutsDstoKKpi");
  cutsDstoKKpi->SetStandardCutsPbPb2010();
  cutsDstoKKpi->SetUseTrackSelectionWithFilterBits(kFALSE);
  cutsDstoKKpi->SetUsePID(kFALSE);
  cutsDstoKKpi->DisableK0starChannel();
  //  Float_t cutsArrayDstoKKpi[20]={0.2,0.3,0.3,0.,0.,0.,0.06,0.02,0.,0.92,0.,100000.,0.02,0.0001,-1.,1.,0.,0.,0.,-1.};
  cutsDstoKKpi->SetPtBins(nptbinsDs+1,ptlimitsDs);
  cutsDstoKKpi->SetCuts(20,nptbinsDs,cutsArrayDstoKKpi);
  cutsDstoKKpi->AddTrackCuts(esdTrackCuts);
  cutsDstoKKpi->SetMinPtCandidate(1.5);
  vHF->SetCutsDstoKKpi(cutsDstoKKpi);
  AliRDHFCutsLctopKpi *cutsLctopKpi = new AliRDHFCutsLctopKpi("CutsLctopKpi");
  cutsLctopKpi->SetStandardCutsPbPb2010();
  cutsLctopKpi->SetUseTrackSelectionWithFilterBits(kFALSE);
  cutsLctopKpi->SetUsePID(kFALSE);
  Float_t cutsArrayLctopKpi[13]={0.13,0.5,0.5,0.,0.,0.,0.06,0.,0.,0.,0.,0.05,0.5};
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


 // D* pt dependent cuts ------------------------------------------

  AliRDHFCutsDStartoKpipi *cutsDStartoKpipi = new AliRDHFCutsDStartoKpipi("CutsDStartoKpipi");
  cutsDStartoKpipi->SetUseTrackSelectionWithFilterBits(kFALSE);
  cutsDStartoKpipi->SetUsePID(kFALSE);
  const Int_t nvars=16;
  const Int_t nptbinsDst=2;
  
  Float_t* ptbins;
  ptbins=new Float_t[nptbinsDst+1];
  ptbins[0]=0.;
  ptbins[1]=5.;
  ptbins[2]=999.;
  
  cutsDStartoKpipi->SetPtBins(nptbinsDst+1,ptbins);
  
  Float_t** rdcutsvalmine;
  rdcutsvalmine=new Float_t*[nvars];
  for(Int_t iv=0;iv<nvars;iv++){
    rdcutsvalmine[iv]=new Float_t[nptbinsDst];
  }

  //0-5
  rdcutsvalmine[0][0]=0.095;  //D0 inv mass window
  rdcutsvalmine[1][0]=0.05;   // dca
  rdcutsvalmine[2][0]=0.9;    // thetastar
  rdcutsvalmine[3][0]=0.5;    // pt Pion
  rdcutsvalmine[4][0]=0.5;    // Pt Kaon
  rdcutsvalmine[5][0]=0.1;    // d0K
  rdcutsvalmine[6][0]=0.1;    // d0Pi
  rdcutsvalmine[7][0]=-0.00001; // d0xd0
  rdcutsvalmine[8][0]=0.9;    // costhetapoint
  rdcutsvalmine[9][0]=0.07;   // Dstar inv mass window
  rdcutsvalmine[10][0]=0.02;  // half width of (M_Kpipi-M_D0)
  rdcutsvalmine[11][0]=0.1;   // Pt min of Pi soft
  rdcutsvalmine[12][0]=100.;  // Pt max of pi soft
  rdcutsvalmine[13][0]=9999.; // theta
  rdcutsvalmine[14][0]=0.96;   // |cosThetaPointXY|
  rdcutsvalmine[15][0]=2.5;    // NormDecayLenghtXY
 //5-999
  rdcutsvalmine[0][1]=0.12;   //D0 inv mass window
  rdcutsvalmine[1][1]=0.06;   // dca
  rdcutsvalmine[2][1]=0.9;    // thetastar
  rdcutsvalmine[3][1]=0.5;    // pt Pion
  rdcutsvalmine[4][1]=0.5;    // Pt Kaon
  rdcutsvalmine[5][1]=0.1;    // d0K
  rdcutsvalmine[6][1]=0.1;    // d0Pi
  rdcutsvalmine[7][1]=0.0001; // d0xd0
  rdcutsvalmine[8][1]=0.7;    // costhetapoint
  rdcutsvalmine[9][1]=0.2;   // Dstar inv mass window
  rdcutsvalmine[10][1]=0.02;  // half width of (M_Kpipi-M_D0)
  rdcutsvalmine[11][1]=0.1;   // Pt min of Pi soft
  rdcutsvalmine[12][1]=100.;  // Pt max of pi soft
  rdcutsvalmine[13][1]=9999.; // theta
  rdcutsvalmine[14][1]=0.8;   // |cosThetaPointXY|
  rdcutsvalmine[15][1]=0.;    // NormDecayLenghtXY

  cutsDStartoKpipi->SetCuts(nvars,nptbinsDst,rdcutsvalmine);
 
  cutsDStartoKpipi->AddTrackCuts(esdTrackCuts);
  cutsDStartoKpipi->AddTrackCutsSoftPi(esdTrackCutsSoftPi);
  cutsDStartoKpipi->SetMinPtCandidate(1.);
  vHF->SetCutsDStartoKpipi(cutsDStartoKpipi);

  //--------------------------------------------------------

  AliRDHFCutsLctoV0 *cutsLctoV0 = new AliRDHFCutsLctoV0("CutsLctoV0");
  Float_t cutsArrayLctoV0[21]={0.2,0.,0.05,0.05,0.5,0.0,0.0,1000.,1000.,0.99,3.,1000.,0.,0.,0.,0.5,9999.,-9999.,-9999.,-9999.,1};
  cutsLctoV0->SetUseTrackSelectionWithFilterBits(kFALSE);
  cutsLctoV0->SetMinPtCandidate(1.);
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


