AliESDtrackCuts* SetupESDcuts();
void AddCFContainers(AliTRDcheckESD* checkESD);

void AddTRDcheckESD()
{
  //AliLog::SetClassDebugLevel("AliTRDcheckESD", 5);
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr) return;
  AliTRDcheckESD *checkESD = new AliTRDcheckESD((char*)"TRDcheckESD");
  checkESD->SetRefTrackFilter(SetupESDcuts());
  mgr->AddTask(checkESD);
  Bool_t mc = mgr->GetMCtruthEventHandler();
  checkESD->SetMC(mc);
  checkESD->SetCollision(/*kFALSE*/);
  checkESD->SetDebugLevel(0);
  checkESD->AddUserTrigger("WU");
  checkESD->AddUserTrigger("QU");
  checkESD->AddUserTrigger("SE");
  checkESD->AddUserTrigger("JT");
  
  AddCFContainers(checkESD);
  
  mgr->ConnectInput(checkESD,  0, mgr->GetCommonInputContainer());  
  mgr->ConnectOutput(checkESD, 1, mgr->CreateContainer(checkESD->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("%s:TRD_Performance",mgr->GetCommonFileName())));
}

AliESDtrackCuts* SetupESDcuts() {
  // Setup ESD cuts for the TPC reference tracks
  AliESDtrackCuts* esdCuts = new AliESDtrackCuts;
  esdCuts->SetPtRange(0.5, 100.0);
  esdCuts->SetEtaRange(-0.9, +0.9);
  esdCuts->SetRequireTPCRefit(kTRUE);
  esdCuts->SetAcceptKinkDaughters(kFALSE);
  esdCuts->SetMaxDCAToVertexXY(1.);
  esdCuts->SetMaxDCAToVertexZ(3.);
  esdCuts->SetMinNClustersTPC(70);
  esdCuts->SetRequireITSRefit(kTRUE);
  esdCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  return esdCuts;
}


void AddCFContainers(AliTRDcheckESD* checkESD) {
  //
  // configure output CF containers
  //
  Double_t bcLimits[3501]; for(Int_t i=0; i<=3500; ++i) bcLimits[i] = -0.5+Double_t(i);
  //Double_t phiLimits[181]; for(Int_t i=0; i<=180; ++i) phiLimits[i] = 0.0+i*2.0*TMath::Pi()/180.0;
  //Double_t phiParamLimits[181]; for(Int_t i=0; i<=180; ++i) phiParamLimits[i] = -1.0*TMath::Pi()+i*2.0*TMath::Pi()/180.0;
  Double_t phiParamLimits[73]; for(Int_t i=0; i<=72; ++i) phiParamLimits[i] = -1.0*TMath::Pi()+i*2.0*TMath::Pi()/72.0;
  Double_t ptLimits[19] = {0.0, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 
                           3.5, 4.0, 4.5,  5.0, 6.0,  7.0, 8.0,  9.0, 10.0};
  Double_t chargeLimits[3] = {-1.5,0.0,1.5};
  //Double_t etaLimits[91]; for(Int_t i=0;i<=90;++i) etaLimits[i] = -0.9+i*1.8/90.;
  Double_t etaLimits[43]; for(Int_t i=0;i<=42;++i) etaLimits[i] = -0.9+i*1.8/42.;
  Double_t trdTrkltsLimits[8] = {-0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5};
  Double_t trigLimits[AliTRDcheckESD::kNMaxAssignedTriggers+1]; for(Int_t i=0;i<=AliTRDcheckESD::kNMaxAssignedTriggers;++i) trigLimits[i] = 0.5+Double_t(i);
  Double_t evMultLimits[6] = {0.0, 700., 1400., 2100., 2800., 3500.};
  Double_t pLimits[18] = {0.0, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.7, 2.0, 2.5, 
                          3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 9.0, 12.0};
  //Double_t qtotLimits[101]; for(Int_t i=0;i<=100;++i) qtotLimits[i] = 0.0+i*10000./100.;
  Double_t qtotLimits[51]; for(Int_t i=0;i<=50;++i) qtotLimits[i] = 0.0+i*10000./50.;
  Double_t layerLimits[7] = {-0.5,0.5,1.5,2.5,3.5,4.5,5.5};
  Double_t sliceLimits[9] = {-0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5};
  //Double_t pLossLimits[101]; for(Int_t i=0;i<=100;++i) pLossLimits[i] = -500.+i*1000.0/100.;
  Double_t pLossLimits[51]; for(Int_t i=0;i<=50;++i) pLossLimits[i] = -500.+i*1000.0/50.;
  Double_t trdClustersLimits[201]; for(Int_t i=0;i<=200;++i) trdClustersLimits[i] = -0.5+i*200.0/200.;
  //Double_t trdQualityLimits[101]; for(Int_t i=0;i<=100;++i) trdQualityLimits[i] = 0.0+i*2.0/100.;
  //Double_t trdBudgetLimits[101]; for(Int_t i=0;i<=100;++i) trdBudgetLimits[i] = 0.0+i*100.0/100.;
  Double_t trdBudgetLimits[51]; for(Int_t i=0;i<=50;++i) trdBudgetLimits[i] = 0.0+i*100.0/50.;
  //Double_t trdChi2Limits[101]; for(Int_t i=0;i<=100;++i) trdChi2Limits[i] = 0.0+i*10.0/100.;
  Double_t trdChi2Limits[26]; for(Int_t i=0;i<=25;++i) trdChi2Limits[i] = 0.0+i*10.0/25.;
  Double_t trackletClustersLimits[31]; for(Int_t i=0;i<=30;++i) trackletClustersLimits[i] = 0.0+i*1.0;
  Double_t trackletClsRowsLimits[31]; for(Int_t i=0;i<=30;++i) trackletClsRowsLimits[i] = 0.0+Double_t(i)*1.0/30.0;
  //Double_t outerParamRadLimits[51]; for(Int_t i=0;i<=50;++i) outerParamRadLimits[i] = 0.0+Double_t(i)*500.0/50.0;
    
  // BC container
  printf("Container 1 :: BunchCrossingsCF\n");
  Int_t bcSteps[3] = {AliTRDcheckESD::kTPCreference, AliTRDcheckESD::kTRD, AliTRDcheckESD::kTOF};
  UInt_t bcVars[3] = {AliTRDcheckESD::kEventBC, AliTRDcheckESD::kTrackPhiTRD, AliTRDcheckESD::kTrackPt};
  TArrayD bcBinning[3] = {TArrayD(3501, bcLimits), TArrayD(73,phiParamLimits), TArrayD(19,ptLimits)};
  checkESD->AddCFContainer("BunchCrossingsCF", "Matching vs bunch crossings", 3, bcSteps, 3, bcVars, bcBinning);
  
  // (phi,eta) matching container
  printf("Container 2 :: MatchingPhiEta\n");
  Int_t matchPhiEtaSteps[3] = {AliTRDcheckESD::kTPCreference, AliTRDcheckESD::kTRD, AliTRDcheckESD::kTOF};
  UInt_t matchPhiEtaVars[5] = {AliTRDcheckESD::kTrackCharge,       AliTRDcheckESD::kTrackPhiTRD,    AliTRDcheckESD::kTrackEtaTRD,
                               AliTRDcheckESD::kTrackTrdTracklets, AliTRDcheckESD::kEventTrigger};
  TArrayD matchPhiEtaBinning[5] = {TArrayD(3,chargeLimits),    TArrayD(73, phiParamLimits), TArrayD(43, etaLimits),
                                   TArrayD(8,trdTrkltsLimits), TArrayD(AliTRDcheckESD::kNMaxAssignedTriggers+1, trigLimits)};
  checkESD->AddCFContainer("MatchingPhiEta", "(phi,eta) matching CF", 3, matchPhiEtaSteps, 5, matchPhiEtaVars, matchPhiEtaBinning);
  
  // pt matching container
  printf("Container 3 :: MatchingPt\n");
  Int_t matchPtSteps[3] = {AliTRDcheckESD::kTPCreference, AliTRDcheckESD::kTRD, AliTRDcheckESD::kTOF};
  UInt_t matchPtVars[6] = {AliTRDcheckESD::kEventMult, AliTRDcheckESD::kTrackCharge,       AliTRDcheckESD::kTrackPhiTRD,
                           AliTRDcheckESD::kTrackPt,   AliTRDcheckESD::kTrackTrdTracklets, AliTRDcheckESD::kEventTrigger};
  TArrayD matchPtBinning[6] = {TArrayD(6,evMultLimits), TArrayD(3, chargeLimits),    TArrayD(73, phiParamLimits),
                               TArrayD(19, ptLimits),   TArrayD(8, trdTrkltsLimits), TArrayD(AliTRDcheckESD::kNMaxAssignedTriggers+1, trigLimits)};
  checkESD->AddCFContainer("MatchingPt", "(pt,phi,multiplicity) matching", 3, matchPtSteps, 6, matchPtVars, matchPtBinning);
  
  // qtot, vs event multiplicity and p
  printf("Container 4 :: CentralityCF\n");
  Int_t qtotCentSteps[1] = {AliTRDcheckESD::kTRD};
  UInt_t qtotCentVars[3] = {AliTRDcheckESD::kEventMult, AliTRDcheckESD::kTrackP,
                            AliTRDcheckESD::kTrackletQtot};
  TArrayD qtotCentBinning[3] = {TArrayD(6,evMultLimits),
                                TArrayD(18, pLimits), TArrayD(51, qtotLimits)};
  checkESD->AddCFContainer("CentralityCF", "qtot vs event multiplicity and p", 1, qtotCentSteps, 3, qtotCentVars, qtotCentBinning);
  
  // qtot, vs event multiplicity and p
  printf("Container 5 :: ClustersCF\n");
  Int_t clsCentSteps[1] = {AliTRDcheckESD::kTRD};
  UInt_t clsCentVars[3] = {AliTRDcheckESD::kEventMult, AliTRDcheckESD::kTrackP,
                           AliTRDcheckESD::kTrackTrdClusters};
  TArrayD clsCentBinning[3] = {TArrayD(6,evMultLimits),
                               TArrayD(18, pLimits), TArrayD(201, trdClustersLimits)};
  checkESD->AddCFContainer("ClustersCF", "clusters/track vs event multiplicity and p", 1, clsCentSteps, 3, clsCentVars, clsCentBinning);
  
  // qtot vs (phi,eta) and layer
  printf("Container 6 :: QtotCF\n");
  Int_t qtotSteps[1] = {AliTRDcheckESD::kTRD};
  UInt_t qtotVars[4] = {AliTRDcheckESD::kTrackPhiTRD, AliTRDcheckESD::kTrackEtaTRD, 
                        AliTRDcheckESD::kTrackletQtot, AliTRDcheckESD::kTrackletLayer};
  TArrayD qtotBinning[4] = {TArrayD(73, phiParamLimits), TArrayD(43, etaLimits),
                            TArrayD(51, qtotLimits), TArrayD(7, layerLimits)};
  checkESD->AddCFContainer("QtotCF", "qtot vs (phi,eta) and layer", 1, qtotSteps, 4, qtotVars, qtotBinning);
  
  // pulse height vs p and slice
  printf("Container 7 :: PulseHeightCF\n");
  Int_t phSteps[1] = {AliTRDcheckESD::kTRD};
  UInt_t phVars[4] = {AliTRDcheckESD::kEventTrigger, AliTRDcheckESD::kTrackletP, 
                      AliTRDcheckESD::kTrackletPHslice, AliTRDcheckESD::kTrackletSlice};
  TArrayD phBinning[4] = {TArrayD(AliTRDcheckESD::kNMaxAssignedTriggers+1, trigLimits), TArrayD(18,pLimits),
                          TArrayD(51, qtotLimits), TArrayD(9,sliceLimits)};
  checkESD->AddCFContainer("PulseHeightCF", "PH vs p and slice", 1, phSteps, 4, phVars, phBinning);  
    
  // TRD quality
  /*printf("Container 8 :: TRD quality\n");
  Int_t trdQualitySteps[1] = {AliTRDcheckESD::kTRD};
  UInt_t trdQualityVars[5] = {AliTRDcheckESD::kTrackP, AliTRDcheckESD::kTrackPhiTRD, AliTRDcheckESD::kTrackEtaTRD,
                              AliTRDcheckESD::kTrackTrdQuality, AliTRDcheckESD::kTrackCharge};
  TArrayD trdQualityBinning[5] = {TArrayD(18, pLimits), TArrayD(181, phiParamLimits), TArrayD(91, etaLimits),
                                  TArrayD(101, trdQualityLimits), TArrayD(3, chargeLimits)};
  checkESD->AddCFContainer("trdQuality", "TRD quality vs (p,phi,eta,charge)", 1, trdQualitySteps, 
			   5, trdQualityVars, trdQualityBinning);
  */
  // TRD chi2
  printf("Container 9 :: TRD chi2\n");
  Int_t trdChi2Steps[1] = {AliTRDcheckESD::kTRD};
  UInt_t trdChi2Vars[5] = {AliTRDcheckESD::kTrackPhiTRD, AliTRDcheckESD::kTrackEtaTRD,
                           AliTRDcheckESD::kTrackTrdChi2,
                           AliTRDcheckESD::kTrackPt, AliTRDcheckESD::kTrackCharge};
  TArrayD trdChi2Binning[5] = {TArrayD(73, phiParamLimits), TArrayD(43, etaLimits),
                               TArrayD(26, trdChi2Limits),
                               TArrayD(19, ptLimits), TArrayD(3, chargeLimits)};
  checkESD->AddCFContainer("trdChi2", "TRD chi2 vs (phi,eta,pt,charge)", 1, trdChi2Steps, 
			   5, trdChi2Vars, trdChi2Binning);
  
  // TRD material budget
  printf("Container 10 :: TRD material budget\n");
  Int_t trdBudgetSteps[1] = {AliTRDcheckESD::kTRD};
  UInt_t trdBudgetVars[3] = {AliTRDcheckESD::kTrackPhiTRD, AliTRDcheckESD::kTrackEtaTRD,
                             AliTRDcheckESD::kTrackTRDBudget};
  TArrayD trdBudgetBinning[3] = {TArrayD(73, phiParamLimits), TArrayD(43, etaLimits),
                                 TArrayD(51, trdBudgetLimits)};
  checkESD->AddCFContainer("trdBudget", "TRD budget vs (phi,eta)", 1, trdBudgetSteps, 
			   3, trdBudgetVars, trdBudgetBinning);
  
  // ploss vs (phi,eta) and layer
  printf("Container 11 :: Momentum loss\n");
  Int_t plossSteps[1] = {AliTRDcheckESD::kTRD};
  UInt_t plossVars[5] = {AliTRDcheckESD::kTrackPhiTRD, AliTRDcheckESD::kTrackEtaTRD, 
                         AliTRDcheckESD::kTrackPlossTRDlayer, AliTRDcheckESD::kTrackletLayer,
                         AliTRDcheckESD::kTrackCharge};
  TArrayD plossBinning[5] = {TArrayD(73, phiParamLimits), TArrayD(43, etaLimits),
                             TArrayD(51, pLossLimits), TArrayD(7, layerLimits),
                             TArrayD(3, chargeLimits)};
  checkESD->AddCFContainer("Ploss", "p loss vs (phi,eta,layer,charge)", 1, plossSteps, 5, plossVars, plossBinning);

  // clusters per tracklet  
  printf("Container 12 :: Clusters per tracklet\n");
  Int_t clustersSteps[1] = {AliTRDcheckESD::kTRD};
  UInt_t clustersVars[5] = {AliTRDcheckESD::kTrackPhiTRD, AliTRDcheckESD::kTrackEtaTRD, 
                            AliTRDcheckESD::kTrackletLayer, AliTRDcheckESD::kTrackletClusters,
                            AliTRDcheckESD::kTrackCharge};
  TArrayD clustersBinning[5] = {TArrayD(73, phiParamLimits), TArrayD(43, etaLimits),
                                TArrayD(7, layerLimits), TArrayD(31, trackletClustersLimits),
                                TArrayD(3,chargeLimits)};
  checkESD->AddCFContainer("clustersPerTracklet", "tracklet clusters vs (phi,eta,layer,charge)", 1, clustersSteps, 5, clustersVars, clustersBinning);

  // clusters/crossed rows  
  printf("Container 13 :: Clusters/crossed rows\n");
  Int_t clsRowsSteps[1] = {AliTRDcheckESD::kTRD};
  UInt_t clsRowsVars[5] = {AliTRDcheckESD::kTrackPhiTRD, AliTRDcheckESD::kTrackEtaTRD, 
                            AliTRDcheckESD::kTrackletLayer, AliTRDcheckESD::kTrackletClustersVsRows,
                            AliTRDcheckESD::kTrackCharge};
  TArrayD clsRowsBinning[5] = {TArrayD(73, phiParamLimits), TArrayD(43, etaLimits),
                                TArrayD(7, layerLimits), TArrayD(31, trackletClsRowsLimits),
                                TArrayD(3,chargeLimits)};
  checkESD->AddCFContainer("clustersVsRows", "tracklet/rows vs (phi,eta,layer,charge)", 1, clsRowsSteps, 5, clsRowsVars, clsRowsBinning);
/*  
  // outer param radius  
  printf("Container 13 :: Outer param\n");
  Int_t outerParamSteps[1] = {AliTRDcheckESD::kTPCreference};
  UInt_t outerParamVars[5] = {AliTRDcheckESD::kTrackPhiTRD, AliTRDcheckESD::kTrackEtaTRD, 
                              AliTRDcheckESD::kTrackOuterParamRadius, AliTRDcheckESD::kTrackTrdTracklets,
                              AliTRDcheckESD::kTrackPt};
  TArrayD outerParamBinning[5] = {TArrayD(181, phiParamLimits), TArrayD(91, etaLimits),
                                  TArrayD(51, outerParamRadLimits), TArrayD(8,trdTrkltsLimits),
                                  TArrayD(19, ptLimits)};
  checkESD->AddCFContainer("outerParam", "outer param radius vs (phi,eta,pt,n-tracklets)", 1, outerParamSteps, 5, outerParamVars, outerParamBinning);
*/
}

