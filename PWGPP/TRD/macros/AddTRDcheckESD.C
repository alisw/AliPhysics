
#if ! defined (__CINT__) || defined (__MAKECINT__)
#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisCuts.h"
#include "AliESDtrackCuts.h"
#include "PWGPP/TRD/AliTRDcheckESD.h"
#endif

AliESDtrackCuts* SetupESDcuts();
void AddCFContainers(AliTRDcheckESD* checkESD);

void AddTRDcheckESD(AliAnalysisManager *mgr)
{
  //AliLog::SetClassDebugLevel("AliTRDcheckESD", 5);
  //  AliInfo("aaaaaa6666666666");
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
  Double_t phiLimits[181]; for(Int_t i=0; i<=180; ++i) phiLimits[i] = 0.0+i*2.0*TMath::Pi()/180.0;
  Double_t phiParamLimits[181]; for(Int_t i=0; i<=180; ++i) phiParamLimits[i] = -1.0*TMath::Pi()+i*2.0*TMath::Pi()/180.0;
  Double_t ptLimits[19] = {0.0, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 
                           3.5, 4.0, 4.5,  5.0, 6.0,  7.0, 8.0,  9.0, 10.0};
  Double_t chargeLimits[3] = {-1.5,0.0,1.5};
  Double_t etaLimits[91]; for(Int_t i=0;i<=90;++i) etaLimits[i] = -0.9+i*1.8/90.;
  Double_t trdTrkltsLimits[8] = {-0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5};
  Double_t trigLimits[AliTRDcheckESD::kNMaxAssignedTriggers+1]; for(Int_t i=0;i<=AliTRDcheckESD::kNMaxAssignedTriggers;++i) trigLimits[i] = 0.5+Double_t(i);
  Double_t evMultLimits[6] = {0.0, 700., 1400., 2100., 2800., 3500.};
  Double_t pLimits[18] = {0.0, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.7, 2.0, 2.5, 
                          3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 9.0, 12.0};
  Double_t qtotLimits[101]; for(Int_t i=0;i<=100;++i) qtotLimits[i] = 0.0+i*10000./100.;
  Double_t layerLimits[7] = {-0.5,0.5,1.5,2.5,3.5,4.5,5.5};
  Double_t sliceLimits[9] = {-0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5};
  Double_t pLossLimits[101]; for(Int_t i=0;i<=100;++i) pLossLimits[i] = -0.5+i*1.0/100.;
  Double_t trdClustersLimits[201]; for(Int_t i=0;i<=200;++i) trdClustersLimits[i] = -0.5+i*200.0/200.;
  Double_t trdQualityLimits[101]; for(Int_t i=0;i<=100;++i) trdQualityLimits[i] = 0.0+i*2.0/100.;
  Double_t trdBudgetLimits[101]; for(Int_t i=0;i<=100;++i) trdBudgetLimits[i] = -100.0+i*200.0/100.;
  Double_t trdChi2Limits[101]; for(Int_t i=0;i<=100;++i) trdChi2Limits[i] = 0.0+i*10.0/100.;
  Double_t trackletClustersLimits[31]; for(Int_t i=0;i<=30;++i) trackletClustersLimits[i] = 0.0+i*1.0;
  Double_t trackletClsRowsLimits[31]; for(Int_t i=0;i<=30;++i) trackletClsRowsLimits[i] = 0.0+Double_t(i)*1.0/30.0;
  Double_t outerParamRadLimits[51]; for(Int_t i=0;i<=50;++i) outerParamRadLimits[i] = 0.0+Double_t(i)*500.0/50.0;
  
  // BC container
  cout << "Container 1 :: BunchCrossingsCF" << endl;
  Int_t bcSteps[3] = {AliTRDcheckESD::kTPCreference, AliTRDcheckESD::kTRD, AliTRDcheckESD::kTOF};
  UInt_t bcVars[3] = {AliTRDcheckESD::kEventBC, AliTRDcheckESD::kTrackPhiTRD, AliTRDcheckESD::kTrackPt};
  TArrayD bcBinning[3] = {TArrayD(3501, bcLimits), TArrayD(181,phiParamLimits), TArrayD(19,ptLimits)};
  checkESD->AddCFContainer("BunchCrossingsCF", "Matching vs bunch crossings", 3, bcSteps, 3, bcVars, bcBinning);
  
  // (phi,eta) matching container
  cout << "Container 2 :: MatchingPhiEta" << endl;
  Int_t matchPhiEtaSteps[3] = {AliTRDcheckESD::kTPCreference, AliTRDcheckESD::kTRD, AliTRDcheckESD::kTOF};
  UInt_t matchPhiEtaVars[5] = {AliTRDcheckESD::kTrackCharge,       AliTRDcheckESD::kTrackPhiTRD,    AliTRDcheckESD::kTrackEtaTRD,
                               AliTRDcheckESD::kTrackTrdTracklets, AliTRDcheckESD::kEventTrigger};
  TArrayD matchPhiEtaBinning[5] = {TArrayD(3,chargeLimits),    TArrayD(181, phiParamLimits), TArrayD(91, etaLimits),
                                   TArrayD(8,trdTrkltsLimits), TArrayD(AliTRDcheckESD::kNMaxAssignedTriggers+1, trigLimits)};
  checkESD->AddCFContainer("MatchingPhiEta", "(phi,eta) matching CF", 3, matchPhiEtaSteps, 5, matchPhiEtaVars, matchPhiEtaBinning);
  
  // pt matching container
  cout << "Container 3 :: MatchingPt" << endl;
  Int_t matchPtSteps[3] = {AliTRDcheckESD::kTPCreference, AliTRDcheckESD::kTRD, AliTRDcheckESD::kTOF};
  UInt_t matchPtVars[6] = {AliTRDcheckESD::kEventMult, AliTRDcheckESD::kTrackCharge,       AliTRDcheckESD::kTrackPhiTRD,
                           AliTRDcheckESD::kTrackPt,   AliTRDcheckESD::kTrackTrdTracklets, AliTRDcheckESD::kEventTrigger};
  TArrayD matchPtBinning[6] = {TArrayD(6,evMultLimits), TArrayD(3, chargeLimits),    TArrayD(181, phiParamLimits),
                               TArrayD(19, ptLimits),   TArrayD(8, trdTrkltsLimits), TArrayD(AliTRDcheckESD::kNMaxAssignedTriggers+1, trigLimits)};
  checkESD->AddCFContainer("MatchingPt", "(pt,phi,multiplicity) matching", 3, matchPtSteps, 6, matchPtVars, matchPtBinning);
  
  // qtot, vs event multiplicity and p
  cout << "Container 4 :: CentralityCF" << endl;
  Int_t qtotCentSteps[1] = {AliTRDcheckESD::kTRD};
  UInt_t qtotCentVars[4] = {AliTRDcheckESD::kEventMult, AliTRDcheckESD::kTrackP,
                            AliTRDcheckESD::kTrackletQtot, AliTRDcheckESD::kTrackTrdClusters};
  TArrayD qtotCentBinning[4] = {TArrayD(6,evMultLimits),
                                TArrayD(18, pLimits), TArrayD(101, qtotLimits), TArrayD(201, trdClustersLimits)};
  checkESD->AddCFContainer("CentralityCF", "qtot vs event multiplicity and p", 1, qtotCentSteps, 4, qtotCentVars, qtotCentBinning);
  
  // qtot vs (phi,eta) and layer
  cout << "Container 5 :: QtotCF" << endl;
  Int_t qtotSteps[1] = {AliTRDcheckESD::kTRD};
  UInt_t qtotVars[4] = {AliTRDcheckESD::kTrackPhiTRD, AliTRDcheckESD::kTrackEtaTRD, 
                        AliTRDcheckESD::kTrackletQtot, AliTRDcheckESD::kTrackletLayer};
  TArrayD qtotBinning[4] = {TArrayD(181, phiParamLimits), TArrayD(91, etaLimits),
                            TArrayD(101, qtotLimits), TArrayD(7, layerLimits)};
  checkESD->AddCFContainer("QtotCF", "qtot vs (phi,eta) and layer", 1, qtotSteps, 4, qtotVars, qtotBinning);
  
  // pulse height vs p and slice
  cout << "Container 6 :: PulseHeightCF" << endl;
  Int_t phSteps[1] = {AliTRDcheckESD::kTRD};
  UInt_t phVars[4] = {AliTRDcheckESD::kEventTrigger, AliTRDcheckESD::kTrackletP, 
                      AliTRDcheckESD::kTrackletPHslice, AliTRDcheckESD::kTrackletSlice};
  TArrayD phBinning[4] = {TArrayD(AliTRDcheckESD::kNMaxAssignedTriggers+1, trigLimits), TArrayD(18,pLimits),
                          TArrayD(101, qtotLimits), TArrayD(9,sliceLimits)};
  checkESD->AddCFContainer("PulseHeightCF", "PH vs p and slice", 1, phSteps, 4, phVars, phBinning);  
    
  // TRD quality
  cout << "Container 7 :: TRD quality" << endl;
  Int_t trdQualitySteps[1] = {AliTRDcheckESD::kTRD};
  UInt_t trdQualityVars[5] = {AliTRDcheckESD::kTrackP, AliTRDcheckESD::kTrackPhiTRD, AliTRDcheckESD::kTrackEtaTRD,
                              AliTRDcheckESD::kTrackTrdQuality, AliTRDcheckESD::kTrackCharge};
  TArrayD trdQualityBinning[5] = {TArrayD(18, pLimits), TArrayD(181, phiParamLimits), TArrayD(91, etaLimits),
                                  TArrayD(101, trdQualityLimits), TArrayD(3, chargeLimits)};
  checkESD->AddCFContainer("trdQuality", "TRD quality vs (p,phi,eta,charge)", 1, trdQualitySteps, 
			   5, trdQualityVars, trdQualityBinning);
  
  // TRD chi2
  cout << "Container 8 :: TRD chi2" << endl;
  Int_t trdChi2Steps[1] = {AliTRDcheckESD::kTRD};
  UInt_t trdChi2Vars[5] = {AliTRDcheckESD::kTrackPhiTRD, AliTRDcheckESD::kTrackEtaTRD,
                           AliTRDcheckESD::kTrackTrdChi2,
                           AliTRDcheckESD::kTrackPt, AliTRDcheckESD::kTrackCharge};
  TArrayD trdChi2Binning[5] = {TArrayD(181, phiParamLimits), TArrayD(91, etaLimits),
                               TArrayD(101, trdChi2Limits),
                               TArrayD(19, ptLimits), TArrayD(3, chargeLimits)};
  checkESD->AddCFContainer("trdChi2", "TRD chi2 vs (phi,eta,pt,charge)", 1, trdChi2Steps, 
			   5, trdChi2Vars, trdChi2Binning);
  
  // TRD material budget
  cout << "Container 9 :: TRD material budget" << endl;
  Int_t trdBudgetSteps[1] = {AliTRDcheckESD::kTRD};
  UInt_t trdBudgetVars[3] = {AliTRDcheckESD::kTrackPhiTRD, AliTRDcheckESD::kTrackEtaTRD,
                             AliTRDcheckESD::kTrackTRDBudget};
  TArrayD trdBudgetBinning[3] = {TArrayD(181, phiParamLimits), TArrayD(91, etaLimits),
                                 TArrayD(101, trdBudgetLimits)};
  checkESD->AddCFContainer("trdBudget", "TRD budget vs (phi,eta)", 1, trdBudgetSteps, 
			   3, trdBudgetVars, trdBudgetBinning);
  
  // ploss vs (phi,eta) and layer
  cout << "Container 10 :: Momentum loss" << endl;
  Int_t plossSteps[1] = {AliTRDcheckESD::kTRD};
  UInt_t plossVars[5] = {AliTRDcheckESD::kTrackPhiTRD, AliTRDcheckESD::kTrackEtaTRD, 
                         AliTRDcheckESD::kTrackPlossTRDlayer, AliTRDcheckESD::kTrackletLayer,
                         AliTRDcheckESD::kTrackCharge};
  TArrayD plossBinning[5] = {TArrayD(181, phiParamLimits), TArrayD(91, etaLimits),
                             TArrayD(101, pLossLimits), TArrayD(7, layerLimits),
                             TArrayD(3, chargeLimits)};
  checkESD->AddCFContainer("Ploss", "p loss vs (phi,eta,layer,charge)", 1, plossSteps, 5, plossVars, plossBinning);

  // clusters per tracklet  
  cout << "Container 11 :: Clusters per tracklet" << endl;
  Int_t clustersSteps[1] = {AliTRDcheckESD::kTRD};
  UInt_t clustersVars[5] = {AliTRDcheckESD::kTrackPhiTRD, AliTRDcheckESD::kTrackEtaTRD, 
                            AliTRDcheckESD::kTrackletLayer, AliTRDcheckESD::kTrackletClusters,
                            AliTRDcheckESD::kTrackCharge};
  TArrayD clustersBinning[5] = {TArrayD(181, phiParamLimits), TArrayD(91, etaLimits),
                                TArrayD(7, layerLimits), TArrayD(31, trackletClustersLimits),
                                TArrayD(3,chargeLimits)};
  checkESD->AddCFContainer("clustersPerTracklet", "tracklet clusters vs (phi,eta,layer,charge)", 1, clustersSteps, 5, clustersVars, clustersBinning);

  // clusters/crossed rows  
  cout << "Container 12 :: Clusters/crossed rows" << endl;
  Int_t clsRowsSteps[1] = {AliTRDcheckESD::kTRD};
  UInt_t clsRowsVars[5] = {AliTRDcheckESD::kTrackPhiTRD, AliTRDcheckESD::kTrackEtaTRD, 
                            AliTRDcheckESD::kTrackletLayer, AliTRDcheckESD::kTrackletClustersVsRows,
                            AliTRDcheckESD::kTrackCharge};
  TArrayD clsRowsBinning[5] = {TArrayD(181, phiParamLimits), TArrayD(91, etaLimits),
                                TArrayD(7, layerLimits), TArrayD(31, trackletClsRowsLimits),
                                TArrayD(3,chargeLimits)};
  checkESD->AddCFContainer("clustersVsRows", "tracklet/rows vs (phi,eta,layer,charge)", 1, clsRowsSteps, 5, clsRowsVars, clsRowsBinning);
/*  
  // outer param radius  
  cout << "Container 13 :: Outer param" << endl;
  Int_t outerParamSteps[1] = {AliTRDcheckESD::kTPCreference};
  UInt_t outerParamVars[5] = {AliTRDcheckESD::kTrackPhiTRD, AliTRDcheckESD::kTrackEtaTRD, 
                              AliTRDcheckESD::kTrackOuterParamRadius, AliTRDcheckESD::kTrackTrdTracklets,
                              AliTRDcheckESD::kTrackPt};
  TArrayD outerParamBinning[5] = {TArrayD(181, phiParamLimits), TArrayD(91, etaLimits),
                                  TArrayD(51, outerParamRadLimits), TArrayD(8,trdTrkltsLimits),
                                  TArrayD(19, ptLimits)};
  checkESD->AddCFContainer("outerParam", "outer param radius vs (phi,eta,pt,n-tracklets)", 1, outerParamSteps, 5, outerParamVars, outerParamBinning);
*/
  
  // TPC-TOF (phi,eta) matching container
  cout << "Container 14 :: MatchingPhiEta_TPCTOF" << endl;
  Int_t matchPhiEtaTPCTOFSteps[3] = {AliTRDcheckESD::kTPCreference, AliTRDcheckESD::kTOFin, AliTRDcheckESD::kTOFout};
  UInt_t matchPhiEtaTPCTOFVars[4] = {AliTRDcheckESD::kTrackPhiTRD, AliTRDcheckESD::kTrackEtaTRD,
                                     AliTRDcheckESD::kTrackPt, AliTRDcheckESD::kTrackCharge};
  TArrayD matchPhiEtaTPCTOFBinning[4] = {TArrayD(181, phiParamLimits), TArrayD(91, etaLimits),
                                         TArrayD(19, ptLimits), TArrayD(3,chargeLimits)};
  checkESD->AddCFContainer("MatchingPhiEta_TPCTOF", "TPC_TOF (phi,eta) matching", 3, matchPhiEtaTPCTOFSteps, 4, matchPhiEtaTPCTOFVars, matchPhiEtaTPCTOFBinning);
}

