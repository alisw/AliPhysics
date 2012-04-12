class AliAnalysisDataContainer;
class AliFlowEventCuts;
class AliFlowEventTrackCuts;
class AliFlowEventTrackSimpleCuts;
class AliRDHFCutsD0toKpi;
class AliRDHFCutsDStartoKpipi; 

void AddTaskFlowD2H(TString fileNameCuts, TString folderName, Int_t nDmeson, Int_t myHarmonic, 
		    Bool_t bDoQC, Bool_t bDoSPTPC, Bool_t bDoSPVZERO, Bool_t bDoEPTPC, Bool_t bDoEPVZERO, 
		    Int_t ptBinWidth, Double_t gapTPC, Double_t etaVZERO1, Double_t etaVZERO2, Double_t etaVZERO3, Double_t etaVZERO4,
		    Bool_t bOldApproach=kFALSE, Bool_t shrinkSP=kFALSE ) {
  TFile *filecuts = TFile::Open( fileNameCuts.Data() );
  if( (!filecuts) || ( filecuts && !filecuts->IsOpen()) ){
    AliFatal("Could not open cuts file.");
  }

  TString fileName = AliAnalysisManager::GetCommonFileName();
  fileName.ReplaceAll(".root","");

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();

  TString thecuts = folderName;
  TString sgn = DMesonName( nDmeson ) + Form("w%d",ptBinWidth);

  //**********************************************************************
  // FLOW TRACK CUTS
  AliFlowTrackCuts* cutsRFPTPC = new AliFlowTrackCuts( "GlobalRFPTPC" );
  cutsRFPTPC->SetParamType(AliFlowTrackCuts::kGlobal);
  cutsRFPTPC->SetPtRange(0.2,5.);
  cutsRFPTPC->SetEtaRange(-0.8,0.8);
  cutsRFPTPC->SetMinNClustersTPC(70);
  cutsRFPTPC->SetMinChi2PerClusterTPC(0.2);
  cutsRFPTPC->SetMaxChi2PerClusterTPC(4.0);
  cutsRFPTPC->SetAcceptKinkDaughters(kFALSE);
  cutsRFPTPC->SetMinimalTPCdedx(10.);
  cutsRFPTPC->SetAODfilterBit(1);
  //  cutsRFPTPC->SetQA(kTRUE);
  AliFlowTrackCuts* cutsRFPVZE = new AliFlowTrackCuts( "GlobalRFPVZE" );
  cutsRFPVZE->SetParamType(AliFlowTrackCuts::kV0);
  cutsRFPVZE->SetEtaRange(-10,+10);
  cutsRFPVZE->SetPhiMin(0);
  cutsRFPVZE->SetPhiMax(TMath::TwoPi());
  //  cutsRFPVZE->SetQA(kTRUE);
  AliFlowTrackCuts* cutsPOI = new AliFlowTrackCuts( "GlobalPOI" );
  cutsPOI->SetParamType(AliFlowTrackCuts::kGlobal);
  cutsPOI->SetPtRange(0.2,-5.);
  cutsPOI->SetEtaRange(0.8,-0.8);

  //**********************************************************************
  // POI FILTER CUTS
  AliFlowTrackSimpleCuts *filterPOIQC[50]; // MASS BANDS
  AliFlowTrackSimpleCuts *filterPOISP[50][2]; // MASS BANDS || ETA || SUBEVENT GAP
  int myMassBands = MassBands( nDmeson, bOldApproach );
  for(int mb=0; mb!=myMassBands; ++mb) {
    filterPOISP[mb][0] = new AliFlowTrackSimpleCuts( Form("FilterPOISP_MB%d_ETANEG",mb) );
    filterPOISP[mb][0]->SetEtaMin( -0.8 );
    filterPOISP[mb][0]->SetEtaMax( 0. );
    filterPOISP[mb][0]->SetMassMin( MassBandLowEdge(nDmeson,mb,bOldApproach) );
    filterPOISP[mb][0]->SetMassMax( MassBandLowEdge(nDmeson,mb+1,bOldApproach) );
    filterPOISP[mb][1] = new AliFlowTrackSimpleCuts( Form("FilterPOISP_MB%d_ETAPOS",mb) );
    filterPOISP[mb][1]->SetEtaMin( 0. );
    filterPOISP[mb][1]->SetEtaMax( +0.8 );
    filterPOISP[mb][1]->SetMassMin( MassBandLowEdge(nDmeson,mb,bOldApproach) );
    filterPOISP[mb][1]->SetMassMax( MassBandLowEdge(nDmeson,mb+1,bOldApproach) );
    filterPOIQC[mb] = new AliFlowTrackSimpleCuts( Form("FilterPOIQC_MB%d",mb) );
    filterPOIQC[mb]->SetEtaMin( -0.8 );
    filterPOIQC[mb]->SetEtaMax( +0.8 );
    filterPOIQC[mb]->SetMassMin( MassBandLowEdge(nDmeson,mb,bOldApproach) );
    filterPOIQC[mb]->SetMassMax( MassBandLowEdge(nDmeson,mb+1,bOldApproach) );
  }

  // * DMESON SELECTOR *************************************************
  AliAnalysisTaskFlowD2H *taskSel;
  switch (nDmeson) {
  case ( AliRDHFCuts::kD0toKpiCuts ):
    AliRDHFCutsD0toKpi *myCutsD0 = (AliRDHFCutsD0toKpi*)filecuts->Get("D0toKpiCuts");
    if(!myCutsD0) {
      AliFatal("Problems reaching D0toKpiCuts");
    }
    taskSel = new AliAnalysisTaskFlowD2H( Form("TaskD0Selector_%s",thecuts.Data()), 
					  cutsRFPTPC, cutsRFPVZE, myCutsD0, nDmeson );
    break;
  case ( AliRDHFCuts::kDstarCuts ):
    AliRDHFCutsDStartoKpipi *myCutsDStar = (AliRDHFCutsDStartoKpipi*)filecuts->Get("DStartoKpipiCuts");
    if(!myCutsDStar) {
      AliFatal("Problems reaching DStarToKpipiCuts");
    }
    taskSel = new AliAnalysisTaskFlowD2H( Form("TaskDStarSelector_%s",thecuts.Data()),
					  cutsRFPTPC, cutsRFPVZE, myCutsDStar, nDmeson);
    break;
  case (AliRDHFCuts::kDplusCuts):
    AliRDHFCutsDplustoKpipi *myCutsDplus = (AliRDHFCutsDplustoKpipi*)filecuts->Get("AnalysisCuts");
    if(!myCutsDplus) {
      AliFatal("Problems reaching AnalysisCuts");
    }
    taskSel = new AliAnalysisTaskFlowD2H( Form("TaskDplusSelector_%s",thecuts.Data()), 
					  cutsRFPTPC, cutsRFPVZE, myCutsDplus, nDmeson );
    break;
  }
  taskSel->SetCommonConstants( MassBins(nDmeson), MinMass(nDmeson), MaxMass(nDmeson), ptBinWidth );

  //taskSel->SelectCollisionCandidates(trigger);
  //taskSel->SetDebug();
  AliAnalysisDataContainer *coutHist = mgr->CreateContainer( Form("%sSelector_%s",sgn.Data(),thecuts.Data()),
							     TList::Class(),AliAnalysisManager::kOutputContainer,
							     Form("%s.root:FlowD2H_%s",fileName.Data(),thecuts.Data()) );
  AliAnalysisDataContainer *exc_TPC = mgr->CreateContainer( Form("TPCEventWithCandidates_%s_%s",sgn.Data(),thecuts.Data()),
							    AliFlowEventSimple::Class(),
							    AliAnalysisManager::kExchangeContainer );
  AliAnalysisDataContainer *exc_VZE = mgr->CreateContainer( Form("VZEEventWithCandidates_%s_%s",sgn.Data(),thecuts.Data()),
							    AliFlowEventSimple::Class(),
							    AliAnalysisManager::kExchangeContainer );
  mgr->AddTask(taskSel);
  mgr->ConnectOutput(taskSel,1,coutHist);
  mgr->ConnectOutput(taskSel,2,exc_TPC);
  mgr->ConnectOutput(taskSel,3,exc_VZE);
  mgr->ConnectInput (taskSel,0,cinput1);

  // * HANGING ANALYSIS TASKS ********************************************
  int harm = myHarmonic;
  for(int mb=0; mb!=myMassBands; ++mb) {
    if(bDoQC) {
      AddQCmethod( Form("%sQCTPCMB%d",sgn.Data(),mb), fileName.Data(), thecuts.Data(), harm, exc_TPC, filterPOIQC[mb]);
    }
    if(bDoSPVZERO) {
      AddSPmethod( Form("%sSPVZEMB%d",sgn.Data(),mb), fileName.Data(), thecuts.Data(), etaVZERO1, etaVZERO2, etaVZERO3, etaVZERO4,
		   "Qa", harm, exc_VZE, 0, filterPOIQC[mb], NULL, false, shrinkSP );
      AddSPmethod( Form("%sSPVZEMB%d",sgn.Data(),mb), fileName.Data(), thecuts.Data(), etaVZERO1, etaVZERO2, etaVZERO3, etaVZERO4,
		   "Qb", harm, exc_VZE, 0, filterPOIQC[mb], NULL, false, shrinkSP );
    }
    if(bDoEPVZERO) {
      AddSPmethod( Form("%sEPVZEMB%d",sgn.Data(),mb), fileName.Data(), thecuts.Data(), etaVZERO1, etaVZERO2, etaVZERO3, etaVZERO4,
		   "QaQb", harm, exc_VZE, 0, filterPOIQC[mb], NULL, true, shrinkSP );
    }
    if(bDoSPTPC) {
      for(int eg=0; eg!=2; ++eg) {
	AddSPmethod( Form("%sSPTPCMB%d",sgn.Data(),mb), fileName.Data(), thecuts.Data(), -0.8, -gapTPC*eg, +gapTPC*eg, +0.8,
		     "Qa", harm, exc_TPC, eg, filterPOISP[mb][1], NULL, false, shrinkSP );
	AddSPmethod( Form("%sSPTPCMB%d",sgn.Data(),mb), fileName.Data(), thecuts.Data(), -0.8, -gapTPC*eg, +gapTPC*eg, +0.8,
		     "Qb", harm, exc_TPC, eg, filterPOISP[mb][0], NULL, false, shrinkSP );
      }
    }
    if(bDoEPTPC) {
      AddSPmethod( Form("%sEPTPCMB%d",sgn.Data(),mb), fileName.Data(), thecuts.Data(), -0.8, -0.0, +0.0, +0.8,
		   "QaQb", harm, exc_TPC, 0, filterPOIQC[mb], NULL, true, shrinkSP );
    }
  }
}

void AddSPmethod(char *name, char *fileName, char *thecuts,
		 double minEtaA, double maxEtaA, double minEtaB, double maxEtaB,
		 char *Qvector, int harmonic, AliAnalysisDataContainer *flowEvent, int eg,
		 AliFlowTrackSimpleCuts *cutsPOI=NULL, AliFlowTrackSimpleCuts *cutsRFP=NULL,
		 bool bEP, bool shrink=false ) {
  TString myFolder = Form("v%d_%s",harmonic,thecuts);
  TString myNameSP = Form("%sSPv%d%sGAP%d_%s",name,harmonic,Qvector,eg,thecuts);
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisDataContainer *flowEvent2 = mgr->CreateContainer( Form("Filter_%s", myNameSP.Data()),
							       AliFlowEventSimple::Class(),
							       AliAnalysisManager::kExchangeContainer );
  AliAnalysisTaskFilterFE *tskFilter = new AliAnalysisTaskFilterFE( Form("TaskFilter_%s",myNameSP.Data()),
								    cutsRFP, cutsPOI);
  tskFilter->SetSubeventEtaRange( minEtaA, maxEtaA, minEtaB, maxEtaB );
  mgr->AddTask(tskFilter);
  mgr->ConnectInput( tskFilter,0,flowEvent);
  mgr->ConnectOutput(tskFilter,1,flowEvent2);
  AliAnalysisDataContainer *outSP = mgr->CreateContainer( myNameSP.Data(),
							  TList::Class(),AliAnalysisManager::kOutputContainer,
							  Form("%s.root:FlowD2H_SP_%s",fileName,myFolder.Data()) );
  AliAnalysisTaskScalarProduct *tskSP = new AliAnalysisTaskScalarProduct( Form("TaskScalarProduct_%s",
									       myNameSP.Data()),kFALSE);
  tskSP->SetApplyCorrectionForNUA(kTRUE);
  tskSP->SetHarmonic(harmonic);
  tskSP->SetTotalQvector(Qvector);
  if(bEP) tskSP->SetBehaveAsEP();
  if(shrink) tskSP->SetBookOnlyBasicCCH(kTRUE);
  mgr->AddTask(tskSP);
  mgr->ConnectInput( tskSP,0,flowEvent2);
  mgr->ConnectOutput(tskSP,1,outSP);
}

void AddQCmethod(char *name, char *fileName, char *thecuts,
		 int harmonic, AliAnalysisDataContainer *flowEvent,
		 AliFlowTrackSimpleCuts *cutsPOI=NULL, AliFlowTrackSimpleCuts *cutsRFP=NULL) {
  TString myFolder = Form("v%d_%s",harmonic,thecuts);
  TString myName = Form("%sv%d_%s",name,harmonic,thecuts);
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisDataContainer *flowEvent2 = mgr->CreateContainer( Form("Filter_%s", myName.Data()),
							       AliFlowEventSimple::Class(),
							       AliAnalysisManager::kExchangeContainer );
  AliAnalysisTaskFilterFE *tskFilter = new AliAnalysisTaskFilterFE( Form("TaskFilter_%s",myName.Data()),
								    cutsRFP, cutsPOI);
  mgr->AddTask(tskFilter);
  mgr->ConnectInput( tskFilter,0,flowEvent);
  mgr->ConnectOutput(tskFilter,1,flowEvent2);
  AliAnalysisDataContainer *outQC = mgr->CreateContainer( myName.Data(),
							  TList::Class(),AliAnalysisManager::kOutputContainer,
							  Form("%s.root:FlowD2H_QC_%s",fileName,myFolder.Data()) );
  AliAnalysisTaskQCumulants *tskQC = new AliAnalysisTaskQCumulants( Form("TaskQCumulants_%s",
									 myName.Data()),kFALSE);
  tskQC->SetApplyCorrectionForNUA(kTRUE);
  tskQC->SetHarmonic(harmonic);
  tskQC->SetBookOnlyBasicCCH(kTRUE);
  mgr->AddTask(tskQC);
  mgr->ConnectInput( tskQC,0,flowEvent2);
  mgr->ConnectOutput(tskQC,1,outQC);
}

int MassBands( int nDmeson, bool bOldApproach=false ) {
  switch (nDmeson) {
  case ( AliRDHFCuts::kDplusCuts ):
  case ( AliRDHFCuts::kD0toKpiCuts ):
    if(bOldApproach) return 5;
    else return 26;
  case ( AliRDHFCuts::kDstarCuts ):
    return 20;
  }
}

double MassBandLowEdge( int nDmeson, int mb, bool bOldApproach=false ) {
  switch (nDmeson) {
  case ( AliRDHFCuts::kDplusCuts ):
  case ( AliRDHFCuts::kD0toKpiCuts ): // 2 + 20 + 4
    double lowEdgeMinimal[5+1] = {1.75,1.80,1.83,1.90,1.93,2.03};
    double lowEdge[26+1] = { 1.66, 1.71, 1.76, 1.77, 1.78, 1.79, 1.80, 1.81, 1.82, 1.83,
			     1.84, 1.85, 1.86, 1.87, 1.88, 1.89, 1.90, 1.91, 1.92, 1.93,
			     1.94, 1.95, 1.96, 2.01, 2.06, 2.11, 2.16 };
    if(bOldApproach) return lowEdgeMinimal[mb];
    else return lowEdge[mb];
  case ( AliRDHFCuts::kDstarCuts ): // 2 + 10 + 3
    double lowEdge[20+1] = {0.1380, 0.1396, 0.1412, 0.1420, 0.1428, 0.1436, 0.1444, 0.1452, 0.1460, 0.1468, 0.1476, 0.1484, 0.1492, 0.1500,
			    0.1508, 0.1516, 0.1524, 0.1532, 0.1548, 0.1564, 0.1580};
    return lowEdge[mb];
  }
}

double MinMass( int nDmeson ) {
  switch (nDmeson) {
  case ( AliRDHFCuts::kDplusCuts ):
  case ( AliRDHFCuts::kD0toKpiCuts ):
    return 1.66;
  case ( AliRDHFCuts::kDstarCuts ):
    return 0.138;
  }
}

double MaxMass( int nDmeson ) {
  switch (nDmeson) {
  case ( AliRDHFCuts::kDplusCuts ):
  case ( AliRDHFCuts::kD0toKpiCuts ):
    return 2.16;
  case ( AliRDHFCuts::kDstarCuts ):
    return 0.158;
  }
}

int MassBins( int nDmeson ) {
  switch (nDmeson) {
  case ( AliRDHFCuts::kDplusCuts ):
  case ( AliRDHFCuts::kD0toKpiCuts ):
    return 50;
  case ( AliRDHFCuts::kDstarCuts ):
    return 25;
  }
}

TString DMesonName( int nDmeson ) {
  TString toReturn;
  switch (nDmeson) {
  case ( AliRDHFCuts::kDplusCuts ):
    toReturn = "DPlus"; break;
  case ( AliRDHFCuts::kD0toKpiCuts ):
    toReturn = "D0"; break;
  case ( AliRDHFCuts::kDstarCuts ):
    toReturn = "DStar"; break;
  }
  return toReturn;
}
