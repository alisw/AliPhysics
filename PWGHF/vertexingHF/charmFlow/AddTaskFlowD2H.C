class AliAnalysisDataContainer;
class AliFlowEventCuts;
class AliFlowEventTrackCuts;
class AliRDHFCutsD0toKpi;
class AliRDHFCutsDStartoKpipi; 

void AddTaskFlowD2H(TString fileNameCuts, TString folderName, Int_t nDmeson, Int_t minCent, Int_t maxCent, Int_t myHarmonic, 
		    Bool_t bDoQC, Bool_t bDoSPTPC, Bool_t bDoSPVZERO, Bool_t bDoEPAlso, Int_t ptBinWidth,
		    Double_t gapTPC, Double_t etaVZERO1, Double_t etaVZERO2, Double_t etaVZERO3, Double_t etaVZERO4,
		    Int_t trigger, Bool_t testOnly=kFALSE ) {
  TFile *filecuts = TFile::Open( fileNameCuts.Data() );
  if( (!filecuts) || ( filecuts && !filecuts->IsOpen()) ){
    AliFatal("Could not open cuts file.");
  }
  filecuts->ls();
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
  int myMassBands = NumberOfMassBins( nDmeson );
  if(testOnly)
    myMassBands = 1;
  for(int mb=0; mb!=myMassBands; ++mb) {
    filterPOISP[mb][0] = new AliFlowTrackSimpleCuts( Form("FilterPOISP_MB%d_ETANEG",mb) );
    filterPOISP[mb][0]->SetEtaMin( -0.8 );
    filterPOISP[mb][0]->SetEtaMax( 0. );
    filterPOISP[mb][0]->SetMassMin( LowMassBin(nDmeson,mb) );
    filterPOISP[mb][0]->SetMassMax( LowMassBin(nDmeson,mb+1) );
    filterPOISP[mb][1] = new AliFlowTrackSimpleCuts( Form("FilterPOISP_MB%d_ETAPOS",mb) );
    filterPOISP[mb][1]->SetEtaMin( 0. );
    filterPOISP[mb][1]->SetEtaMax( +0.8 );
    filterPOISP[mb][1]->SetMassMin( LowMassBin(nDmeson,mb) );
    filterPOISP[mb][1]->SetMassMax( LowMassBin(nDmeson,mb+1) );
    filterPOIQC[mb] = new AliFlowTrackSimpleCuts( Form("FilterPOIQC_MB%d",mb) );
    filterPOIQC[mb]->SetEtaMin( -0.8 );
    filterPOIQC[mb]->SetEtaMax( +0.8 );
    filterPOIQC[mb]->SetMassMin( LowMassBin(nDmeson,mb) );
    filterPOIQC[mb]->SetMassMax( LowMassBin(nDmeson,mb+1) );
  }

  int minCC = minCent;
  int maxCC = maxCent;
  TString ccName = Form("cc%d%d",minCC,maxCC);
  // * EVENT CUTS ******************************************************
  AliFlowEventCuts* cutsEvent = new AliFlowEventCuts( Form("EventCuts_%s",ccName.Data()) );
  cutsEvent->SetCentralityPercentileRange(minCC,maxCC);
  cutsEvent->SetCentralityPercentileMethod(AliFlowEventCuts::kV0);
  cutsEvent->SetNContributorsRange(2);
  cutsEvent->SetPrimaryVertexZrange(-9.,9.);
  cutsEvent->SetQA(kTRUE);
  cutsEvent->SetRefMultMethod(AliESDtrackCuts::kTrackletsITSTPC);
  // * DMESON SELECTOR *************************************************
  AliAnalysisTaskFlowD2H *taskSel;
  switch (nDmeson) {
  case ( AliRDHFCuts::kD0toKpiCuts ):
    AliRDHFCutsD0toKpi *myCutsD0;
    myCutsD0 = (AliRDHFCutsD0toKpi*)filecuts->Get("D0toKpiCuts");
    if(!myCutsD0) {
      AliFatal("Problems reaching D0toKpiCuts");
    }
    taskSel = new AliAnalysisTaskFlowD2H( Form("TaskD0Selector_%s",ccName.Data()), 
					  cutsEvent, myCutsD0, nDmeson );
    break;
  case ( AliRDHFCuts::kDstarCuts ):
    AliRDHFCutsDStartoKpipi *myCutsDStar;
    myCutsDStar = (AliRDHFCutsDStartoKpipi*)filecuts->Get("DStartoKpipiCuts");
    if(!myCutsDStar) {
      AliFatal("Problems reaching DStarToKpipiCuts");
    }
    taskSel = new AliAnalysisTaskFlowD2H( Form("TaskDStarSelector_%s",ccName.Data()),
					  cutsEvent, myCutsDStar, nDmeson);
    break;
  case (AliRDHFCuts::kDplusCuts):
    AliRDHFCutsDplustoKpipi *myCutsDplus;
    myCutsDplus = (AliRDHFCutsDplustoKpipi*)filecuts->Get("AnalysisCuts");
    if(!myCutsDplus) {
      AliFatal("Problems reaching AnalysisCuts");
    }
    taskSel = new AliAnalysisTaskFlowD2H( Form("TaskDplusSelector_%s",ccName.Data()), 
					  cutsEvent, myCutsDplus, nDmeson );
    break;
  }
  taskSel->SelectCollisionCandidates(trigger);
  //taskSel->SetDebug();
  AliAnalysisDataContainer *coutHist = mgr->CreateContainer( Form("%sSelector_%s",sgn.Data(),ccName.Data()),
							     TList::Class(),AliAnalysisManager::kOutputContainer,
							     Form("%s.root:FlowD2H_%s_%s",fileName.Data(),ccName.Data(),thecuts.Data()) );
  AliAnalysisDataContainer *coutArr = mgr->CreateContainer( Form("TaskSelectorCandidates_%s_%s",ccName.Data(),sgn.Data()),
							    TObjArray::Class(),AliAnalysisManager::kExchangeContainer);
  mgr->AddTask(taskSel);
  mgr->ConnectOutput(taskSel,1,coutHist);
  mgr->ConnectOutput(taskSel,2,coutArr);
  mgr->ConnectInput (taskSel,0,cinput1);
  // * RAW FLOWEVENT TPC *************************************************
  AliAnalysisDataContainer *exc_TPC = MakeFlowEvent("TPC", ccName.Data(), fileName.Data(), thecuts.Data(), ptBinWidth, cinput1,
						    coutArr, cutsEvent, cutsRFPTPC, cutsPOI, true, nDmeson, 1.0);
  // * RAW FLOWEVENT VZERO ***********************************************
  AliAnalysisDataContainer *exc_VZE = MakeFlowEvent("VZE", ccName.Data(), fileName.Data(), thecuts.Data(), ptBinWidth, cinput1,
						    coutArr, cutsEvent, cutsRFPVZE, cutsPOI, false, nDmeson, 5.0);
  // * HANGING ANALYSIS TASKS ********************************************
  int harm = myHarmonic;
  for(int mb=0; mb!=myMassBands; ++mb) {
    if(bDoQC) {
      AddQCmethod( Form("%sQCTPCMB%d",sgn.Data(),mb), ccName.Data(), fileName.Data(), thecuts.Data(), harm, exc_TPC, filterPOIQC[mb]);
    }
    if(bDoSPVZERO) {
      AddSPmethod( Form("%sSPVZEMB%d",sgn.Data(),mb), ccName.Data(), fileName.Data(), thecuts.Data(), etaVZERO1, etaVZERO2, etaVZERO3, etaVZERO4,
		   "Qa", harm, exc_VZE, 0, filterPOIQC[mb], NULL, bDoEPAlso );
      AddSPmethod( Form("%sSPVZEMB%d",sgn.Data(),mb), ccName.Data(), fileName.Data(), thecuts.Data(), etaVZERO1, etaVZERO2, etaVZERO3, etaVZERO4,
		   "Qb", harm, exc_VZE, 0, filterPOIQC[mb], NULL, bDoEPAlso );
    }
    if(bDoSPTPC) {
      for(int eg=0; eg!=2; ++eg) {
	AddSPmethod( Form("%sSPTPCMB%d",sgn.Data(),mb), ccName.Data(), fileName.Data(), thecuts.Data(), -0.8, -gapTPC*eg, +gapTPC*eg, +0.8,
		     "Qa", harm, exc_TPC, eg, filterPOISP[mb][1], NULL, bDoEPAlso );
	AddSPmethod( Form("%sSPTPCMB%d",sgn.Data(),mb), ccName.Data(), fileName.Data(), thecuts.Data(), -0.8, -gapTPC*eg, +gapTPC*eg, +0.8,
		     "Qb", harm, exc_TPC, eg, filterPOISP[mb][0], NULL, bDoEPAlso );
      }
    }
  }
}

void AddSPmethod(char *name, char *ccName, char *fileName, char *thecuts,
		 double minEtaA, double maxEtaA, double minEtaB, double maxEtaB,
		 char *Qvector, int harmonic, AliAnalysisDataContainer *flowEvent, int eg,
		 AliFlowTrackSimpleCuts *cutsPOI=NULL, AliFlowTrackSimpleCuts *cutsRFP=NULL,
		 bool bEP ) {
  TString myFolder = Form("%sv%d_%s",ccName,harmonic,thecuts);
  TString myNameSP = Form("%s%sSPv%d%sGAP%d",name,ccName,harmonic,Qvector,eg);
  TString myNameEP = Form("%s%sEPv%d%sGAP%d",name,ccName,harmonic,Qvector,eg);
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
  //SP
  AliAnalysisDataContainer *outSP = mgr->CreateContainer( myNameSP.Data(),
							  TList::Class(),AliAnalysisManager::kOutputContainer,
							  Form("%s.root:FlowD2H_SP_%s",fileName,myFolder.Data()) );
  AliAnalysisTaskScalarProduct *tskSP = new AliAnalysisTaskScalarProduct( Form("TaskScalarProduct_%s",
									       myNameSP.Data()),kFALSE);
  tskSP->SetApplyCorrectionForNUA(kTRUE);
  tskSP->SetHarmonic(harmonic);
  tskSP->SetTotalQvector(Qvector);
  mgr->AddTask(tskSP);
  mgr->ConnectInput( tskSP,0,flowEvent2);
  mgr->ConnectOutput(tskSP,1,outSP);

  //EP
  if(bEP) {
    AliAnalysisDataContainer *outEP = mgr->CreateContainer( myNameEP.Data(),
							    TList::Class(),AliAnalysisManager::kOutputContainer,
							    Form("%s.root:FlowD2H_EP_%s",fileName,myFolder.Data()) );
    AliAnalysisTaskScalarProduct *tskEP = new AliAnalysisTaskScalarProduct( Form("TaskEventsPlane_%s",
										 myNameEP.Data()),kFALSE);
    tskEP->SetApplyCorrectionForNUA(kTRUE);
    tskEP->SetHarmonic(harmonic);
    tskEP->SetTotalQvector(Qvector);
    tskEP->SetBehaveAsEP();
    mgr->AddTask(tskEP);
    mgr->ConnectInput( tskEP,0,flowEvent2);
    mgr->ConnectOutput(tskEP,1,outEP);
  }
}

void AddQCmethod(char *name, char *ccName, char *fileName, char *thecuts,
		 int harmonic, AliAnalysisDataContainer *flowEvent,
		 AliFlowTrackSimpleCuts *cutsPOI=NULL, AliFlowTrackSimpleCuts *cutsRFP=NULL) {
  TString myFolder = Form("%sv%d_%s",ccName,harmonic,thecuts);
  TString myName = Form("%s%sv%d",name,ccName,harmonic);
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
  mgr->AddTask(tskQC);
  mgr->ConnectInput( tskQC,0,flowEvent2);
  mgr->ConnectOutput(tskQC,1,outQC);
}

AliAnalysisDataContainer* MakeFlowEvent(char* RFPName, char* ccName, char *fileName, char *thecuts, Int_t ptbins,
					AliAnalysisDataContainer *input, 
					AliAnalysisDataContainer *array,
					AliFlowEventCuts *eventCuts,
					AliFlowTrackCuts *rfpCuts, AliFlowTrackCuts *poiCuts, 
					bool turnOnQA, int nDmeson, double etaMax=1.0) {
  TString meson = DMesonName(nDmeson);
  TString myName = Form("%sFlowEventQA%s%s",meson.Data(),ccName,RFPName);
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisDataContainer *exc = mgr->CreateContainer( Form("exchange_%s", myName.Data()),
							AliFlowEventSimple::Class(), 
							AliAnalysisManager::kExchangeContainer );
  AliAnalysisDataContainer *out = mgr->CreateContainer( myName.Data(),
							TList::Class(),AliAnalysisManager::kOutputContainer,
							Form("%s.root:FlowD2H_%s_%s",fileName,ccName,thecuts) );
  AliAnalysisTaskFlowEvent *fe = new AliAnalysisTaskFlowEvent( Form("Task_%s",myName.Data()),"", kFALSE, 666, kTRUE );
  fe->SetCutsEvent(eventCuts);
  fe->SetCutsRP(rfpCuts);
  fe->SetCutsPOI(poiCuts);
  fe->SetQAOn(turnOnQA);
  fe->SelectCollisionCandidates(AliVEvent::kCentral|AliVEvent::kSemiCentral|AliVEvent::kMB);

  fe->SetNbinsEta(20); fe->SetEtaMin(-etaMax); fe->SetEtaMax(+etaMax);

  int bins = NumberOfMassBins(nDmeson);
  fe->SetNbinsMass( bins );
  fe->SetMassMin( LowMassBin(nDmeson,0) );
  fe->SetMassMax( LowMassBin(nDmeson,bins) );

  fe->SetNbinsPt(24/ptbins); fe->SetPtMin(0); fe->SetPtMax(24);
  fe->SetNbinsMult(2); fe->SetMultMin(2); fe->SetMultMax(10);
  fe->SetNbinsPhi(20); fe->SetPhiMin(0); fe->SetPhiMax(TMath::TwoPi());
  fe->SetNbinsQ(10); fe->SetQMin( 0 ); fe->SetQMax( 3 );

  mgr->AddTask(fe);
  mgr->ConnectInput (fe,0,input);
  mgr->ConnectInput (fe,1,array);
  mgr->ConnectOutput(fe,1,exc);
  mgr->ConnectOutput(fe,2,out);
  return exc;
}

double LowMassBin( int nDmeson, int mb) {
  double min,max;
  int bins = NumberOfMassBins( nDmeson );
  if(mb>bins) return 0.0;
  switch (nDmeson) {
  case ( AliRDHFCuts::kD0toKpiCuts ):
    min = 1.700;
    max = 2.200;
    break;
  case ( AliRDHFCuts::kDstarCuts ):
    min = 0.138;
    max = 0.158;
    break;
  }
  return min+mb*(max-min)/bins;
}

int NumberOfMassBins( int nDmeson ) {
  switch (nDmeson) {
  case ( AliRDHFCuts::kD0toKpiCuts ):
    return 50;
  case ( AliRDHFCuts::kDstarCuts ):
    return 25;
  }
}

TString DMesonName( int nDmeson ) {
  TString toReturn;
  switch (nDmeson) {
  case ( AliRDHFCuts::kD0toKpiCuts ):
    toReturn = "D0"; break;
  case ( AliRDHFCuts::kDstarCuts ):
    toReturn = "DStar"; break;
  }
  return toReturn;
}
