void AddTaskFlowCascade(int trigger, int centrMin, int centrMax, int cut=1, int harmonic=2) {
  AddTaskFlowCascade(trigger, centrMin, centrMax, "Xi", 
		     Form("Xi%d%dcut%d", centrMin, centrMax, cut), 0, cut, "V0M", "VZESP", harmonic);

  AddTaskFlowCascade(trigger, centrMin, centrMax, "Omega", 
		     Form("Omega%d%dcut%d", centrMin,centrMax, cut), 1, cut, "V0M", "VZESP", harmonic);
}

void AddTaskFlowCascade(int trigger, float centrMin, float centrMax, 
			TString folderName="myFolder", TString suffixName="mySuffix", 
			int specie=0, int cuts=1, char* MULT="V0M", 
			TString method="VZESP QC TPCSP", int harmonic=2) {
  TString fileName = AliAnalysisManager::GetCommonFileName();
  //  fileName.ReplaceAll(".root","");

  //-E-V-E-N-T- -c-u-t-s--------------------------------------------------------
  AliFlowEventCuts* cutsEvent = new AliFlowEventCuts(Form("event_cuts_%s",suffixName.Data()));
  cutsEvent->SetCentralityPercentileRange(centrMin,centrMax);
  if(MULT=="V0M")
    cutsEvent->SetCentralityPercentileMethod(AliFlowEventCuts::kV0);
  else
    cutsEvent->SetCentralityPercentileMethod(AliFlowEventCuts::kTPConly);
  cutsEvent->SetRefMultMethod(AliFlowEventCuts::kSPDtracklets);
  cutsEvent->SetNContributorsRange(2);
  cutsEvent->SetPrimaryVertexZrange(-10.,+10.);

  //-R-P---c-u-t-s--------------------------------------------------------------
  AliFlowTrackCuts *cutsRPTPC 
    = AliFlowTrackCuts::GetStandardTPCStandaloneTrackCuts();
  cutsRPTPC->SetParamType( AliFlowTrackCuts::kGlobal );
  cutsRPTPC->SetAODfilterBit(1); // for AOD compatibility
  
  AliFlowTrackCuts *cutsRPVZE = AliFlowTrackCuts::GetStandardVZEROOnlyTrackCuts();
  //  AliFlowTrackCuts* trkCuts = new AliFlowTrackCuts();
  //AliFlowTrackCuts *cutsRPVZE = trkCuts->GetStandardVZEROOnlyTrackCuts();
  cutsRPVZE->SetApplyRecentering(kTRUE);
  
  //-D-A-U-G-H-T-E-R-S---c-u-t-s------------------------------------------------
  //  AliESDtrackCuts* cutsDaughter 
  //  = new AliESDtrackCuts(Form("daughter_cuts_%s",suffixName.Data()) );
  AliFlowTrackCuts * cutsDaughter 
    = new AliFlowTrackCuts(Form("daughter_cuts_%s",suffixName.Data()));
  cutsDaughter->SetPtRange(0.15,10.0);
  cutsDaughter->SetEtaRange(-0.8, 0.8 );
  cutsDaughter->SetMinNClustersTPC(70);
  cutsDaughter->SetMaxChi2PerClusterTPC(4.0);
  cutsDaughter->SetRequireTPCRefit(kTRUE);
  cutsDaughter->SetAcceptKinkDaughters(kFALSE);

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();

  //-----------------Cascade TASK----------------------------
  AliAnalysisTaskFlowCascade *taskSel 
    = new AliAnalysisTaskFlowCascade(Form("Cascade_%s",suffixName.Data()),
				     cutsEvent, cutsRPTPC, cutsRPVZE,
				     cutsDaughter );
  taskSel->SelectCollisionCandidates(trigger);
  //taskSel->SetDebug();
  taskSel->SetCuts2010(cuts);
  taskSel->SetSpecie(specie);
  //printf( "CMM %d %f %f\n", MassBins(specie), MinMass(specie), MSFT_MaxMass(specie) );
  taskSel->SetCommonConstants( MSFT_MassBins(specie), MSFT_MinMass(specie), MSFT_MaxMass(specie) );
  AliAnalysisDataContainer *cOutHist 
    = mgr->CreateContainer(Form("OutHistos_%s",suffixName.Data()),
			   TList::Class(),
			   AliAnalysisManager::kOutputContainer,
			   Form("%s:Selector_%s",fileName.Data(),
				folderName.Data()));
  AliAnalysisDataContainer *exc_TPC 
    = mgr->CreateContainer( Form("TPCEventWithCandidates_%s",suffixName.Data()),
			    AliFlowEventSimple::Class(),
			    AliAnalysisManager::kExchangeContainer );
  AliAnalysisDataContainer *exc_VZE 
    = mgr->CreateContainer( Form("VZEEventWithCandidates_%s",suffixName.Data()),
			    AliFlowEventSimple::Class(),
			    AliAnalysisManager::kExchangeContainer );
  mgr->AddTask(taskSel);
  mgr->ConnectInput (taskSel, 0, cinput1);
  mgr->ConnectOutput(taskSel, 1, exc_TPC);
  mgr->ConnectOutput(taskSel, 2, exc_VZE);
  mgr->ConnectOutput(taskSel, 3, cOutHist);

  //-------------------FLOW TASKS----------------------------
  AliFlowTrackSimpleCuts *filter[15], *filterhf[15][2]; // MASS BANDS
  for(int mb=0; mb!=MSFT_MassBands(0); ++mb) {
    filter[mb] = new AliFlowTrackSimpleCuts( Form("Filter_MB%d",mb) );
    filter[mb]->SetEtaMin( -0.8 ); 
    filter[mb]->SetEtaMax( +0.8 );
    filter[mb]->SetMassMin( MSFT_MassBandLowEdge(specie, mb) ); 
    filter[mb]->SetMassMax( MSFT_MassBandLowEdge(specie, mb+1) );

    filterhf[mb][0] = new AliFlowTrackSimpleCuts( Form("Filterhf0_MB%d",mb) );
    filterhf[mb][0]->SetEtaMin( 0.0 ); 
    filterhf[mb][0]->SetEtaMax( +0.8 );
    filterhf[mb][0]->SetMassMin( MSFT_MassBandLowEdge(specie, mb) ); 
    filterhf[mb][0]->SetMassMax( MSFT_MassBandLowEdge(specie, mb+1) );

    filterhf[mb][1] = new AliFlowTrackSimpleCuts( Form("Filterhf1_MB%d",mb) );
    filterhf[mb][1]->SetEtaMin( -0.8 ); 
    filterhf[mb][1]->SetEtaMax( 0.0 );
    filterhf[mb][1]->SetMassMin( MSFT_MassBandLowEdge(specie, mb) ); 
    filterhf[mb][1]->SetMassMax( MSFT_MassBandLowEdge(specie, mb+1) );

    if(method.Contains("QC"))
      MSFT_AddQCmethod( Form("QCTPCMB%d", mb), folderName.Data(), suffixName.Data(), harmonic, 
		   exc_TPC, filter[mb]); // QC TPC
    if(method.Contains("TPCSP")){
      MSFT_AddSPmethod( Form("SPTPCMB%d", mb), folderName.Data(), suffixName.Data(), harmonic, 
		   exc_TPC, filterhf[mb][0], "Qa" ); // SP TPC Qa
      MSFT_AddSPmethod( Form("SPTPCMB%d", mb), folderName.Data(), suffixName.Data(), harmonic, 
		 exc_TPC, filterhf[mb][1], "Qb" ); // SP TPC Qb
    }
    if(method.Contains("VZESP")){
      MSFT_AddSPmethod( Form("SPVZEMB%d", mb), folderName.Data(), suffixName.Data(), harmonic, 
		   exc_VZE, filter[mb], "Qa" ); // SP VZE Qa
      MSFT_AddSPmethod( Form("SPVZEMB%d",mb), folderName.Data(), suffixName.Data(), harmonic, 
		   exc_VZE, filter[mb], "Qb" ); // SP VZE Qa
      MSFT_AddSPmethod( Form("SPVZEMB%d",mb), folderName.Data(), suffixName.Data(), harmonic,
		   exc_VZE, filter[mb], "QaQb" ); // SP VZE QaQb
    }
  }
}

void MSFT_AddQCmethod(char *name, TString myFolder, char *thecuts, int harmonic, 
		      AliAnalysisDataContainer *flowEvent, 
		      AliFlowTrackSimpleCuts *cutsPOI=NULL) {
  TString fileName = AliAnalysisManager::GetCommonFileName();
  myFolder.Append( Form("v%d",harmonic) );
  TString myName = Form("%sv%d_%s", name, harmonic, thecuts);
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisDataContainer *flowEvent2 
    = mgr->CreateContainer( Form("Filter_%s", myName.Data()),
			    AliFlowEventSimple::Class(),
			    AliAnalysisManager::kExchangeContainer );
  AliAnalysisTaskFilterFE *tskFilter 
    = new AliAnalysisTaskFilterFE( Form("TaskFilter_%s",myName.Data()),
				   NULL, cutsPOI);
  mgr->AddTask(tskFilter);
  mgr->ConnectInput( tskFilter, 0, flowEvent);
  mgr->ConnectOutput(tskFilter, 1, flowEvent2);

  AliAnalysisDataContainer *outQC 
    = mgr->CreateContainer( myName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
			    Form("%s:FlowCascade_QC_%s", fileName.Data(), myFolder.Data()) );
  AliAnalysisTaskQCumulants *tskQC 
    = new AliAnalysisTaskQCumulants( Form("TaskQCumulants_%s",myName.Data()),kFALSE );
  tskQC->SetApplyCorrectionForNUA(kTRUE);
  tskQC->SetHarmonic(harmonic);
  tskQC->SetBookOnlyBasicCCH(kTRUE);
  mgr->AddTask(tskQC);
  mgr->ConnectInput( tskQC, 0, flowEvent2);
  mgr->ConnectOutput(tskQC, 1, outQC);
}


void MSFT_AddSPmethod(char *name, TString myFolder, char *thecuts, int harmonic,
		      AliAnalysisDataContainer *flowEvent, AliFlowTrackSimpleCuts *cutsPOI=NULL,
		      char *Qvector) {
  TString fileName = AliAnalysisManager::GetCommonFileName();
  myFolder.Append( Form("v%d",harmonic) );
  TString myNameSP = Form("%sv%d%s_%s", name, harmonic, Qvector, thecuts);
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisDataContainer *flowEvent2 
    = mgr->CreateContainer( Form("Filter_%s", myNameSP.Data()),
			    AliFlowEventSimple::Class(),
			    AliAnalysisManager::kExchangeContainer );
  AliAnalysisTaskFilterFE *tskFilter 
    = new AliAnalysisTaskFilterFE( Form("TaskFilter_%s",myNameSP.Data()),
				   NULL, cutsPOI);
  tskFilter->SetSubeventEtaRange( -5.0, -1.0, 1.0, +5.0 );
  mgr->AddTask(tskFilter);
  mgr->ConnectInput( tskFilter, 0, flowEvent);
  mgr->ConnectOutput(tskFilter, 1, flowEvent2);
  AliAnalysisDataContainer *outSP 
    = mgr->CreateContainer( myNameSP.Data(),TList::Class(),AliAnalysisManager::kOutputContainer,
			    Form("%s:FlowCascade_SP_%s",fileName.Data(),myFolder.Data()) );
  AliAnalysisTaskScalarProduct *tskSP 
    = new AliAnalysisTaskScalarProduct( Form("TaskScalarProduct_%s",myNameSP.Data()),kFALSE);
  tskSP->SetApplyCorrectionForNUA(kTRUE);
  tskSP->SetHarmonic(harmonic);
  tskSP->SetTotalQvector(Qvector);
  tskSP->SetBookOnlyBasicCCH(kTRUE);
  mgr->AddTask(tskSP);
  mgr->ConnectInput( tskSP,0,flowEvent2);
  mgr->ConnectOutput(tskSP,1,outSP);
}

double MSFT_MassBandLowEdge( int nv0, int mb ) {
  if(nv0==0) {
    double lowEdge[15+1]={ 1.282, 1.292, 1.302, 1.312, 1.317, 1.319, 1.321, 1.322, 
			   1.323, 1.325, 1.327, 1.332, 1.342, 1.352, 1.362, 1.372};
  } if(nv0==1) {
    double lowEdge[15+1]={ 1.632, 1.642, 1.652, 1.662, 1.667, 1.669, 1.671, 1.672, 
			   1.673, 1.675, 1.677, 1.682, 1.692, 1.702, 1.712, 1.722 };
  }
  return lowEdge[mb];
}

int MSFT_MassBands( int nv0 ) {
  if(nv0==0) {
    return 15;
  } else if(nv0==1) {
    return 15;
  }
}

int MSFT_MassBins( int nv0 ) {
  if(nv0==0) {
    return 90;
  } else if(nv0==1) {
    return 90;
  }
}

double MSFT_MinMass( int nv0 ) {
  return MSFT_MassBandLowEdge( nv0, 0 );
}

double MSFT_MaxMass( int nv0 ) {
  return MSFT_MassBandLowEdge( nv0, MSFT_MassBands(nv0) );
}
