void AddTaskFlowStrange(int trigger, int centrMin, int centrMax, int harmonic=2) {
  AddTaskFlowStrange(trigger,centrMin,centrMax,"K0",Form("K%d%dc0",centrMin,centrMax),   0,0,"V0M",harmonic);
  AddTaskFlowStrange(trigger,centrMin,centrMax,"K0",Form("K%d%dc1",centrMin,centrMax),   0,1,"V0M",harmonic);
  AddTaskFlowStrange(trigger,centrMin,centrMax,"K0",Form("K%d%dc2",centrMin,centrMax),   0,2,"V0M",harmonic);
  AddTaskFlowStrange(trigger,centrMin,centrMax,"K0",Form("K%d%dc3",centrMin,centrMax),   0,3,"V0M",harmonic);
  AddTaskFlowStrange(trigger,centrMin,centrMax,"K0",Form("K%d%dc2TRK",centrMin,centrMax),0,2,"TRK",harmonic);

  AddTaskFlowStrange(trigger,centrMin,centrMax,"L0",Form("L%d%dc0",centrMin,centrMax),   1,0,"V0M",harmonic);
  AddTaskFlowStrange(trigger,centrMin,centrMax,"L0",Form("L%d%dc1",centrMin,centrMax),   1,1,"V0M",harmonic);
  AddTaskFlowStrange(trigger,centrMin,centrMax,"L0",Form("L%d%dc2",centrMin,centrMax),   1,2,"V0M",harmonic);
  AddTaskFlowStrange(trigger,centrMin,centrMax,"L0",Form("L%d%dc3",centrMin,centrMax),   1,3,"V0M",harmonic);
  AddTaskFlowStrange(trigger,centrMin,centrMax,"L0",Form("L%d%dc2TRK",centrMin,centrMax),1,2,"TRK",harmonic);
}
void AddTaskFlowStrange(int trigger, float centrMin, float centrMax, TString folderName="myFolder", TString suffixName="mySuffix", 
			int specie=0, int cuts=1, char* MULT="V0M", int harmonic=2) {
  TString fileName = AliAnalysisManager::GetCommonFileName();
  fileName.ReplaceAll(".root","");

  //-E-V-E-N-T- -c-u-t-s--------------------------------------------------------
  AliFlowEventCuts* cutsEvent = new AliFlowEventCuts(Form("event_cuts_%s",suffixName.Data()));
  cutsEvent->SetCentralityPercentileRange(centrMin,centrMax);
  if(MULT=="V0M")
    cutsEvent->SetCentralityPercentileMethod(AliFlowEventCuts::kV0);
  else
    cutsEvent->SetCentralityPercentileMethod(AliFlowEventCuts::kTPConly);
  cutsEvent->SetRefMultMethod(AliFlowEventCuts::kSPDtracklets);
  cutsEvent->SetNContributorsRange(2);
  cutsEvent->SetPrimaryVertexZrange(-9.5,+9.5);

  //-R-P---c-u-t-s--------------------------------------------------------------
  AliFlowTrackCuts *cutsRPTPC = AliFlowTrackCuts::GetStandardTPCStandaloneTrackCuts();
  cutsRPTPC->SetParamType( AliFlowTrackCuts::kGlobal );
  cutsRPTPC->SetAODfilterBit(1); // for AOD compatibility
  AliFlowTrackCuts *cutsRPVZE = AliFlowTrackCuts::GetStandardVZEROOnlyTrackCuts();

  //-D-A-U-G-H-T-E-R-S---c-u-t-s------------------------------------------------
  AliESDtrackCuts* cutsDaughter = new AliESDtrackCuts(Form("daughter_cuts_%s",suffixName.Data()) );
  cutsDaughter->SetPtRange(0.2,10.0);
  cutsDaughter->SetEtaRange(-0.8, 0.8 );
  cutsDaughter->SetMinNClustersTPC(80);
  cutsDaughter->SetMaxChi2PerClusterTPC(4.0);
  cutsDaughter->SetRequireTPCRefit(kTRUE);
  cutsDaughter->SetAcceptKinkDaughters(kFALSE);

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();

  //-----------------STRANGE TASK----------------------------
  AliAnalysisTaskFlowStrange *taskSel = new AliAnalysisTaskFlowStrange(Form("Strange_%s",suffixName.Data()),
								       cutsEvent, cutsRPTPC, cutsRPVZE,
								       cutsDaughter );
  taskSel->SelectCollisionCandidates(trigger);
  //taskSel->SetDebug();
  taskSel->SetCuts2010(cuts);
  taskSel->SetK0L0(specie);
  //printf( "CMM %d %f %f\n", MassBins(specie), MinMass(specie), MaxMass(specie) );
  taskSel->SetCommonConstants( MassBins(specie), MinMass(specie), MaxMass(specie) );
  AliAnalysisDataContainer *cOutHist = mgr->CreateContainer(Form("OutHistos_%s",suffixName.Data()),
							    TList::Class(),
							    AliAnalysisManager::kOutputContainer,
							    Form("%s.root:Selector_%s",fileName.Data(),
								 folderName.Data()));
  AliAnalysisDataContainer *exc_TPC = mgr->CreateContainer( Form("TPCEventWithCandidates_%s",suffixName.Data()),
                                                            AliFlowEventSimple::Class(),
                                                            AliAnalysisManager::kExchangeContainer );
  AliAnalysisDataContainer *exc_VZE = mgr->CreateContainer( Form("VZEEventWithCandidates_%s",suffixName.Data()),
                                                            AliFlowEventSimple::Class(),
                                                            AliAnalysisManager::kExchangeContainer );
  mgr->AddTask(taskSel);
  mgr->ConnectInput (taskSel,0,cinput1);
  mgr->ConnectOutput(taskSel,1,exc_TPC);
  mgr->ConnectOutput(taskSel,2,exc_VZE);
  mgr->ConnectOutput(taskSel,3,cOutHist);

  //-------------------FLOW TASKS----------------------------
  AliFlowTrackSimpleCuts *filter[15], *filterhf[15][2]; // MASS BANDS
  for(int mb=0; mb!=MassBands(0); ++mb) {
    filter[mb] = new AliFlowTrackSimpleCuts( Form("Filter_MB%d",mb) );
    filter[mb]->SetEtaMin( -0.8 ); filter[mb]->SetEtaMax( +0.8 );
    filter[mb]->SetMassMin( MassBandLowEdge(specie,mb) ); filter[mb]->SetMassMax( MassBandLowEdge(specie,mb+1) );

    filterhf[mb][0] = new AliFlowTrackSimpleCuts( Form("Filterhf0_MB%d",mb) );
    filterhf[mb][0]->SetEtaMin( 0.0 ); filterhf[mb][0]->SetEtaMax( +0.8 );
    filterhf[mb][0]->SetMassMin( MassBandLowEdge(specie,mb) ); filterhf[mb][0]->SetMassMax( MassBandLowEdge(specie,mb+1) );

    filterhf[mb][1] = new AliFlowTrackSimpleCuts( Form("Filterhf1_MB%d",mb) );
    filterhf[mb][1]->SetEtaMin( -0.8 ); filterhf[mb][1]->SetEtaMax( 0.0 );
    filterhf[mb][1]->SetMassMin( MassBandLowEdge(specie,mb) ); filterhf[mb][1]->SetMassMax( MassBandLowEdge(specie,mb+1) );

    AddQCmethod( Form("QCTPCMB%d",mb), folderName.Data(), suffixName.Data(), harmonic, exc_TPC, filter[mb]); // QC TPC
    AddSPmethod( Form("SPTPCMB%d",mb), folderName.Data(), suffixName.Data(), harmonic, exc_TPC, filterhf[mb][0], "Qa" ); // SP TPC Qa
    AddSPmethod( Form("SPTPCMB%d",mb), folderName.Data(), suffixName.Data(), harmonic, exc_TPC, filterhf[mb][1], "Qb" ); // SP TPC Qb
    AddSPmethod( Form("SPVZEMB%d",mb), folderName.Data(), suffixName.Data(), harmonic, exc_VZE, filter[mb], "Qa" ); // SP VZE Qa
    AddSPmethod( Form("SPVZEMB%d",mb), folderName.Data(), suffixName.Data(), harmonic, exc_VZE, filter[mb], "Qb" ); // SP VZE Qa
  }
}

void AddQCmethod(char *name, TString myFolder, char *thecuts, int harmonic, 
		 AliAnalysisDataContainer *flowEvent, AliFlowTrackSimpleCuts *cutsPOI=NULL) {
  TString fileName = AliAnalysisManager::GetCommonFileName();
  myFolder.Append( Form("v%d",harmonic) );
  TString myName = Form("%sv%d_%s",name,harmonic,thecuts);
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisDataContainer *flowEvent2 = mgr->CreateContainer( Form("Filter_%s", myName.Data()),
                                                               AliFlowEventSimple::Class(),
                                                               AliAnalysisManager::kExchangeContainer );
  AliAnalysisTaskFilterFE *tskFilter = new AliAnalysisTaskFilterFE( Form("TaskFilter_%s",myName.Data()),
                                                                    NULL, cutsPOI);
  mgr->AddTask(tskFilter);
  mgr->ConnectInput( tskFilter,0,flowEvent);
  mgr->ConnectOutput(tskFilter,1,flowEvent2);
  AliAnalysisDataContainer *outQC = mgr->CreateContainer( myName.Data(),TList::Class(),AliAnalysisManager::kOutputContainer,
                                                          Form("%s:FlowStrange_QC_%s",fileName.Data(),myFolder.Data()) );
  AliAnalysisTaskQCumulants *tskQC = new AliAnalysisTaskQCumulants( Form("TaskQCumulants_%s",myName.Data()),kFALSE );
  tskQC->SetApplyCorrectionForNUA(kTRUE);
  tskQC->SetHarmonic(harmonic);
  tskQC->SetBookOnlyBasicCCH(kTRUE);
  mgr->AddTask(tskQC);
  mgr->ConnectInput( tskQC,0,flowEvent2);
  mgr->ConnectOutput(tskQC,1,outQC);
}


void AddSPmethod(char *name, TString myFolder, char *thecuts, int harmonic,
		 AliAnalysisDataContainer *flowEvent, AliFlowTrackSimpleCuts *cutsPOI=NULL,
                 char *Qvector) {
  TString fileName = AliAnalysisManager::GetCommonFileName();
  myFolder.Append( Form("v%d",harmonic) );
  TString myNameSP = Form("%sv%d%s_%s",name,harmonic,Qvector,thecuts);
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisDataContainer *flowEvent2 = mgr->CreateContainer( Form("Filter_%s", myNameSP.Data()),
                                                               AliFlowEventSimple::Class(),
                                                               AliAnalysisManager::kExchangeContainer );
  AliAnalysisTaskFilterFE *tskFilter = new AliAnalysisTaskFilterFE( Form("TaskFilter_%s",myNameSP.Data()),
                                                                    NULL, cutsPOI);
  tskFilter->SetSubeventEtaRange( -5.0, 0.0, 0.0, +5.0 );
  mgr->AddTask(tskFilter);
  mgr->ConnectInput( tskFilter,0,flowEvent);
  mgr->ConnectOutput(tskFilter,1,flowEvent2);
  AliAnalysisDataContainer *outSP = mgr->CreateContainer( myNameSP.Data(),TList::Class(),AliAnalysisManager::kOutputContainer,
                                                          Form("%s:FlowStrange_SP_%s",fileName.Data(),myFolder.Data()) );
  AliAnalysisTaskScalarProduct *tskSP = new AliAnalysisTaskScalarProduct( Form("TaskScalarProduct_%s",myNameSP.Data()),kFALSE);
  tskSP->SetApplyCorrectionForNUA(kTRUE);
  tskSP->SetHarmonic(harmonic);
  tskSP->SetTotalQvector(Qvector);
  tskSP->SetBookOnlyBasicCCH(kTRUE);
  mgr->AddTask(tskSP);
  mgr->ConnectInput( tskSP,0,flowEvent2);
  mgr->ConnectOutput(tskSP,1,outSP);
}

double MassBandLowEdge( int nv0, int mb ) {
  if(nv0==0) {
    double lowEdge[15+1]={ 0.452, 0.462, 0.472, 0.482, 0.489, 0.492, 0.495, 0.498,
			   0.501, 0.504, 0.507, 0.514, 0.524, 0.534, 0.544, 0.554 };
  } if(nv0==1) {
    double lowEdge[15+1]={ 1.094, 1.099, 1.104, 1.109, 1.112, 1.114, 1.115, 1.116, 
			   1.117, 1.118, 1.120, 1.123, 1.128, 1.133, 1.138, 1.143 };
  }
  return lowEdge[mb];
}

int MassBands( int nv0 ) {
  if(nv0==0) {
    return 15;
  } else if(nv0==1) {
    return 15;
  }
}

int MassBins( int nv0 ) {
  if(nv0==0) {
    return 102;
  } else if(nv0==1) {
    return 49;
  }
}

double MinMass( int nv0 ) {
  return MassBandLowEdge( nv0, 0 );
}

double MaxMass( int nv0 ) {
  return MassBandLowEdge( nv0, MassBands(nv0) );
}
