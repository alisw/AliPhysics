// ONE PRECONF WAGON
void AddTaskFlowStrange(int trigger, float centrMin, float centrMax, int set, TString folderName="myFolder", TString suffixName="mySuffix", 
			int specie=0, char* MULT="V0M", int harmonic=2, int matchMC=-1, bool doQC=true, bool doSPTPC=true, bool doSPVZE=true,
			bool doQA=false, bool useFlowEventCuts=true) {
  Double_t c[11];
  switch(set) {
  case(0): // Filtering cuts //
    c[0]=0; c[1]=+1; c[2]=+0.500; c[3]=0;    c[4]=1e+6; c[5]=0;    c[6]=-1;   c[7]=+1;   c[8]=1e+6; c[9]=1e+6; c[10]=0;
    break;
  case(1): // Topological cuts // 
    c[0]=0; c[1]=+1; c[2]=+0.998; c[3]=+0.1; c[4]=0;    c[5]=+0.2; c[6]=-0.8; c[7]=+0.8; c[8]=+3.0; c[9]=+3.0; c[10]=5;
    break;
  default:
    return;
  }
  AddTaskFlowStrange(trigger, centrMin, centrMax, folderName, suffixName, specie,
		     c[0],c[1],c[2],c[3],c[4],c[5],c[6],c[7],c[8],c[9],c[10],
		     MULT, harmonic, matchMC, doQC, doSPTPC, doSPVZE, doQA, useFlowEventCuts);
}
// ALL CENTRALITIES
void AddTaskFlowStrange(int trigger, TString folderName="myFolder", TString suffix="mySuffix", int specie=0,
                        double cut0, double cut1, double cut2, double cut3, double cut4, double cut5, double cut6, double cut7, double cut8,
                        double cut9, double cut10, char* MULT="V0M", int harmonic=2, int matchMC=-1,
                        bool doQC=true, bool doSPTPC=true, bool doSPVZE=true, bool doQA=false, bool useFlowEventCuts=true ) {
  int centMin[8] = {00,05,10,20,30,40,50,60};
  int centMax[8] = {05,10,20,30,40,50,60,80};
  TString particle="none";
  if(specie==0) particle="kze";
  if(specie==1) particle="lda";
  if(specie==90) particle="kch";
  TString name;
  for(int i=0; i!=8; ++i) {
    name = Form("%s%02d%02d%s",particle.Data(),centMin[i],centMax[i],suffix.Data());
    AddTaskFlowStrange(trigger, centMin[i], centMax[i], folderName, name, specie,
		       cut0, cut1, cut2, cut3, cut4, cut5, cut6, cut7, cut8, cut9, cut10, MULT, harmonic, matchMC,
		       doQC, doSPTPC, doSPVZE, doQA, useFlowEventCuts);
  }
}
// ONE WAGON
void AddTaskFlowStrange(int trigger, float centrMin, float centrMax, TString folderName="myFolder", TString suffixName="mySuffix", int specie=0, 
			double cut0, double cut1, double cut2, double cut3, double cut4, double cut5, double cut6, double cut7, double cut8,
			double cut9, double cut10, char* MULT="V0M", int harmonic=2, int matchMC=-1, 
			bool doQC=true, bool doSPTPC=true, bool doSPVZE=true, bool doQA=false, bool useFlowEventCuts=true) {
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
  cutsEvent->SetPrimaryVertexZrange(-10.0,+10.0);

  //-R-P---c-u-t-s--------------------------------------------------------------
  AliFlowTrackCuts *cutsRPTPC = AliFlowTrackCuts::GetStandardTPCStandaloneTrackCuts();
  cutsRPTPC->SetParamType( AliFlowTrackCuts::kGlobal );
  cutsRPTPC->SetAODfilterBit(1); // for AOD compatibility
  AliFlowTrackCuts *cutsRPVZE = AliFlowTrackCuts::GetStandardVZEROOnlyTrackCuts();

  //-D-A-U-G-H-T-E-R-S---c-u-t-s------------------------------------------------
  AliESDtrackCuts* cutsDaughter = new AliESDtrackCuts(Form("daughter_cuts_%s",suffixName.Data()) );
  cutsDaughter->SetPtRange(0.2,12.0);
  cutsDaughter->SetEtaRange(-0.8, 0.8 );
  cutsDaughter->SetMinNClustersTPC(70);
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
  Double_t cuts[11];
  cuts[0]=cut0; cuts[1]=cut1; cuts[2]=cut2; cuts[3]=cut3; cuts[4]=cut4; cuts[5]=cut5; cuts[6]=cut6; cuts[7]=cut7; cuts[8]=cut8; cuts[9]=cut9; cuts[10]=cut10;
  taskSel->SetCuts(cuts);
  taskSel->SetQA(doQA);
  taskSel->SetK0L0(specie);
  taskSel->SetMCmatch(matchMC);
  taskSel->SetUseEventSelection(useFlowEventCuts);
  taskSel->SetCommonConstants( SFT_MassBins(specie), SFT_MinMass(specie), SFT_MaxMass(specie) );
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

  if( (!doQC) && (!doSPVZE) && (!doSPVZE) ) return;

  //-------------------FLOW TASKS----------------------------
  AliFlowTrackSimpleCuts *filter[12], *filterhf[12][2]; // MASS BANDS
  for(int mb=0; mb!=12; ++mb) {
    filter[mb] = new AliFlowTrackSimpleCuts( Form("Filter_MB%d",mb) );
    filter[mb]->SetEtaMin( -0.8 ); filter[mb]->SetEtaMax( +0.8 );
    filter[mb]->SetMassMin( SFT_MassBandLowEdge(specie,mb) ); filter[mb]->SetMassMax( SFT_MassBandLowEdge(specie,mb+1) );

    filterhf[mb][0] = new AliFlowTrackSimpleCuts( Form("Filterhf0_MB%d",mb) );
    filterhf[mb][0]->SetEtaMin( 0.0 ); filterhf[mb][0]->SetEtaMax( +0.8 );
    filterhf[mb][0]->SetMassMin( SFT_MassBandLowEdge(specie,mb) ); filterhf[mb][0]->SetMassMax( SFT_MassBandLowEdge(specie,mb+1) );

    filterhf[mb][1] = new AliFlowTrackSimpleCuts( Form("Filterhf1_MB%d",mb) );
    filterhf[mb][1]->SetEtaMin( -0.8 ); filterhf[mb][1]->SetEtaMax( 0.0 );
    filterhf[mb][1]->SetMassMin( SFT_MassBandLowEdge(specie,mb) ); filterhf[mb][1]->SetMassMax( SFT_MassBandLowEdge(specie,mb+1) );

    if(doQC) {
      SFT_AddQCmethod( Form("QCTPCMB%d",mb), folderName.Data(), suffixName.Data(), harmonic, exc_TPC, filter[mb]); // QC TPC
    }
    if(doSPTPC) {
      SFT_AddSPmethod( Form("SPTPCMB%d",mb), folderName.Data(), suffixName.Data(), harmonic, exc_TPC, filterhf[mb][0], "Qa" ); // SP TPC Qa
      SFT_AddSPmethod( Form("SPTPCMB%d",mb), folderName.Data(), suffixName.Data(), harmonic, exc_TPC, filterhf[mb][1], "Qb" ); // SP TPC Qb
    }
    if(doSPVZE) {
      SFT_AddSPmethod( Form("SPVZEMB%d",mb), folderName.Data(), suffixName.Data(), harmonic, exc_VZE, filter[mb], "Qa" ); // SP VZE Qa
      SFT_AddSPmethod( Form("SPVZEMB%d",mb), folderName.Data(), suffixName.Data(), harmonic, exc_VZE, filter[mb], "Qb" ); // SP VZE Qa
    }
  }
}
// ADDING QC METHOD
void SFT_AddQCmethod(char *name, TString myFolder, char *thecuts, int harmonic, 
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
// ADDING SP METHOD
void SFT_AddSPmethod(char *name, TString myFolder, char *thecuts, int harmonic,
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
// MASSBAND HELPER
double SFT_MassBandLowEdge( int nv0, int mb ) {
  switch(nv0) {
  case(0):
    double lowEdge[13]={ 0.412, 0.440, 0.468, 0.484, 0.492, 0.496, 0.500, 0.504, 0.508, 0.516, 0.532, 0.560, 0.588 };
    break;
  case(1):
    double lowEdge[13]={ 1.075, 1.090, 1.100, 1.108, 1.112, 1.114, 1.116, 1.118, 1.120, 1.124, 1.132, 1.142, 1.167 };
    break;
  default:
    double lowEdge[13]={ 0.000, 0.100, 0.200, 0.300, 0.400, 0.500, 0.600, 0.700, 0.800, 0.840, 0.860, 0.900, 1.000 };
  }
  return lowEdge[mb];
}
// MASSBAND HELPER
int SFT_MassBands( int nv0 ) {
  return 12;
}
// MASSBAND HELPER
int SFT_MassBins( int nv0 ) {
  int bins=88;
  switch(nv0) {
  case(0):
    bins=88;
    break;
  case(1):
    bins=92;
    break;
  default:
    bins=13;
    break;
  }
  return bins;
}
// MASSBAND HELPER
double SFT_MinMass( int nv0 ) {
  return SFT_MassBandLowEdge( nv0, 0 );
}
// MASSBAND HELPER
double SFT_MaxMass( int nv0 ) {
  return SFT_MassBandLowEdge( nv0, SFT_MassBands(nv0) );
}
