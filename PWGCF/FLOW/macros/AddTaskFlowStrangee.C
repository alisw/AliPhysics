UInt_t SFT_gbTrigger;
Bool_t SFT_gbReadESD;
Bool_t SFT_gbReadMC;
Int_t SFT_gbMatchMC;
Bool_t SFT_gbAvoidExec;
Bool_t SFT_gbSkipCentrality;
Bool_t SFT_gbExtraEventCut;
TString SFT_gbCentMethod;
Int_t SFT_gbCentPerMin,SFT_gbCentPerMax;
Bool_t SFT_gbRunPP;
Bool_t SFT_gbRunPA;
Int_t SFT_gbSpecie;
Bool_t SFT_gbHomemade;
Bool_t SFT_gbOnline;
Int_t SFT_gbMinNClsTPC;
Int_t SFT_gbMinNClsITS;
Int_t SFT_gbMinXRows;
Double_t SFT_gbMaxChi2PerNClsTPC;
Double_t SFT_gbMinXRowsOverNClsFTPC;
Double_t SFT_gbMinEta;
Double_t SFT_gbMaxEta;
Double_t SFT_gbMinPt;
Double_t SFT_gbMinImpactParameterXY;
Double_t SFT_gbMaxNSigmaPID;
Double_t SFT_gbMaxRapidity;
Double_t SFT_gbMaxDCAdaughters;
Double_t SFT_gbMinCosinePointingAngleXY;
Double_t SFT_gbPIDPt;
Double_t SFT_gbMinQt;
Bool_t   SFT_gbQtPie;
Double_t SFT_gbMinRadXY;
Double_t SFT_gbMaxDecayLength;
Double_t SFT_gbMaxProductIPXY;
Int_t SFT_gbDebug;
Int_t SFT_gbQA;
TString SFT_gbFolder;
TString SFT_gbSuffix;
TString SFT_gbStamp;
Int_t SFT_gbRFPFilterBit;
Double_t SFT_gbRFPminPt;
Double_t SFT_gbRFPmaxPt;
Double_t SFT_gbRFPAminEta;
Double_t SFT_gbRFPAmaxEta;
Double_t SFT_gbRFPCminEta;
Double_t SFT_gbRFPCmaxEta;
Double_t SFT_gbRFPTPCsignal;
Double_t SFT_gbRFPmaxIPxy;
Double_t SFT_gbRFPmaxIPz;
Int_t SFT_gbRFPTPCncls;
Bool_t SFT_gbAddPitoMCRP;
Bool_t SFT_gbAllCC;
Bool_t SFT_gbSkipSelection;
Bool_t SFT_gbSkipVn;
Int_t SFT_gbWhichPsi;
Bool_t SFT_gbFlowPackage;
Bool_t SFT_gbShrinkFP;
Bool_t SFT_gbSPVZE;
Bool_t SFT_gbSPTPC;
Bool_t SFT_gbSPVZEhalf;
Bool_t SFT_gbQCTPC;
Bool_t SFT_gbMCEP;
Int_t SFT_gbHarmonic;
TString SFT_gbVZEload;
Bool_t SFT_gbVZEsave;
Bool_t SFT_gbVZEmb;
Bool_t SFT_gbVZEpdisk;
Int_t SFT_gbV0CRingMin;
Int_t SFT_gbV0CRingMax;
Int_t SFT_gbV0ARingMin;
Int_t SFT_gbV0ARingMax;
Int_t SFT_gbDauITS0On;
Int_t SFT_gbDauITS1On;
Int_t SFT_gbDauITS2On;
Int_t SFT_gbDauITS3On;
Int_t SFT_gbDauITS4On;
Int_t SFT_gbDauITS5On;
Bool_t SFT_gbDauSPDany;
Bool_t SFT_gbDauITSrefit;
Int_t SFT_gbmaxITSclusterShared;
Double_t SFT_gbmaxITSChi2;

Bool_t SFT_gbUntagDaughter;
Int_t SFT_gbPostMatched;
Double_t SFT_gbVertexZcut;


void AddTaskFlowStrangee(TString configFile, TString alienaddress, Bool_t skipTerminate=kFALSE) {
  Int_t ret = gSystem->Exec( Form("alien_cp %s/%s .",alienaddress.Data(),configFile.Data()) );
  printf("FlowStrange copying from grid %d\n",ret);
  AddTaskFlowStrangee(configFile,skipTerminate);
}
void AddTaskFlowStrangee(TString configFile, Bool_t skipTerminate=kFALSE) {
  SFT_ReadConfig(configFile);
  if(SFT_gbAllCC) {
    int centMin[8] = {00,05,10,20,30,40,50,60};
    int ncent=7;
    if(SFT_gbRunPP) {
      ncent=3;
      centMin[0]=00;
      centMin[1]=20;
      centMin[2]=40;
      centMin[3]=60;
    } else if(SFT_gbRunPA) {
      ncent=3;
      centMin[0]=00;
      centMin[1]=20;
      centMin[2]=40;
      centMin[3]=60;
    }
    TString antSuffix = SFT_gbSuffix;
    for(int cc=0; cc!=ncent; ++cc) {
      SFT_gbCentPerMin = centMin[cc];
      SFT_gbCentPerMax = centMin[cc+1];
      SFT_gbSuffix = Form("%s%02d%02d",antSuffix.Data(),SFT_gbCentPerMin,SFT_gbCentPerMax);
      AddTaskFlowStrangee(skipTerminate);
    }
  } else {
    AddTaskFlowStrangee(skipTerminate);
  }
}
void AddTaskFlowStrangee(Bool_t skipTerminate) {
  SFT_PrintConfig();

  TString fileName = AliAnalysisManager::GetCommonFileName();
  fileName.ReplaceAll(".root","");
  SFT_gbStamp = SFT_gbFolder + SFT_gbSuffix;

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();

  //-----------------STRANGE TASK----------------------------
  AliAnalysisTaskFlowStrangee *taskSel = new AliAnalysisTaskFlowStrangee(Form("FS_%s",SFT_gbStamp.Data()) );
  taskSel->SelectCollisionCandidates(SFT_gbTrigger);
  taskSel->SetReadESD(SFT_gbReadESD);
  taskSel->SetPostMatched(SFT_gbPostMatched);
  taskSel->SetReadMC(SFT_gbReadMC);
  taskSel->SetAvoidExec(SFT_gbAvoidExec);
  taskSel->SetSkipCentralitySelection(SFT_gbSkipCentrality);
  taskSel->SetSkipSelection(SFT_gbSkipSelection);
  taskSel->SetSkipVn(SFT_gbSkipVn);
  taskSel->SetExtraEventRejection(SFT_gbExtraEventCut);
  taskSel->SetCentralityRange(SFT_gbCentMethod,SFT_gbCentPerMin,SFT_gbCentPerMax);
  taskSel->SetSkipTerminate(skipTerminate);
  taskSel->SetHarmonic(SFT_gbHarmonic);
  Int_t npt = 24;
  Double_t ptbin[25] = { 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8,
                         2.0, 2.2, 2.4, 2.8, 3.2, 3.6, 4.0, 4.5, 5.0, 6.0,
                         8.0, 10., 12., 16., 20.};
  taskSel->SetPtEdges(npt,ptbin);
  if(SFT_gbRunPP) taskSel->Setpp();
  if(SFT_gbRunPA) taskSel->SetpA();
  taskSel->SetDebug(SFT_gbDebug);
  taskSel->SetK0L0(SFT_gbSpecie);
  taskSel->SetOnline( SFT_gbOnline );
  taskSel->SetMass( SFT_MassBins(SFT_gbSpecie),
                    SFT_MinMass(SFT_gbSpecie),
                    SFT_MaxMass(SFT_gbSpecie) );

  taskSel->SetWhichPsi(SFT_gbWhichPsi);
  taskSel->SetRFPFilterBit(SFT_gbRFPFilterBit);
  taskSel->SetRFPMinPt(SFT_gbRFPminPt);
  taskSel->SetRFPMaxPt(SFT_gbRFPmaxPt);
  taskSel->SetRFPAMinEta(SFT_gbRFPAminEta);
  taskSel->SetRFPAMaxEta(SFT_gbRFPAmaxEta);
  taskSel->SetRFPCMinEta(SFT_gbRFPCminEta);
  taskSel->SetRFPCMaxEta(SFT_gbRFPCmaxEta);
  taskSel->SetRFPTPCSignal(SFT_gbRFPTPCsignal);
  taskSel->SetRFPMaxIPxy(SFT_gbRFPmaxIPxy);
  taskSel->SetRFPMaxIPz(SFT_gbRFPmaxIPz);
  taskSel->SetRFPMinTPCCls(SFT_gbRFPTPCncls);

  taskSel->SetAddPiToMCReactionPlane(SFT_gbAddPitoMCRP);
  taskSel->SetDauUnTagProcedure(SFT_gbUntagDaughter);
  taskSel->SetVertexZcut(SFT_gbVertexZcut);

  taskSel->SetDauMinNClsTPC(SFT_gbMinNClsTPC);
  taskSel->SetDauMinNClsITS(SFT_gbMinNClsITS);
  taskSel->SetDauMinXRows(SFT_gbMinXRows);
  taskSel->SetDauMaxChi2PerNClsTPC(SFT_gbMaxChi2PerNClsTPC);
  taskSel->SetDauMinXRowsOverNClsFTPC(SFT_gbMinXRowsOverNClsFTPC);
  taskSel->SetDauMinEta(SFT_gbMinEta);
  taskSel->SetDauMaxEta(SFT_gbMaxEta);
  taskSel->SetDauMinPt(SFT_gbMinPt);
  taskSel->SetDauMinImpactParameterXY(SFT_gbMinImpactParameterXY);
  taskSel->SetDauMaxNSigmaPID(SFT_gbMaxNSigmaPID);
  taskSel->SetDauITSLayer(0,SFT_gbDauITS0On);
  taskSel->SetDauITSLayer(1,SFT_gbDauITS1On);
  taskSel->SetDauITSLayer(2,SFT_gbDauITS2On);
  taskSel->SetDauITSLayer(3,SFT_gbDauITS3On);
  taskSel->SetDauITSLayer(4,SFT_gbDauITS4On);
  taskSel->SetDauITSLayer(5,SFT_gbDauITS5On);
  taskSel->SetDauSPDRequireAny(SFT_gbDauSPDany);
  taskSel->SetDauITSrefit(SFT_gbDauITSrefit);

  //newITScuts
  taskSel->SetMaxSharedITSCluster(SFT_gbmaxITSclusterShared);
  taskSel->SetMaxChi2perITSCluster(SFT_gbmaxITSChi2);
  
  taskSel->SetMaxRapidity(SFT_gbMaxRapidity);
  taskSel->SetMaxDCAdaughters(SFT_gbMaxDCAdaughters);
  taskSel->SetMinCosinePointingAngleXY(SFT_gbMinCosinePointingAngleXY);
  taskSel->SetMinQt(SFT_gbMinQt,SFT_gbQtPie);
  taskSel->SetStopPIDAtPt(SFT_gbPIDPt);
  taskSel->SetMinRadXY(SFT_gbMinRadXY);
  taskSel->SetMaxDecayLength(SFT_gbMaxDecayLength);
  taskSel->SetMaxProductIPXY(SFT_gbMaxProductIPXY);
  taskSel->SetMinEta(SFT_gbMinEta);
  taskSel->SetMaxEta(SFT_gbMaxEta);
  taskSel->SetMinPt(SFT_gbMinPt);
  taskSel->SetUseFlowPackage(SFT_gbFlowPackage);

  taskSel->SetQAlevel(SFT_gbQA);
  /*
  AliFlowTrackCuts* cutsRP = new AliFlowTrackCuts(Form("RFPcuts%s",SFT_gbSuffix));
  if(!cutsRP) {
    //if(debug) cout << " Fatal error: no RP cuts found, could be a library problem! " << endl;
    return 0x0;
  }
  cutsRP = cutsRP->GetStandardVZEROOnlyTrackCuts(); // select vzero tracks
  cutsRP->SetVZEROgainEqualizationPerRing(kFALSE);
  cutsRP->SetApplyRecentering(kTRUE);
  taskSel->SetRPCuts(cutsRP);
  */
  
  if(SFT_gbVZEload.Length()>4) {
    TFile *ocalib = TFile::Open(SFT_gbVZEload);
    if(ocalib->IsOpen()) {
      TList *vzero = ocalib->Get("VZECALIB");
      taskSel->LoadVZEResponse(vzero,SFT_gbVZEmb,SFT_gbVZEpdisk);
    } else {
      printf("ADDTASKFLOWSTRANGEE COULD NOT OPEN %s. NO VZE CALIBRATION LOADED!\n",SFT_gbVZEload.Data());
    }
  }
  printf("Loading %d %d %d %d as VZE configuration\n",SFT_gbV0CRingMin, SFT_gbV0CRingMax, SFT_gbV0ARingMin, SFT_gbV0ARingMax);
  taskSel->SetRFPVZERingRange( SFT_gbV0CRingMin, SFT_gbV0CRingMax, SFT_gbV0ARingMin, SFT_gbV0ARingMax );
  taskSel->SetStoreVZEResponse(SFT_gbVZEsave);
  

  AliAnalysisDataContainer *cOutHist = mgr->CreateContainer(Form("FS_OH_%s",SFT_gbStamp.Data()),
							    TList::Class(),
							    AliAnalysisManager::kOutputContainer,
							    Form("%s.root:Selector_%s",fileName.Data(),
								 SFT_gbFolder.Data()));
  AliAnalysisDataContainer *exc_TPC = mgr->CreateContainer( Form("FS_TPC_%s",SFT_gbStamp.Data()),
                                                            AliFlowEventSimple::Class(),
                                                            AliAnalysisManager::kExchangeContainer );
  AliAnalysisDataContainer *exc_VZE = mgr->CreateContainer( Form("FS_VZE_%s",SFT_gbStamp.Data()),
                                                            AliFlowEventSimple::Class(),
                                                            AliAnalysisManager::kExchangeContainer );
  mgr->AddTask(taskSel);
  mgr->ConnectInput (taskSel,0,cinput1);
  mgr->ConnectOutput(taskSel,1,cOutHist);
  mgr->ConnectOutput(taskSel,2,exc_TPC);
  mgr->ConnectOutput(taskSel,3,exc_VZE);

  if(!SFT_gbFlowPackage) return;
  if( (!SFT_gbQCTPC) && (!SFT_gbSPVZE) && (!SFT_gbSPTPC) ) return;
  //-------------------FLOWPACKAGE TASKS----------------------------
  AliFlowTrackSimpleCuts *filter[20], *filterhf[20][2]; // MASS BANDS
  int mbs = SFT_MassBands(SFT_gbSpecie);
  if(SFT_gbPostMatched) mbs = 1;
  for(int mb=0; mb!=mbs; ++mb) {
    filter[mb] = new AliFlowTrackSimpleCuts( Form("Filter_MB%d",mb) );
    filter[mb]->SetEtaMin( -0.8 ); filter[mb]->SetEtaMax( +0.8 );
    Double_t minmass = SFT_MassBandLowEdge(SFT_gbSpecie,mb);
    Double_t maxmass = SFT_MassBandLowEdge(SFT_gbSpecie,mb+1);
    if(SFT_gbPostMatched) maxmass = SFT_MassBandLowEdge(SFT_gbSpecie,SFT_MassBands(SFT_gbSpecie));
    filter[mb]->SetMassMin( minmass ); filter[mb]->SetMassMax( maxmass );
    //half window for POIs
    filterhf[mb][0] = new AliFlowTrackSimpleCuts( Form("Filterhf0_MB%d",mb) );
    filterhf[mb][0]->SetEtaMin( +0.0 ); filterhf[mb][0]->SetEtaMax( +0.8 );
    filterhf[mb][0]->SetMassMin( minmass ); filterhf[mb][0]->SetMassMax( maxmass );
    filterhf[mb][1] = new AliFlowTrackSimpleCuts( Form("Filterhf1_MB%d",mb) );
    filterhf[mb][1]->SetEtaMin( -0.8 ); filterhf[mb][1]->SetEtaMax( -0.0 );
    filterhf[mb][1]->SetMassMin( minmass ); filterhf[mb][1]->SetMassMax( maxmass );
    if(SFT_gbQCTPC) {
      SFT_AddQCmethod( Form("QCTPCMB%d",mb), exc_TPC, filter[mb]); // QC TPC
    }
    if(SFT_gbSPTPC) {
      SFT_AddSPmethod( Form("SPTPCMB%d",mb), exc_TPC, filterhf[mb][0], "Qa" ); // SP TPC Qa
      SFT_AddSPmethod( Form("SPTPCMB%d",mb), exc_TPC, filterhf[mb][1], "Qb" ); // SP TPC Qb
      SFT_AddSPmethod( Form("SPTPC2MB%d",mb), exc_TPC, filterhf[mb][0], "Qa", 0.2 ); // SP TPC Qa
      SFT_AddSPmethod( Form("SPTPC2MB%d",mb), exc_TPC, filterhf[mb][1], "Qb", 0.2 ); // SP TPC Qb
      SFT_AddSPmethod( Form("SPTPC4MB%d",mb), exc_TPC, filterhf[mb][0], "Qa", 0.4 ); // SP TPC Qa
      SFT_AddSPmethod( Form("SPTPC4MB%d",mb), exc_TPC, filterhf[mb][1], "Qb", 0.4 ); // SP TPC Qb
    }
    if(SFT_gbMCEP) {
      SFT_AddMCEPmethod( Form("MCEPMB%d",mb), exc_TPC, filter[mb]); // MCEP TPC
    }
    if(SFT_gbSPVZE) {
      if(SFT_gbSPVZEhalf) {
	SFT_AddSPmethod( Form("SPVZEMB%d",mb), exc_VZE, filterhf[mb][0], "Qa", 1.0 ); // SP VZE Qa
	SFT_AddSPmethod( Form("SPVZEMB%d",mb), exc_VZE, filterhf[mb][1], "Qb", 1.0 ); // SP VZE Qa
      } else {
	SFT_AddSPmethod( Form("SPVZEMB%d",mb), exc_VZE, filter[mb], "Qa", 1.0 ); // SP VZE Qa
	SFT_AddSPmethod( Form("SPVZEMB%d",mb), exc_VZE, filter[mb], "Qb", 1.0 ); // SP VZE Qa
      }
    }
  }
}
void SFT_AddMCEPmethod(char *name, AliAnalysisDataContainer *flowEvent, AliFlowTrackSimpleCuts *cutsPOI=NULL) {
  TString fileName = AliAnalysisManager::GetCommonFileName();
  TString myFolder = Form("%sv%d",SFT_gbFolder.Data(),SFT_gbHarmonic);
  TString myName = Form("%sv%d_%s",name,SFT_gbHarmonic,SFT_gbSuffix.Data());
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
                                                          Form("%s:FlowStrange_MCEP_%s",fileName.Data(),myFolder.Data()) );
  AliAnalysisTaskMCEventPlane *tskQC = new AliAnalysisTaskMCEventPlane( Form("TaskMCEP_%s",myName.Data()) );
  tskQC->SetHarmonic(SFT_gbHarmonic);
  mgr->AddTask(tskQC);
  mgr->ConnectInput( tskQC,0,flowEvent2);
  mgr->ConnectOutput(tskQC,1,outQC);
}
void SFT_AddQCmethod(char *name, AliAnalysisDataContainer *flowEvent, AliFlowTrackSimpleCuts *cutsPOI=NULL) {
  TString fileName = AliAnalysisManager::GetCommonFileName();
  TString myFolder = Form("%sv%d",SFT_gbFolder.Data(),SFT_gbHarmonic);
  TString myName = Form("%sv%d_%s",name,SFT_gbHarmonic,SFT_gbSuffix.Data());
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
  tskQC->SetHarmonic(SFT_gbHarmonic);
  tskQC->SetBookOnlyBasicCCH(SFT_gbShrinkFP);
  mgr->AddTask(tskQC);
  mgr->ConnectInput( tskQC,0,flowEvent2);
  mgr->ConnectOutput(tskQC,1,outQC);
}
void SFT_AddSPmethod(char *name, AliAnalysisDataContainer *flowEvent, AliFlowTrackSimpleCuts *cutsPOI=NULL, char *Qvector, Double_t gap=0.0) {
  TString fileName = AliAnalysisManager::GetCommonFileName();
  TString myFolder = Form("%sv%d",SFT_gbFolder.Data(),SFT_gbHarmonic);
  TString myNameSP = Form("%sv%d%s_%s",name,SFT_gbHarmonic,Qvector,SFT_gbSuffix.Data());
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisDataContainer *flowEvent2 = mgr->CreateContainer( Form("Filter_%s", myNameSP.Data()),
                                                               AliFlowEventSimple::Class(),
                                                               AliAnalysisManager::kExchangeContainer );
  AliAnalysisTaskFilterFE *tskFilter = new AliAnalysisTaskFilterFE( Form("TaskFilter_%s",myNameSP.Data()),
                                                                    NULL, cutsPOI);
  tskFilter->SetSubeventEtaRange( -5.0, -gap, +gap, +5.0 );
  mgr->AddTask(tskFilter);
  mgr->ConnectInput( tskFilter,0,flowEvent);
  mgr->ConnectOutput(tskFilter,1,flowEvent2);
  AliAnalysisDataContainer *outSP = mgr->CreateContainer( myNameSP.Data(),TList::Class(),AliAnalysisManager::kOutputContainer,
                                                          Form("%s:FlowStrange_SP_%s",fileName.Data(),myFolder.Data()) );
  AliAnalysisTaskScalarProduct *tskSP = new AliAnalysisTaskScalarProduct( Form("TaskScalarProduct_%s",myNameSP.Data()),kFALSE);
  tskSP->SetApplyCorrectionForNUA(kTRUE);
  tskSP->SetHarmonic(SFT_gbHarmonic);
  tskSP->SetTotalQvector(Qvector);
  tskSP->SetBookOnlyBasicCCH(SFT_gbShrinkFP);
  mgr->AddTask(tskSP);
  mgr->ConnectInput( tskSP,0,flowEvent2);
  mgr->ConnectOutput(tskSP,1,outSP);
}
double SFT_MassBandLowEdge( int nv0, int mb ) {
  if(nv0>10&&mb==0) return -5;
  if(nv0>10&&mb==1) return +5;
  switch(nv0) {
  case(0):
    double lowEdge[14]={0.398, 0.420, 0.444, 0.468, 0.486,
                         0.490, 0.494, 0.498, 0.502, 0.506, 
                         0.524, 0.548, 0.572, 0.598};
    break;
  default:
    double lowEdge[13]={1.084, 1.094, 1.104, 1.110, 1.114,
			1.116, 1.118, 1.122, 1.128, 1.138,
			1.148, 1.158, 1.168};
    break;
  }
  return lowEdge[mb];
}
int SFT_MassBands( int nv0 ) {
  int bands=1;
  switch(nv0) {
  case(0):
    bands = 13;
    break;
  default:
    bands = 12;
  }
  if(nv0>10) bands=1;
  return bands;
}
int SFT_MassBins( int nv0 ) {
  int bins=100;
  switch(nv0) {
  case(0)://kZERO
    bins=100;
    break;
  default://LAMBDA
    bins=84;
    break;
  }
  if(nv0>10) bins=100;//CHARGED
  return bins;
}
double SFT_MinMass( int nv0 ) {
  return SFT_MassBandLowEdge( nv0, 0 );
}
double SFT_MaxMass( int nv0 ) {
  return SFT_MassBandLowEdge( nv0, SFT_MassBands(nv0) );
}
void SFT_PrintConfig() {
  printf("***********************************\n");
  printf("* STRANGE FLOW TASK CONFIGURATION *\n");
  printf("* SUFFIX  %s *\n", SFT_gbSuffix.Data() );
  printf("* TRIGGER  %d *\n", SFT_gbTrigger );
  printf("* RUNPP  %d *\n", SFT_gbRunPP );
  printf("* RUNPA  %d *\n", SFT_gbRunPA );
  printf("* AVOIDEXEC  %d *\n", SFT_gbAvoidExec );
  printf("* SKIPCENTRALITY  %d *\n", SFT_gbSkipCentrality );
  printf("* ESD  %d *\n", SFT_gbReadESD );
  printf("* MC  %d *\n", SFT_gbReadMC );
  printf("* HARMONIC  %d *\n", SFT_gbHarmonic );
  printf("* ADDPITOMCEP  %d *\n", SFT_gbAddPitoMCRP );
  printf("* POSTMATCHED  %d *\n", SFT_gbPostMatched );
  printf("* EXTRAEVENTCUT  %d *\n", SFT_gbExtraEventCut );
  printf("* CENTMETHOD  %s *\n", SFT_gbCentMethod.Data() );
  printf("* CENTPERMIN  %d *\n", SFT_gbCentPerMin );
  printf("* CENTPERMAX  %d *\n", SFT_gbCentPerMax );
  printf("* VERTEXZ  %f *\n", SFT_gbVertexZcut );
  printf("* SPECIE  %d *\n", SFT_gbSpecie );
  printf("* HOMEMADE  %d *\n", SFT_gbHomemade );
  printf("* ONLINE  %d *\n", SFT_gbOnline );
  printf("* MINNCLSTTPC  %d *\n", SFT_gbMinNClsTPC );
  printf("* MINNCLSITS  %d *\n", SFT_gbMinNClsITS );
  printf("* MINXROWS  %d *\n", SFT_gbMinXRows );
  printf("* MAXCHI2NCLSTPC  %f *\n", SFT_gbMaxChi2PerNClsTPC );
  printf("* MINXROWSNFCLSTPC  %f *\n", SFT_gbMinXRowsOverNClsFTPC );
  printf("* MINETA  %f *\n", SFT_gbMinEta );
  printf("* MAXETA  %f *\n", SFT_gbMaxEta );
  printf("* MINPT  %f *\n", SFT_gbMinPt );
  printf("* UNTAG  %f *\n", SFT_gbUntagDaughter );
  printf("* MIND0XY  %f *\n", SFT_gbMinImpactParameterXY );
  printf("* MAXSIGMAPID  %f *\n", SFT_gbMaxNSigmaPID );
  printf("* MAXY  %f *\n", SFT_gbMaxRapidity );
  printf("* MAXDCA  %f *\n", SFT_gbMaxDCAdaughters );
  printf("* MINCTP  %f *\n", SFT_gbMinCosinePointingAngleXY );
  printf("* MINQT  %f *\n", SFT_gbMinQt );
  printf("* QTPIE  %f *\n", SFT_gbQtPie );
  printf("* STOPPID  %f *\n", SFT_gbPIDPt );
  printf("* MINRADXY  %f *\n", SFT_gbMinRadXY );
  printf("* MAXDL  %f *\n", SFT_gbMaxDecayLength );
  printf("* D0D0XY  %f *\n", SFT_gbMaxProductIPXY );
  printf("* DEBUG  %d *\n", SFT_gbDebug );
  printf("* QA  %d *\n", SFT_gbQA );
  printf("* SKIPSELECTION  %d *\n", SFT_gbSkipSelection );
  printf("* SKIPVN  %d *\n", SFT_gbSkipVn );
  printf("* USEFP  %d *\n", SFT_gbFlowPackage );
  printf("* SPVZE  %d *\n", SFT_gbSPVZE );
  printf("* SPVZEHALF  %d *\n", SFT_gbSPVZEhalf );
  printf("* SPTPC  %d *\n", SFT_gbSPTPC );
  printf("* QCTPC  %d *\n", SFT_gbQCTPC );
  printf("* MCEP  %d *\n", SFT_gbMCEP );
  printf("* SHRINKFP  %d *\n", SFT_gbShrinkFP );
  printf("* RFFILTERBIT  %d *\n", SFT_gbRFPFilterBit );
  printf("* RFMINPT  %f *\n", SFT_gbRFPminPt );
  printf("* RFMAXPT  %f *\n", SFT_gbRFPmaxPt );
  printf("* RFCMINETA  %f *\n", SFT_gbRFPCminEta );
  printf("* RFCMAXETA  %f *\n", SFT_gbRFPCmaxEta );
  printf("* RFAMINETA  %f *\n", SFT_gbRFPAminEta );
  printf("* RFAMAXETA  %f *\n", SFT_gbRFPAmaxEta );
  printf("* RFTPCSIGNAL  %f *\n", SFT_gbRFPTPCsignal );
  printf("* RFMAXIPXY  %f *\n", SFT_gbRFPmaxIPxy );
  printf("* RFMAXIPZ  %f *\n", SFT_gbRFPmaxIPz );
  printf("* RFTPCNCLS  %d *\n", SFT_gbRFPTPCncls );
  printf("* RFVZEC_RingMin  %d *\n", SFT_gbV0CRingMin );
  printf("* RFVZEC_RingMax  %d *\n", SFT_gbV0CRingMax );
  printf("* RFVZEA_RingMin  %d *\n", SFT_gbV0ARingMin );
  printf("* RFVZEA_RingMax  %d *\n", SFT_gbV0ARingMax );
  printf("* WHICHPSI  %d *\n", SFT_gbWhichPsi );
  printf("* VZELOAD  %s *\n", SFT_gbVZEload.Data() );
  printf("* VZELINEAR  %d *\n", SFT_gbVZEmb );
  printf("* VZEPERDISK  %d *\n", SFT_gbVZEpdisk );
  printf("* VZESAVE  %d *\n", SFT_gbVZEsave );
  printf("* DAUITS0  %d *\n", SFT_gbDauITS0On );
  printf("* DAUITS1  %d *\n", SFT_gbDauITS1On );
  printf("* DAUITS2  %d *\n", SFT_gbDauITS2On );
  printf("* DAUITS3  %d *\n", SFT_gbDauITS3On );
  printf("* DAUITS4  %d *\n", SFT_gbDauITS4On );
  printf("* DAUITS5  %d *\n", SFT_gbDauITS5On );
  printf("* DAUSPDANY  %d *\n", SFT_gbDauSPDany );
  printf("* DAUITSrefit  %d *\n", SFT_gbDauITSrefit );

  printf("* DAUmaxITSclusterShared  %d *\n", SFT_gbmaxITSclusterShared );
  printf("* DAUmaxITSChi2  %f *\n", SFT_gbmaxITSChi2 );
  printf("***********************************\n");
}
void SFT_ReadConfig(TString ipf) {
  SFT_ResetVars();
  printf("Reading %s\n",ipf.Data());
  ifstream input(ipf.Data());
  TString varname;
  Double_t vardouble;
  Int_t varint;
  UInt_t varuint;
  Bool_t varbool;
  for(;input.good();) {
    input >> varname;
    if(!input.good()) {
      break;
    } else if(!varname.CompareTo("SUFFIX")) {
      input >> SFT_gbSuffix;
    } else if(!varname.CompareTo("TRIGGER")) {
      input >> SFT_gbTrigger;
    } else if(!varname.CompareTo("RUNPP")) {
      input >> SFT_gbRunPP;
    } else if(!varname.CompareTo("RUNPA")) {
      input >> SFT_gbRunPA;
    } else if(!varname.CompareTo("AVOIDEXEC")) {
      input >> SFT_gbAvoidExec;
    } else if(!varname.CompareTo("SKIPCENTRALITY")) {
      input >> SFT_gbSkipCentrality;
    } else if(!varname.CompareTo("ESD")) {
      input >> SFT_gbReadESD;
    } else if(!varname.CompareTo("MC")) {
      input >> SFT_gbReadMC;
    } else if(!varname.CompareTo("EXTRAEVENTCUT")) {
      input >> SFT_gbExtraEventCut;
    } else if(!varname.CompareTo("CENTMETHOD")) {
      input >> SFT_gbCentMethod;
    } else if(!varname.CompareTo("CENTPERMIN")) {
      input >> SFT_gbCentPerMin;
    } else if(!varname.CompareTo("CENTPERMAX")) {
      input >> SFT_gbCentPerMax;
    } else if(!varname.CompareTo("SPECIE")) {
      input >> SFT_gbSpecie;
    } else if(!varname.CompareTo("HOMEMADE")) {
      input >> SFT_gbHomemade;
    } else if(!varname.CompareTo("ONLINE")) {
      input >> SFT_gbOnline;
    } else if(!varname.CompareTo("MINNCLSTTPC")) {
      input >> SFT_gbMinNClsTPC;
    } else if(!varname.CompareTo("MINNCLSITS")) {
      input >> SFT_gbMinNClsITS;
    } else if(!varname.CompareTo("MINXROWS")) {
      input >> SFT_gbMinXRows;
    } else if(!varname.CompareTo("MAXCHI2NCLSTPC")) {
      input >> SFT_gbMaxChi2PerNClsTPC;
    } else if(!varname.CompareTo("MINXROWSNFCLSTPC")) {
      input >> SFT_gbMinXRowsOverNClsFTPC;
    } else if(!varname.CompareTo("MINETA")) {
      input >> SFT_gbMinEta;
    } else if(!varname.CompareTo("MAXETA")) {
      input >> SFT_gbMaxEta;
    } else if(!varname.CompareTo("MINPT")) {
      input >> SFT_gbMinPt;
    } else if(!varname.CompareTo("MIND0XY")) {
      input >> SFT_gbMinImpactParameterXY;
    } else if(!varname.CompareTo("MAXSIGMAPID")) {
      input >> SFT_gbMaxNSigmaPID;
    } else if(!varname.CompareTo("MAXY")) {
      input >> SFT_gbMaxRapidity;
    } else if(!varname.CompareTo("MAXDCA")) {
      input >> SFT_gbMaxDCAdaughters;
    } else if(!varname.CompareTo("MINCTP")) {
      input >> SFT_gbMinCosinePointingAngleXY;
    } else if(!varname.CompareTo("MINQT")) {
      input >> SFT_gbMinQt;
    } else if(!varname.CompareTo("QTPIE")) {
      input >> SFT_gbQtPie;
    } else if(!varname.CompareTo("STOPPID")) {
      input >> SFT_gbPIDPt;
    } else if(!varname.CompareTo("MINRADXY")) {
      input >> SFT_gbMinRadXY;
    } else if(!varname.CompareTo("MAXDL")) {
      input >> SFT_gbMaxDecayLength;
    } else if(!varname.CompareTo("D0D0XY")) {
      input >> SFT_gbMaxProductIPXY;
    } else if(!varname.CompareTo("DEBUG")) {
      input >> SFT_gbDebug;
    } else if(!varname.CompareTo("QA")) {
      input >> SFT_gbQA;
    } else if(!varname.CompareTo("SKIPSELECTION")) {
      input >> SFT_gbSkipSelection;
    } else if(!varname.CompareTo("SKIPVN")) {
      input >> SFT_gbSkipVn;
    } else if(!varname.CompareTo("USEFP")) {
      input >> SFT_gbFlowPackage;
    } else if(!varname.CompareTo("SPVZE")) {
      input >> SFT_gbSPVZE;
    } else if(!varname.CompareTo("SPTPC")) {
      input >> SFT_gbSPTPC;
    } else if(!varname.CompareTo("SPVZEHALF")) {      
      input >> SFT_gbSPVZEhalf;
    } else if(!varname.CompareTo("QCTPC")) {
      input >> SFT_gbQCTPC;
    } else if(!varname.CompareTo("MCEP")) {
      input >> SFT_gbMCEP;
    } else if(!varname.CompareTo("ADDPITOMCEP")) {
      input >> SFT_gbAddPitoMCRP;
    } else if(!varname.CompareTo("SHRINKFP")) {
      input >> SFT_gbShrinkFP;
    } else if(!varname.CompareTo("RFFILTERBIT")) {
      input >> SFT_gbRFPFilterBit;
    } else if(!varname.CompareTo("RFMINPT")) {
      input >> SFT_gbRFPminPt;
    } else if(!varname.CompareTo("RFMAXPT")) {
      input >> SFT_gbRFPmaxPt;
    } else if(!varname.CompareTo("RFCMINETA")) {
      input >> SFT_gbRFPCminEta;
    } else if(!varname.CompareTo("RFCMAXETA")) {
      input >> SFT_gbRFPCmaxEta;
    } else if(!varname.CompareTo("RFAMINETA")) {
      input >> SFT_gbRFPAminEta;
    } else if(!varname.CompareTo("RFAMAXETA")) {
      input >> SFT_gbRFPAmaxEta;
    } else if(!varname.CompareTo("RFTPCSIGNAL")) {
      input >> SFT_gbRFPTPCsignal;
    } else if(!varname.CompareTo("RFMAXIPXY")) {
      input >> SFT_gbRFPmaxIPxy;
    } else if(!varname.CompareTo("RFMAXIPZ")) {
      input >> SFT_gbRFPmaxIPz;
    } else if(!varname.CompareTo("RFTPCNCLS")) {
      input >> SFT_gbRFPTPCncls;
    } else if(!varname.CompareTo("VZELOAD")) {
      input >> SFT_gbVZEload;
    } else if(!varname.CompareTo("VZELINEAR")) {
      input >> SFT_gbVZEmb;
    } else if(!varname.CompareTo("VZEPERDISK")) {
      input >> SFT_gbVZEpdisk;
    } else if(!varname.CompareTo("VZESAVE")) {
      input >> SFT_gbVZEsave;
    } else if(!varname.CompareTo("ALLCC")) {
      input >> SFT_gbAllCC;
    } else if(!varname.CompareTo("UNTAG")) {
      input >> SFT_gbUntagDaughter;
    } else if(!varname.CompareTo("POSTMATCHED")) {
      input >> SFT_gbPostMatched;
    } else if(!varname.CompareTo("VERTEXZ")) {
      input >> SFT_gbVertexZcut;
    } else if(!varname.CompareTo("WHICHPSI")) {
      input >> SFT_gbWhichPsi;
    } else if(!varname.CompareTo("DAUITS0")) {
      input >> SFT_gbDauITS0On;
    } else if(!varname.CompareTo("DAUITS1")) {
      input >> SFT_gbDauITS1On;
    } else if(!varname.CompareTo("DAUITS2")) {
      input >> SFT_gbDauITS2On;
    } else if(!varname.CompareTo("DAUITS3")) {
      input >> SFT_gbDauITS3On;
    } else if(!varname.CompareTo("DAUITS4")) {
      input >> SFT_gbDauITS4On;
    } else if(!varname.CompareTo("DAUITS5")) {
      input >> SFT_gbDauITS5On;
    } else if(!varname.CompareTo("DAUSPDANY")) {
      input >> SFT_gbDauSPDany;
    } else if(!varname.CompareTo("DAUITSREFIT")) {
      input >> SFT_gbDauITSrefit;
    } else if(!varname.CompareTo("maxITSclusterShared")) {
      input >> SFT_gbmaxITSclusterShared;
    } else if(!varname.CompareTo("maxITSChi2")) {
      input >> SFT_gbmaxITSChi2;
    } else if(!varname.CompareTo("RFVZEC_RingMin")) {
      input >> SFT_gbV0CRingMin;
    } else if(!varname.CompareTo("RFVZEC_RingMax")) {
      input >> SFT_gbV0CRingMax;
    } else if(!varname.CompareTo("RFVZEA_RingMin")) {
      input >> SFT_gbV0ARingMin;
    } else if(!varname.CompareTo("RFVZEA_RingMax")) {
      input >> SFT_gbV0ARingMax;
    } else if(!varname.CompareTo("HARMONIC")) {
      input >> SFT_gbHarmonic;
    } else {
      printf("I dont understand %s\n",varname.Data());
    }
  }
}
void SFT_ResetVars() {
  SFT_gbTrigger=1;
  SFT_gbReadESD=0;
  SFT_gbReadMC=0;
  SFT_gbMatchMC=0;
  SFT_gbAvoidExec=0;
  SFT_gbExtraEventCut=0;
  SFT_gbCentMethod="V0MTRK";
  SFT_gbCentPerMin=0;
  SFT_gbCentPerMax=100;
  SFT_gbRunPP=0;
  SFT_gbRunPA=0;
  SFT_gbSpecie=0;
  SFT_gbHomemade=0;
  SFT_gbOnline=0;
  SFT_gbMinNClsTPC=70;
  SFT_gbMinNClsITS=-1;
  SFT_gbMinXRows=70;
  SFT_gbMaxChi2PerNClsTPC=4.0;
  SFT_gbMinXRowsOverNClsFTPC=0.8;
  SFT_gbMinEta=-0.8;
  SFT_gbMaxEta=+0.8;
  SFT_gbMinPt=0.1;
  SFT_gbMinImpactParameterXY=0.1;
  SFT_gbMaxNSigmaPID=3.0;
  SFT_gbMaxRapidity=0.5;
  SFT_gbMaxDCAdaughters=1.0;
  SFT_gbMinCosinePointingAngleXY=0.998;
  SFT_gbMinQt=0.2;
  SFT_gbQtPie=kTRUE;
  SFT_gbPIDPt=3.0;
  SFT_gbMinRadXY=5.0;
  SFT_gbMaxDecayLength=3.0;
  SFT_gbMaxProductIPXY=0.0;
  SFT_gbDebug=0;
  SFT_gbQA=0;
  SFT_gbFolder="FlowStrange";
  SFT_gbSuffix="NOTFOUND";
  SFT_gbRFPFilterBit=1;
  SFT_gbRFPminPt=0.2;
  SFT_gbRFPmaxPt=5.0;
  SFT_gbRFPCminEta=-0.8;
  SFT_gbRFPCmaxEta=0.0;
  SFT_gbRFPAminEta=0.0;
  SFT_gbRFPAmaxEta=+0.8;
  SFT_gbRFPTPCsignal=10;
  SFT_gbRFPmaxIPxy=2.4;
  SFT_gbRFPmaxIPz=3.2;
  SFT_gbRFPTPCncls=70;
  SFT_gbV0CRingMin=0;
  SFT_gbV0CRingMax=3;
  SFT_gbV0ARingMin=0;
  SFT_gbV0ARingMax=3;
  SFT_gbAllCC=0;
  SFT_gbSkipSelection=0;
  SFT_gbSkipVn=0;
  SFT_gbWhichPsi=1;
  SFT_gbFlowPackage=0;
  SFT_gbSPVZE=0;
  SFT_gbSPVZEhalf=kFALSE;
  SFT_gbSPTPC=0;
  SFT_gbQCTPC=0;
  SFT_gbMCEP=0;
  SFT_gbHarmonic=2;
  SFT_gbVZEload="no";
  SFT_gbVZEsave=0;
  SFT_gbVZEmb=0;
  SFT_gbVZEpdisk=0;
  SFT_gbUntagDaughter=kTRUE;
  SFT_gbPostMatched=kFALSE;
  SFT_gbShrinkFP=kTRUE;
  SFT_gbVertexZcut=10.0;
  SFT_gbDauITS0On=-1;
  SFT_gbDauITS1On=-1;
  SFT_gbDauITS2On=-1;
  SFT_gbDauITS3On=-1;
  SFT_gbDauITS4On=-1;
  SFT_gbDauITS5On=-1;
  SFT_gbSkipCentrality=kFALSE;
  SFT_gbAddPitoMCRP=kFALSE;
  SFT_gbDauSPDany=kFALSE;
  SFT_gbDauITSrefit=kFALSE;
  SFT_gbmaxITSclusterShared=6;
  SFT_gbmaxITSChi2=36.;
}
