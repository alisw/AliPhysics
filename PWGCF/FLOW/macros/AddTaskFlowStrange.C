UInt_t gbTrigger;
Bool_t gbReadESD;
Bool_t gbReadMC;
Int_t gbMatchMC;
Bool_t gbAvoidExec;
Bool_t gbExtraEventCut=kFALSE;
TString gbCentMethod;
Int_t gbCentPerMin,gbCentPerMax;
Bool_t gbRunPP;
Bool_t gbRunPA;
Int_t gbSpecie;
Bool_t gbHomemade;
Bool_t gbOnline;
Int_t gbMinNClsTPC;
Int_t gbMinXRows;
Double_t gbMaxChi2PerNClsTPC;
Double_t gbMinXRowsOverNClsFTPC;
Double_t gbMinEta;
Double_t gbMaxEta;
Double_t gbMinPt;
Double_t gbMinImpactParameterXY;
Double_t gbMaxNSigmaPID;
Double_t gbMaxRapidity;
Double_t gbMaxDCAdaughters;
Double_t gbMinCosinePointingAngleXY;
Double_t gbMinQt;
Bool_t   gbQtPie=kTRUE;
Double_t gbMinRadXY;
Double_t gbMaxDecayLength;
Double_t gbMaxProductIPXY;
Int_t gbDebug=0;
Int_t gbQA;
TString gbFolder="FlowStrange";
TString gbSuffix="kze";
TString gbStamp;
Int_t gbRFPFilterBit=1;
Double_t gbRFPminPt=0.2;
Double_t gbRFPmaxPt=5.0;
Double_t gbRFPminEta=-0.8;
Double_t gbRFPmaxEta=+0.8;
Double_t gbRFPTPCsignal=10;
Double_t gbRFPmaxIPxy=2.4;
Double_t gbRFPmaxIPz=3.2;
Int_t gbRFPTPCncls=70;

Bool_t gbAllCC=0;
Bool_t gbSkipSelection=0;
Bool_t gbSkipFlow=0;
Int_t gbWhichPsi=2;
Bool_t gbFlowPackage=1;
Bool_t gbSPVZE=1;//1
Bool_t gbSPTPC=1;
Bool_t gbQCTPC=1;//1
Int_t gbHarmonic=2;
TString gbVZEload="";
Bool_t gbVZEsave=0;


void AddTaskFlowStrange(TString configFile) {
  SFT_ReadConfig(configFile);
  if(gbAllCC) {
    int centMin[9] = {00,05,10,20,30,40,50,60,70};
    int centMax[9] = {05,10,20,30,40,50,60,70,80};
    int ncent=9;
    if(gbRunPP) {
      ncent=3;
      centMin[0]=10; centMax[0]=30;
      centMin[1]=30; centMax[1]=50;
      centMin[2]=0; centMax[2]=100;
    } else if(gbRunPA) {
      ncent=4;
      centMin[0]=00; centMax[0]=20;
      centMin[1]=20; centMax[1]=40;
      centMin[2]=40; centMax[2]=60;
      centMin[3]=60; centMax[3]=80;
    } else {
    }
    TString antSuffix = gbSuffix;
    for(int cc=0; cc!=ncent; ++cc) {
      gbCentPerMin = centMin[cc];
      gbCentPerMax = centMax[cc];
      gbSuffix = Form("%s%d%d",antSuffix.Data(),gbCentPerMin,gbCentPerMax);
      AddTaskFlowStrange();
    }
  } else {
    AddTaskFlowStrange();
  }
}
void AddTaskFlowStrange() {
  SFT_PrintConfig();

  TString fileName = AliAnalysisManager::GetCommonFileName();
  fileName.ReplaceAll(".root","");
  gbStamp = gbFolder + gbSuffix;

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();

  //-----------------STRANGE TASK----------------------------
  AliAnalysisTaskFlowStrange *taskSel = new AliAnalysisTaskFlowStrange(Form("FS_%s",gbStamp.Data()) );
  taskSel->SelectCollisionCandidates(gbTrigger);
  taskSel->SetReadESD(gbReadESD);
  taskSel->SetReadMC(gbReadMC);
  taskSel->SetAvoidExec(gbAvoidExec);
  taskSel->SetSkipSelection(gbSkipSelection);
  taskSel->SetSkipFlow(gbSkipFlow);
  taskSel->SetExtraEventRejection(gbExtraEventCut);
  taskSel->SetCentralityRange(gbCentMethod,gbCentPerMin,gbCentPerMax);
  if(gbRunPP) taskSel->Setpp();
  if(gbRunPA) taskSel->SetpA();
  taskSel->SetDebug(gbDebug);
  taskSel->SetK0L0(gbSpecie);
  taskSel->SetOnline( gbOnline );
  taskSel->SetMass( SFT_MassBins(gbSpecie),
                    SFT_MinMass(gbSpecie),
                    SFT_MaxMass(gbSpecie) );

  taskSel->SetWhichPsi(gbWhichPsi);
  taskSel->SetRFPFilterBit(gbRFPFilterBit);
  taskSel->SetRFPMinPt(gbRFPminPt);
  taskSel->SetRFPMaxPt(gbRFPmaxPt);
  taskSel->SetRFPMinEta(gbRFPminEta);
  taskSel->SetRFPMaxEta(gbRFPmaxEta);
  taskSel->SetRFPTPCSignal(gbRFPTPCsignal);
  taskSel->SetRFPMaxIPxy(gbRFPmaxIPxy);
  taskSel->SetRFPMaxIPz(gbRFPmaxIPz);
  taskSel->SetRFPMinTPCCls(gbRFPTPCncls);

  taskSel->SetDauMinNClsTPC(gbMinNClsTPC);
  taskSel->SetDauMaxChi2PerNClsTPC(gbMaxChi2PerNClsTPC);
  taskSel->SetDauMinXRowsOverNClsFTPC(gbMinXRowsOverNClsFTPC);
  taskSel->SetDauMinEta(gbMinEta);
  taskSel->SetDauMaxEta(gbMaxEta);
  taskSel->SetDauMinPt(gbMinPt);
  taskSel->SetDauMinImpactParameterXY(gbMinImpactParameterXY);
  taskSel->SetDauMaxNSigmaPID(gbMaxNSigmaPID);

  taskSel->SetMaxRapidity(gbMaxRapidity);
  taskSel->SetMaxDCAdaughters(gbMaxDCAdaughters);
  taskSel->SetMinCosinePointingAngleXY(gbMinCosinePointingAngleXY);
  taskSel->SetMinQt(gbMinQt,gbQtPie);
  taskSel->SetMinRadXY(gbMinRadXY);
  taskSel->SetMaxDecayLength(gbMaxDecayLength);
  taskSel->SetMaxProductIPXY(gbMaxProductIPXY);
  taskSel->SetMinEta(gbMinEta);
  taskSel->SetMaxEta(gbMaxEta);
  taskSel->SetMinPt(gbMinPt);
  taskSel->SetUseFlowPackage(gbFlowPackage);

  taskSel->SetQAlevel(gbQA);
  if(gbVZEload.Length()>4) {
    TFile *ocalib = TFile::Open(gbVZEload);
    if(ocalib->IsOpen()) {
      TList *vzero = ocalib->Get("VZECALIB");
      taskSel->LoadVZEResponse(vzero);
    } else {
      printf("ADDTASKFLOWSTRANGE COULD NOT OPEN %s. NO VZE CALIBRATION LOADED!\n",gbVZEload.Data());
    }
  }

  taskSel->SetStoreVZEResponse(gbVZEsave);
  AliAnalysisDataContainer *cOutHist = mgr->CreateContainer(Form("FS_OH_%s",gbStamp.Data()),
							    TList::Class(),
							    AliAnalysisManager::kOutputContainer,
							    Form("%s.root:Selector_%s",fileName.Data(),
								 gbFolder.Data()));
  AliAnalysisDataContainer *exc_TPC = mgr->CreateContainer( Form("FS_TPC_%s",gbStamp.Data()),
                                                            AliFlowEventSimple::Class(),
                                                            AliAnalysisManager::kExchangeContainer );
  AliAnalysisDataContainer *exc_VZE = mgr->CreateContainer( Form("FS_VZE_%s",gbStamp.Data()),
                                                            AliFlowEventSimple::Class(),
                                                            AliAnalysisManager::kExchangeContainer );
  mgr->AddTask(taskSel);
  mgr->ConnectInput (taskSel,0,cinput1);
  mgr->ConnectOutput(taskSel,1,cOutHist);
  mgr->ConnectOutput(taskSel,2,exc_TPC);
  mgr->ConnectOutput(taskSel,3,exc_VZE);

  if(!gbFlowPackage) return;
  if( (!gbQCTPC) && (!gbSPVZE) && (!gbSPTPC) ) return;
  //-------------------FLOWPACKAGE TASKS----------------------------
  AliFlowTrackSimpleCuts *filter[20], *filterhf[20][2]; // MASS BANDS
  for(int mb=0; mb!=SFT_MassBands(gbSpecie); ++mb) {
    filter[mb] = new AliFlowTrackSimpleCuts( Form("Filter_MB%d",mb) );
    filter[mb]->SetEtaMin( -0.8 ); filter[mb]->SetEtaMax( +0.8 );
    filter[mb]->SetMassMin( SFT_MassBandLowEdge(gbSpecie,mb) ); filter[mb]->SetMassMax( SFT_MassBandLowEdge(gbSpecie,mb+1) );
    //half window for POIs
    filterhf[mb][0] = new AliFlowTrackSimpleCuts( Form("Filterhf0_MB%d",mb) );
    filterhf[mb][0]->SetEtaMin( +0.0 ); filterhf[mb][0]->SetEtaMax( +0.8 );
    filterhf[mb][0]->SetMassMin( SFT_MassBandLowEdge(gbSpecie,mb) ); filterhf[mb][0]->SetMassMax( SFT_MassBandLowEdge(gbSpecie,mb+1) );
    filterhf[mb][1] = new AliFlowTrackSimpleCuts( Form("Filterhf1_MB%d",mb) );
    filterhf[mb][1]->SetEtaMin( -0.8 ); filterhf[mb][1]->SetEtaMax( -0.0 );
    filterhf[mb][1]->SetMassMin( SFT_MassBandLowEdge(gbSpecie,mb) ); filterhf[mb][1]->SetMassMax( SFT_MassBandLowEdge(gbSpecie,mb+1) );
    if(gbQCTPC) {
      SFT_AddQCmethod( Form("QCTPCMB%d",mb), exc_TPC, filter[mb]); // QC TPC
    }
    if(gbSPTPC) {
      SFT_AddSPmethod( Form("SPTPCMB%d",mb), exc_TPC, filterhf[mb][0], "Qa" ); // SP TPC Qa
      SFT_AddSPmethod( Form("SPTPCMB%d",mb), exc_TPC, filterhf[mb][1], "Qb" ); // SP TPC Qb
      SFT_AddSPmethod( Form("SPTPC2MB%d",mb), exc_TPC, filterhf[mb][0], "Qa", 0.2 ); // SP TPC Qa
      SFT_AddSPmethod( Form("SPTPC2MB%d",mb), exc_TPC, filterhf[mb][1], "Qb", 0.2 ); // SP TPC Qb
      SFT_AddSPmethod( Form("SPTPC4MB%d",mb), exc_TPC, filterhf[mb][0], "Qa", 0.4 ); // SP TPC Qa
      SFT_AddSPmethod( Form("SPTPC4MB%d",mb), exc_TPC, filterhf[mb][1], "Qb", 0.4 ); // SP TPC Qb
    }
    if(gbSPVZE) {
      SFT_AddSPmethod( Form("SPVZEMB%d",mb), exc_VZE, filter[mb], "Qa" ); // SP VZE Qa
      SFT_AddSPmethod( Form("SPVZEMB%d",mb), exc_VZE, filter[mb], "Qb" ); // SP VZE Qa
    }
  }
}
// ADDING QC METHOD
void SFT_AddQCmethod(char *name, AliAnalysisDataContainer *flowEvent, AliFlowTrackSimpleCuts *cutsPOI=NULL) {
  TString fileName = AliAnalysisManager::GetCommonFileName();
  TString myFolder = Form("%sv%d",gbFolder.Data(),gbHarmonic);
  TString myName = Form("%sv%d_%s",name,gbHarmonic,gbSuffix.Data());
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
  tskQC->SetHarmonic(gbHarmonic);
  tskQC->SetBookOnlyBasicCCH(kTRUE);
  mgr->AddTask(tskQC);
  mgr->ConnectInput( tskQC,0,flowEvent2);
  mgr->ConnectOutput(tskQC,1,outQC);
}
// ADDING SP METHOD
void SFT_AddSPmethod(char *name, AliAnalysisDataContainer *flowEvent, AliFlowTrackSimpleCuts *cutsPOI=NULL, char *Qvector, Double_t gap=0.0) {
  TString fileName = AliAnalysisManager::GetCommonFileName();
  TString myFolder = Form("%sv%d",gbFolder.Data(),gbHarmonic);
  TString myNameSP = Form("%sv%d%s_%s",name,gbHarmonic,Qvector,gbSuffix.Data());
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
  tskSP->SetHarmonic(gbHarmonic);
  tskSP->SetTotalQvector(Qvector);
  //tskSP->SetBookOnlyBasicCCH(kTRUE);
  tskSP->SetBookOnlyBasicCCH(kFALSE);
  mgr->AddTask(tskSP);
  mgr->ConnectInput( tskSP,0,flowEvent2);
  mgr->ConnectOutput(tskSP,1,outSP);
}
// MASSBAND HELPER
double SFT_MassBandLowEdge( int nv0, int mb ) {
  switch(nv0) {
  case(0):
    double lowEdge[14]={0.398, 0.420, 0.444, 0.468, 0.486,
                         0.490, 0.494, 0.498, 0.502, 0.506, 
                         0.524, 0.548, 0.572, 0.598};
    //0.492, 0.502, 0.526, 0.550, 0.574, 0.598};
    break;
  default:
    double lowEdge[10]={1.084, 1.094, 1.104, 1.114, 1.118, 1.128, 1.138, 1.148, 1.158, 1.168};
    break;
  }
  if(nv0>10&&mb==0) return -5;
  if(nv0>10&&mb==1) return +5;
  return lowEdge[mb];
}
// MASSBAND HELPER
int SFT_MassBands( int nv0 ) {
  int bands=1;
  switch(nv0) {
  case(0):
    bands = 13;
    break;
  default:
    bands = 9;
  }
  if(nv0>10) bands=1;
  return bands;
}
// MASSBAND HELPER
int SFT_MassBins( int nv0 ) {
  int bins=100;
  switch(nv0) {
  case(0)://kZERO
    bins=100;
    break;
  default://LAMBDA
    bins=112;
    break;
  }
  if(nv0>10) bins=100;//CHARGED
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
void SFT_PrintConfig() {
  printf("***********************************\n");
  printf("* STRANGE FLOW TASK CONFIGURATION *\n");
  printf("* SUFFIX  %8s                *\n", gbSuffix.Data() );
  printf("* TRIGGER  %8d               *\n", gbTrigger );
  printf("* RUNPP  %3d                      *\n", gbRunPP );
  printf("* RUNPA  %3d                      *\n", gbRunPA );
  printf("* AVOIDEXEC  %3d                  *\n", gbAvoidExec );
  printf("* ESD  %3d                        *\n", gbReadESD );
  printf("* MC  %3d                         *\n", gbReadMC );
  printf("* EXTRAEVENTCUT  %3d              *\n", gbExtraEventCut );
  printf("* CENTMETHOD  %8s                 *\n", gbCentMethod.Data() );
  printf("* CENTPERMIN  %3d                 *\n", gbCentPerMin );
  printf("* CENTPERMAX  %3d                 *\n", gbCentPerMax );
  printf("* SPECIE  %3d                     *\n", gbSpecie );
  printf("* HOMEMADE  %3d                   *\n", gbHomemade );
  printf("* ONLINE  %3d                     *\n", gbOnline );
  printf("* MINNCLSTTPC  %3d                *\n", gbMinNClsTPC );
  printf("* MINXROWS  %3d                   *\n", gbMinXRows );
  printf("* MAXCHI2NCLSTPC  %+9.6f      *\n", gbMaxChi2PerNClsTPC );
  printf("* MINXROWSNFCLSTPC  %+9.6f    *\n", gbMinXRowsOverNClsFTPC );
  printf("* MINETA  %+9.6f              *\n", gbMinEta );
  printf("* MAXETA  %+9.6f              *\n", gbMaxEta );
  printf("* MINPT  %+9.6f               *\n", gbMinPt );
  printf("* MIND0XY  %+9.6f             *\n", gbMinImpactParameterXY );
  printf("* MAXSIGMAPID  %+9.6f         *\n", gbMaxNSigmaPID );
  printf("* MAXY  %+9.6f                *\n", gbMaxRapidity );
  printf("* MAXDCA  %+9.6f              *\n", gbMaxDCAdaughters );
  printf("* MINCTP  %+9.6f              *\n", gbMinCosinePointingAngleXY );
  printf("* MINQT  %+9.6f               *\n", gbMinQt );
  printf("* QTPIE  %+9.6f               *\n", gbQtPie );
  printf("* MINRADXY  %+9.6f            *\n", gbMinRadXY );
  printf("* MAXDL  %+9.6f               *\n", gbMaxDecayLength );
  printf("* D0D0XY  %+9.6f              *\n", gbMaxProductIPXY );
  printf("* DEBUG  %3d                      *\n", gbDebug );
  printf("* QA  %3d                         *\n", gbQA );
  printf("* SKIPSELECTION  %3d              *\n", gbSkipSelection );
  printf("* SKIPFLOW  %3d                   *\n", gbSkipFlow );
  printf("* USEFP  %3d                      *\n", gbFlowPackage );
  printf("* SPVZE  %3d                      *\n", gbSPVZE );
  printf("* SPTPC  %3d                      *\n", gbSPTPC );
  printf("* QCTPC  %3d                      *\n", gbQCTPC );
  printf("* RFFILTERBIT  %3d                *\n", gbRFPFilterBit );
  printf("* RFMINPT  %+9.6f             *\n", gbRFPminPt );
  printf("* RFMAXPT  %+9.6f             *\n", gbRFPmaxPt );
  printf("* RFMINETA  %+9.6f            *\n", gbRFPminEta );
  printf("* RFMAXETA  %+9.6f            *\n", gbRFPmaxEta );
  printf("* RFTPCSIGNAL  %+9.6f         *\n", gbRFPTPCsignal );
  printf("* RFMAXIPXY  %+9.6f           *\n", gbRFPmaxIPxy );
  printf("* RFMAXIPZ  %+9.6f            *\n", gbRFPmaxIPz );
  printf("* RFTPCNCLS  %3d                  *\n", gbRFPTPCncls );
  printf("* VZELOAD  %8s            *\n", gbVZEload.Data() );
  printf("* VZESAVE  %3d                    *\n", gbVZEsave );
  printf("***********************************\n");
}
void SFT_ReadConfig(TString ipf) {
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
      input >> gbSuffix;
    } else if(!varname.CompareTo("TRIGGER")) {
      input >> gbTrigger;
    } else if(!varname.CompareTo("RUNPP")) {
      input >> gbRunPP;
    } else if(!varname.CompareTo("RUNPA")) {
      input >> gbRunPA;
    } else if(!varname.CompareTo("AVOIDEXEC")) {
      input >> gbAvoidExec;
    } else if(!varname.CompareTo("ESD")) {
      input >> gbReadESD;
    } else if(!varname.CompareTo("MC")) {
      input >> gbReadMC;
    } else if(!varname.CompareTo("EXTRAEVENTCUT")) {
      input >> gbExtraEventCut;
    } else if(!varname.CompareTo("CENTMETHOD")) {
      input >> gbCentMethod;
    } else if(!varname.CompareTo("CENTPERMIN")) {
      input >> gbCentPerMin;
    } else if(!varname.CompareTo("CENTPERMAX")) {
      input >> gbCentPerMax;
    } else if(!varname.CompareTo("SPECIE")) {
      input >> gbSpecie;
    } else if(!varname.CompareTo("HOMEMADE")) {
      input >> gbHomemade;
    } else if(!varname.CompareTo("ONLINE")) {
      input >> gbOnline;
    } else if(!varname.CompareTo("MINNCLSTTPC")) {
      input >> gbMinNClsTPC;
    } else if(!varname.CompareTo("MINXROWS")) {
      input >> gbMinXRows;
    } else if(!varname.CompareTo("MAXCHI2NCLSTPC")) {
      input >> gbMaxChi2PerNClsTPC;
    } else if(!varname.CompareTo("MINXROWSNFCLSTPC")) {
      input >> gbMinXRowsOverNClsFTPC;
    } else if(!varname.CompareTo("MINETA")) {
      input >> gbMinEta;
    } else if(!varname.CompareTo("MAXETA")) {
      input >> gbMaxEta;
    } else if(!varname.CompareTo("MINPT")) {
      input >> gbMinPt;
    } else if(!varname.CompareTo("MIND0XY")) {
      input >> gbMinImpactParameterXY;
    } else if(!varname.CompareTo("MAXSIGMAPID")) {
      input >> gbMaxNSigmaPID;
    } else if(!varname.CompareTo("MAXY")) {
      input >> gbMaxRapidity;
    } else if(!varname.CompareTo("MAXDCA")) {
      input >> gbMaxDCAdaughters;
    } else if(!varname.CompareTo("MINCTP")) {
      input >> gbMinCosinePointingAngleXY;
    } else if(!varname.CompareTo("MINQT")) {
      input >> gbMinQt;
    } else if(!varname.CompareTo("QTPIE")) {
      input >> gbQtPie;
    } else if(!varname.CompareTo("MINRADXY")) {
      input >> gbMinRadXY;
    } else if(!varname.CompareTo("MAXDL")) {
      input >> gbMaxDecayLength;
    } else if(!varname.CompareTo("D0D0XY")) {
      input >> gbMaxProductIPXY;
    } else if(!varname.CompareTo("DEBUG")) {
      input >> gbDebug;
    } else if(!varname.CompareTo("QA")) {
      input >> gbQA;
    } else if(!varname.CompareTo("SKIPSELECTION")) {
      input >> gbSkipSelection;
    } else if(!varname.CompareTo("SKIPFLOW")) {
      input >> gbSkipFlow;
    } else if(!varname.CompareTo("USEFP")) {
      input >> gbFlowPackage;
    } else if(!varname.CompareTo("SPVZE")) {
      input >> gbSPVZE;
    } else if(!varname.CompareTo("SPTPC")) {
      input >> gbSPTPC;
    } else if(!varname.CompareTo("QCTPC")) {
      input >> gbQCTPC;
    } else if(!varname.CompareTo("RFFILTERBIT")) {
      input >> gbRFPFilterBit;
    } else if(!varname.CompareTo("RFMINPT")) {
      input >> gbRFPminPt;
    } else if(!varname.CompareTo("RFMAXPT")) {
      input >> gbRFPmaxPt;
    } else if(!varname.CompareTo("RFMINETA")) {
      input >> gbRFPminEta;
    } else if(!varname.CompareTo("RFMAXETA")) {
      input >> gbRFPmaxEta;
    } else if(!varname.CompareTo("RFTPCSIGNAL")) {
      input >> gbRFPTPCsignal;
    } else if(!varname.CompareTo("RFMAXIPXY")) {
      input >> gbRFPmaxIPxy;
    } else if(!varname.CompareTo("RFMAXIPZ")) {
      input >> gbRFPmaxIPz;
    } else if(!varname.CompareTo("RFTPCNCLS")) {
      input >> gbRFPTPCncls;
    } else if(!varname.CompareTo("VZELOAD")) {
      input >> gbVZEload;
    } else if(!varname.CompareTo("VZESAVE")) {
      input >> gbVZEsave;
    } else if(!varname.CompareTo("ALLCC")) {
      input >> gbAllCC;
    } else {
      printf("I dont understand %s\n",varname.Data());
    }
  }
}
