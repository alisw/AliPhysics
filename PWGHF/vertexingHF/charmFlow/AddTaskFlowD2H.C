void AddTaskFlowD2H( int nDmeson, TString fileName, int RPsource,
                     int EtaGapSP=3, double aMax=-0.2, double bMin=+0.2,
                     int minCent=30, int maxCent=50, TString suffixName="C" )
{
  suffixName += Form("%d%d_%s",minCent,maxCent,fileName.Data());
  //-R-P---c-u-t-s--------------------------------------------------------------
  AliFlowTrackCuts* cutsRP;
  double aMin, aMax, bMin, bMax;    // For SP method
  switch (RPsource) {
    case (0):
      cutsRP = (AliFlowTrackCuts*) 
               AliFlowTrackCuts::GetStandardTPCStandaloneTrackCuts2010();
      cutsRP->SetEtaRange( -0.8, +0.8 );
      cutsRP->SetPtRange( 0.2, 5.0 );
      cutsRP->SetAODfilterBit(1);
      cutsRP->SetName( Form("rp_cuts_%s",suffixName.Data()) );
      aMin=-0.8, bMax=0.8;
      break;
    case (1):
      cutsRP = (AliFlowTrackCuts*) 
               AliFlowTrackCuts::GetStandardVZEROOnlyTrackCuts();
      aMin = -3.7; aMax = -2.2 ; bMin = 3.4; bMax = +5.1;
      break;
  }
  fileName.Append(".root");
  //-P-O-I---c-u-t-s------------------------------------------------------------
  double flowBands[5][2];
  switch (nDmeson) { //thanks to ROOT for allowing this
    case (AliRDHFCuts::kD0toKpiCuts):
      AliRDHFCutsD0toKpi *cutsPOI = new AliRDHFCutsD0toKpi( Form("D0toKpi_%s",suffixName.Data()) );
      flowBands[0][0]=1.700; flowBands[0][1]=1.750;
      flowBands[1][0]=1.750; flowBands[1][1]=1.800;
      flowBands[2][0]=1.830; flowBands[2][1]=1.900;
      flowBands[3][0]=1.930; flowBands[3][1]=1.970;
      flowBands[4][0]=1.970; flowBands[4][1]=2.030;
      break;
    case (AliRDHFCuts::kDstarCuts):
      AliRDHFCutsDStartoKpipi *cutsPOI = new AliRDHFCutsDStartoKpipi( Form("DStartoKpipi_%s",suffixName.Data()) );
      flowBands[0][0]=0.140; flowBands[0][1]=0.142;
      flowBands[1][0]=0.144; flowBands[1][1]=0.147;
      flowBands[2][0]=0.150; flowBands[2][1]=0.155;
      flowBands[3][0]=0.155; flowBands[3][1]=0.160;
      flowBands[4][0]=0.160; flowBands[4][1]=0.165;
      break;
  }
  cutsPOI->SetStandardCutsPbPb2010();
  cutsPOI->SetUseCentrality(AliRDHFCuts::kCentV0M);
  cutsPOI->SetMinCentrality(minCent); cutsPOI->SetMaxCentrality(maxCent);
  cutsPOI->SetUseAOD049();

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();

  //-----------------FLOWD2H TASK----------------------------
  AliAnalysisTaskFlowD2H *taskSel = new AliAnalysisTaskFlowD2H(
                                    Form("FlowD2H_%s",suffixName.Data()),
                                    cutsRP, cutsPOI, nDmeson);
  taskSel->SelectCollisionCandidates(AliVEvent::kMB);
  taskSel->SetDebug();
  taskSel->SetFlowEtaRangeAB( aMin, aMax, bMin, bMax );
  TString Qvector;
  switch (EtaGapSP) {
    case 0:
      taskSel->SetPOIEtaRange( -1.5, aMax );
      Qvector = "Qb";
      break;
    case 1:
      taskSel->SetPOIEtaRange( bMin, +1.5 );
      Qvector = "Qa";
      break;
    case 3:
      taskSel->SetPOIEtaRange( -1.5, 1.5 );
      Qvector = "QaQb"; // dummy
      break;
  }
  taskSel->SetFlowPtRange( 2, 16 );
  for(int i=0; i!=5; ++i)
    taskSel->SetFlowBandRange( i, flowBands[i][0], flowBands[i][1] );
  mgr->AddTask(taskSel);
  AliAnalysisDataContainer *coutputCandidatesQA = mgr->CreateContainer(
     Form("QA_%s",suffixName.Data()),TList::Class(),
     AliAnalysisManager::kOutputContainer, fileName+":QACandidates");
  mgr->ConnectOutput(taskSel,1,coutputCandidatesQA);
  mgr->ConnectInput (taskSel,0,cinput1);
  AliAnalysisDataContainer *coutputCandidates[5];
  for(int r=0; r!=5; ++r) {
    coutputCandidates[r] = mgr->CreateContainer(
       Form("Flow_MassBand%d_%s",r,suffixName.Data()),
       AliFlowEventSimple::Class(), AliAnalysisManager::kExchangeContainer);
    mgr->ConnectOutput(taskSel,2+r,coutputCandidates[r]);
  }
  // Scalar Product
  AliAnalysisTaskScalarProduct *taskSP[5];
  AliAnalysisDataContainer *coutputSP[5];
  for(int r=0; r!=5; ++r) {
    taskSP[r] = new AliAnalysisTaskScalarProduct(
             Form("SP_MassBand%d_%s",r,suffixName.Data()),kFALSE);
    taskSP[r]->SetRelDiffMsub(1.0);
    taskSP[r]->SetApplyCorrectionForNUA(kTRUE);
    taskSP[r]->SetTotalQvector( Qvector.Data() );
    mgr->AddTask( taskSP[r] );
    TString outputSP = fileName + Form(":outputSP_MassBand%d",r);
    coutputSP[r] = mgr->CreateContainer(
         Form("cobjSP_MassBand%d_%s",r,suffixName.Data()),TList::Class(),
         AliAnalysisManager::kOutputContainer,outputSP.Data());
    mgr->ConnectInput(taskSP[r],0,coutputCandidates[r]);
    mgr->ConnectOutput(taskSP[r],1,coutputSP[r]);
  }
  if(RPsource!=0)
    return;
  if(EtaGapSP!=3)
    return;
  cout << "HERE!\n";
  // Event Plane
  AliAnalysisTaskEventPlane *taskEP[5];
  AliAnalysisDataContainer *coutputEP[5];
  for(int r=0; r!=5; ++r) {
    taskEP[r] = new AliAnalysisTaskEventPlane(
             Form("EP_MassBand%d_%s",r,suffixName.Data()),kFALSE);
    taskEP[r]->SetApplyCorrectionForNUA(kTRUE);
    taskEP[r]->SetNormalizationType( 0 );
    mgr->AddTask( taskEP[r] );
    TString outputEP = fileName + Form(":outputEP_MassBand%d",r);
    coutputEP[r] = mgr->CreateContainer(
         Form("cobjEP_MassBand%d_%s",r,suffixName.Data()),TList::Class(),
         AliAnalysisManager::kOutputContainer,outputEP.Data());
    mgr->ConnectInput(taskEP[r],0,coutputCandidates[r]);
    mgr->ConnectOutput(taskEP[r],1,coutputEP[r]);
  }
  // Event Plane 2
  AliAnalysisTaskEventPlane *taskEP2[5];
  AliAnalysisDataContainer *coutputEP2[5];
  for(int r=0; r!=5; ++r) {
    taskEP2[r] = new AliAnalysisTaskEventPlane(
             Form("EP2_MassBand%d_%s",r,suffixName.Data()),kFALSE);
    taskEP2[r]->SetApplyCorrectionForNUA(kTRUE);
    taskEP2[r]->SetNormalizationType( 0 );
    mgr->AddTask( taskEP2[r] );
    TString outputEP2 = fileName + Form(":outputEP2_MassBand%d",r);
    coutputEP2[r] = mgr->CreateContainer(
         Form("cobjEP2_MassBand%d_%s",r,suffixName.Data()),TList::Class(),
         AliAnalysisManager::kOutputContainer,outputEP2.Data());
    mgr->ConnectInput(taskEP2[r],0,coutputCandidates[r]);
    mgr->ConnectOutput(taskEP2[r],1,coutputEP2[r]);
  }
  // Q-Cumulants
  AliAnalysisTaskQCumulants *taskQC[5];
  AliAnalysisDataContainer *coutputQC[5];
  for(int r=0; r!=5; ++r) {
    taskQC[r] = new AliAnalysisTaskQCumulants(
             Form("QC_MassBand%d_%s",r,suffixName.Data()),kFALSE);
    taskQC[r]->SetCalculateCumulantsVsM(kFALSE);
    taskQC[r]->SetnBinsMult(10000);
    taskQC[r]->SetMinMult(0.);
    taskQC[r]->SetMaxMult(10000.);
    taskQC[r]->SetApplyCorrectionForNUA(kTRUE);
    taskQC[r]->SetFillMultipleControlHistograms(kFALSE);
    mgr->AddTask( taskQC[r] );
    TString outputQC = fileName + Form(":outputQC_MassBand%d",r);
    coutputQC[r] = mgr->CreateContainer(
         Form("cobjQC_MassBand%d_%s",r,suffixName.Data()),TList::Class(),
         AliAnalysisManager::kOutputContainer,outputQC.Data());
    mgr->ConnectInput(taskQC[r],0,coutputCandidates[r]);
    mgr->ConnectOutput(taskQC[r],1,coutputQC[r]);
  }
}

