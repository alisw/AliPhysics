void AddTaskCMEAnalysis(Bool_t isPbPb = kTRUE,
			Double_t whichData = 2011,
			Bool_t isData = kTRUE,
			Bool_t checkPileup = kTRUE,
			AliFlowEventCuts::refMultMethod gCentralityEstimator = AliFlowEventCuts::kVZERO,
			TString qVector = "Qa",
			Double_t gEtaGap = 0.0,
			Bool_t isVZERO = kFALSE,
			Int_t gAODfilterBit = 768,
			Double_t gDCAvtxXY = -1.,
			Double_t gDCAvtxZ = -1.,
			Double_t gChargeRP = 0,
			Double_t gChargePOI = 0,
			Bool_t gRunSP = kTRUE,
			Bool_t gRunQC = kTRUE,
			Bool_t gRunMHLS = kTRUE,
			Bool_t gRunMHUS = kTRUE,
			Bool_t doQA = kFALSE,
			Bool_t doHigherHarmonic = kFALSE,
			Int_t vnHarmonic=2, TString sLabel = "FB96") {
  //Macro to be used for studies of CME for charged particles
  //The macro uses as an input a configuration macro that
  //creates the AliFlowEventCuts and AliFlowTrackCuts objects
  //Author:Panos.Christakoglou@nikhef.nl
  //=============================================//

  //Define trigger class based on year
  UInt_t triggerSelectionString = AliVEvent::kINT7;
  if(isPbPb) {
    switch(whichData) {
    case 2010:
      triggerSelectionString = AliVEvent::kMB;
      break;
    case 2011:
      triggerSelectionString = AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral;
      break;
    case 2015:
      triggerSelectionString = AliVEvent::kINT7;
      break;
    default:
      break;
    }
  }
  else{
      triggerSelectionString = AliVEvent::kINT7;
  }
  //=============================================//
      
  //=============================================//
  //Define centralities
  const int nCentralities = 10;

  Double_t gCentrality[nCentralities];
  

  gCentrality[0] = 0; gCentrality[1] = 5;

  for(Int_t i = 2; i < 10; i++)
     gCentrality[i] = (i - 1)*10;
    


  //======================================================================//
  // FLOW cut objects config
  AliFlowEventCuts* cutsEvent_wQA[nCentralities];//with QA on
  AliFlowTrackCuts* cutsRP_wQA[nCentralities];   //with QA on
  AliFlowEventCuts* cutsEvent[nCentralities];    //without QA
  AliFlowTrackCuts* cutsRP[nCentralities];       //without QA

  AliFlowTrackCuts* cutsPOI_PP[nCentralities];
  AliFlowTrackCuts* cutsPOI_NN[nCentralities];
  AliFlowTrackCuts* cutsPOI_PN[nCentralities];
  AliFlowTrackCuts* cutsPOI_PNwQA[nCentralities]; //with QA on


  TString suffixNamePP[nCentralities];
  TString suffixNameNN[nCentralities];
  TString suffixNamePN[nCentralities];

  TString outputSlotPP_NameSPv2[nCentralities];
  TString outputSlotPP_NameSPv4[nCentralities];
  TString outputSlotPP_NameQCv2[nCentralities];
  TString outputSlotPP_NameQCv4[nCentralities];
  TString outputSlotPP_NameMHLS[nCentralities];
  TString outputSlotPP_NameMHLS2[nCentralities];
  TString outputSlotPP_NameMHUS[nCentralities];
  TString outputSlotPP_NameMHUS2[nCentralities];

  TString outputSlotNN_NameSPv2[nCentralities];
  TString outputSlotNN_NameSPv4[nCentralities];
  TString outputSlotNN_NameQCv2[nCentralities];
  TString outputSlotNN_NameQCv4[nCentralities];
  TString outputSlotNN_NameMHLS[nCentralities];
  TString outputSlotNN_NameMHLS2[nCentralities];
  TString outputSlotNN_NameMHUS[nCentralities];
  TString outputSlotNN_NameMHUS2[nCentralities];

  TString outputSlotPN_NameSPv2[nCentralities];
  TString outputSlotPN_NameSPv4[nCentralities];
  TString outputSlotPN_NameQCv2[nCentralities];
  TString outputSlotPN_NameQCv4[nCentralities];
  TString outputSlotPN_NameMHLS[nCentralities];
  TString outputSlotPN_NameMHLS2[nCentralities];
  TString outputSlotPN_NameMHUS[nCentralities];
  TString outputSlotPN_NameMHUS2[nCentralities];


  
  for(Int_t iCentralityBin = 0; iCentralityBin < nCentralities - 1; iCentralityBin++) {
    //Create the event cut object
    cutsEvent_wQA[iCentralityBin] = createFlowEventCutObject(gCentrality[iCentralityBin],gCentrality[iCentralityBin+1],isPbPb,whichData,gCentralityEstimator,doQA,checkPileup);
    cutsEvent[iCentralityBin] = createFlowEventCutObject(gCentrality[iCentralityBin],gCentrality[iCentralityBin+1],isPbPb,whichData,gCentralityEstimator,kFALSE,checkPileup);
        
    //Create the RP cut object
    cutsRP_wQA[iCentralityBin] = createFlowRPCutObject(gCentrality[iCentralityBin],gCentrality[iCentralityBin+1],isVZERO,gAODfilterBit,whichData,gChargeRP,doQA);
    cutsRP[iCentralityBin] = createFlowRPCutObject(gCentrality[iCentralityBin],gCentrality[iCentralityBin+1],isVZERO,gAODfilterBit,whichData,gChargeRP,kFALSE);
    
    
    //Create the POI cut object for ++,-- and +-
    cutsPOI_PP[iCentralityBin] = createFlowPOICutObject(gCentrality[iCentralityBin],gCentrality[iCentralityBin+1],qVector,gEtaGap,isVZERO,isPbPb,gAODfilterBit,whichData,gDCAvtxXY,gDCAvtxZ,1,kFALSE);
    cutsPOI_NN[iCentralityBin] = createFlowPOICutObject(gCentrality[iCentralityBin],gCentrality[iCentralityBin+1],qVector,gEtaGap,isVZERO,isPbPb,gAODfilterBit,whichData,gDCAvtxXY,gDCAvtxZ,-1,kFALSE);
    cutsPOI_PN[iCentralityBin] = createFlowPOICutObject(gCentrality[iCentralityBin],gCentrality[iCentralityBin+1],qVector,gEtaGap,isVZERO,isPbPb,gAODfilterBit,whichData,gDCAvtxXY,gDCAvtxZ,0,kFALSE);
    cutsPOI_PNwQA[iCentralityBin] = createFlowPOICutObject(gCentrality[iCentralityBin],gCentrality[iCentralityBin+1],qVector,gEtaGap,isVZERO,isPbPb,gAODfilterBit,whichData,gDCAvtxXY,gDCAvtxZ,0,doQA);

        
    //change container names, Add FB in extension. So that wagons with different FB can run in same train. 
    suffixNamePP[iCentralityBin] += "Centrality";
    suffixNamePP[iCentralityBin] += gCentrality[iCentralityBin];
    suffixNamePP[iCentralityBin] += "To";
    suffixNamePP[iCentralityBin] += gCentrality[iCentralityBin+1];
    suffixNamePP[iCentralityBin] += sLabel;
  
    suffixNameNN[iCentralityBin] += "Centrality";
    suffixNameNN[iCentralityBin] += gCentrality[iCentralityBin];
    suffixNameNN[iCentralityBin] += "To";
    suffixNameNN[iCentralityBin] += gCentrality[iCentralityBin+1];
    suffixNameNN[iCentralityBin] += sLabel;

    suffixNamePN[iCentralityBin] += "Centrality";
    suffixNamePN[iCentralityBin] += gCentrality[iCentralityBin];
    suffixNamePN[iCentralityBin] += "To";
    suffixNamePN[iCentralityBin] += gCentrality[iCentralityBin+1];
    suffixNamePN[iCentralityBin] += sLabel;

    /*
      suffixNamePP[iCentralityBin] += "Centrality";
      suffixNamePP[iCentralityBin] += gCentrality[iCentralityBin];
      suffixNamePP[iCentralityBin] += "To";
      suffixNamePP[iCentralityBin] += gCentrality[iCentralityBin+1];

      suffixNameNN[iCentralityBin] += "Centrality";
      suffixNameNN[iCentralityBin] += gCentrality[iCentralityBin];
      suffixNameNN[iCentralityBin] += "To";
      suffixNameNN[iCentralityBin] += gCentrality[iCentralityBin+1];

      suffixNamePN[iCentralityBin] += "Centrality";
      suffixNamePN[iCentralityBin] += gCentrality[iCentralityBin];
      suffixNamePN[iCentralityBin] += "To";
      suffixNamePN[iCentralityBin] += gCentrality[iCentralityBin+1];
    */

    /*if(gChargePOI == 1) 
	suffixName[iCentralityBin] += "PlusPlus";
      else if(gChargePOI == -1) 
	suffixName[iCentralityBin] += "MinusMinus";*/


    //Giving names to output slots
    suffixNamePP[iCentralityBin] += "PlusPlus";    
    outputSlotPP_NameSPv2[iCentralityBin] = Form("%s_SP_v2",suffixNamePP[iCentralityBin].Data());
    outputSlotPP_NameSPv4[iCentralityBin] = Form("%s_SP_v4",suffixNamePP[iCentralityBin].Data());
    outputSlotPP_NameQCv2[iCentralityBin] = Form("%s_QC_v2",suffixNamePP[iCentralityBin].Data());
    outputSlotPP_NameQCv4[iCentralityBin] = Form("%s_QC_v4",suffixNamePP[iCentralityBin].Data());
    outputSlotPP_NameMHLS[iCentralityBin] = Form("%s_MHLS_v1",suffixNamePP[iCentralityBin].Data());
    outputSlotPP_NameMHLS2[iCentralityBin] = Form("%s_MHLS_v2",suffixNamePP[iCentralityBin].Data());
    outputSlotPP_NameMHUS[iCentralityBin]  = Form("%s_MHUS_v1",suffixNamePP[iCentralityBin].Data());
    outputSlotPP_NameMHUS2[iCentralityBin] = Form("%s_MHUS_v2",suffixNamePP[iCentralityBin].Data());

    suffixNameNN[iCentralityBin] += "MinusMinus";    
    outputSlotNN_NameSPv2[iCentralityBin] = Form("%s_SP_v2",suffixNameNN[iCentralityBin].Data());
    outputSlotNN_NameSPv4[iCentralityBin] = Form("%s_SP_v4",suffixNameNN[iCentralityBin].Data());
    outputSlotNN_NameQCv2[iCentralityBin] = Form("%s_QC_v2",suffixNameNN[iCentralityBin].Data());
    outputSlotNN_NameQCv4[iCentralityBin] = Form("%s_QC_v4",suffixNameNN[iCentralityBin].Data());
    outputSlotNN_NameMHLS[iCentralityBin] = Form("%s_MHLS_v1",suffixNameNN[iCentralityBin].Data());
    outputSlotNN_NameMHLS2[iCentralityBin] = Form("%s_MHLS_v2",suffixNameNN[iCentralityBin].Data());
    outputSlotNN_NameMHUS[iCentralityBin]  = Form("%s_MHUS_v1",suffixNameNN[iCentralityBin].Data());
    outputSlotNN_NameMHUS2[iCentralityBin] = Form("%s_MHUS_v2",suffixNameNN[iCentralityBin].Data());

    suffixNamePN[iCentralityBin] += "PlusMinus";    
    outputSlotPN_NameSPv2[iCentralityBin] = Form("%s_SP_v2",suffixNamePN[iCentralityBin].Data());
    outputSlotPN_NameSPv4[iCentralityBin] = Form("%s_SP_v4",suffixNamePN[iCentralityBin].Data());
    outputSlotPN_NameQCv2[iCentralityBin] = Form("%s_QC_v2",suffixNamePN[iCentralityBin].Data());
    outputSlotPN_NameQCv4[iCentralityBin] = Form("%s_QC_v4",suffixNamePN[iCentralityBin].Data());
    outputSlotPN_NameMHLS[iCentralityBin] = Form("%s_MHLS_v1",suffixNamePN[iCentralityBin].Data());
    outputSlotPN_NameMHLS2[iCentralityBin] = Form("%s_MHLS_v2",suffixNamePN[iCentralityBin].Data());
    outputSlotPN_NameMHUS[iCentralityBin]  = Form("%s_MHUS_v1",suffixNamePN[iCentralityBin].Data());
    outputSlotPN_NameMHUS2[iCentralityBin] = Form("%s_MHUS_v2",suffixNamePN[iCentralityBin].Data());


  }//loop over centrality classes

















  TString fileName = "AnalysisResults";
  fileName.Append(".root");


  //========================FLOWPACKAGE TASKS=========================//
  AliAnalysisDataContainer *cinput1_PP[nCentralities];
  AliAnalysisDataContainer *coutputFE_PP[nCentralities];
  AliAnalysisDataContainer *coutputFE_PPQA[nCentralities];
  AliAnalysisTaskFlowEvent *taskFE_PP[nCentralities];
  AliAnalysisDataContainer *flowEvent_PP[nCentralities];
  AliAnalysisTaskFilterFE  *taskFilter_PP[nCentralities];


  AliAnalysisDataContainer *coutputSPv2_PP[nCentralities];
  AliAnalysisTaskScalarProduct *taskSPv2_PP[nCentralities];
  AliAnalysisDataContainer *coutputQCv2_PP[nCentralities];
  AliAnalysisTaskQCumulants *taskQCv2_PP[nCentralities];
  AliAnalysisDataContainer *coutputMHLS_PP[nCentralities];
  AliAnalysisTaskMixedHarmonics *taskMHLS_PP[nCentralities];

  AliAnalysisDataContainer *coutputSPv4_PP[nCentralities];
  AliAnalysisTaskScalarProduct *taskSPv4_PP[nCentralities];
  AliAnalysisDataContainer *coutputQCv4_PP[nCentralities];
  AliAnalysisTaskQCumulants *taskQCv4_PP[nCentralities];
  AliAnalysisDataContainer *coutputMHLS_PP2[nCentralities];
  AliAnalysisTaskMixedHarmonics *taskMHLS_PP2[nCentralities];
  

  //----NN----
  AliAnalysisDataContainer *cinput1_NN[nCentralities];
  AliAnalysisDataContainer *coutputFE_NN[nCentralities];
  AliAnalysisDataContainer *coutputFE_NNQA[nCentralities];
  AliAnalysisTaskFlowEvent *taskFE_NN[nCentralities];
  AliAnalysisDataContainer *flowEvent_NN[nCentralities];
  AliAnalysisTaskFilterFE *taskFilter_NN[nCentralities];

  AliAnalysisDataContainer *coutputSPv2_NN[nCentralities];
  AliAnalysisTaskScalarProduct *taskSPv2_NN[nCentralities];
  AliAnalysisDataContainer *coutputQCv2_NN[nCentralities];
  AliAnalysisTaskQCumulants *taskQCv2_NN[nCentralities];
  AliAnalysisDataContainer *coutputMHLS_NN[nCentralities];
  AliAnalysisTaskMixedHarmonics *taskMHLS_NN[nCentralities];

  AliAnalysisDataContainer *coutputSPv4_NN[nCentralities];
  AliAnalysisTaskScalarProduct *taskSPv4_NN[nCentralities];
  AliAnalysisDataContainer *coutputQCv4_NN[nCentralities];
  AliAnalysisTaskQCumulants *taskQCv4_NN[nCentralities];
  AliAnalysisDataContainer *coutputMHLS_NN2[nCentralities];
  AliAnalysisTaskMixedHarmonics *taskMHLS_NN2[nCentralities];


 //----PN----
  AliAnalysisDataContainer *cinput1_PN[nCentralities];
  AliAnalysisDataContainer *coutputFE_PN[nCentralities];
  AliAnalysisDataContainer *coutputFE_PNQA[nCentralities];
  AliAnalysisTaskFlowEvent *taskFE_PN[nCentralities];
  AliAnalysisDataContainer *flowEvent_PN[nCentralities];
  AliAnalysisTaskFilterFE *taskFilter_PN[nCentralities];

  AliAnalysisDataContainer *coutputSPv2_PN[nCentralities];
  AliAnalysisTaskScalarProduct *taskSPv2_PN[nCentralities];
  AliAnalysisDataContainer *coutputQCv2_PN[nCentralities];
  AliAnalysisTaskQCumulants *taskQCv2_PN[nCentralities];
  AliAnalysisDataContainer *coutputMHLS_PN[nCentralities];
  AliAnalysisTaskMixedHarmonics *taskMHLS_PN[nCentralities];
  //unlike sign:
  AliAnalysisDataContainer *coutputMHUS[nCentralities];
  AliAnalysisTaskMixedHarmonics *taskMHUS[nCentralities];
  
  AliAnalysisDataContainer *coutputSPv4_PN[nCentralities];
  AliAnalysisTaskScalarProduct *taskSPv4_PN[nCentralities];
  AliAnalysisDataContainer *coutputQCv4_PN[nCentralities];
  AliAnalysisTaskQCumulants *taskQCv4_PN[nCentralities];
  AliAnalysisDataContainer *coutputMHLS_PN2[nCentralities];
  AliAnalysisTaskMixedHarmonics *taskMHLS_PN2[nCentralities];
  //unlike sign:
  AliAnalysisDataContainer *coutputMHUS2[nCentralities];
  AliAnalysisTaskMixedHarmonics *taskMHUS2[nCentralities];

  TString outputQA_PP[nCentralities];
  TString outputQA_NN[nCentralities];
  TString outputQA_PN[nCentralities];

  TString myNameSPv2[nCentralities];
  TString myNameSPv4[nCentralities];
  TString myNameQCv2[nCentralities];
  TString myNameQCv4[nCentralities];
  TString myNameMHLS[nCentralities];
  TString myNameMHLS2[nCentralities];
  TString myNameMHUS[nCentralities];
  TString myNameMHUS2[nCentralities];
  TString        slot[nCentralities];

  
  for(Int_t iCentralityBin = 0; iCentralityBin < nCentralities - 1; iCentralityBin++) {

    //======================== Initiate v2 Tasks ==========================

    //------- PP --------------
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
      Error("AddTaskFlowEvent", "No analysis manager to connect to.");
      return NULL;
    }
    
    //Check the analysis type using the event handlers connected to the analysis
    // manager. The availability of MC handler can also be checked here.
    if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskFlowEvent", "This task requires an input event handler");
      return NULL;
    }
    
    taskFE_PP[iCentralityBin] = new AliAnalysisTaskFlowEvent(Form("TaskFlowEvent_%s",suffixNamePP[iCentralityBin].Data()),"",kFALSE);
    taskFE_PP[iCentralityBin]->SelectCollisionCandidates(triggerSelectionString);
    
    //Sub events
    Double_t minA = -0.8;//
    Double_t maxA = -0.5*gEtaGap;//
    Double_t minB = +0.5*gEtaGap;//
    Double_t maxB = 0.8;//
    
    if(!isVZERO)
      taskFE_PP[iCentralityBin]->SetSubeventEtaRange(minA, maxA, minB, maxB);
    if(isVZERO)
      taskFE_PP[iCentralityBin]->SetSubeventEtaRange(-5,-1.5,+1.5,5);
    mgr->AddTask(taskFE_PP[iCentralityBin]);
    
    // Pass cuts for RPs and POIs to the task:
    taskFE_PP[iCentralityBin]->SetCutsEvent(cutsEvent[iCentralityBin]);
    taskFE_PP[iCentralityBin]->SetCutsRP(cutsRP[iCentralityBin]);

    taskFE_PP[iCentralityBin]->SetCutsPOI(cutsPOI_PP[iCentralityBin]);

    if (cutsRP[iCentralityBin]->GetParamType()==AliFlowTrackCuts::kVZERO) {
      taskFE_PP[iCentralityBin]->SetHistWeightvsPhiMin(0.);
      taskFE_PP[iCentralityBin]->SetHistWeightvsPhiMax(200.);
    }
    
    if(!mgr){
      AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    }
    cinput1_PP[iCentralityBin] = mgr->GetCommonInputContainer();
    
    coutputFE_PP[iCentralityBin] = mgr->CreateContainer(Form("FlowEvent_%s",suffixNamePP[iCentralityBin].Data()),AliFlowEventSimple::Class(),AliAnalysisManager::kExchangeContainer);
    
    mgr->ConnectInput(taskFE_PP[iCentralityBin],0,cinput1_PP[iCentralityBin]);
    mgr->ConnectOutput(taskFE_PP[iCentralityBin],1,coutputFE_PP[iCentralityBin]);
    //==========================================================
    
    //======================================================================//
    flowEvent_PP[iCentralityBin] = mgr->CreateContainer(Form("Filter_%s",outputSlotPP_NameSPv2[iCentralityBin].Data()),AliFlowEventSimple::Class(),AliAnalysisManager::kExchangeContainer);
    taskFilter_PP[iCentralityBin] = new AliAnalysisTaskFilterFE(Form("TaskFilter_%s",outputSlotPP_NameSPv2[iCentralityBin].Data()),cutsRP[iCentralityBin], NULL);
    
    if(!isVZERO) 
      taskFilter_PP[iCentralityBin]->SetSubeventEtaRange(minA, maxA, minB, maxB);
    if(isVZERO) 
      taskFilter_PP[iCentralityBin]->SetSubeventEtaRange(-5,-1.5,+1.5,5);
    mgr->AddTask(taskFilter_PP[iCentralityBin]);
    mgr->ConnectInput(taskFilter_PP[iCentralityBin],0,coutputFE_PP[iCentralityBin]);
    mgr->ConnectOutput(taskFilter_PP[iCentralityBin],1,flowEvent_PP[iCentralityBin]);
    //======================================================================//
  
    //======================================================================//
    if(gRunSP) {
      TString outputSPv2_PP = fileName;
      outputSPv2_PP += ":outputSPv2analysis";
      outputSPv2_PP += "PlusPlus";

      taskSPv2_PP[iCentralityBin] = new AliAnalysisTaskScalarProduct(Form("TaskScalarProduct_%s",outputSlotPP_NameSPv2[iCentralityBin].Data()),kFALSE);
      taskSPv2_PP[iCentralityBin]->SetHarmonic(vnHarmonic);
      taskSPv2_PP[iCentralityBin]->SelectCollisionCandidates(triggerSelectionString);
      taskSPv2_PP[iCentralityBin]->SetRelDiffMsub(1.0);
      taskSPv2_PP[iCentralityBin]->SetTotalQvector(qVector);
      taskSPv2_PP[iCentralityBin]->SetApplyCorrectionForNUA(kTRUE);
      
      coutputSPv2_PP[iCentralityBin] = mgr->CreateContainer(Form("%s",outputSlotPP_NameSPv2[iCentralityBin].Data()),TList::Class(),AliAnalysisManager::kOutputContainer,outputSPv2_PP);
      mgr->AddTask(taskSPv2_PP[iCentralityBin]);
      mgr->ConnectInput(taskSPv2_PP[iCentralityBin],0,flowEvent_PP[iCentralityBin]);
      mgr->ConnectInput(taskSPv2_PP[iCentralityBin],0,coutputFE_PP[iCentralityBin]);
      mgr->ConnectOutput(taskSPv2_PP[iCentralityBin],1,coutputSPv2_PP[iCentralityBin]);
      
	/*
        else if(doHigherHarmonic){
	TString outputSPv4_PP = fileName;
	outputSPv4_PP += ":outputSPv4analysis";
	outputSPv4_PP += "PlusPlus";  
    
	taskSPv4_PP[iCentralityBin] = new AliAnalysisTaskScalarProduct(Form("TaskScalarProduct_%s",outputSlotPP_NameSPv4[iCentralityBin].Data()),kFALSE);
	taskSPv4_PP[iCentralityBin]->SetHarmonic(4);
	taskSPv4_PP[iCentralityBin]->SelectCollisionCandidates(triggerSelectionString);
	taskSPv4_PP[iCentralityBin]->SetRelDiffMsub(1.0);
	taskSPv4_PP[iCentralityBin]->SetTotalQvector(qVector);
	taskSPv4_PP[iCentralityBin]->SetApplyCorrectionForNUA(kTRUE);
      
	coutputSPv4_PP[iCentralityBin] = mgr->CreateContainer(Form("%s",outputSlotPP_NameSPv4[iCentralityBin].Data()),TList::Class(),AliAnalysisManager::kOutputContainer,outputSPv4_PP);
	mgr->AddTask(taskSPv4_PP[iCentralityBin]);
	mgr->ConnectInput(taskSPv4_PP[iCentralityBin],0,flowEvent_PP[iCentralityBin]);
	mgr->ConnectInput(taskSPv4_PP[iCentralityBin],0,coutputFE_PP[iCentralityBin]);
	mgr->ConnectOutput(taskSPv4_PP[iCentralityBin],1,coutputSPv4_PP[iCentralityBin]);
      }*/
    }
    //======================================================================//
   

    
    //======================================================================//
    if(gRunQC) {

      TString outputQCv2_PP = fileName;
      outputQCv2_PP += ":outputQCv2analysis";
      outputQCv2_PP += "PlusPlus";
      
      taskQCv2_PP[iCentralityBin] = new AliAnalysisTaskQCumulants(Form("TaskQCumulant_%s",outputSlotPP_NameQCv2[iCentralityBin].Data()),kFALSE);
      taskQCv2_PP[iCentralityBin]->SetHarmonic(vnHarmonic);
      taskQCv2_PP[iCentralityBin]->SelectCollisionCandidates(triggerSelectionString);
      taskQCv2_PP[iCentralityBin]->SetUsePhiWeights(kFALSE); 
      taskQCv2_PP[iCentralityBin]->SetUsePtWeights(kFALSE);
      taskQCv2_PP[iCentralityBin]->SetUseEtaWeights(kFALSE); 
      taskQCv2_PP[iCentralityBin]->SetCalculateCumulantsVsM(kFALSE);
      taskQCv2_PP[iCentralityBin]->SetnBinsMult(10000);
      taskQCv2_PP[iCentralityBin]->SetMinMult(0.);
      taskQCv2_PP[iCentralityBin]->SetMaxMult(10000.);
      taskQCv2_PP[iCentralityBin]->SetApplyCorrectionForNUA(kTRUE);
      taskQCv2_PP[iCentralityBin]->SetFillMultipleControlHistograms(kFALSE);     
      
      coutputQCv2_PP[iCentralityBin] = mgr->CreateContainer(Form("%s",outputSlotPP_NameQCv2[iCentralityBin].Data()),TList::Class(),AliAnalysisManager::kOutputContainer,outputQCv2_PP);
      mgr->AddTask(taskQCv2_PP[iCentralityBin]);
      mgr->ConnectInput(taskQCv2_PP[iCentralityBin],0,flowEvent_PP[iCentralityBin]);
      mgr->ConnectInput(taskQCv2_PP[iCentralityBin],0,coutputFE_PP[iCentralityBin]);
      mgr->ConnectOutput(taskQCv2_PP[iCentralityBin],1,coutputQCv2_PP[iCentralityBin]);
      
 
      /*else if(doHigherHarmonic){
	TString outputQCv4_PP = fileName;
	outputQCv4_PP += ":outputQCv4analysis";
	outputQCv4_PP += "PlusPlus";

	taskQCv4_PP[iCentralityBin] = new AliAnalysisTaskQCumulants(Form("TaskQCumulant_%s",outputSlotPP_NameQCv4[iCentralityBin].Data()),kFALSE);
	taskQCv4_PP[iCentralityBin]->SetHarmonic(4);
	taskQCv4_PP[iCentralityBin]->SelectCollisionCandidates(triggerSelectionString);
	taskQCv4_PP[iCentralityBin]->SetUsePhiWeights(kFALSE); 
	taskQCv4_PP[iCentralityBin]->SetUsePtWeights(kFALSE);
	taskQCv4_PP[iCentralityBin]->SetUseEtaWeights(kFALSE); 
	taskQCv4_PP[iCentralityBin]->SetCalculateCumulantsVsM(kFALSE);
	taskQCv4_PP[iCentralityBin]->SetnBinsMult(10000);
	taskQCv4_PP[iCentralityBin]->SetMinMult(0.);
	taskQCv4_PP[iCentralityBin]->SetMaxMult(10000.);
	taskQCv4_PP[iCentralityBin]->SetApplyCorrectionForNUA(kTRUE);
	taskQCv4_PP[iCentralityBin]->SetFillMultipleControlHistograms(kFALSE);     
      
	coutputQCv4_PP[iCentralityBin] = mgr->CreateContainer(Form("%s",outputSlotPP_NameQCv4[iCentralityBin].Data()),TList::Class(),AliAnalysisManager::kOutputContainer,outputQCv4_PP);
	mgr->AddTask(taskQCv4_PP[iCentralityBin]);
	mgr->ConnectInput(taskQCv4_PP[iCentralityBin],0,flowEvent_PP[iCentralityBin]);
	mgr->ConnectInput(taskQCv4_PP[iCentralityBin],0,coutputFE_PP[iCentralityBin]);
	mgr->ConnectOutput(taskQCv4_PP[iCentralityBin],1,coutputQCv4_PP[iCentralityBin]);
	}*/
    }



   //-------------- NN --------------
    if(!mgr){
      AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    }

    
    //Check the analysis type using the event handlers connected to the analysis
    // manager. The availability of MC handler can also be checked here.
    if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskFlowEvent", "This task requires an input event handler");
      return NULL;
    }
    
    taskFE_NN[iCentralityBin] = new AliAnalysisTaskFlowEvent(Form("TaskFlowEvent_%s",suffixNameNN[iCentralityBin].Data()),"",kFALSE);
    taskFE_NN[iCentralityBin]->SelectCollisionCandidates(triggerSelectionString);
    
    //Sub events
    Double_t minA = -0.8;//
    Double_t maxA = -0.5*gEtaGap;//
    Double_t minB = +0.5*gEtaGap;//
    Double_t maxB = 0.8;//
    
    if(!isVZERO)
      taskFE_NN[iCentralityBin]->SetSubeventEtaRange(minA, maxA, minB, maxB);
    if(isVZERO)
      taskFE_NN[iCentralityBin]->SetSubeventEtaRange(-5,-1.5,+1.5,5);
    mgr->AddTask(taskFE_NN[iCentralityBin]);
    
    // Pass cuts for RPs and POIs to the task:
    taskFE_NN[iCentralityBin]->SetCutsEvent(cutsEvent[iCentralityBin]);
    taskFE_NN[iCentralityBin]->SetCutsRP(cutsRP[iCentralityBin]);

    taskFE_NN[iCentralityBin]->SetCutsPOI(cutsPOI_NN[iCentralityBin]);

    if (cutsRP[iCentralityBin]->GetParamType()==AliFlowTrackCuts::kVZERO) {
      taskFE_NN[iCentralityBin]->SetHistWeightvsPhiMin(0.);
      taskFE_NN[iCentralityBin]->SetHistWeightvsPhiMax(200.);
    }

    if(!mgr){
      AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    }
    cinput1_NN[iCentralityBin] = mgr->GetCommonInputContainer();
    
    coutputFE_NN[iCentralityBin] = mgr->CreateContainer(Form("FlowEvent_%s",suffixNameNN[iCentralityBin].Data()),AliFlowEventSimple::Class(),AliAnalysisManager::kExchangeContainer);
    
    mgr->ConnectInput(taskFE_NN[iCentralityBin],0,cinput1_NN[iCentralityBin]);
    mgr->ConnectOutput(taskFE_NN[iCentralityBin],1,coutputFE_NN[iCentralityBin]);
    //==========================================================
    
    //======================================================================//
    flowEvent_NN[iCentralityBin] = mgr->CreateContainer(Form("Filter_%s",outputSlotNN_NameSPv2[iCentralityBin].Data()),AliFlowEventSimple::Class(),AliAnalysisManager::kExchangeContainer);
    taskFilter_NN[iCentralityBin] = new AliAnalysisTaskFilterFE(Form("TaskFilter_%s",outputSlotNN_NameSPv2[iCentralityBin].Data()),cutsRP[iCentralityBin], NULL);
    
    if(!isVZERO) 
      taskFilter_NN[iCentralityBin]->SetSubeventEtaRange(minA, maxA, minB, maxB);
    if(isVZERO) 
      taskFilter_NN[iCentralityBin]->SetSubeventEtaRange(-5,-1.5,+1.5,5);
    mgr->AddTask(taskFilter_NN[iCentralityBin]);
    mgr->ConnectInput(taskFilter_NN[iCentralityBin],0,coutputFE_NN[iCentralityBin]);
    mgr->ConnectOutput(taskFilter_NN[iCentralityBin],1,flowEvent_NN[iCentralityBin]);
    //======================================================================//
      
    //======================================================================//
    if(gRunSP) {

      TString outputSPv2_NN = fileName;
      outputSPv2_NN += ":outputSPv2analysis";
      outputSPv2_NN += "MinusMinus";

      taskSPv2_NN[iCentralityBin] = new AliAnalysisTaskScalarProduct(Form("TaskScalarProduct_%s",outputSlotNN_NameSPv2[iCentralityBin].Data()),kFALSE);
      taskSPv2_NN[iCentralityBin]->SetHarmonic(vnHarmonic);
      taskSPv2_NN[iCentralityBin]->SelectCollisionCandidates(triggerSelectionString);
      taskSPv2_NN[iCentralityBin]->SetRelDiffMsub(1.0);
      taskSPv2_NN[iCentralityBin]->SetTotalQvector(qVector);
      taskSPv2_NN[iCentralityBin]->SetApplyCorrectionForNUA(kTRUE);
      
      coutputSPv2_NN[iCentralityBin] = mgr->CreateContainer(Form("%s",outputSlotNN_NameSPv2[iCentralityBin].Data()),TList::Class(),AliAnalysisManager::kOutputContainer,outputSPv2_NN);
      mgr->AddTask(taskSPv2_NN[iCentralityBin]);
      mgr->ConnectInput(taskSPv2_NN[iCentralityBin],0,flowEvent_NN[iCentralityBin]);
      mgr->ConnectInput(taskSPv2_NN[iCentralityBin],0,coutputFE_NN[iCentralityBin]);
      mgr->ConnectOutput(taskSPv2_NN[iCentralityBin],1,coutputSPv2_NN[iCentralityBin]);
      
      /*else if(doHigherHarmonic){
	TString outputSPv4_NN = fileName;
	outputSPv4_NN += ":outputSPv4analysis";
	outputSPv4_NN += "MinusMinus";
      
	taskSPv4_NN[iCentralityBin] = new AliAnalysisTaskScalarProduct(Form("TaskScalarProduct_%s",outputSlotNN_NameSPv4[iCentralityBin].Data()),kFALSE);
	taskSPv4_NN[iCentralityBin]->SetHarmonic(4);
	taskSPv4_NN[iCentralityBin]->SelectCollisionCandidates(triggerSelectionString);
	taskSPv4_NN[iCentralityBin]->SetRelDiffMsub(1.0);
	taskSPv4_NN[iCentralityBin]->SetTotalQvector(qVector);
	taskSPv4_NN[iCentralityBin]->SetApplyCorrectionForNUA(kTRUE);
      
	coutputSPv4_NN[iCentralityBin] = mgr->CreateContainer(Form("%s",outputSlotNN_NameSPv4[iCentralityBin].Data()),TList::Class(),AliAnalysisManager::kOutputContainer,outputSPv4_NN);
	mgr->AddTask(taskSPv4_NN[iCentralityBin]);
	mgr->ConnectInput(taskSPv4_NN[iCentralityBin],0,flowEvent_NN[iCentralityBin]);
	mgr->ConnectInput(taskSPv4_NN[iCentralityBin],0,coutputFE_NN[iCentralityBin]);
	mgr->ConnectOutput(taskSPv4_NN[iCentralityBin],1,coutputSPv4_NN[iCentralityBin]);
	}*/
    }
    //======================================================================//

    //======================================================================//
    if(gRunQC) {

      TString outputQCv2_NN = fileName;
      outputQCv2_NN += ":outputQCv2analysis";
      outputQCv2_NN += "MinusMinus";
      
      taskQCv2_NN[iCentralityBin] = new AliAnalysisTaskQCumulants(Form("TaskQCumulant_%s",outputSlotNN_NameQCv2[iCentralityBin].Data()),kFALSE);
      taskQCv2_NN[iCentralityBin]->SetHarmonic(vnHarmonic);
      taskQCv2_NN[iCentralityBin]->SelectCollisionCandidates(triggerSelectionString);
      taskQCv2_NN[iCentralityBin]->SetUsePhiWeights(kFALSE); 
      taskQCv2_NN[iCentralityBin]->SetUsePtWeights(kFALSE);
      taskQCv2_NN[iCentralityBin]->SetUseEtaWeights(kFALSE); 
      taskQCv2_NN[iCentralityBin]->SetCalculateCumulantsVsM(kFALSE);
      taskQCv2_NN[iCentralityBin]->SetnBinsMult(10000);
      taskQCv2_NN[iCentralityBin]->SetMinMult(0.);
      taskQCv2_NN[iCentralityBin]->SetMaxMult(10000.);
      taskQCv2_NN[iCentralityBin]->SetApplyCorrectionForNUA(kTRUE);
      taskQCv2_NN[iCentralityBin]->SetFillMultipleControlHistograms(kFALSE);     
      
      coutputQCv2_NN[iCentralityBin] = mgr->CreateContainer(Form("%s",outputSlotNN_NameQCv2[iCentralityBin].Data()),TList::Class(),AliAnalysisManager::kOutputContainer,outputQCv2_NN);
      mgr->AddTask(taskQCv2_NN[iCentralityBin]);
      mgr->ConnectInput(taskQCv2_NN[iCentralityBin],0,flowEvent_NN[iCentralityBin]);
      mgr->ConnectInput(taskQCv2_NN[iCentralityBin],0,coutputFE_NN[iCentralityBin]);
      mgr->ConnectOutput(taskQCv2_NN[iCentralityBin],1,coutputQCv2_NN[iCentralityBin]);
      
      /*else if(doHigherHarmonic){
	TString outputQCv4_NN = fileName;
	outputQCv4_NN += ":outputQCv4analysis";
	outputQCv4_NN += "MinusMinus";
     
	taskQCv4_NN[iCentralityBin] = new AliAnalysisTaskQCumulants(Form("TaskQCumulant_%s",outputSlotNN_NameQCv4[iCentralityBin].Data()),kFALSE);
	taskQCv4_NN[iCentralityBin]->SetHarmonic(4);
	taskQCv4_NN[iCentralityBin]->SelectCollisionCandidates(triggerSelectionString);
	taskQCv4_NN[iCentralityBin]->SetUsePhiWeights(kFALSE); 
	taskQCv4_NN[iCentralityBin]->SetUsePtWeights(kFALSE);
	taskQCv4_NN[iCentralityBin]->SetUseEtaWeights(kFALSE); 
	taskQCv4_NN[iCentralityBin]->SetCalculateCumulantsVsM(kFALSE);
	taskQCv4_NN[iCentralityBin]->SetnBinsMult(10000);
	taskQCv4_NN[iCentralityBin]->SetMinMult(0.);
	taskQCv4_NN[iCentralityBin]->SetMaxMult(10000.);
	taskQCv4_NN[iCentralityBin]->SetApplyCorrectionForNUA(kTRUE);
	taskQCv4_NN[iCentralityBin]->SetFillMultipleControlHistograms(kFALSE);     
      
	coutputQCv4_NN[iCentralityBin] = mgr->CreateContainer(Form("%s",outputSlotNN_NameQCv4[iCentralityBin].Data()),TList::Class(),AliAnalysisManager::kOutputContainer,outputQCv4_NN);
	mgr->AddTask(taskQCv4_NN[iCentralityBin]);
	mgr->ConnectInput(taskQCv4_NN[iCentralityBin],0,flowEvent_NN[iCentralityBin]);
	mgr->ConnectInput(taskQCv4_NN[iCentralityBin],0,coutputFE_NN[iCentralityBin]);
	mgr->ConnectOutput(taskQCv4_NN[iCentralityBin],1,coutputQCv4_NN[iCentralityBin]);
	}*/
    }

   //------------------- PN --------------
    if(!mgr){
      AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    }
    
    //Check the analysis type using the event handlers connected to the analysis
    // manager. The availability of MC handler can also be checked here.
    if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskFlowEvent", "This task requires an input event handler");
      return NULL;
    }
    
    taskFE_PN[iCentralityBin] = new AliAnalysisTaskFlowEvent(Form("TaskFlowEvent_%s",suffixNamePN[iCentralityBin].Data()),"",doQA);
    taskFE_PN[iCentralityBin]->SelectCollisionCandidates(triggerSelectionString);
    
    //Sub events
    Double_t minA = -0.8;//
    Double_t maxA = -0.5*gEtaGap;//
    Double_t minB = +0.5*gEtaGap;//
    Double_t maxB = 0.8;//
    
    if(!isVZERO)
      taskFE_PN[iCentralityBin]->SetSubeventEtaRange(minA, maxA, minB, maxB);
    if(isVZERO)
      taskFE_PN[iCentralityBin]->SetSubeventEtaRange(-5,-1.5,+1.5,5);
    mgr->AddTask(taskFE_PN[iCentralityBin]);
    
    // Pass cuts for RPs and POIs to the task:
    taskFE_PN[iCentralityBin]->SetCutsEvent(cutsEvent_wQA[iCentralityBin]);
    taskFE_PN[iCentralityBin]->SetCutsRP(cutsRP_wQA[iCentralityBin]);
    taskFE_PN[iCentralityBin]->SetCutsPOI(cutsPOI_PNwQA[iCentralityBin]);

    if (cutsRP[iCentralityBin]->GetParamType()==AliFlowTrackCuts::kVZERO) {
      taskFE_PN[iCentralityBin]->SetHistWeightvsPhiMin(0.);
      taskFE_PN[iCentralityBin]->SetHistWeightvsPhiMax(200.);
    }
    
    if(!mgr){
      AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    }
    cinput1_PN[iCentralityBin] = mgr->GetCommonInputContainer();
    
    coutputFE_PN[iCentralityBin] = mgr->CreateContainer(Form("FlowEvent_%s",suffixNamePN[iCentralityBin].Data()),AliFlowEventSimple::Class(),AliAnalysisManager::kExchangeContainer);
    
    mgr->ConnectInput(taskFE_PN[iCentralityBin],0,cinput1_PN[iCentralityBin]);
    mgr->ConnectOutput(taskFE_PN[iCentralityBin],1,coutputFE_PN[iCentralityBin]);
    //==========================================================
    
    //======================================================================//
    flowEvent_PN[iCentralityBin] = mgr->CreateContainer(Form("Filter_%s",outputSlotPN_NameSPv2[iCentralityBin].Data()),AliFlowEventSimple::Class(),AliAnalysisManager::kExchangeContainer);
    taskFilter_PN[iCentralityBin] = new AliAnalysisTaskFilterFE(Form("TaskFilter_%s",outputSlotPN_NameSPv2[iCentralityBin].Data()),cutsRP[iCentralityBin], NULL);
    
    if(!isVZERO) 
      taskFilter_PN[iCentralityBin]->SetSubeventEtaRange(minA, maxA, minB, maxB);
    if(isVZERO) 
      taskFilter_PN[iCentralityBin]->SetSubeventEtaRange(-5,-1.5,+1.5,5);
    mgr->AddTask(taskFilter_PN[iCentralityBin]);
    mgr->ConnectInput(taskFilter_PN[iCentralityBin],0,coutputFE_PN[iCentralityBin]);
    mgr->ConnectOutput(taskFilter_PN[iCentralityBin],1,flowEvent_PN[iCentralityBin]);
    //======================================================================//
      
    //======================================================================//
    if(gRunSP) {

      TString outputSPv2_PN = fileName;
      outputSPv2_PN += ":outputSPv2analysis";
      outputSPv2_PN += "PlusMinus";

      taskSPv2_PN[iCentralityBin] = new AliAnalysisTaskScalarProduct(Form("TaskScalarProduct_%s",outputSlotPN_NameSPv2[iCentralityBin].Data()),kFALSE);
      taskSPv2_PN[iCentralityBin]->SetHarmonic(vnHarmonic);
      taskSPv2_PN[iCentralityBin]->SelectCollisionCandidates(triggerSelectionString);
      taskSPv2_PN[iCentralityBin]->SetRelDiffMsub(1.0);
      taskSPv2_PN[iCentralityBin]->SetTotalQvector(qVector);
      taskSPv2_PN[iCentralityBin]->SetApplyCorrectionForNUA(kTRUE);
      
      coutputSPv2_PN[iCentralityBin] = mgr->CreateContainer(Form("%s",outputSlotPN_NameSPv2[iCentralityBin].Data()),TList::Class(),AliAnalysisManager::kOutputContainer,outputSPv2_PN);
      mgr->AddTask(taskSPv2_PN[iCentralityBin]);
      mgr->ConnectInput(taskSPv2_PN[iCentralityBin],0,flowEvent_PN[iCentralityBin]);
      mgr->ConnectInput(taskSPv2_PN[iCentralityBin],0,coutputFE_PN[iCentralityBin]);
      mgr->ConnectOutput(taskSPv2_PN[iCentralityBin],1,coutputSPv2_PN[iCentralityBin]);

      /*else if(doHigherHarmonic){
	TString outputSPv4_PN = fileName;
	outputSPv4_PN += ":outputSPv4analysis";
	outputSPv4_PN += "PlusMinus";
      
	taskSPv4_PN[iCentralityBin] = new AliAnalysisTaskScalarProduct(Form("TaskScalarProduct_%s",outputSlotPN_NameSPv4[iCentralityBin].Data()),kFALSE);
	taskSPv4_PN[iCentralityBin]->SetHarmonic(4);
	taskSPv4_PN[iCentralityBin]->SelectCollisionCandidates(triggerSelectionString);
	taskSPv4_PN[iCentralityBin]->SetRelDiffMsub(1.0);
	taskSPv4_PN[iCentralityBin]->SetTotalQvector(qVector);
	taskSPv4_PN[iCentralityBin]->SetApplyCorrectionForNUA(kTRUE);
      
	coutputSPv4_PN[iCentralityBin] = mgr->CreateContainer(Form("%s",outputSlotPN_NameSPv4[iCentralityBin].Data()),TList::Class(),AliAnalysisManager::kOutputContainer,outputSPv4_PN);
	mgr->AddTask(taskSPv4_PN[iCentralityBin]);
	mgr->ConnectInput(taskSPv4_PN[iCentralityBin],0,flowEvent_PN[iCentralityBin]);
	mgr->ConnectInput(taskSPv4_PN[iCentralityBin],0,coutputFE_PN[iCentralityBin]);
	mgr->ConnectOutput(taskSPv4_PN[iCentralityBin],1,coutputSPv4_PN[iCentralityBin]);
	}*/
    }
    //======================================================================//

    //======================================================================//
    if(gRunQC) {

      TString outputQCv2_PN = fileName;
      outputQCv2_PN += ":outputQCv2analysis";
      outputQCv2_PN += "PlusMinus";

      taskQCv2_PN[iCentralityBin] = new AliAnalysisTaskQCumulants(Form("TaskQCumulant_%s",outputSlotPN_NameQCv2[iCentralityBin].Data()),kFALSE);
      taskQCv2_PN[iCentralityBin]->SetHarmonic(vnHarmonic);
      taskQCv2_PN[iCentralityBin]->SelectCollisionCandidates(triggerSelectionString);
      taskQCv2_PN[iCentralityBin]->SetUsePhiWeights(kFALSE); 
      taskQCv2_PN[iCentralityBin]->SetUsePtWeights(kFALSE);
      taskQCv2_PN[iCentralityBin]->SetUseEtaWeights(kFALSE); 
      taskQCv2_PN[iCentralityBin]->SetCalculateCumulantsVsM(kFALSE);
      taskQCv2_PN[iCentralityBin]->SetnBinsMult(10000);
      taskQCv2_PN[iCentralityBin]->SetMinMult(0.);
      taskQCv2_PN[iCentralityBin]->SetMaxMult(10000.);
      taskQCv2_PN[iCentralityBin]->SetApplyCorrectionForNUA(kTRUE);
      taskQCv2_PN[iCentralityBin]->SetFillMultipleControlHistograms(kFALSE);     
      
      coutputQCv2_PN[iCentralityBin] = mgr->CreateContainer(Form("%s",outputSlotPN_NameQCv2[iCentralityBin].Data()),TList::Class(),AliAnalysisManager::kOutputContainer,outputQCv2_PN);
      mgr->AddTask(taskQCv2_PN[iCentralityBin]);
      mgr->ConnectInput(taskQCv2_PN[iCentralityBin],0,flowEvent_PN[iCentralityBin]);
      mgr->ConnectInput(taskQCv2_PN[iCentralityBin],0,coutputFE_PN[iCentralityBin]);
      mgr->ConnectOutput(taskQCv2_PN[iCentralityBin],1,coutputQCv2_PN[iCentralityBin]);

      /*else if(doHigherHarmonic){
	TString outputQCv4_PN = fileName;
	outputQCv4_PN += ":outputQCv4analysis";
	outputQCv4_PN += "PlusMinus";

	taskQCv4_PN[iCentralityBin] = new AliAnalysisTaskQCumulants(Form("TaskQCumulant_%s",outputSlotPN_NameQCv4[iCentralityBin].Data()),kFALSE);
	taskQCv4_PN[iCentralityBin]->SetHarmonic(4);
	taskQCv4_PN[iCentralityBin]->SelectCollisionCandidates(triggerSelectionString);
	taskQCv4_PN[iCentralityBin]->SetUsePhiWeights(kFALSE); 
	taskQCv4_PN[iCentralityBin]->SetUsePtWeights(kFALSE);
	taskQCv4_PN[iCentralityBin]->SetUseEtaWeights(kFALSE); 
	taskQCv4_PN[iCentralityBin]->SetCalculateCumulantsVsM(kFALSE);
	taskQCv4_PN[iCentralityBin]->SetnBinsMult(10000);
	taskQCv4_PN[iCentralityBin]->SetMinMult(0.);
	taskQCv4_PN[iCentralityBin]->SetMaxMult(10000.);
	taskQCv4_PN[iCentralityBin]->SetApplyCorrectionForNUA(kTRUE);
	taskQCv4_PN[iCentralityBin]->SetFillMultipleControlHistograms(kFALSE);     
      
	coutputQCv4_PN[iCentralityBin] = mgr->CreateContainer(Form("%s",outputSlotPN_NameQCv4[iCentralityBin].Data()),TList::Class(),AliAnalysisManager::kOutputContainer,outputQCv4_PN);
	mgr->AddTask(taskQCv4_PN[iCentralityBin]);
	mgr->ConnectInput(taskQCv4_PN[iCentralityBin],0,flowEvent_PN[iCentralityBin]);
	mgr->ConnectInput(taskQCv4_PN[iCentralityBin],0,coutputFE_PN[iCentralityBin]);
	mgr->ConnectOutput(taskQCv4_PN[iCentralityBin],1,coutputQCv4_PN[iCentralityBin]);
	}*/
    }

    //======================================================================//















    //============================== now CME correlator ================================//

    //---------- PP -------------
      if(!doHigherHarmonic){
	TString outputMHLS_PP = fileName;
	outputMHLS_PP += ":outputMHLSanalysis";
	outputMHLS_PP += "PlusPlus";
	taskMHLS_PP[iCentralityBin] = new AliAnalysisTaskMixedHarmonics(Form("TaskMixedHarmonicsLS_%s",outputSlotPP_NameMHLS[iCentralityBin].Data()),kFALSE);
	taskMHLS_PP[iCentralityBin]->SetHarmonic(1); // n in cos[n(phi1+phi2-2phi3)] and cos[n(psi1+psi2-2phi3)]
	taskMHLS_PP[iCentralityBin]->SetNoOfMultipicityBins(10000);
	taskMHLS_PP[iCentralityBin]->SetMultipicityBinWidth(1.);
	taskMHLS_PP[iCentralityBin]->SetMinMultiplicity(1.);
	taskMHLS_PP[iCentralityBin]->SetCorrectForDetectorEffects(kTRUE);
	taskMHLS_PP[iCentralityBin]->SetEvaluateDifferential3pCorrelator(kTRUE); // evaluate <<cos[n(psi1+psi2-2phi3)]>> (Remark: two nested loops)
	taskMHLS_PP[iCentralityBin]->SetOppositeChargesPOI(kFALSE); //
	if(isPbPb){
	  taskMHLS_PP[iCentralityBin]->SetRejectPileUp(checkPileup);
	  taskMHLS_PP[iCentralityBin]->SetRejectPileUpTight(checkPileup);
	}
	else{
	  taskMHLS_PP[iCentralityBin]->SetRejectPileUp(kFALSE);
	  taskMHLS_PP[iCentralityBin]->SetRejectPileUpTight(kFALSE);
	}

	coutputMHLS_PP[iCentralityBin] = mgr->CreateContainer(Form("%s",outputSlotPP_NameMHLS[iCentralityBin].Data()),TList::Class(),AliAnalysisManager::kOutputContainer,outputMHLS_PP);
	mgr->AddTask(taskMHLS_PP[iCentralityBin]);
	mgr->ConnectInput(taskMHLS_PP[iCentralityBin],0,flowEvent_PP[iCentralityBin]);
	mgr->ConnectInput(taskMHLS_PP[iCentralityBin],0,coutputFE_PP[iCentralityBin]);
	mgr->ConnectOutput(taskMHLS_PP[iCentralityBin],1,coutputMHLS_PP[iCentralityBin]);
      }
      else if(doHigherHarmonic){
	TString outputMHLS_PP2 = fileName;
	outputMHLS_PP2 += ":outputMHLS2analysis";
	outputMHLS_PP2 += "PlusPlus";     
	taskMHLS_PP2[iCentralityBin] = new AliAnalysisTaskMixedHarmonics(Form("TaskMixedHarmonicsLS_%s",outputSlotPP_NameMHLS2[iCentralityBin].Data()),kFALSE);
	taskMHLS_PP2[iCentralityBin]->SetHarmonic(2); // n in cos[n(phi1+phi2-2phi3)] and cos[n(psi1+psi2-2phi3)]
	taskMHLS_PP2[iCentralityBin]->SetNoOfMultipicityBins(10000);
	taskMHLS_PP2[iCentralityBin]->SetMultipicityBinWidth(1.);
	taskMHLS_PP2[iCentralityBin]->SetMinMultiplicity(1.);
	taskMHLS_PP2[iCentralityBin]->SetCorrectForDetectorEffects(kTRUE);
	taskMHLS_PP2[iCentralityBin]->SetEvaluateDifferential3pCorrelator(kTRUE); // evaluate <<cos[n(psi1+psi2-2phi3)]>> (Remark: two nested loops)
	taskMHLS_PP2[iCentralityBin]->SetOppositeChargesPOI(kFALSE); //
	if(isPbPb){
	  taskMHLS_PP2[iCentralityBin]->SetRejectPileUp(checkPileup);
	  taskMHLS_PP2[iCentralityBin]->SetRejectPileUpTight(checkPileup);
	}
	else{
	  taskMHLS_PP2[iCentralityBin]->SetRejectPileUp(kFALSE);
	  taskMHLS_PP2[iCentralityBin]->SetRejectPileUpTight(kFALSE);
	}

	coutputMHLS_PP2[iCentralityBin] = mgr->CreateContainer(Form("%s",outputSlotPP_NameMHLS2[iCentralityBin].Data()),TList::Class(),AliAnalysisManager::kOutputContainer,outputMHLS_PP2);
	mgr->AddTask(taskMHLS_PP2[iCentralityBin]);
	mgr->ConnectInput(taskMHLS_PP2[iCentralityBin],0,flowEvent_PP[iCentralityBin]);
	mgr->ConnectInput(taskMHLS_PP2[iCentralityBin],0,coutputFE_PP[iCentralityBin]);
	mgr->ConnectOutput(taskMHLS_PP2[iCentralityBin],1,coutputMHLS_PP2[iCentralityBin]);
      }

      //-------- NN -----------
      if(!doHigherHarmonic){
	TString outputMHLS_NN = fileName;
	outputMHLS_NN += ":outputMHLSanalysis";
	outputMHLS_NN += "MinusMinus";
	taskMHLS_NN[iCentralityBin] = new AliAnalysisTaskMixedHarmonics(Form("TaskMixedHarmonicsLS_%s",outputSlotNN_NameMHLS[iCentralityBin].Data()),kFALSE);
	taskMHLS_NN[iCentralityBin]->SetHarmonic(1); // n in cos[n(phi1+phi2-2phi3)] and cos[n(psi1+psi2-2phi3)]
	taskMHLS_NN[iCentralityBin]->SetNoOfMultipicityBins(10000);
	taskMHLS_NN[iCentralityBin]->SetMultipicityBinWidth(1.);
	taskMHLS_NN[iCentralityBin]->SetMinMultiplicity(1.);
	taskMHLS_NN[iCentralityBin]->SetCorrectForDetectorEffects(kTRUE);
	taskMHLS_NN[iCentralityBin]->SetEvaluateDifferential3pCorrelator(kTRUE); // evaluate <<cos[n(psi1+psi2-2phi3)]>> (Remark: two nested loops)
	taskMHLS_NN[iCentralityBin]->SetOppositeChargesPOI(kFALSE); //
	if(isPbPb){
	  taskMHLS_NN[iCentralityBin]->SetRejectPileUp(checkPileup);
	  taskMHLS_NN[iCentralityBin]->SetRejectPileUpTight(checkPileup);
	}
	else{
	  taskMHLS_NN[iCentralityBin]->SetRejectPileUp(kFALSE);
	  taskMHLS_NN[iCentralityBin]->SetRejectPileUpTight(kFALSE);
	}

	coutputMHLS_NN[iCentralityBin] = mgr->CreateContainer(Form("%s",outputSlotNN_NameMHLS[iCentralityBin].Data()),TList::Class(),AliAnalysisManager::kOutputContainer,outputMHLS_NN);
	mgr->AddTask(taskMHLS_NN[iCentralityBin]);
	mgr->ConnectInput(taskMHLS_NN[iCentralityBin],0,flowEvent_NN[iCentralityBin]);
	mgr->ConnectInput(taskMHLS_NN[iCentralityBin],0,coutputFE_NN[iCentralityBin]);
	mgr->ConnectOutput(taskMHLS_NN[iCentralityBin],1,coutputMHLS_NN[iCentralityBin]);
      }
      else if(doHigherHarmonic){
	TString outputMHLS_NN2 = fileName;
	outputMHLS_NN2 += ":outputMHLS2analysis";
	outputMHLS_NN2 += "MinusMinus";     
	taskMHLS_NN2[iCentralityBin] = new AliAnalysisTaskMixedHarmonics(Form("TaskMixedHarmonicsLS_%s",outputSlotNN_NameMHLS2[iCentralityBin].Data()),kFALSE);
	taskMHLS_NN2[iCentralityBin]->SetHarmonic(2); // n in cos[n(phi1+phi2-2phi3)] and cos[n(psi1+psi2-2phi3)]
	taskMHLS_NN2[iCentralityBin]->SetNoOfMultipicityBins(10000);
	taskMHLS_NN2[iCentralityBin]->SetMultipicityBinWidth(1.);
	taskMHLS_NN2[iCentralityBin]->SetMinMultiplicity(1.);
	taskMHLS_NN2[iCentralityBin]->SetCorrectForDetectorEffects(kTRUE);
	taskMHLS_NN2[iCentralityBin]->SetEvaluateDifferential3pCorrelator(kTRUE); // evaluate <<cos[n(psi1+psi2-2phi3)]>> (Remark: two nested loops)
	taskMHLS_NN2[iCentralityBin]->SetOppositeChargesPOI(kFALSE); //
	if(isPbPb){
	  taskMHLS_NN2[iCentralityBin]->SetRejectPileUp(checkPileup);
	  taskMHLS_NN2[iCentralityBin]->SetRejectPileUpTight(checkPileup);
	}
	else{
	  taskMHLS_NN2[iCentralityBin]->SetRejectPileUp(kFALSE);
	  taskMHLS_NN2[iCentralityBin]->SetRejectPileUpTight(kFALSE);
	}



	coutputMHLS_NN2[iCentralityBin] = mgr->CreateContainer(Form("%s",outputSlotNN_NameMHLS2[iCentralityBin].Data()),TList::Class(),AliAnalysisManager::kOutputContainer,outputMHLS_NN2);
	mgr->AddTask(taskMHLS_NN2[iCentralityBin]);
	mgr->ConnectInput(taskMHLS_NN2[iCentralityBin],0,flowEvent_NN[iCentralityBin]);
	mgr->ConnectInput(taskMHLS_NN2[iCentralityBin],0,coutputFE_NN[iCentralityBin]);
	mgr->ConnectOutput(taskMHLS_NN2[iCentralityBin],1,coutputMHLS_NN2[iCentralityBin]);
      }
      //--------- PN --------------
      if(!doHigherHarmonic){
	TString outputMHUS = fileName;
	outputMHUS += ":outputMHUSanalysis";

	taskMHUS[iCentralityBin] = new AliAnalysisTaskMixedHarmonics(Form("TaskMixedHarmonicsLS_%s",outputSlotPN_NameMHUS[iCentralityBin].Data()),kFALSE);
	taskMHUS[iCentralityBin]->SetHarmonic(1); // n in cos[n(phi1+phi2-2phi3)] and cos[n(psi1+psi2-2phi3)]
	taskMHUS[iCentralityBin]->SetNoOfMultipicityBins(10000);
	taskMHUS[iCentralityBin]->SetMultipicityBinWidth(1.);
	taskMHUS[iCentralityBin]->SetMinMultiplicity(1.);
	taskMHUS[iCentralityBin]->SetCorrectForDetectorEffects(kTRUE);
	taskMHUS[iCentralityBin]->SetEvaluateDifferential3pCorrelator(kTRUE); // evaluate <<cos[n(psi1+psi2-2phi3)]>> (Remark: two nested loops)
	taskMHUS[iCentralityBin]->SetOppositeChargesPOI(kTRUE); //
 	taskMHUS[iCentralityBin]->SetFillQAHistograms(doQA);
	if(isPbPb){
	  taskMHUS[iCentralityBin]->SetRejectPileUp(checkPileup);
	  taskMHUS[iCentralityBin]->SetRejectPileUpTight(checkPileup);
	}
	else{
	  taskMHUS[iCentralityBin]->SetRejectPileUp(kFALSE);
	  taskMHUS[iCentralityBin]->SetRejectPileUpTight(kFALSE);
	}

      
	coutputMHUS[iCentralityBin] = mgr->CreateContainer(Form("%s",outputSlotPN_NameMHUS[iCentralityBin].Data()),TList::Class(),AliAnalysisManager::kOutputContainer,outputMHUS);
	mgr->AddTask(taskMHUS[iCentralityBin]);
	mgr->ConnectInput(taskMHUS[iCentralityBin],0,flowEvent_PN[iCentralityBin]);
	mgr->ConnectInput(taskMHUS[iCentralityBin],0,coutputFE_PN[iCentralityBin]);
	mgr->ConnectOutput(taskMHUS[iCentralityBin],1,coutputMHUS[iCentralityBin]);
      }
      else if(doHigherHarmonic){
	TString outputMHUS2 = fileName;
	outputMHUS2 += ":outputMHUS2analysis";
      
	taskMHUS2[iCentralityBin] = new AliAnalysisTaskMixedHarmonics(Form("TaskMixedHarmonicsLS_%s",outputSlotPN_NameMHUS2[iCentralityBin].Data()),kFALSE);
	taskMHUS2[iCentralityBin]->SetHarmonic(2); // n in cos[n(phi1+phi2-2phi3)] and cos[n(psi1+psi2-2phi3)]
	taskMHUS2[iCentralityBin]->SetNoOfMultipicityBins(10000);
	taskMHUS2[iCentralityBin]->SetMultipicityBinWidth(1.);
	taskMHUS2[iCentralityBin]->SetMinMultiplicity(1.);
	taskMHUS2[iCentralityBin]->SetCorrectForDetectorEffects(kTRUE);
	taskMHUS2[iCentralityBin]->SetEvaluateDifferential3pCorrelator(kTRUE); // evaluate <<cos[n(psi1+psi2-2phi3)]>> (Remark: two nested loops)
	taskMHUS2[iCentralityBin]->SetOppositeChargesPOI(kTRUE); //
 	taskMHUS2[iCentralityBin]->SetFillQAHistograms(doQA);
	if(isPbPb){
	  taskMHUS2[iCentralityBin]->SetRejectPileUp(checkPileup);
	  taskMHUS2[iCentralityBin]->SetRejectPileUpTight(checkPileup);
	}
	else{
	  taskMHUS2[iCentralityBin]->SetRejectPileUp(kFALSE);
	  taskMHUS2[iCentralityBin]->SetRejectPileUpTight(kFALSE);
	}


	coutputMHUS2[iCentralityBin] = mgr->CreateContainer(Form("%s",outputSlotPN_NameMHUS2[iCentralityBin].Data()),TList::Class(),AliAnalysisManager::kOutputContainer,outputMHUS2);
	mgr->AddTask(taskMHUS2[iCentralityBin]);
	mgr->ConnectInput(taskMHUS2[iCentralityBin],0,flowEvent_PN[iCentralityBin]);
	mgr->ConnectInput(taskMHUS2[iCentralityBin],0,coutputFE_PN[iCentralityBin]);
	mgr->ConnectOutput(taskMHUS2[iCentralityBin],1,coutputMHUS2[iCentralityBin]);
      }

    //======================================================================//
    


    //============================== QAs ====================================//
    if(taskFE_PP[iCentralityBin]->GetQAOn()) {
      outputQA_PP[iCentralityBin] = fileName;
      outputQA_PP[iCentralityBin] += ":QA";
      coutputFE_PPQA[iCentralityBin] = mgr->CreateContainer(Form("QA_%s",suffixNamePP[iCentralityBin].Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputQA_PP[iCentralityBin]);
      mgr->ConnectOutput(taskFE_PP[iCentralityBin],2,coutputFE_PPQA[iCentralityBin]);
    }
    if(taskFE_NN[iCentralityBin]->GetQAOn()) {
      outputQA_NN[iCentralityBin] = fileName;
      outputQA_NN[iCentralityBin] += ":QA";
      coutputFE_NNQA[iCentralityBin] = mgr->CreateContainer(Form("QA_%s",suffixNameNN[iCentralityBin].Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputQA_NN[iCentralityBin]);
      mgr->ConnectOutput(taskFE_NN[iCentralityBin],2,coutputFE_NNQA[iCentralityBin]);
    }
    if(taskFE_PN[iCentralityBin]->GetQAOn()) {
      outputQA_PN[iCentralityBin] = fileName;
      outputQA_PN[iCentralityBin] += ":QA";
      coutputFE_PNQA[iCentralityBin] = mgr->CreateContainer(Form("QA_%s",suffixNamePN[iCentralityBin].Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputQA_PN[iCentralityBin]);
      mgr->ConnectOutput(taskFE_PN[iCentralityBin],2,coutputFE_PNQA[iCentralityBin]);
    }

    //======================================================================//
  


  } //loop over centralities

   

}

//_________________________________________________________//
AliFlowEventCuts *createFlowEventCutObject(Int_t gCentralityMin = -1,
                                           Int_t gCentralityMax = -1,
                                           Bool_t isPbPb = kTRUE,
                                           Double_t whichData = 2011,
                                           AliFlowEventCuts::refMultMethod gCentralityEstimator = AliFlowEventCuts::kVZERO,
                                           Bool_t doQA = kFALSE, Bool_t bPileup=kFALSE) {
    //Part of the code that creates the event cut objects
    Double_t gVertexZmin = -10., gVertexZmax = 10.;
    
    //Create the event cut objects
    AliFlowEventCuts *cutsEvent = new AliFlowEventCuts(Form("eventcutsCentrality%dTo%d",gCentralityMin,gCentralityMax));

    if(isPbPb) {
      //cutsEvent->SetCentralityPercentileRange(gCentralityMin,gCentralityMax);
      cutsEvent->SetCentralityPercentileMethod(gCentralityEstimator);

      cutsEvent->SetCentralityPercentileRange(gCentralityMin,gCentralityMax,kFALSE);

      if(whichData==2011){ 
	cutsEvent->SetLHC11h(kTRUE);
      }
      else if(whichData!=2011){
	cutsEvent->SetLHC11h(kFALSE);
	if(whichData==2015){ 
	  cutsEvent->SetCentralityPercentileRange(gCentralityMin,gCentralityMax,kTRUE);
        //cutsEvent->SetCheckPileup(bPileup); //
	}
      }
      cutsEvent->SetCutTPCmultiplicityOutliersAOD(kFALSE);
    }
    else{
       cutsEvent->SetCentralityPercentileMethod(gCentralityEstimator);
       cutsEvent->SetCentralityPercentileRange(gCentralityMin,gCentralityMax,kTRUE);
       cutsEvent->SetCutTPCmultiplicityOutliersAOD(kFALSE);
    }  
    
    cutsEvent->SetPrimaryVertexZrange(gVertexZmin,gVertexZmax);
    cutsEvent->SetQA(doQA);
    
    //return the object
    return cutsEvent;
}

//_________________________________________________________//
AliFlowTrackCuts *createFlowRPCutObject(Int_t gCentralityMin = -1,
                                        Int_t gCentralityMax = -1,
                                        Bool_t isVZERO = kFALSE,
                                        Int_t gAODfilterBit = 768,
					Double_t whichData = 2011,
                                        Double_t gChargeRP = 0.,
                                        Bool_t doQA = kFALSE) {
    //Part of the code that creates the RP cut objects
    Double_t gEtaMin = -0.8, gEtaMax = 0.8;
    Int_t gMinTPCdedx = 10;
   //  Int_t gMinTPCdedx = -999;
 
    //Create the event cut objects
    AliFlowTrackCuts *cutsRP = new AliFlowTrackCuts(Form("rpCutsCentrality%dTo%d",gCentralityMin,gCentralityMax));
    if(!isVZERO) {
        cutsRP->SetPtRange(0.2,5.);
        cutsRP->SetEtaRange(gEtaMin,gEtaMax);
        cutsRP->SetMinimalTPCdedx(gMinTPCdedx);
        cutsRP->SetAODfilterBit(gAODfilterBit);       
        if (gChargeRP!=0) cutsRP->SetCharge(gChargeRP);
    }
    else if(isVZERO) { // use vzero sub analysis
        if(whichData = 2010)
            cutsRP = AliFlowTrackCuts::GetStandardVZEROOnlyTrackCuts2010();
        if(whichData = 2011)
            cutsRP = AliFlowTrackCuts::GetStandardVZEROOnlyTrackCuts2011();
	else{ cout<<"VZEROonly track cuts for this dataset are automatically set to 2010 settings!"<<endl; cutsRP = AliFlowTrackCuts::GetStandardVZEROOnlyTrackCuts2010();}
   }//VZERO
    
    
    cutsRP->SetQA(doQA);
    
    //return the object
    return cutsRP;
}

//_________________________________________________________//
AliFlowTrackCuts *createFlowPOICutObject(Int_t gCentralityMin = -1,
                                         Int_t gCentralityMax = -1,
                                         TString gQvector = "Qa",
                                         Double_t gEtaGap = 0.0,
                                         Bool_t isVZERO = kFALSE,
                                         Bool_t isPbPb = kTRUE,
                                         Int_t gAODfilterBit = 768,
					 Double_t whichData = 2011,
					 Double_t gDCAvtxXY = -1.,
                                         Double_t gDCAvtxZ = -1.,
                                         Double_t gChargePOI = 0.,
                                         Bool_t doQA = kFALSE) {
    //Part of the code that creates the POI cut objects
    Double_t gEtaMin = -0.8, gEtaMax = 0.8;
    Int_t gMinTPCdedx = 10;

    AliFlowTrackCuts *cutsPOI = new AliFlowTrackCuts(Form("poiCutsCentrality%dTo%d",gCentralityMin,gCentralityMax));
        
    /*if(!isVZERO) {
        if(gQvector.Contains("Qa")) {
            //cout<<"(Panos): "<<gQvector.Data()<<endl;
            cutsPOI->SetEtaRange(+0.5*gEtaGap,gEtaMax );
        }
        else if(gQvector.Contains("Qb")) {
            //cout<<"(Panos): "<<gQvector.Data()<<endl;
            cutsPOI->SetEtaRange(gEtaMin,-0.5*gEtaGap);
        }
    }
    else if(isVZERO)
    cutsPOI->SetEtaRange(gEtaMin,gEtaMax);*/
    
    
    cutsPOI->SetAcceptKinkDaughters(kFALSE);
    cutsPOI->SetPtRange(0.2,5.);
    cutsPOI->SetEtaRange(gEtaMin,gEtaMax);
    cutsPOI->SetMinimalTPCdedx(gMinTPCdedx);
    cutsPOI->SetAODfilterBit(gAODfilterBit);       
    if (gChargePOI!=0) 
      cutsPOI->SetCharge(gChargePOI);
    if(gDCAvtxXY > 0)
        cutsPOI->SetMaxDCAToVertexXY(gDCAvtxXY);
    if(gDCAvtxZ > 0)
        cutsPOI->SetMaxDCAToVertexZ(gDCAvtxZ);
    cutsPOI->SetQA(doQA);
    
    //return the object
    return cutsPOI;
}
