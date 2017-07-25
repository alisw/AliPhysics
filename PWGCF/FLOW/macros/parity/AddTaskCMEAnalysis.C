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
			Bool_t doQA = kFALSE) {
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
  //=============================================//
      
  //=============================================//
  //Define centralities
  static const Int_t nCentralities = 10;
  Double_t gCentrality[nCentralities];
  
  if(isPbPb) {
    if(!isData) {
      for(Int_t i = 0; i < 9; i++)
	gCentrality[i] = 0 + i*10;
    }
    else {
      gCentrality[0] = 0; gCentrality[1] = 5;
      for(Int_t i = 2; i < 10; i++)
	gCentrality[i] = (i - 1)*10;
    }
  }
  
  //======================================================================//
  // FLOW cut objects config
  AliFlowEventCuts* cutsEvent[nCentralities];
  AliFlowTrackCuts* cutsRP[nCentralities];
  AliFlowTrackCuts* cutsPOI[nCentralities];
  TString suffixName[nCentralities];
  TString outputSlotNameSPv2[nCentralities];
  TString outputSlotNameSPv4[nCentralities];
  TString outputSlotNameQCv2[nCentralities];
  TString outputSlotNameQCv4[nCentralities];
  TString outputSlotNameMHLS[nCentralities];
  TString outputSlotNameMHLS2[nCentralities];
  TString outputSlotNameMHUS[nCentralities];
  TString outputSlotNameMHUS2[nCentralities];
  
  for(Int_t iCentralityBin = 0; iCentralityBin < nCentralities - 1; iCentralityBin++) {
    //Create the event cut object
    cutsEvent[iCentralityBin] = createFlowEventCutObject(gCentrality[iCentralityBin],gCentrality[iCentralityBin+1],isPbPb,whichData,gCentralityEstimator,doQA);
    
    //Create the RP cut object
    cutsRP[iCentralityBin] = createFlowRPCutObject(gCentrality[iCentralityBin],gCentrality[iCentralityBin+1],isVZERO,gAODfilterBit,whichData,gChargeRP,doQA);
    
    //Create the POI cut object
    cutsPOI[iCentralityBin] = createFlowPOICutObject(gCentrality[iCentralityBin],gCentrality[iCentralityBin+1],qVector,gEtaGap,isVZERO,isPbPb,gAODfilterBit,whichData,gDCAvtxXY,gDCAvtxZ,gChargePOI,doQA);
    
    if(isPbPb) {
      suffixName[iCentralityBin] += "Centrality";
      suffixName[iCentralityBin] += gCentrality[iCentralityBin];
      suffixName[iCentralityBin] += "To";
      suffixName[iCentralityBin] += gCentrality[iCentralityBin+1];
      if(gChargePOI == 1) 
	suffixName[iCentralityBin] += "PlusPlus";
      else if(gChargePOI == -1) 
	suffixName[iCentralityBin] += "MinusMinus";
    }
    
    //Giving names to output slots
    outputSlotNameSPv2[iCentralityBin] = Form("%s_SP_v2",suffixName[iCentralityBin].Data());
    outputSlotNameSPv4[iCentralityBin] = Form("%s_SP_v4",suffixName[iCentralityBin].Data());
    outputSlotNameQCv2[iCentralityBin] = Form("%s_QC_v2",suffixName[iCentralityBin].Data());
    outputSlotNameQCv4[iCentralityBin] = Form("%s_QC_v4",suffixName[iCentralityBin].Data());
    outputSlotNameMHLS[iCentralityBin] = Form("%s_MHLS_v1",suffixName[iCentralityBin].Data());
    outputSlotNameMHLS2[iCentralityBin] = Form("%s_MHLS_v2",suffixName[iCentralityBin].Data());
    outputSlotNameMHUS[iCentralityBin] = Form("%s_MHUS_v1",suffixName[iCentralityBin].Data());
    outputSlotNameMHUS2[iCentralityBin] = Form("%s_MHUS_v2",suffixName[iCentralityBin].Data());
  }//loop over centrality classes

  TString fileName = "AnalysisResults";
  fileName.Append(".root");


  //========================FLOWPACKAGE TASKS=========================//
  AliAnalysisDataContainer *cinput1[nCentralities];
  AliAnalysisDataContainer *coutputFE[nCentralities];
  AliAnalysisDataContainer *coutputFEQA[nCentralities];
  AliAnalysisTaskFlowEvent *taskFE[nCentralities];

  AliAnalysisDataContainer *flowEvent[nCentralities];
  AliAnalysisTaskFilterFE *taskFilter[nCentralities];

  AliAnalysisDataContainer *coutputSPv2[nCentralities];
  AliAnalysisTaskScalarProduct *taskSPv2[nCentralities];
  AliAnalysisDataContainer *coutputSPv4[nCentralities];
  AliAnalysisTaskScalarProduct *taskSPv4[nCentralities];
  AliAnalysisDataContainer *coutputQCv2[nCentralities];
  AliAnalysisTaskQCumulants *taskQCv2[nCentralities];
  AliAnalysisDataContainer *coutputQCv4[nCentralities];
  AliAnalysisTaskQCumulants *taskQCv4[nCentralities];
  AliAnalysisDataContainer *coutputMHLS[nCentralities];
  AliAnalysisTaskMixedHarmonics *taskMHLS[nCentralities];
  AliAnalysisDataContainer *coutputMHLS2[nCentralities];
  AliAnalysisTaskMixedHarmonics *taskMHLS2[nCentralities];
  AliAnalysisDataContainer *coutputMHUS[nCentralities];
  AliAnalysisTaskMixedHarmonics *taskMHUS[nCentralities];
  AliAnalysisDataContainer *coutputMHUS2[nCentralities];
  AliAnalysisTaskMixedHarmonics *taskMHUS2[nCentralities];

  TString outputQA[nCentralities];
  TString myNameSPv2[nCentralities];
  TString myNameSPv4[nCentralities];
  TString myNameQCv2[nCentralities];
  TString myNameQCv4[nCentralities];
  TString myNameMHLS[nCentralities];
  TString myNameMHLS2[nCentralities];
  TString myNameMHUS[nCentralities];
  TString myNameMHUS2[nCentralities];
  TString slot[nCentralities];
  
  for(Int_t iCentralityBin = 0; iCentralityBin < nCentralities - 1; iCentralityBin++) {
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
    
    taskFE[iCentralityBin] = new AliAnalysisTaskFlowEvent(Form("TaskFlowEvent_%s",suffixName[iCentralityBin].Data()),"",doQA);
    taskFE[iCentralityBin]->SelectCollisionCandidates(triggerSelectionString);
    
    //Sub events
    Double_t minA = -0.8;//
    Double_t maxA = -0.5*gEtaGap;//
    Double_t minB = +0.5*gEtaGap;//
    Double_t maxB = 0.8;//
    
    if(!isVZERO)
      taskFE[iCentralityBin]->SetSubeventEtaRange(minA, maxA, minB, maxB);
    if(isVZERO)
      taskFE[iCentralityBin]->SetSubeventEtaRange(-5,-1.5,+1.5,5);
    mgr->AddTask(taskFE[iCentralityBin]);
    
    // Pass cuts for RPs and POIs to the task:
    taskFE[iCentralityBin]->SetCutsEvent(cutsEvent[iCentralityBin]);
    taskFE[iCentralityBin]->SetCutsRP(cutsRP[iCentralityBin]);
    taskFE[iCentralityBin]->SetCutsPOI(cutsPOI[iCentralityBin]);
    if (cutsRP[iCentralityBin]->GetParamType()==AliFlowTrackCuts::kVZERO) {
      taskFE[iCentralityBin]->SetHistWeightvsPhiMin(0.);
      taskFE[iCentralityBin]->SetHistWeightvsPhiMax(200.);
    }
    
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    cinput1[iCentralityBin] = mgr->GetCommonInputContainer();
    
    coutputFE[iCentralityBin] = mgr->CreateContainer(Form("FlowEvent_%s",suffixName[iCentralityBin].Data()),AliFlowEventSimple::Class(),AliAnalysisManager::kExchangeContainer);
    
    mgr->ConnectInput(taskFE[iCentralityBin],0,cinput1[iCentralityBin]);
    mgr->ConnectOutput(taskFE[iCentralityBin],1,coutputFE[iCentralityBin]);
    //==========================================================
    
    //======================================================================//
    flowEvent[iCentralityBin] = mgr->CreateContainer(Form("Filter_%s",outputSlotNameSPv2[iCentralityBin].Data()),AliFlowEventSimple::Class(),AliAnalysisManager::kExchangeContainer);
    taskFilter[iCentralityBin] = new AliAnalysisTaskFilterFE(Form("TaskFilter_%s",outputSlotNameSPv2[iCentralityBin].Data()),cutsRP[iCentralityBin], NULL);
    
    if(!isVZERO) 
      taskFilter[iCentralityBin]->SetSubeventEtaRange(minA, maxA, minB, maxB);
    if(isVZERO) 
      taskFilter[iCentralityBin]->SetSubeventEtaRange(-5,-1.5,+1.5,5);
    mgr->AddTask(taskFilter[iCentralityBin]);
    mgr->ConnectInput(taskFilter[iCentralityBin],0,coutputFE[iCentralityBin]);
    mgr->ConnectOutput(taskFilter[iCentralityBin],1,flowEvent[iCentralityBin]);
    //======================================================================//
      
    //======================================================================//
    if(gRunSP) {
      TString outputSPv2 = fileName;
      outputSPv2 += ":outputSPv2analysis";
      if(gChargePOI == 1) outputSPv2 += "PlusPlus";
      else if(gChargePOI == -1) outputSPv2 += "MinusMinus";

      taskSPv2[iCentralityBin] = new AliAnalysisTaskScalarProduct(Form("TaskScalarProduct_%s",outputSlotNameSPv2[iCentralityBin].Data()),kFALSE);
      taskSPv2[iCentralityBin]->SetHarmonic(2);
      taskSPv2[iCentralityBin]->SelectCollisionCandidates(triggerSelectionString);
      taskSPv2[iCentralityBin]->SetRelDiffMsub(1.0);
      taskSPv2[iCentralityBin]->SetTotalQvector(qVector);
      taskSPv2[iCentralityBin]->SetApplyCorrectionForNUA(kTRUE);
      
      coutputSPv2[iCentralityBin] = mgr->CreateContainer(Form("%s",outputSlotNameSPv2[iCentralityBin].Data()),TList::Class(),AliAnalysisManager::kOutputContainer,outputSPv2);
      mgr->AddTask(taskSPv2[iCentralityBin]);
      mgr->ConnectInput(taskSPv2[iCentralityBin],0,flowEvent[iCentralityBin]);
      mgr->ConnectInput(taskSPv2[iCentralityBin],0,coutputFE[iCentralityBin]);
      mgr->ConnectOutput(taskSPv2[iCentralityBin],1,coutputSPv2[iCentralityBin]);
      
      TString outputSPv4 = fileName;
      outputSPv4 += ":outputSPv4analysis";
      if(gChargePOI == 1) outputSPv4 += "PlusPlus";
      else if(gChargePOI == -1) outputSPv4 += "MinusMinus";
      
      taskSPv4[iCentralityBin] = new AliAnalysisTaskScalarProduct(Form("TaskScalarProduct_%s",outputSlotNameSPv4[iCentralityBin].Data()),kFALSE);
      taskSPv4[iCentralityBin]->SetHarmonic(4);
      taskSPv4[iCentralityBin]->SelectCollisionCandidates(triggerSelectionString);
      taskSPv4[iCentralityBin]->SetRelDiffMsub(1.0);
      taskSPv4[iCentralityBin]->SetTotalQvector(qVector);
      taskSPv4[iCentralityBin]->SetApplyCorrectionForNUA(kTRUE);
      
      coutputSPv4[iCentralityBin] = mgr->CreateContainer(Form("%s",outputSlotNameSPv4[iCentralityBin].Data()),TList::Class(),AliAnalysisManager::kOutputContainer,outputSPv4);
      mgr->AddTask(taskSPv4[iCentralityBin]);
      mgr->ConnectInput(taskSPv4[iCentralityBin],0,flowEvent[iCentralityBin]);
      mgr->ConnectInput(taskSPv4[iCentralityBin],0,coutputFE[iCentralityBin]);
      mgr->ConnectOutput(taskSPv4[iCentralityBin],1,coutputSPv4[iCentralityBin]);
    }
    //======================================================================//

    //======================================================================//
    if(gRunQC) {
      TString outputQCv2 = fileName;
      outputQCv2 += ":outputQCv2analysis";
      if(gChargePOI == 1) outputQCv2 += "PlusPlus";
      else if(gChargePOI == -1) outputQCv2 += "MinusMinus";
      
      taskQCv2[iCentralityBin] = new AliAnalysisTaskQCumulants(Form("TaskQCumulant_%s",outputSlotNameQCv2[iCentralityBin].Data()),kFALSE);
      taskQCv2[iCentralityBin]->SetHarmonic(2);
      taskQCv2[iCentralityBin]->SelectCollisionCandidates(triggerSelectionString);
      taskQCv2[iCentralityBin]->SetUsePhiWeights(kFALSE); 
      taskQCv2[iCentralityBin]->SetUsePtWeights(kFALSE);
      taskQCv2[iCentralityBin]->SetUseEtaWeights(kFALSE); 
      taskQCv2[iCentralityBin]->SetCalculateCumulantsVsM(kFALSE);
      taskQCv2[iCentralityBin]->SetnBinsMult(10000);
      taskQCv2[iCentralityBin]->SetMinMult(0.);
      taskQCv2[iCentralityBin]->SetMaxMult(10000.);
      taskQCv2[iCentralityBin]->SetApplyCorrectionForNUA(kTRUE);
      taskQCv2[iCentralityBin]->SetFillMultipleControlHistograms(kFALSE);     
      
      coutputQCv2[iCentralityBin] = mgr->CreateContainer(Form("%s",outputSlotNameQCv2[iCentralityBin].Data()),TList::Class(),AliAnalysisManager::kOutputContainer,outputQCv2);
      mgr->AddTask(taskQCv2[iCentralityBin]);
      mgr->ConnectInput(taskQCv2[iCentralityBin],0,flowEvent[iCentralityBin]);
      mgr->ConnectInput(taskQCv2[iCentralityBin],0,coutputFE[iCentralityBin]);
      mgr->ConnectOutput(taskQCv2[iCentralityBin],1,coutputQCv2[iCentralityBin]);
      
      TString outputQCv4 = fileName;
      outputQCv4 += ":outputQCv4analysis";
      if(gChargePOI == 1) outputQCv4 += "PlusPlus";
      else if(gChargePOI == -1) outputQCv4 += "MinusMinus";
     
      taskQCv4[iCentralityBin] = new AliAnalysisTaskQCumulants(Form("TaskQCumulant_%s",outputSlotNameQCv4[iCentralityBin].Data()),kFALSE);
      taskQCv4[iCentralityBin]->SetHarmonic(4);
      taskQCv4[iCentralityBin]->SelectCollisionCandidates(triggerSelectionString);
      taskQCv4[iCentralityBin]->SetUsePhiWeights(kFALSE); 
      taskQCv4[iCentralityBin]->SetUsePtWeights(kFALSE);
      taskQCv4[iCentralityBin]->SetUseEtaWeights(kFALSE); 
      taskQCv4[iCentralityBin]->SetCalculateCumulantsVsM(kFALSE);
      taskQCv4[iCentralityBin]->SetnBinsMult(10000);
      taskQCv4[iCentralityBin]->SetMinMult(0.);
      taskQCv4[iCentralityBin]->SetMaxMult(10000.);
      taskQCv4[iCentralityBin]->SetApplyCorrectionForNUA(kTRUE);
      taskQCv4[iCentralityBin]->SetFillMultipleControlHistograms(kFALSE);     
      
      coutputQCv4[iCentralityBin] = mgr->CreateContainer(Form("%s",outputSlotNameQCv4[iCentralityBin].Data()),TList::Class(),AliAnalysisManager::kOutputContainer,outputQCv4);
      mgr->AddTask(taskQCv4[iCentralityBin]);
      mgr->ConnectInput(taskQCv4[iCentralityBin],0,flowEvent[iCentralityBin]);
      mgr->ConnectInput(taskQCv4[iCentralityBin],0,coutputFE[iCentralityBin]);
      mgr->ConnectOutput(taskQCv4[iCentralityBin],1,coutputQCv4[iCentralityBin]);
    }
    //======================================================================//

    //======================================================================//
    if(gRunMHLS) {
      TString outputMHLS = fileName;
      outputMHLS += ":outputMHLSanalysis";
      if(gChargePOI == 1) outputMHLS += "PlusPlus";
      else if(gChargePOI == -1) outputMHLS += "MinusMinus";
      
      taskMHLS[iCentralityBin] = new AliAnalysisTaskMixedHarmonics(Form("TaskMixedHarmonicsLS_%s",outputSlotNameMHLS[iCentralityBin].Data()),kFALSE);
      taskMHLS[iCentralityBin]->SetHarmonic(1); // n in cos[n(phi1+phi2-2phi3)] and cos[n(psi1+psi2-2phi3)]
      taskMHLS[iCentralityBin]->SetNoOfMultipicityBins(10000);
      taskMHLS[iCentralityBin]->SetMultipicityBinWidth(1.);
      taskMHLS[iCentralityBin]->SetMinMultiplicity(1.);
      taskMHLS[iCentralityBin]->SetCorrectForDetectorEffects(kTRUE);
      taskMHLS[iCentralityBin]->SetEvaluateDifferential3pCorrelator(kTRUE); // evaluate <<cos[n(psi1+psi2-2phi3)]>> (Remark: two nested loops)
      taskMHLS[iCentralityBin]->SetOppositeChargesPOI(kFALSE); //
      
      coutputMHLS[iCentralityBin] = mgr->CreateContainer(Form("%s",outputSlotNameMHLS[iCentralityBin].Data()),TList::Class(),AliAnalysisManager::kOutputContainer,outputMHLS);
      mgr->AddTask(taskMHLS[iCentralityBin]);
      mgr->ConnectInput(taskMHLS[iCentralityBin],0,flowEvent[iCentralityBin]);
      mgr->ConnectInput(taskMHLS[iCentralityBin],0,coutputFE[iCentralityBin]);
      mgr->ConnectOutput(taskMHLS[iCentralityBin],1,coutputMHLS[iCentralityBin]);
      
      TString outputMHLS2 = fileName;
      outputMHLS2 += ":outputMHLS2analysis";
      if(gChargePOI == 1) outputMHLS2 += "PlusPlus";
      else if(gChargePOI == -1) outputMHLS2 += "MinusMinus";
      
      taskMHLS2[iCentralityBin] = new AliAnalysisTaskMixedHarmonics(Form("TaskMixedHarmonicsLS_%s",outputSlotNameMHLS2[iCentralityBin].Data()),kFALSE);
      taskMHLS2[iCentralityBin]->SetHarmonic(2); // n in cos[n(phi1+phi2-2phi3)] and cos[n(psi1+psi2-2phi3)]
      taskMHLS2[iCentralityBin]->SetNoOfMultipicityBins(10000);
      taskMHLS2[iCentralityBin]->SetMultipicityBinWidth(1.);
      taskMHLS2[iCentralityBin]->SetMinMultiplicity(1.);
      taskMHLS2[iCentralityBin]->SetCorrectForDetectorEffects(kTRUE);
      taskMHLS2[iCentralityBin]->SetEvaluateDifferential3pCorrelator(kTRUE); // evaluate <<cos[n(psi1+psi2-2phi3)]>> (Remark: two nested loops)
      taskMHLS2[iCentralityBin]->SetOppositeChargesPOI(kFALSE); //
      
      coutputMHLS2[iCentralityBin] = mgr->CreateContainer(Form("%s",outputSlotNameMHLS2[iCentralityBin].Data()),TList::Class(),AliAnalysisManager::kOutputContainer,outputMHLS2);
      mgr->AddTask(taskMHLS2[iCentralityBin]);
      mgr->ConnectInput(taskMHLS2[iCentralityBin],0,flowEvent[iCentralityBin]);
      mgr->ConnectInput(taskMHLS2[iCentralityBin],0,coutputFE[iCentralityBin]);
      mgr->ConnectOutput(taskMHLS2[iCentralityBin],1,coutputMHLS2[iCentralityBin]);
    }
    //======================================================================//

    //======================================================================//
    if(gRunMHUS) {
      TString outputMHUS = fileName;
      outputMHUS += ":outputMHUSanalysis";
      
      taskMHUS[iCentralityBin] = new AliAnalysisTaskMixedHarmonics(Form("TaskMixedHarmonicsLS_%s",outputSlotNameMHUS[iCentralityBin].Data()),kFALSE);
      taskMHUS[iCentralityBin]->SetHarmonic(1); // n in cos[n(phi1+phi2-2phi3)] and cos[n(psi1+psi2-2phi3)]
      taskMHUS[iCentralityBin]->SetNoOfMultipicityBins(10000);
      taskMHUS[iCentralityBin]->SetMultipicityBinWidth(1.);
      taskMHUS[iCentralityBin]->SetMinMultiplicity(1.);
      taskMHUS[iCentralityBin]->SetCorrectForDetectorEffects(kTRUE);
      taskMHUS[iCentralityBin]->SetEvaluateDifferential3pCorrelator(kTRUE); // evaluate <<cos[n(psi1+psi2-2phi3)]>> (Remark: two nested loops)
      taskMHUS[iCentralityBin]->SetOppositeChargesPOI(kTRUE); //
      
      coutputMHUS[iCentralityBin] = mgr->CreateContainer(Form("%s",outputSlotNameMHUS[iCentralityBin].Data()),TList::Class(),AliAnalysisManager::kOutputContainer,outputMHUS);
      mgr->AddTask(taskMHUS[iCentralityBin]);
      mgr->ConnectInput(taskMHUS[iCentralityBin],0,flowEvent[iCentralityBin]);
      mgr->ConnectInput(taskMHUS[iCentralityBin],0,coutputFE[iCentralityBin]);
      mgr->ConnectOutput(taskMHUS[iCentralityBin],1,coutputMHUS[iCentralityBin]);
      
      TString outputMHUS2 = fileName;
      outputMHUS2 += ":outputMHUS2analysis";
      
      taskMHUS2[iCentralityBin] = new AliAnalysisTaskMixedHarmonics(Form("TaskMixedHarmonicsLS_%s",outputSlotNameMHUS2[iCentralityBin].Data()),kFALSE);
      taskMHUS2[iCentralityBin]->SetHarmonic(2); // n in cos[n(phi1+phi2-2phi3)] and cos[n(psi1+psi2-2phi3)]
      taskMHUS2[iCentralityBin]->SetNoOfMultipicityBins(10000);
      taskMHUS2[iCentralityBin]->SetMultipicityBinWidth(1.);
      taskMHUS2[iCentralityBin]->SetMinMultiplicity(1.);
      taskMHUS2[iCentralityBin]->SetCorrectForDetectorEffects(kTRUE);
      taskMHUS2[iCentralityBin]->SetEvaluateDifferential3pCorrelator(kTRUE); // evaluate <<cos[n(psi1+psi2-2phi3)]>> (Remark: two nested loops)
      taskMHUS2[iCentralityBin]->SetOppositeChargesPOI(kTRUE); //
      
      coutputMHUS2[iCentralityBin] = mgr->CreateContainer(Form("%s",outputSlotNameMHUS2[iCentralityBin].Data()),TList::Class(),AliAnalysisManager::kOutputContainer,outputMHUS2);
      mgr->AddTask(taskMHUS2[iCentralityBin]);
      mgr->ConnectInput(taskMHUS2[iCentralityBin],0,flowEvent[iCentralityBin]);
      mgr->ConnectInput(taskMHUS2[iCentralityBin],0,coutputFE[iCentralityBin]);
      mgr->ConnectOutput(taskMHUS2[iCentralityBin],1,coutputMHUS2[iCentralityBin]);
    }
    //======================================================================//
    
    //======================================================================//
    if (taskFE[iCentralityBin]->GetQAOn()) {
      outputQA[iCentralityBin] = fileName;
      outputQA[iCentralityBin] += ":QA";
      coutputFEQA[iCentralityBin] = mgr->CreateContainer(Form("QA_%s",suffixName[iCentralityBin].Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputQA[iCentralityBin]);
      mgr->ConnectOutput(taskFE[iCentralityBin],2,coutputFEQA[iCentralityBin]);
    }
    //======================================================================//
  }//loop over centralities
}

//_________________________________________________________//
AliFlowEventCuts *createFlowEventCutObject(Int_t gCentralityMin = -1,
                                           Int_t gCentralityMax = -1,
                                           Bool_t isPbPb = kTRUE,
                                           Double_t whichData = 2011,
                                           AliFlowEventCuts::refMultMethod gCentralityEstimator = AliFlowEventCuts::kVZERO,
                                           Bool_t doQA = kFALSE) {
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
	if(whichData==2015) 
	  cutsEvent->SetCentralityPercentileRange(gCentralityMin,gCentralityMax,kTRUE);
      }
      cutsEvent->SetCutTPCmultiplicityOutliersAOD(kTRUE);
    }
    //    else
    //cutsEvent->SetCheckPileup(checkPileup);
    
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
