void AddTaskCMEAnalysis(Bool_t isPbPb,
			Double_t whichData,
			Bool_t isData,
			Bool_t checkPileup,
			AliFlowEventCuts::refMultMethod gCentralityEstimator,
			TString qVector,
			Double_t gEtaGap,
			Bool_t isVZERO,
			Int_t gAODfilterBit,
			Double_t gDCAvtxXY,
			Double_t gDCAvtxZ,
			Bool_t gRunSP,
			Bool_t gRunQC,
			Bool_t gRunMHLS,
			Bool_t gRunMHUS,
			Bool_t doQA,
			Int_t harmonicMH,
			Int_t vnHarmonic,
			Bool_t kZDCStudy,
			TString sLabel);
AliFlowEventCuts *createFlowEventCutObject(Int_t gCentralityMin,
                                           Int_t gCentralityMax,
                                           Bool_t isPbPb,
                                           Double_t whichData,
                                           AliFlowEventCuts::refMultMethod gCentralityEstimator,
                                           Bool_t doQA, 
					   Bool_t bPileup);
AliFlowTrackCuts *createFlowRPCutObject(Int_t gCentralityMin,
                                        Int_t gCentralityMax,
                                        Bool_t isVZERO,
                                        Int_t gAODfilterBit,
					Double_t whichData,
                                        Double_t gChargeRP,
                                        Bool_t doQA);
AliFlowTrackCuts *createFlowPOICutObject(Int_t gCentralityMin,
                                         Int_t gCentralityMax,
                                         TString gQvector,
                                         Double_t gEtaGap,
                                         Bool_t isVZERO,
                                         Bool_t isPbPb,
                                         Int_t gAODfilterBit,
					 Double_t whichData,
					 Double_t gDCAvtxXY,
                                         Double_t gDCAvtxZ,
                                         Double_t gChargePOI,
                                         Bool_t doQA);

void AddTaskCMEAnalysis(Bool_t isPbPb = kTRUE,
			Int_t whichData = 2010,
			Bool_t isData = kTRUE,
			Bool_t checkPileup = kTRUE,
			AliFlowEventCuts::refMultMethod gCentralityEstimator = AliFlowEventCuts::kVZERO,
			TString qVector = "Qa",
			Double_t gEtaGap = 0.0,
			Bool_t isVZERO = kFALSE,
			Int_t gAODfilterBit = 768,
			Double_t gDCAvtxXY = -1.,
			Double_t gDCAvtxZ = -1.,
			Bool_t gRunSP = kFALSE,
			Bool_t gRunQC = kFALSE,
			Bool_t gRunMHLS = kTRUE,
			Bool_t gRunMHUS = kTRUE,
			Bool_t doQA = kFALSE,
			Int_t harmonicMH = 1,
			Int_t vnHarmonic = 2,
			Bool_t kZDCStudy = kTRUE,
                        Bool_t kUseGainEqualisationMap = kTRUE,
  			const char* lhcPeriod = "LHC10h",
                        const char* fileForZDCGainMap = "gainZDC.LHC10h.root",
  			TString sLabel = "FB768") {
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
    case 2018:
      triggerSelectionString = AliVEvent::kINT7 | AliVEvent::kCentral | AliVEvent::kSemiCentral;
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
    
  //=================ZDC calibration====================//
  TH3F *gHistZNChannel = 0x0;
  TH3F *gHistZPChannel = 0x0;
  if(kUseGainEqualisationMap) {
    //============Get the equilization map============//
    TFile *calibrationFile = TFile::Open(fileForZDCGainMap);
    if((!calibrationFile)||(!calibrationFile->IsOpen())) {
      cout<<"No calibration file found!!!"<<endl;
      return;
    }

    TList *list = dynamic_cast<TList *>(calibrationFile->Get(lhcPeriod));
    if(!list) {
      cout<<"Calibration TList not found!!!"<<endl;
      return;
    }

    //Get the objects from the list
    gHistZNChannel = dynamic_cast<TH3F *>(list->FindObject("gHistZNChannelGainEqualizationMap"));
    if(!gHistZNChannel) {
      cout<<"ZN channel calibration object not found!!!"<<endl;
      return;
    }

    TH3F *gHistZPChannel = dynamic_cast<TH3F *>(list->FindObject("gHistZPChannelGainEqualizationMap"));
    if(!gHistZPChannel) {
      cout<<"ZP channel calibration object not found!!!"<<endl;
      return;
    }
  }//Use gain equalisation map                   

  //======================================================================//
  // FLOW cut objects config
  AliFlowEventCuts* cutsEvent_wQA[nCentralities];//with QA on
  AliFlowTrackCuts* cutsRP_wQA[nCentralities];   //with QA on
  AliFlowEventCuts* cutsEvent[nCentralities];    //without QA
  AliFlowTrackCuts* cutsRP[nCentralities];       //without QA

  AliFlowTrackCuts* cutsPOI_PP[nCentralities];
  AliFlowTrackCuts* cutsPOI_NN[nCentralities];
  AliFlowTrackCuts* cutsPOI_PN[nCentralities];
  AliFlowTrackCuts* cutsPOI_All[nCentralities];
  AliFlowTrackCuts* cutsPOI_PNwQA[nCentralities]; //with QA on

  TString suffixNamePP[nCentralities];
  TString suffixNameNN[nCentralities];
  TString suffixNamePN[nCentralities];
  TString suffixNameAll[nCentralities];

  TString outputSlotAll_NameSPv2[nCentralities];
  TString outputSlotAll_NameQCv2[nCentralities];

  TString outputSlotPP_NameMHLS[nCentralities];
  TString outputSlotNN_NameMHLS[nCentralities];
  TString outputSlotPN_NameMHUS[nCentralities];
  
  for(Int_t iCentralityBin = 0; iCentralityBin < nCentralities - 1; iCentralityBin++) {
    //Create the event cut object
    cutsEvent_wQA[iCentralityBin] = createFlowEventCutObject(gCentrality[iCentralityBin],gCentrality[iCentralityBin+1],isPbPb,whichData,gCentralityEstimator,doQA,checkPileup);
    cutsEvent[iCentralityBin] = createFlowEventCutObject(gCentrality[iCentralityBin],gCentrality[iCentralityBin+1],isPbPb,whichData,gCentralityEstimator,kFALSE,checkPileup);
        
    //Create the RP cut object
    cutsRP_wQA[iCentralityBin] = createFlowRPCutObject(gCentrality[iCentralityBin],gCentrality[iCentralityBin+1],isVZERO,gAODfilterBit,whichData,0,doQA);
    cutsRP[iCentralityBin] = createFlowRPCutObject(gCentrality[iCentralityBin],gCentrality[iCentralityBin+1],isVZERO,gAODfilterBit,whichData,0,kFALSE);
        
    //Create the POI cut object for ++,-- and +-
    cutsPOI_PP[iCentralityBin] = createFlowPOICutObject(gCentrality[iCentralityBin],gCentrality[iCentralityBin+1],qVector,gEtaGap,isVZERO,isPbPb,gAODfilterBit,whichData,gDCAvtxXY,gDCAvtxZ,1,kFALSE);
    cutsPOI_NN[iCentralityBin] = createFlowPOICutObject(gCentrality[iCentralityBin],gCentrality[iCentralityBin+1],qVector,gEtaGap,isVZERO,isPbPb,gAODfilterBit,whichData,gDCAvtxXY,gDCAvtxZ,-1,kFALSE);
    cutsPOI_PN[iCentralityBin] = createFlowPOICutObject(gCentrality[iCentralityBin],gCentrality[iCentralityBin+1],qVector,gEtaGap,isVZERO,isPbPb,gAODfilterBit,whichData,gDCAvtxXY,gDCAvtxZ,0,kFALSE);
    cutsPOI_All[iCentralityBin] = createFlowPOICutObject(gCentrality[iCentralityBin],gCentrality[iCentralityBin+1],qVector,gEtaGap,isVZERO,isPbPb,gAODfilterBit,whichData,gDCAvtxXY,gDCAvtxZ,0,kFALSE);
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

    suffixNameAll[iCentralityBin] += "Centrality";
    suffixNameAll[iCentralityBin] += gCentrality[iCentralityBin];
    suffixNameAll[iCentralityBin] += "To";
    suffixNameAll[iCentralityBin] += gCentrality[iCentralityBin+1];
    suffixNameAll[iCentralityBin] += sLabel;

    //Giving names to output slots
    suffixNameAll[iCentralityBin] += "All";    
    outputSlotAll_NameSPv2[iCentralityBin] = Form("%s_SP_v2",suffixNamePP[iCentralityBin].Data());
    outputSlotAll_NameQCv2[iCentralityBin] = Form("%s_QC_v2",suffixNamePP[iCentralityBin].Data());

    suffixNamePP[iCentralityBin] += "PlusPlus";    
    outputSlotPP_NameMHLS[iCentralityBin] = Form("%s_MHLS_v1",suffixNamePP[iCentralityBin].Data());

    suffixNameNN[iCentralityBin] += "MinusMinus";    
    outputSlotNN_NameMHLS[iCentralityBin] = Form("%s_MHLS_v1",suffixNameNN[iCentralityBin].Data());

    suffixNamePN[iCentralityBin] += "PlusMinus";    
    outputSlotPN_NameMHUS[iCentralityBin]  = Form("%s_MHUS_v1",suffixNamePN[iCentralityBin].Data());
  }//loop over centrality classes

  TString fileName = "AnalysisResults";
  fileName.Append(".root");

  //========================FLOWPACKAGE TASKS=========================//
  //----PP----
  AliAnalysisDataContainer *cinput1_PP[nCentralities];
  AliAnalysisDataContainer *coutputFE_PP[nCentralities];
  AliAnalysisDataContainer *coutputFE_PPQA[nCentralities];
  AliAnalysisTaskFlowEvent *taskFE_PP[nCentralities];
  AliAnalysisDataContainer *flowEvent_PP[nCentralities];
  AliAnalysisTaskFilterFE  *taskFilter_PP[nCentralities];
  AliAnalysisDataContainer *coutputMHLS_PP[nCentralities];
  AliAnalysisTaskMixedHarmonics *taskMHLS_PP[nCentralities];
    
  //----NN----
  AliAnalysisDataContainer *cinput1_NN[nCentralities];
  AliAnalysisDataContainer *coutputFE_NN[nCentralities];
  AliAnalysisDataContainer *coutputFE_NNQA[nCentralities];
  AliAnalysisTaskFlowEvent *taskFE_NN[nCentralities];
  AliAnalysisDataContainer *flowEvent_NN[nCentralities];
  AliAnalysisTaskFilterFE *taskFilter_NN[nCentralities];
  AliAnalysisDataContainer *coutputMHLS_NN[nCentralities];
  AliAnalysisTaskMixedHarmonics *taskMHLS_NN[nCentralities];

  //----PN----
  AliAnalysisDataContainer *cinput1_PN[nCentralities];
  AliAnalysisDataContainer *coutputFE_PN[nCentralities];
  AliAnalysisDataContainer *coutputFE_PNQA[nCentralities];
  AliAnalysisTaskFlowEvent *taskFE_PN[nCentralities];
  AliAnalysisDataContainer *flowEvent_PN[nCentralities];
  AliAnalysisTaskFilterFE *taskFilter_PN[nCentralities];
  AliAnalysisDataContainer *coutputMHUS_PN[nCentralities];
  AliAnalysisTaskMixedHarmonics *taskMHUS_PN[nCentralities];

  //----All----
  AliAnalysisDataContainer *cinput1_All[nCentralities];
  AliAnalysisDataContainer *coutputFE_All[nCentralities];
  AliAnalysisDataContainer *coutputFE_AllQA[nCentralities];
  AliAnalysisTaskFlowEvent *taskFE_All[nCentralities];
  AliAnalysisDataContainer *flowEvent_All[nCentralities];
  AliAnalysisTaskFilterFE *taskFilter_All[nCentralities];
  AliAnalysisDataContainer *coutputSPv2[nCentralities];
  AliAnalysisTaskScalarProduct *taskSPv2[nCentralities];
  AliAnalysisDataContainer *coutputQCv2[nCentralities];
  AliAnalysisTaskQCumulants *taskQCv2[nCentralities];
  
  TString outputQA_PP[nCentralities];
  TString outputQA_NN[nCentralities];
  TString outputQA_PN[nCentralities];
  TString outputQA_All[nCentralities];

  TString myNameSPv2[nCentralities];
  TString myNameQCv2[nCentralities];
  TString myNameMHLS[nCentralities];
  TString myNameMHUS[nCentralities];
  TString        slot[nCentralities];

  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  for(Int_t iCentralityBin = 0; iCentralityBin < nCentralities - 1; iCentralityBin++) {
    //======================== Initiate v2 Tasks ==========================
    //------- ALL --------------
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
    
    taskFE_All[iCentralityBin] = new AliAnalysisTaskFlowEvent(Form("TaskFlowEvent_%s",suffixNameAll[iCentralityBin].Data()),"",kFALSE);
    taskFE_All[iCentralityBin]->SelectCollisionCandidates(triggerSelectionString);

    //Sub events
    Double_t minA = -0.8;//
    Double_t maxA = -0.5*gEtaGap;//
    Double_t minB = +0.5*gEtaGap;//
    Double_t maxB = 0.8;//
    
    if(!isVZERO)
      taskFE_All[iCentralityBin]->SetSubeventEtaRange(minA, maxA, minB, maxB);
    if(isVZERO)
      taskFE_All[iCentralityBin]->SetSubeventEtaRange(-5,-1.5,+1.5,5);
    mgr->AddTask(taskFE_All[iCentralityBin]);
    
    //Add the gain map top the task                                               
    if(kUseGainEqualisationMap) {
      if(gHistZNChannel)
	taskFE_All[iCentralityBin]->SetZNChannelGainEqualizationMap(gHistZNChannel);
      if(gHistZPChannel)
	taskFE_All[iCentralityBin]->SetZPChannelGainEqualizationMap(gHistZPChannel);    
      cout<<"(Panos All)"<<endl;
    }

    // Pass cuts for RPs and POIs to the task:
    taskFE_All[iCentralityBin]->SetCutsEvent(cutsEvent[iCentralityBin]);
    taskFE_All[iCentralityBin]->SetCutsRP(cutsRP[iCentralityBin]);
    taskFE_All[iCentralityBin]->SetCutsPOI(cutsPOI_All[iCentralityBin]);

    if (cutsRP[iCentralityBin]->GetParamType()==AliFlowTrackCuts::kVZERO) {
      taskFE_All[iCentralityBin]->SetHistWeightvsPhiMin(0.);
      taskFE_All[iCentralityBin]->SetHistWeightvsPhiMax(200.);
    }
    
    if(!mgr){
      mgr = AliAnalysisManager::GetAnalysisManager();
    }
    cinput1_All[iCentralityBin] = mgr->GetCommonInputContainer();
    
    coutputFE_All[iCentralityBin] = mgr->CreateContainer(Form("FlowEvent_%s",suffixNameAll[iCentralityBin].Data()),AliFlowEventSimple::Class(),AliAnalysisManager::kExchangeContainer);
    
    mgr->ConnectInput(taskFE_All[iCentralityBin],0,cinput1_All[iCentralityBin]);
    mgr->ConnectOutput(taskFE_All[iCentralityBin],1,coutputFE_All[iCentralityBin]);
    //==========================================================
    
    //======================================================================//
    flowEvent_All[iCentralityBin] = mgr->CreateContainer(Form("Filter_%s",outputSlotAll_NameSPv2[iCentralityBin].Data()),AliFlowEventSimple::Class(),AliAnalysisManager::kExchangeContainer);
    taskFilter_All[iCentralityBin] = new AliAnalysisTaskFilterFE(Form("TaskFilter_%s",outputSlotAll_NameSPv2[iCentralityBin].Data()),cutsRP[iCentralityBin], NULL);
    
    if(!isVZERO) 
      taskFilter_All[iCentralityBin]->SetSubeventEtaRange(minA, maxA, minB, maxB);
    if(isVZERO) 
      taskFilter_All[iCentralityBin]->SetSubeventEtaRange(-5,-1.5,+1.5,5);
    mgr->AddTask(taskFilter_All[iCentralityBin]);
    mgr->ConnectInput(taskFilter_All[iCentralityBin],0,coutputFE_All[iCentralityBin]);
    mgr->ConnectOutput(taskFilter_All[iCentralityBin],1,flowEvent_All[iCentralityBin]);
    //======================================================================//
  
    //======================================================================//
    if(gRunSP) {
      TString outputSPv2 = fileName;
      outputSPv2 += ":outputSPv2analysis";

      taskSPv2[iCentralityBin] = new AliAnalysisTaskScalarProduct(Form("TaskScalarProduct_%s",outputSlotAll_NameSPv2[iCentralityBin].Data()),kFALSE);
      taskSPv2[iCentralityBin]->SetHarmonic(vnHarmonic);
      taskSPv2[iCentralityBin]->SelectCollisionCandidates(triggerSelectionString);
      taskSPv2[iCentralityBin]->SetRelDiffMsub(1.0);
      taskSPv2[iCentralityBin]->SetTotalQvector(qVector);
      taskSPv2[iCentralityBin]->SetApplyCorrectionForNUA(kTRUE);
      
      coutputSPv2[iCentralityBin] = mgr->CreateContainer(Form("%s",outputSlotAll_NameSPv2[iCentralityBin].Data()),TList::Class(),AliAnalysisManager::kOutputContainer,outputSPv2);
      mgr->AddTask(taskSPv2[iCentralityBin]);
      mgr->ConnectInput(taskSPv2[iCentralityBin],0,flowEvent_All[iCentralityBin]);
      mgr->ConnectInput(taskSPv2[iCentralityBin],0,coutputFE_All[iCentralityBin]);
      mgr->ConnectOutput(taskSPv2[iCentralityBin],1,coutputSPv2[iCentralityBin]);
    }
    //======================================================================//
       
    //======================================================================//
    if(gRunQC) {
      TString outputQCv2 = fileName;
      outputQCv2 += ":outputQCv2analysis";
       
      taskQCv2[iCentralityBin] = new AliAnalysisTaskQCumulants(Form("TaskQCumulant_%s",outputSlotAll_NameQCv2[iCentralityBin].Data()),kFALSE);
      taskQCv2[iCentralityBin]->SetHarmonic(vnHarmonic);
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
      
      coutputQCv2[iCentralityBin] = mgr->CreateContainer(Form("%s",outputSlotAll_NameQCv2[iCentralityBin].Data()),TList::Class(),AliAnalysisManager::kOutputContainer,outputQCv2);
      mgr->AddTask(taskQCv2[iCentralityBin]);
      mgr->ConnectInput(taskQCv2[iCentralityBin],0,flowEvent_All[iCentralityBin]);
      mgr->ConnectInput(taskQCv2[iCentralityBin],0,coutputFE_All[iCentralityBin]);
      mgr->ConnectOutput(taskQCv2[iCentralityBin],1,coutputQCv2[iCentralityBin]);
    }
    //======================================================================//

    //======================================================================//
    //------- PP --------------
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
    //Add the gain map top the task                                               
    if(kUseGainEqualisationMap) {
      if(gHistZNChannel)
	taskFE_PP[iCentralityBin]->SetZNChannelGainEqualizationMap(gHistZNChannel);
      if(gHistZPChannel)
	taskFE_PP[iCentralityBin]->SetZPChannelGainEqualizationMap(gHistZPChannel);    
      cout<<"(Panos PP)"<<endl;
    }

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
      mgr = AliAnalysisManager::GetAnalysisManager();
    }
    cinput1_PP[iCentralityBin] = mgr->GetCommonInputContainer();
    
    coutputFE_PP[iCentralityBin] = mgr->CreateContainer(Form("FlowEvent_%s",suffixNamePP[iCentralityBin].Data()),AliFlowEventSimple::Class(),AliAnalysisManager::kExchangeContainer);
    
    mgr->ConnectInput(taskFE_PP[iCentralityBin],0,cinput1_PP[iCentralityBin]);
    mgr->ConnectOutput(taskFE_PP[iCentralityBin],1,coutputFE_PP[iCentralityBin]);
    //==========================================================
    
    //======================================================================//
    flowEvent_PP[iCentralityBin] = mgr->CreateContainer(Form("Filter_%s",outputSlotPP_NameMHLS[iCentralityBin].Data()),AliFlowEventSimple::Class(),AliAnalysisManager::kExchangeContainer);
    taskFilter_PP[iCentralityBin] = new AliAnalysisTaskFilterFE(Form("TaskFilter_%s",outputSlotPP_NameMHLS[iCentralityBin].Data()),cutsRP[iCentralityBin], NULL);
    
    if(!isVZERO) 
      taskFilter_PP[iCentralityBin]->SetSubeventEtaRange(minA, maxA, minB, maxB);
    if(isVZERO) 
      taskFilter_PP[iCentralityBin]->SetSubeventEtaRange(-5,-1.5,+1.5,5);
    mgr->AddTask(taskFilter_PP[iCentralityBin]);
    mgr->ConnectInput(taskFilter_PP[iCentralityBin],0,coutputFE_PP[iCentralityBin]);
    mgr->ConnectOutput(taskFilter_PP[iCentralityBin],1,flowEvent_PP[iCentralityBin]);
    //======================================================================//
  
    //======================================================================//
    if(gRunMHLS) {
      TString outputMHLS_PP = fileName;
      outputMHLS_PP += ":outputMHLSanalysis";
      outputMHLS_PP += "PlusPlus";
      taskMHLS_PP[iCentralityBin] = new AliAnalysisTaskMixedHarmonics(Form("TaskMixedHarmonicsLS_%s",outputSlotPP_NameMHLS[iCentralityBin].Data()),kFALSE);
      taskMHLS_PP[iCentralityBin]->SelectCollisionCandidates(triggerSelectionString);
      taskMHLS_PP[iCentralityBin]->SetHarmonic(harmonicMH);// n in cos[n(phi1+phi2-2phi3)] and cos[n(psi1+psi2-2phi3)]                   
      taskMHLS_PP[iCentralityBin]->SetNoOfMultipicityBins(10000);
      taskMHLS_PP[iCentralityBin]->SetMultipicityBinWidth(1.);
      taskMHLS_PP[iCentralityBin]->SetMinMultiplicity(1.);
      taskMHLS_PP[iCentralityBin]->SetCorrectForDetectorEffects(kTRUE);
      taskMHLS_PP[iCentralityBin]->SetEvaluateDifferential3pCorrelator(kTRUE); // evaluate <<cos[n(psi1+psi2-2phi3)]>> (Remark: two nested loops)         
      taskMHLS_PP[iCentralityBin]->SetOppositeChargesPOI(kFALSE); //                                                                                      
      if(kZDCStudy) taskMHLS_PP[iCentralityBin]->SetCalculateVsZDC(kTRUE);
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
    //======================================================================//

    //======================================================================//
    //------- NN --------------
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
    
    taskFE_NN[iCentralityBin] = new AliAnalysisTaskFlowEvent(Form("TaskFlowEvent_%s",suffixNameNN[iCentralityBin].Data()),"",kFALSE);
    taskFE_NN[iCentralityBin]->SelectCollisionCandidates(triggerSelectionString);
    //Add the gain map top the task                                               
    if(kUseGainEqualisationMap) {
      if(gHistZNChannel)
	taskFE_NN[iCentralityBin]->SetZNChannelGainEqualizationMap(gHistZNChannel);
      if(gHistZPChannel)
	taskFE_NN[iCentralityBin]->SetZPChannelGainEqualizationMap(gHistZPChannel);    
      cout<<"(Panos NN)"<<endl;
    }

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
      mgr = AliAnalysisManager::GetAnalysisManager();
    }
    cinput1_NN[iCentralityBin] = mgr->GetCommonInputContainer();
    
    coutputFE_NN[iCentralityBin] = mgr->CreateContainer(Form("FlowEvent_%s",suffixNameNN[iCentralityBin].Data()),AliFlowEventSimple::Class(),AliAnalysisManager::kExchangeContainer);
    
    mgr->ConnectInput(taskFE_NN[iCentralityBin],0,cinput1_NN[iCentralityBin]);
    mgr->ConnectOutput(taskFE_NN[iCentralityBin],1,coutputFE_NN[iCentralityBin]);
    //==========================================================
    
    //======================================================================//
    flowEvent_NN[iCentralityBin] = mgr->CreateContainer(Form("Filter_%s",outputSlotNN_NameMHLS[iCentralityBin].Data()),AliFlowEventSimple::Class(),AliAnalysisManager::kExchangeContainer);
    taskFilter_NN[iCentralityBin] = new AliAnalysisTaskFilterFE(Form("TaskFilter_%s",outputSlotNN_NameMHLS[iCentralityBin].Data()),cutsRP[iCentralityBin], NULL);
    
    if(!isVZERO) 
      taskFilter_NN[iCentralityBin]->SetSubeventEtaRange(minA, maxA, minB, maxB);
    if(isVZERO) 
      taskFilter_NN[iCentralityBin]->SetSubeventEtaRange(-5,-1.5,+1.5,5);
    mgr->AddTask(taskFilter_NN[iCentralityBin]);
    mgr->ConnectInput(taskFilter_NN[iCentralityBin],0,coutputFE_NN[iCentralityBin]);
    mgr->ConnectOutput(taskFilter_NN[iCentralityBin],1,flowEvent_NN[iCentralityBin]);
    //======================================================================//
  
    //======================================================================//
    if(gRunMHLS) {
      TString outputMHLS_NN = fileName;
      outputMHLS_NN += ":outputMHLSanalysis";
      outputMHLS_NN += "MinusMinus";
      taskMHLS_NN[iCentralityBin] = new AliAnalysisTaskMixedHarmonics(Form("TaskMixedHarmonicsLS_%s",outputSlotNN_NameMHLS[iCentralityBin].Data()),kFALSE);
      taskMHLS_NN[iCentralityBin]->SelectCollisionCandidates(triggerSelectionString);
      taskMHLS_NN[iCentralityBin]->SetHarmonic(harmonicMH); // n in cos[n(phi1+phi2-2phi3)] and cos[n(psi1+psi2-2phi3)]                     
      taskMHLS_NN[iCentralityBin]->SetNoOfMultipicityBins(10000);
      taskMHLS_NN[iCentralityBin]->SetMultipicityBinWidth(1.);
      taskMHLS_NN[iCentralityBin]->SetMinMultiplicity(1.);
      taskMHLS_NN[iCentralityBin]->SetCorrectForDetectorEffects(kTRUE);
      taskMHLS_NN[iCentralityBin]->SetEvaluateDifferential3pCorrelator(kTRUE); // evaluate <<cos[n(psi1+psi2-2phi3)]>> (Remark: two nested loops)         
      taskMHLS_NN[iCentralityBin]->SetOppositeChargesPOI(kFALSE); //                                                                                      
      if(kZDCStudy) taskMHLS_NN[iCentralityBin]->SetCalculateVsZDC(kTRUE);
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
    //======================================================================//

    //======================================================================//
    //------- PN --------------
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
    
    taskFE_PN[iCentralityBin] = new AliAnalysisTaskFlowEvent(Form("TaskFlowEvent_%s",suffixNamePN[iCentralityBin].Data()),"",kFALSE);
    taskFE_PN[iCentralityBin]->SelectCollisionCandidates(triggerSelectionString);
    //Add the gain map top the task                                               
    if(kUseGainEqualisationMap) {
      if(gHistZNChannel)
	taskFE_PN[iCentralityBin]->SetZNChannelGainEqualizationMap(gHistZNChannel);
      if(gHistZPChannel)
	taskFE_PN[iCentralityBin]->SetZPChannelGainEqualizationMap(gHistZPChannel);    
      cout<<"(Panos PN)"<<endl;
    }

    if(!isVZERO)
      taskFE_PN[iCentralityBin]->SetSubeventEtaRange(minA, maxA, minB, maxB);
    if(isVZERO)
      taskFE_PN[iCentralityBin]->SetSubeventEtaRange(-5,-1.5,+1.5,5);
    mgr->AddTask(taskFE_PN[iCentralityBin]);
    
    // Pass cuts for RPs and POIs to the task:
    taskFE_PN[iCentralityBin]->SetCutsEvent(cutsEvent[iCentralityBin]);
    taskFE_PN[iCentralityBin]->SetCutsRP(cutsRP[iCentralityBin]);
    taskFE_PN[iCentralityBin]->SetCutsPOI(cutsPOI_PN[iCentralityBin]);

    if (cutsRP[iCentralityBin]->GetParamType()==AliFlowTrackCuts::kVZERO) {
      taskFE_PN[iCentralityBin]->SetHistWeightvsPhiMin(0.);
      taskFE_PN[iCentralityBin]->SetHistWeightvsPhiMax(200.);
    }
    
    if(!mgr){
      mgr = AliAnalysisManager::GetAnalysisManager();
    }
    cinput1_PN[iCentralityBin] = mgr->GetCommonInputContainer();
    
    coutputFE_PN[iCentralityBin] = mgr->CreateContainer(Form("FlowEvent_%s",suffixNamePN[iCentralityBin].Data()),AliFlowEventSimple::Class(),AliAnalysisManager::kExchangeContainer);
    
    mgr->ConnectInput(taskFE_PN[iCentralityBin],0,cinput1_PN[iCentralityBin]);
    mgr->ConnectOutput(taskFE_PN[iCentralityBin],1,coutputFE_PN[iCentralityBin]);
    //==========================================================
    
    //======================================================================//
    flowEvent_PN[iCentralityBin] = mgr->CreateContainer(Form("Filter_%s",outputSlotPN_NameMHUS[iCentralityBin].Data()),AliFlowEventSimple::Class(),AliAnalysisManager::kExchangeContainer);
    taskFilter_PN[iCentralityBin] = new AliAnalysisTaskFilterFE(Form("TaskFilter_%s",outputSlotPN_NameMHUS[iCentralityBin].Data()),cutsRP[iCentralityBin], NULL);
    
    if(!isVZERO) 
      taskFilter_PN[iCentralityBin]->SetSubeventEtaRange(minA, maxA, minB, maxB);
    if(isVZERO) 
      taskFilter_PN[iCentralityBin]->SetSubeventEtaRange(-5,-1.5,+1.5,5);
    mgr->AddTask(taskFilter_PN[iCentralityBin]);
    mgr->ConnectInput(taskFilter_PN[iCentralityBin],0,coutputFE_PN[iCentralityBin]);
    mgr->ConnectOutput(taskFilter_PN[iCentralityBin],1,flowEvent_PN[iCentralityBin]);
    //======================================================================//
  
    //======================================================================//
    if(gRunMHUS) {
      TString outputMHUS_PN = fileName;
      outputMHUS_PN += ":outputMHUSanalysis";
      outputMHUS_PN += "PlusMinus";
      taskMHUS_PN[iCentralityBin] = new AliAnalysisTaskMixedHarmonics(Form("TaskMixedHarmonicsUS_%s",outputSlotPN_NameMHUS[iCentralityBin].Data()),kFALSE);
      taskMHUS_PN[iCentralityBin]->SelectCollisionCandidates(triggerSelectionString);
      taskMHUS_PN[iCentralityBin]->SetHarmonic(harmonicMH); // n in cos[n(phi1+phi2-2phi3)] and cos[n(psi1+psi2-2phi3)]                     
      taskMHUS_PN[iCentralityBin]->SetNoOfMultipicityBins(10000);
      taskMHUS_PN[iCentralityBin]->SetMultipicityBinWidth(1.);
      taskMHUS_PN[iCentralityBin]->SetMinMultiplicity(1.);
      taskMHUS_PN[iCentralityBin]->SetCorrectForDetectorEffects(kTRUE);
      taskMHUS_PN[iCentralityBin]->SetEvaluateDifferential3pCorrelator(kTRUE); // evaluate <<cos[n(psi1+psi2-2phi3)]>> (Remark: two nested loops)         
      taskMHUS_PN[iCentralityBin]->SetOppositeChargesPOI(kTRUE); //                                                                                      
      if(kZDCStudy) taskMHUS_PN[iCentralityBin]->SetCalculateVsZDC(kTRUE);
      if(isPbPb){
	taskMHUS_PN[iCentralityBin]->SetRejectPileUp(checkPileup);
	taskMHUS_PN[iCentralityBin]->SetRejectPileUpTight(checkPileup);
      }
      else{
	taskMHUS_PN[iCentralityBin]->SetRejectPileUp(kFALSE);
	taskMHUS_PN[iCentralityBin]->SetRejectPileUpTight(kFALSE);
      }
      
      coutputMHUS_PN[iCentralityBin] = mgr->CreateContainer(Form("%s",outputSlotPN_NameMHUS[iCentralityBin].Data()),TList::Class(),AliAnalysisManager::kOutputContainer,outputMHUS_PN);
      mgr->AddTask(taskMHUS_PN[iCentralityBin]);
      mgr->ConnectInput(taskMHUS_PN[iCentralityBin],0,flowEvent_PN[iCentralityBin]);
      mgr->ConnectInput(taskMHUS_PN[iCentralityBin],0,coutputFE_PN[iCentralityBin]);
      mgr->ConnectOutput(taskMHUS_PN[iCentralityBin],1,coutputMHUS_PN[iCentralityBin]);
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
                                           Bool_t doQA = kFALSE, 
					   Bool_t bPileup = kFALSE) 
{
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
	if((whichData==2015)||(whichData==2018)) { 
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
        if(whichData == 2010)
            cutsRP = AliFlowTrackCuts::GetStandardVZEROOnlyTrackCuts2010();
        if(whichData == 2011)
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
