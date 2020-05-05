void AddTaskESEFlowQC(TString particleSpecies = "",
                      Int_t uptoWhichHarmonics = 5,
                      Bool_t isMC = kFALSE,
		      AliFlowEventCuts::refMultMethod gMultiplicityEstimator = AliFlowEventCuts::kVZERO,
                      UInt_t triggerSelectionString = AliVEvent::kINT7,
		      Bool_t gCutOnSphericity = kFALSE,
		      Double_t gMinSphericity = 0.0,
		      Double_t gMaxSphericity = 1.0,
		      Int_t gAODfilterBit = 768,
                      Double_t gProbPID = 0.9,
		      Double_t gCharge = 0,
                      Bool_t doQA = kTRUE) {
  //======================================================================//
  // FLOW cut objects config  
  const int nHarmonics = uptoWhichHarmonics - 1;
  AliFlowEventCuts* cutsEvent;
  AliFlowTrackCuts* cutsRP;
  AliFlowTrackCuts* cutsPOI;
  TString suffixName;
  TString outputSlotName[nHarmonics];
  
  //Create the event cut object
  cutsEvent = createFlowEventCutObject(isMC,gMultiplicityEstimator,gCutOnSphericity,gMinSphericity,gMaxSphericity,doQA);
    
  //Create the RP cut object
  cutsRP = createFlowRPCutObject(isMC,gAODfilterBit,gCharge,doQA);
    
  //Create the POI cut object
  cutsPOI = createFlowPOICutObject(isMC,particleSpecies,gAODfilterBit,gProbPID,gCharge,doQA);
        
  suffixName += "QC";
  if(particleSpecies != "") {
    suffixName += "_";
    suffixName += particleSpecies.Data();
    suffixName += "_";
  }
  else {
    suffixName += "_AllCharged";
    suffixName = Form("%s_%i",suffixName.Data(),gAODfilterBit);
    if(gCutOnSphericity) 
      suffixName = Form("%s_%.1f_%.1f",suffixName.Data(),gMinSphericity,gMaxSphericity);
  }
        
  //=====================================================================
  for(int nHarmonic = 2; nHarmonic <= uptoWhichHarmonics; nHarmonic++)
    outputSlotName[nHarmonic - 2] = Form("%s_v%d",suffixName.Data(),nHarmonic);
    
  TString fileName = "AnalysisResults";
  fileName.Append(".root");
  
  //========================FLOWPACKAGE TASKS=========================//
  AliAnalysisDataContainer *cinput1;
  AliAnalysisDataContainer *coutputFE;
  AliAnalysisDataContainer *coutputFEQA;
  AliAnalysisTaskFlowEvent *taskFE;
  
  AliAnalysisDataContainer *flowEvent[nHarmonics];
  AliAnalysisTaskFilterFE *tskFilter[nHarmonics];
  
  AliAnalysisDataContainer *coutputQC[nHarmonics];
  AliAnalysisTaskQCumulants *taskQC[nHarmonics];
  
  TString outputQA;
  TString slot[nHarmonics];
  
  //Create the analysis manager
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
  
  taskFE = new AliAnalysisTaskFlowEvent(Form("TaskFlowEvent_%s",suffixName.Data()),"",doQA);
  taskFE->SelectCollisionCandidates(triggerSelectionString);
    
  mgr->AddTask(taskFE);
    
  // Pass cuts for RPs and POIs to the task:
  taskFE->SetCutsEvent(cutsEvent);
  taskFE->SetCutsRP(cutsRP);
  taskFE->SetCutsPOI(cutsPOI);
  if (cutsRP->GetParamType()==AliFlowTrackCuts::kVZERO) {
    taskFE->SetHistWeightvsPhiMin(0.);
    taskFE->SetHistWeightvsPhiMax(200.);
  }
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  cinput1 = mgr->GetCommonInputContainer();
  
  coutputFE = mgr->CreateContainer(Form("FlowEvent_%s",suffixName.Data()),AliFlowEventSimple::Class(),AliAnalysisManager::kExchangeContainer);
  
  mgr->ConnectInput(taskFE,0,cinput1);
  mgr->ConnectOutput(taskFE,1,coutputFE);
  //==========================================================
  
  for(int nHarmonic = 2; nHarmonic <= uptoWhichHarmonics; nHarmonic++){
    flowEvent[nHarmonic-2] = mgr->CreateContainer(Form("Filter_%s",outputSlotName[nHarmonic-2].Data()),AliFlowEventSimple::Class(),AliAnalysisManager::kExchangeContainer);
    tskFilter[nHarmonic-2] = new AliAnalysisTaskFilterFE(Form("TaskFilter_%s",outputSlotName[nHarmonic-2].Data()),cutsRP, NULL);
    
    mgr->AddTask(tskFilter[nHarmonic-2]);
    mgr->ConnectInput( tskFilter[nHarmonic-2],0,coutputFE);
    mgr->ConnectOutput(tskFilter[nHarmonic-2],1,flowEvent[nHarmonic-2]);
    
    //QC task
    taskQC[nHarmonic-2] = new AliAnalysisTaskQCumulants(Form("TaskQCumulants_%s",outputSlotName[nHarmonic-2].Data()),kFALSE);
    taskQC[nHarmonic-2]->SetHarmonic(nHarmonic);
    taskQC[nHarmonic-2]->SelectCollisionCandidates(triggerSelectionString);
    taskQC[nHarmonic-2]->SetApplyCorrectionForNUA(kTRUE);
    taskQC[nHarmonic-2]->SetMaxCommonResultsHistogram(8);
    taskQC[nHarmonic-2]->SetFillMultipleControlHistograms(kTRUE);
    taskQC[nHarmonic-2]->SetCalculateCumulantsVsM(kTRUE);
    taskQC[nHarmonic-2]->SetCalculateAllCorrelationsVsM(kTRUE);
    taskQC[nHarmonic-2]->SetFillProfilesVsMUsingWeights(kTRUE);
    taskQC[nHarmonic-2]->SetMinMult(0);
    taskQC[nHarmonic-2]->SetMaxMult(500);
    taskQC[nHarmonic-2]->SetnBinsMult(500);
    
    TString outputQC = fileName;
    outputQC += ":outputQCanalysis";
    
    coutputQC[nHarmonic-2] = mgr->CreateContainer(Form("%s",outputSlotName[nHarmonic-2].Data()),TList::Class(),AliAnalysisManager::kOutputContainer,outputQC);
    
    mgr->AddTask(taskQC[nHarmonic-2]);
    
    mgr->ConnectInput(taskQC[nHarmonic-2],0,flowEvent[nHarmonic-2]);
    mgr->ConnectInput(taskQC[nHarmonic-2],0,coutputFE);
    mgr->ConnectOutput(taskQC[nHarmonic-2],1,coutputQC[nHarmonic-2]);
  }//loop over the harmonics
  
  if (taskFE->GetQAOn()) {
    outputQA = fileName;
    outputQA += ":QA";
    coutputFEQA = mgr->CreateContainer(Form("QA_%s",suffixName.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputQA);
    mgr->ConnectOutput(taskFE,2,coutputFEQA);
  }
}

//_________________________________________________________//
AliFlowEventCuts *createFlowEventCutObject(Bool_t isMC, AliFlowEventCuts::refMultMethod gMultiplicityEstimator = AliFlowEventCuts::kVZERO, Bool_t gCutOnSphericity = kFALSE, Double_t gMinSphericity = 0.0, Double_t gMaxSphericity = 1.0, Bool_t doQA = kFALSE) {
  //Part of the code that creates the event cut objects
  Double_t gVertexZmin = -10., gVertexZmax = 10.;
  
  //Create the event cut objects
  AliFlowEventCuts *cutsEvent = new AliFlowEventCuts("eventcuts");
  cutsEvent->SetPrimaryVertexZrange(gVertexZmin,gVertexZmax);
  if(gCutOnSphericity)
    cutsEvent->SetSphericityCut(gMinSphericity,gMaxSphericity);

  if(!isMC) {
    cutsEvent->SetCentralityPercentileMethod(gMultiplicityEstimator);
    //cutsEvent->SetCutTPCmultiplicityOutliersAOD(kTRUE);
    cutsEvent->SetCheckPileup(kTRUE);
  }

  if(isMC) {
    AliFlowTrackCuts* cutsRefMult = new AliFlowTrackCuts("MCRefMult");
    cutsRefMult->SetParamType(AliFlowTrackCuts::kMC);
    cutsRefMult->SetParamMix(AliFlowTrackCuts::kTrackWithMCkine);
    cutsRefMult->SetMCisPrimary(kTRUE);
    cutsRefMult->SetPtRange(0.2,10.);
    cutsRefMult->SetEtaRange(-0.8,0.8);

    cutsRefMult->SetQA(kFALSE);
    cutsEvent->SetRefMultCuts(cutsRefMult);
  }

  cutsEvent->SetQA(doQA);
  
  //return the object
  return cutsEvent;
}

//_________________________________________________________//
AliFlowTrackCuts *createFlowRPCutObject(Bool_t isMC,
					Int_t gAODfilterBit = 768,
                                        Double_t gCharge = 0.,
                                        Bool_t doQA = kFALSE) {
  //Part of the code that creates the RP cut objects
  Double_t gEtaMin = -0.8, gEtaMax = 0.8;
  Int_t gMinTPCdedx = 10;
  
  //Create the event cut objects
  AliFlowTrackCuts *cutsRP = new AliFlowTrackCuts("rpCuts");
  
  if (gCharge!=0) cutsRP->SetCharge(gCharge);
    cutsRP->SetPtRange(0.2,10.);
  cutsRP->SetEtaRange(gEtaMin,gEtaMax);

  if(isMC) {
    cutsRP->SetParamType(AliFlowTrackCuts::kMC);
    cutsRP->SetParamMix(AliFlowTrackCuts::kTrackWithMCkine);
    cutsRP->SetMCisPrimary(kTRUE);
  }
  else {  
    cutsRP->SetMinNClustersTPC(70);
    cutsRP->SetMinimalTPCdedx(gMinTPCdedx);
    cutsRP->SetAODfilterBit(gAODfilterBit);
  }
  cutsRP->SetQA(doQA);
  
  //return the object
  return cutsRP;
}

//_________________________________________________________//
AliFlowTrackCuts *createFlowPOICutObject(TString particleSpecies = "Pion",
					 Int_t gAODfilterBit = 768,
                                         Double_t gProbPID = 0.9,
					 Double_t gCharge = 0.,
                                         Bool_t doQA = kFALSE) {
  //Part of the code that creates the POI cut objects
  Double_t gEtaMin = -0.8, gEtaMax = 0.8;
  Int_t gMinTPCdedx = 10;
  
  AliFlowTrackCuts *cutsPOI = new AliFlowTrackCuts("poiCuts");
  
  if (gCharge!=0) cutsPOI->SetCharge(gCharge);
  cutsPOI->SetPtRange(0.2,10.);
  cutsPOI->SetEtaRange(gEtaMin,gEtaMax);

  if(isMC) {
    cutsPOI->SetParamType(AliFlowTrackCuts::kMC);
    cutsPOI->SetParamMix(AliFlowTrackCuts::kTrackWithMCkine);
    cutsPOI->SetMCisPrimary(kTRUE);
  }
  else {
    cutsPOI->SetMinNClustersTPC(70);
    cutsPOI->SetMinimalTPCdedx(gMinTPCdedx);
    cutsPOI->SetAODfilterBit(gAODfilterBit);
    
    Bool_t isPID = kFALSE;
    if(particleSpecies != "")
      isPID = kTRUE;
    
    if(isPID) {
      AliFlowTrackCuts::PIDsource sourcePID=AliFlowTrackCuts::kTOFbayesian;
      AliPID::EParticleType particleType = AliPID::kPion;
      
      if(particleSpecies.Contains("Pion")) {
	particleType = AliPID::kPion;
      }
      else if(particleSpecies.Contains("Kaon")) {
	particleType = AliPID::kKaon;
      }
      else if(particleSpecies.Contains("Proton")) {
	particleType = AliPID::kProton;
      }
      
      cutsPOI->SetPID(particleType, sourcePID, gProbPID);
    }
  }

  cutsPOI->SetQA(doQA);
    
  //return the object
  return cutsPOI;
}
