void AddTaskPIDFlowSP(TString particleSpecies = "Pion",
		      TString qVector = "Qa",
		      Int_t uptoWhichHarmonics = 5,
		      Double_t gEtaGap = 0.0,
		      Bool_t isVZERO = kFALSE,
		      Bool_t isPbPb = kTRUE,
		      Bool_t is2011 = kTRUE,
		      AliFlowEventCuts::refMultMethod gCentralityEstimator = AliFlowEventCuts::kVZERO,
		      Int_t gAODfilterBit = 768,
		      Double_t gProbPID = 0.9,
		      UInt_t triggerSelectionString = AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral,
		      Double_t gDCAvtxXY = -1.,
		      Double_t gDCAvtxZ = -1.,
		      Bool_t doQA = kTRUE) {
  //Macro to be used for studies of vn for identified particles
  //The macro uses as an input a configuration macro that 
  //creates the AliFlowEventCuts and AliFlowTrackCuts objects
  //Author:Panos.Christakoglou@nikhef.nl
  //Define centralities  
  static const Int_t nCentralitiesMax = 24;
  Double_t gCentrality[nCentralitiesMax] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,30,40,50};    
  //static const Int_t nCentralitiesMax = 3;
  //Double_t gCentrality[nCentralitiesMax] = {0,10,20};
  
  //======================================================================//
  // FLOW cut objects config
  //gROOT->LoadMacro(configFile);
  //gROOT->Macro(configFile);

  const int nHarmonics = uptoWhichHarmonics - 1;
  AliFlowEventCuts* cutsEvent[nCentralitiesMax];
  AliFlowTrackCuts* cutsRP[nCentralitiesMax];
  AliFlowTrackCuts* cutsPOI[nCentralitiesMax];
  TString suffixName[nCentralitiesMax];
  TString outputSlotName[nCentralitiesMax][nHarmonics];
  for(Int_t iCentralityBin = 0; iCentralityBin < nCentralitiesMax - 1; iCentralityBin++) {
    //Create the event cut object
    cutsEvent[iCentralityBin] = createFlowEventCutObject(gCentrality[iCentralityBin],gCentrality[iCentralityBin+1],isPbPb,is2011,gCentralityEstimator,doQA);

    //Create the RP cut object
    cutsRP[iCentralityBin] = createFlowRPCutObject(gCentrality[iCentralityBin],gCentrality[iCentralityBin+1],isVZERO,is2011,gAODfilterBit,doQA);        

    //Create the POI cut object
    cutsPOI[iCentralityBin] = createFlowPOICutObject(gCentrality[iCentralityBin],gCentrality[iCentralityBin+1],particleSpecies,qVector,gEtaGap,isVZERO,is2011,isPbPb,gAODfilterBit,gProbPID,gDCAvtxXY,gDCAvtxZ,doQA); 
    
    suffixName[iCentralityBin] = particleSpecies.Data();
    if(isPbPb) {
      suffixName[iCentralityBin] += "_";
      suffixName[iCentralityBin] += qVector.Data();
      suffixName[iCentralityBin] += "_Centrality";
      suffixName[iCentralityBin] += gCentrality[iCentralityBin];
      suffixName[iCentralityBin] += "To";
      suffixName[iCentralityBin] += gCentrality[iCentralityBin+1];
    }

    //=====================================================================
    for(int nHarmonic = 2; nHarmonic <= uptoWhichHarmonics; nHarmonic++) 
      outputSlotName[iCentralityBin][nHarmonic - 2] = Form("%s_v%d",suffixName[iCentralityBin].Data(),nHarmonic);
  }//loop over centrality classes
    
  TString fileName = "AnalysisResults";
  fileName.Append(".root");
  
  //========================FLOWPACKAGE TASKS=========================//
  AliAnalysisDataContainer *cinput1[nCentralitiesMax];
  AliAnalysisDataContainer *coutputFE[nCentralitiesMax];
  AliAnalysisDataContainer *coutputFEQA[nCentralitiesMax];
  AliAnalysisTaskFlowEvent *taskFE[nCentralitiesMax];
  
  AliAnalysisDataContainer *flowEvent[nCentralitiesMax][nHarmonics];
  AliAnalysisTaskFilterFE *tskFilter[nCentralitiesMax][nHarmonics];
  
  AliAnalysisDataContainer *coutputSP[nCentralitiesMax][nHarmonics];
  AliAnalysisTaskScalarProduct *taskSP[nCentralitiesMax][nHarmonics];
  
  TString outputQA[nCentralitiesMax];
  TString myNameSP[nCentralitiesMax][nHarmonics];
  TString slot[nCentralitiesMax][nHarmonics];
  
  for(Int_t iCentralityBin = 0; iCentralityBin < nCentralitiesMax - 1; iCentralityBin++) {
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

    for(int nHarmonic = 2; nHarmonic <= uptoWhichHarmonics; nHarmonic++){
      flowEvent[iCentralityBin][nHarmonic-2] = mgr->CreateContainer(Form("Filter_%s",outputSlotName[iCentralityBin][nHarmonic-2].Data()),AliFlowEventSimple::Class(),AliAnalysisManager::kExchangeContainer);
      tskFilter[iCentralityBin][nHarmonic-2] = new AliAnalysisTaskFilterFE(Form("TaskFilter_%s",outputSlotName[iCentralityBin][nHarmonic-2].Data()),cutsRP[iCentralityBin], NULL);
      if(!isVZERO){
	tskFilter[iCentralityBin][nHarmonic-2]->SetSubeventEtaRange(minA, maxA, minB, maxB);
      }
      if(isVZERO) tskFilter[iCentralityBin][nHarmonic-2]->SetSubeventEtaRange(-5,-1.5,+1.5,5);
      mgr->AddTask(tskFilter[iCentralityBin][nHarmonic-2]);
      mgr->ConnectInput( tskFilter[iCentralityBin][nHarmonic-2],0,coutputFE[iCentralityBin]);
      mgr->ConnectOutput(tskFilter[iCentralityBin][nHarmonic-2],1,flowEvent[iCentralityBin][nHarmonic-2]);
      
      taskSP[iCentralityBin][nHarmonic-2] = new AliAnalysisTaskScalarProduct(Form("TaskScalarProduct_%s",outputSlotName[iCentralityBin][nHarmonic-2].Data()),kFALSE);
      taskSP[iCentralityBin][nHarmonic-2]->SetHarmonic(nHarmonic);
      taskSP[iCentralityBin][nHarmonic-2]->SelectCollisionCandidates(triggerSelectionString);
      taskSP[iCentralityBin][nHarmonic-2]->SetRelDiffMsub(1.0);
      taskSP[iCentralityBin][nHarmonic-2]->SetTotalQvector(qVector);
      taskSP[iCentralityBin][nHarmonic-2]->SetApplyCorrectionForNUA(kTRUE);
      
      TString outputSP = fileName;
      outputSP += ":outputSPanalysis";
            
      coutputSP[iCentralityBin][nHarmonic-2] = mgr->CreateContainer(Form("%s",outputSlotName[iCentralityBin][nHarmonic-2].Data()),TList::Class(),AliAnalysisManager::kOutputContainer,outputSP);
      mgr->AddTask(taskSP[iCentralityBin][nHarmonic-2]);
      mgr->ConnectInput(taskSP[iCentralityBin][nHarmonic-2],0,flowEvent[iCentralityBin][nHarmonic-2]);
      mgr->ConnectInput(taskSP[iCentralityBin][nHarmonic-2],0,coutputFE[iCentralityBin]);
      mgr->ConnectOutput(taskSP[iCentralityBin][nHarmonic-2],1,coutputSP[iCentralityBin][nHarmonic-2]);
    }//loop over the harmonics
  
    if (taskFE[iCentralityBin]->GetQAOn()) {
      outputQA[iCentralityBin] = fileName;
      outputQA[iCentralityBin] += ":QA";
      coutputFEQA[iCentralityBin] = mgr->CreateContainer(Form("QA_%s",suffixName[iCentralityBin].Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputQA[iCentralityBin]);
      mgr->ConnectOutput(taskFE[iCentralityBin],2,coutputFEQA[iCentralityBin]);
    }
  }//loop over centralities  
}

//_________________________________________________________//
AliFlowEventCuts *createFlowEventCutObject(Int_t gCentralityMin = -1,
					   Int_t gCentralityMax = -1,
					   Bool_t isPbPb = kTRUE,
					   Bool_t is2011 = kTRUE,
					   AliFlowEventCuts::refMultMethod gCentralityEstimator = AliFlowEventCuts::kVZERO,
					   Bool_t doQA = kFALSE) {
  //Part of the code that creates the event cut objects
  Double_t gVertexZmin = -10., gVertexZmax = 10.;
  
  //Create the event cut objects
  AliFlowEventCuts *cutsEvent = new AliFlowEventCuts(Form("eventcutsCentrality%dTo%d",gCentralityMin,gCentralityMax));
  if(isPbPb) {
    cutsEvent->SetCentralityPercentileRange(gCentralityMin,gCentralityMax);
    cutsEvent->SetCentralityPercentileMethod(gCentralityEstimator);
    cutsEvent->SetLHC11h(is2011);
    cutsEvent->SetCutTPCmultiplicityOutliersAOD(kTRUE);
  }
  else 
    cutsEvent->SetCheckPileup(kTRUE);
             
  cutsEvent->SetPrimaryVertexZrange(gVertexZmin,gVertexZmax);
  cutsEvent->SetQA(doQA);

  //return the object
  return cutsEvent;
}

//_________________________________________________________//
AliFlowTrackCuts *createFlowRPCutObject(Int_t gCentralityMin = -1,
					Int_t gCentralityMax = -1,
					Bool_t isVZERO = kFALSE,
					Bool_t is2011 = kTRUE,
					Int_t gAODfilterBit = 768,
					Bool_t doQA = kFALSE) {
  //Part of the code that creates the RP cut objects
  Double_t gEtaMin = -0.8, gEtaMax = 0.8;
  Int_t gMinTPCdedx = 10;

  //Create the event cut objects
  AliFlowTrackCuts *cutsRP = new AliFlowTrackCuts(Form("rpCutsCentrality%dTo%d",gCentralityMin,gCentralityMax));
  if(!isVZERO) {
    cutsRP->SetPtRange(0.2,5.);
    cutsRP->SetEtaRange(gEtaMin,gEtaMax);
    cutsRP->SetMinimalTPCdedx(gMinTPCdedx);
    cutsRP->SetAODfilterBit(gAODfilterBit);
  }
  else if(isVZERO) { // use vzero sub analysis
    if(!is2011) 
      cutsRP = AliFlowTrackCuts::GetStandardVZEROOnlyTrackCuts2010(); 
    if(is2011)  
      cutsRP = AliFlowTrackCuts::GetStandardVZEROOnlyTrackCuts2011(); 
  }//VZERO
  
  cutsRP->SetQA(doQA);
  
  //return the object
  return cutsRP;
}

//_________________________________________________________//
AliFlowTrackCuts *createFlowPOICutObject(Int_t gCentralityMin = -1,
					 Int_t gCentralityMax = -1,
					 TString particleSpecies = "Pion",
					 TString gQvector = "Qa",
					 Double_t gEtaGap = 0.0,
					 Bool_t isVZERO = kFALSE,
					 Bool_t is2011 = kTRUE,
					 Bool_t isPbPb = kTRUE,
					 Int_t gAODfilterBit = 768,
					 Double_t gProbPID = 0.9,
					 Double_t gDCAvtxXY = -1.,
					 Double_t gDCAvtxZ = -1.,
					 Bool_t doQA = kFALSE) {
  //Part of the code that creates the POI cut objects
  Double_t gEtaMin = -0.8, gEtaMax = 0.8;
  Int_t gMinTPCdedx = 10;
  Bool_t isPID = kTRUE;
  Int_t gCharge = 0;  
  AliFlowTrackCuts::PIDsource sourcePID=AliFlowTrackCuts::kTOFbayesian;
  AliPID::EParticleType particleType = AliPID::kPion;

  if(particleSpecies.Contains("Pion")) {
    //cout<<"(Panos) Pions"<<endl;
    particleType = AliPID::kPion;
  }
  else if(particleSpecies.Contains("Kaon")) {
    //cout<<"(Panos) Kaons"<<endl;
    particleType = AliPID::kKaon;
  }
  else if(particleSpecies.Contains("Proton")) {
    //cout<<"(Panos) Protons"<<endl;
    particleType = AliPID::kProton;
  }
    
  AliFlowTrackCuts *cutsPOI = new AliFlowTrackCuts(Form("poiCutsCentrality%dTo%d",gCentralityMin,gCentralityMax));

  //for 2010 data to use old TPC PID Response instead of the official one
  if(!is2011) 
    cutsPOI->GetBayesianResponse()->ForceOldDedx(); 
  
  if(!isVZERO) {
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
    cutsPOI->SetEtaRange(gEtaMin,gEtaMax);
  
  cutsPOI->SetAcceptKinkDaughters(kFALSE);
  
  if(isPID) {
    cutsPOI->SetPID(particleType, sourcePID, gProbPID);
    
    if(isPbPb) 
      cutsPOI->SetPriors((gCentralityMin + gCentralityMax)*0.5);
  }
  
  if (gCharge!=0) cutsPOI->SetCharge(gCharge);

  cutsPOI->SetMinNClustersTPC(70);
  cutsPOI->SetPtRange(0.2,6.);
  cutsPOI->SetAODfilterBit(gAODfilterBit);
  if(gDCAvtxXY > 0) 
    cutsPOI->SetMaxDCAToVertexXY(gDCAvtxXY);
  if(gDCAvtxZ > 0)
    cutsPOI->SetMaxDCAToVertexZ(gDCAvtxZ);
  cutsPOI->SetMinimalTPCdedx(gMinTPCdedx);
  cutsPOI->SetQA(doQA);

  //return the object
  return cutsPOI;
}
