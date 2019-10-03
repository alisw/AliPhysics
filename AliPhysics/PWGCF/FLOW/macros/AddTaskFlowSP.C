void AddTaskFlowSP(TString particleSpecies = "",
                      TString qVector = "Qa",
                      Int_t uptoWhichHarmonics = 5,
                      Double_t gEtaGap = 0.0,
                      Bool_t isVZERO = kFALSE,
                      Bool_t isPbPb = kTRUE,
                      Bool_t isHijing = kFALSE,
                      Bool_t isPID = kFALSE,
                      Double_t whichData = 2011,
                      AliFlowEventCuts::refMultMethod gCentralityEstimator = AliFlowEventCuts::kVZERO,
                      Int_t gAODfilterBit = 768,
                      Double_t gProbPID = 0.9,
                      UInt_t triggerSelectionString = AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral,
                      Double_t gDCAvtxXY = -1.,
                      Double_t gDCAvtxZ = -1.,
                      Double_t gCharge = 0,
                      Bool_t doQA = kTRUE) {
  //Macro to be used for studies of vn for identified particles
  //The macro uses as an input a configuration macro that
  //creates the AliFlowEventCuts and AliFlowTrackCuts objects
  //Author:Panos.Christakoglou@nikhef.nl
  //Define centralities
  static const Int_t nCentralitiesMax = 24;
  Double_t gCentrality[nCentralitiesMax];
  
  if(isHijing) {
    gCentrality[0] = 0;  
    gCentrality[1] = 50;  
  }
  else {
    if(whichData==2011) {
      for(Int_t i = 0; i < nCentralitiesMax-3; i++) 
	gCentrality[i] = i;
      gCentrality[21] = 30;
      gCentrality[22] = 40;
      gCentrality[23] = 50;
    }
    else if(whichData!=2011) {
      for(Int_t i = 0; i < 11; i++)
	gCentrality[i] = 0 + i*10;
    }
  }
  
  Int_t nCentralities = 2;
  if((isPbPb)&&(whichData==2011)&&(!isHijing)) nCentralities = nCentralitiesMax;
  else if((isPbPb)&&(whichData!=2011)&&(!isHijing)) nCentralities = 11;
 
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
  //Int_t nCentralities = 2;
  //if((isPbPb)&&(is2011)&&(!isHijing)) nCentralities = nCentralitiesMax;
  //else if((isPbPb)&&(!is2011)&&(!isHijing)) nCentralities = 6;
  for(Int_t iCentralityBin = 0; iCentralityBin < nCentralities - 1; iCentralityBin++) {
  //cout<<"the centrality is "<<iCentralityBin<<endl;     
//Create the event cut object
    cutsEvent[iCentralityBin] = createFlowEventCutObject(gCentrality[iCentralityBin],gCentrality[iCentralityBin+1],isPbPb,whichData,gCentralityEstimator,doQA);
    
    //Create the RP cut object
    cutsRP[iCentralityBin] = createFlowRPCutObject(gCentrality[iCentralityBin],gCentrality[iCentralityBin+1],isVZERO,whichData,gAODfilterBit,gCharge,doQA);
    
    //Create the POI cut object
    cutsPOI[iCentralityBin] = createFlowPOICutObject(gCentrality[iCentralityBin],gCentrality[iCentralityBin+1],particleSpecies,qVector,gEtaGap,isPID,isVZERO,whichData,isPbPb,gAODfilterBit,gProbPID,gDCAvtxXY,gDCAvtxZ,gCharge,doQA);
    
    if(isPbPb) {
      suffixName[iCentralityBin] = "SP_Centrality";
      suffixName[iCentralityBin] += gCentrality[iCentralityBin];
      suffixName[iCentralityBin] += "To";
      suffixName[iCentralityBin] += gCentrality[iCentralityBin+1];
    }
      
    if(isPID){
          suffixName[iCentralityBin] += "_";
          suffixName[iCentralityBin] += particleSpecies.Data();
          suffixName[iCentralityBin] += "_";
          suffixName[iCentralityBin] += qVector.Data();
    }
    else{
          suffixName[iCentralityBin] += "_AllCharged_";
          suffixName[iCentralityBin] += qVector.Data();
          suffixName[iCentralityBin] = Form("%s_%.f",suffixName[iCentralityBin].Data(),10*gEtaGap);
        
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
                                           Double_t whichData = 2011,
                                           AliFlowEventCuts::refMultMethod gCentralityEstimator = AliFlowEventCuts::kVZERO,
                                           Bool_t doQA = kFALSE) {
  //Part of the code that creates the event cut objects
  Double_t gVertexZmin = -10., gVertexZmax = 10.;
  
  //Create the event cut objects
    AliFlowEventCuts *cutsEvent = new AliFlowEventCuts(Form("eventcutsCentrality%dTo%d",gCentralityMin,gCentralityMax));
    if(isPbPb) {
        if(whichData==2015)
            cutsEvent->SetCentralityPercentileRange(gCentralityMin,gCentralityMax,kTRUE);
        else
            cutsEvent->SetCentralityPercentileRange(gCentralityMin,gCentralityMax,kFALSE);
        cutsEvent->SetCentralityPercentileMethod(gCentralityEstimator);
        if(whichData==2010)
            cutsEvent->SetLHC11h(kFALSE);
        else if (whichData==2011)
            cutsEvent->SetLHC11h(kTRUE);
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
                                        Double_t whichData = 2011,
                                        Int_t gAODfilterBit = 768,
                                        Double_t gCharge = 0.,
                                        Bool_t doQA = kFALSE) {
    //Part of the code that creates the RP cut objects
    Double_t gEtaMin = -0.8, gEtaMax = 0.8;
    Int_t gMinTPCdedx = 10;
    
    //Create the event cut objects
    AliFlowTrackCuts *cutsRP = new AliFlowTrackCuts(Form("rpCutsCentrality%dTo%d",gCentralityMin,gCentralityMax));
    
    if(!isVZERO) {
        cutsRP->SetPtRange(0.2,10.);
        cutsRP->SetEtaRange(gEtaMin,gEtaMax);
        cutsRP->SetMinimalTPCdedx(gMinTPCdedx);
        cutsRP->SetAODfilterBit(gAODfilterBit);
        
        if (gCharge!=0) cutsRP->SetCharge(gCharge);
    }

    else if(isVZERO) { // use vzero sub analysis
        if(whichData != 2011)
            cutsRP = AliFlowTrackCuts::GetStandardVZEROOnlyTrackCuts2010();
        if(whichData == 2011)
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
                                         Bool_t isPID = kFALSE,
                                         Bool_t isVZERO = kFALSE,
                                         Double_t whichData = 2011,
                                         Bool_t isPbPb = kTRUE,
                                         Int_t gAODfilterBit = 768,
                                         Double_t gProbPID = 0.9,
                                         Double_t gDCAvtxXY = -1.,
                                         Double_t gDCAvtxZ = -1.,
                                         Double_t gCharge = 0.,
                                         Bool_t doQA = kFALSE) {
    //Part of the code that creates the POI cut objects
    Double_t gEtaMin = -0.8, gEtaMax = 0.8;
    Int_t gMinTPCdedx = 10;
   
    if (particleSpecies != "")
        isPID = kTRUE;
    
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
    if(whichData != 2011)
        cutsPOI->GetBayesianResponse()->ForceOldDedx(); 
    
    if(!isVZERO) {
        if(gQvector.Contains("Qa")) {
            //cout<<"(Panos): "<<gQvector.Data()<<endl;
            cutsPOI->SetEtaRange(+0.5*gEtaGap,gEtaMax);
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
    cutsPOI->SetPtRange(0.2,10.);
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
