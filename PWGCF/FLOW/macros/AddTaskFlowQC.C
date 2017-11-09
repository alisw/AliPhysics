void AddTaskFlowQC(TString particleSpecies = "",
                      Int_t uptoWhichHarmonics = 5,
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
        //Create the event cut object
        cutsEvent[iCentralityBin] = createFlowEventCutObject(gCentrality[iCentralityBin],gCentrality[iCentralityBin+1],isPbPb,whichData,gCentralityEstimator,doQA);
        
        //Create the RP cut object
        cutsRP[iCentralityBin] = createFlowRPCutObject(gCentrality[iCentralityBin],gCentrality[iCentralityBin+1],whichData,gAODfilterBit,gCharge,doQA);
        
        //Create the POI cut object
        cutsPOI[iCentralityBin] = createFlowPOICutObject(gCentrality[iCentralityBin],gCentrality[iCentralityBin+1],particleSpecies,isPID,whichData,isPbPb,gAODfilterBit,gProbPID,gDCAvtxXY,gDCAvtxZ,gCharge,doQA);
        
        suffixName[iCentralityBin] = particleSpecies.Data();
        
        if(isPbPb) {
            suffixName[iCentralityBin] += "QC_Centrality";
            suffixName[iCentralityBin] += gCentrality[iCentralityBin];
            suffixName[iCentralityBin] += "To";
            suffixName[iCentralityBin] += gCentrality[iCentralityBin+1];
        }
        
        if(isPID){
            suffixName[iCentralityBin] += "_";
            suffixName[iCentralityBin] += particleSpecies.Data();
            suffixName[iCentralityBin] += "_";

        }
        else{
            suffixName[iCentralityBin] += "_AllCharged";
            suffixName[iCentralityBin] = Form("%s_%i",suffixName[iCentralityBin].Data(),gAODfilterBit);
            
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
    
    AliAnalysisDataContainer *coutputQC[nCentralitiesMax][nHarmonics];
    AliAnalysisTaskQCumulants *taskQC[nCentralitiesMax][nHarmonics];
    
    TString outputQA[nCentralitiesMax];
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
            
            mgr->AddTask(tskFilter[iCentralityBin][nHarmonic-2]);
            mgr->ConnectInput( tskFilter[iCentralityBin][nHarmonic-2],0,coutputFE[iCentralityBin]);
            mgr->ConnectOutput(tskFilter[iCentralityBin][nHarmonic-2],1,flowEvent[iCentralityBin][nHarmonic-2]);
            
            
            
            taskQC[iCentralityBin][nHarmonic-2] = new AliAnalysisTaskQCumulants(Form("TaskQCumulants_%s",outputSlotName[iCentralityBin][nHarmonic-2].Data()),kFALSE);
            taskQC[iCentralityBin][nHarmonic-2]->SetHarmonic(nHarmonic);
            taskQC[iCentralityBin][nHarmonic-2]->SelectCollisionCandidates(triggerSelectionString);
            taskQC[iCentralityBin][nHarmonic-2]->SetApplyCorrectionForNUA(kTRUE);
            
            
            TString outputQC = fileName;
            outputQC += ":outputQCanalysis";
            
            coutputQC[iCentralityBin][nHarmonic-2] = mgr->CreateContainer(Form("%s",outputSlotName[iCentralityBin][nHarmonic-2].Data()),TList::Class(),AliAnalysisManager::kOutputContainer,outputQC);
            
            mgr->AddTask(taskQC[iCentralityBin][nHarmonic-2]);
            
            
            mgr->ConnectInput(taskQC[iCentralityBin][nHarmonic-2],0,flowEvent[iCentralityBin][nHarmonic-2]);
            mgr->ConnectInput(taskQC[iCentralityBin][nHarmonic-2],0,coutputFE[iCentralityBin]);
            mgr->ConnectOutput(taskQC[iCentralityBin][nHarmonic-2],1,coutputQC[iCentralityBin][nHarmonic-2]);
       	  
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
AliFlowEventCuts *createFlowEventCutObject(Double_t gCentralityMin = -1,
                                           Double_t gCentralityMax = -1,
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
                                        Double_t whichData = 2011,
                                        Int_t gAODfilterBit = 768,
                                        Double_t gCharge = 0.,
                                        Bool_t doQA = kFALSE) {
    //Part of the code that creates the RP cut objects
    Double_t gEtaMin = -0.8, gEtaMax = 0.8;
    Int_t gMinTPCdedx = 10;
    
    //Create the event cut objects
    AliFlowTrackCuts *cutsRP = new AliFlowTrackCuts(Form("rpCutsCentrality%dTo%d",gCentralityMin,gCentralityMax));
    
    if (gCharge!=0) cutsRP->SetCharge(gCharge);
    
        cutsRP->SetMinNClustersTPC(70);
        cutsRP->SetPtRange(0.2,10.);
        cutsRP->SetEtaRange(gEtaMin,gEtaMax);
        cutsRP->SetMinimalTPCdedx(gMinTPCdedx);
        cutsRP->SetAODfilterBit(gAODfilterBit);
        cutsRP->SetQA(doQA);
    
    //return the object
    return cutsRP;
}

//_________________________________________________________//
AliFlowTrackCuts *createFlowPOICutObject(Int_t gCentralityMin = -1,
                                         Int_t gCentralityMax = -1,
                                         TString particleSpecies = "Pion",
					 Bool_t isPID = kFALSE,
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
    
    AliFlowTrackCuts *cutsPOI = new AliFlowTrackCuts(Form("poiCutsCentrality%dTo%d",gCentralityMin,gCentralityMax));
    
    if (gCharge!=0) cutsPOI->SetCharge(gCharge);
    
    cutsPOI->SetMinNClustersTPC(70);
    cutsPOI->SetPtRange(0.2,10.);
    cutsPOI->SetEtaRange(gEtaMin,gEtaMax);
    cutsPOI->SetMinimalTPCdedx(gMinTPCdedx);
    cutsPOI->SetAODfilterBit(gAODfilterBit);
    if(gDCAvtxXY > 0)
        cutsPOI->SetMaxDCAToVertexXY(gDCAvtxXY);
    if(gDCAvtxZ > 0)
        cutsPOI->SetMaxDCAToVertexZ(gDCAvtxZ);
        cutsPOI->SetQA(doQA);
    
    if (particleSpecies != "")
    isPID = kTRUE;
    
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
    
    
    
    //for 2010 data to use old TPC PID Response instead of the official one
    if(whichData!=2011)
        cutsPOI->GetBayesianResponse()->ForceOldDedx();
    
    
    if(isPID) {
        cutsPOI->SetPID(particleType, sourcePID, gProbPID);
        
    if(isPbPb)
            cutsPOI->SetPriors((gCentralityMin + gCentralityMax)*0.5);
    }
    
        //return the object
    return cutsPOI;
}
