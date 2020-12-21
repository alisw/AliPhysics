AliFlowEventCuts *createFlowEventCutObject(Double_t gCentralityMin,
                                           Double_t gCentralityMax,
                                           Bool_t isPbPb,
                                           Double_t whichData,
                                           AliFlowEventCuts::refMultMethod gCentralityEstimator,
                                           Bool_t doQA);
AliFlowTrackCuts *createFlowRPCutObject(Double_t gCentralityMin,
                                        Double_t gCentralityMax,
                                        Double_t gPtMin,
										Double_t gPtMax,
                                        Double_t whichData,
                                        Int_t gAODfilterBit,
                                        Double_t gDCAvtxXY,
                                        Double_t gDCAvtxZ,
                                        Double_t gCharge,
                                        Bool_t doQA);
AliFlowTrackCuts *createFlowPOICutObject(Double_t gCentralityMin,
                                         Double_t gCentralityMax,
                                         Double_t gPtMin,
										 Double_t gPtMax,
                                         TString particleSpecies,
					                     Bool_t isPID,
                                         Double_t whichData,
                                         Bool_t isPbPb,
                                         Int_t gAODfilterBit,
                                         Double_t gProbPID,	
                                         Double_t gDCAvtxXY,
                                         Double_t gDCAvtxZ,
                                         Double_t gCharge,
                                         Bool_t doQA);

void AddTaskFlowQCFineCentBinning(TString particleSpecies = "",
                      Int_t uptoWhichHarmonics = 2,
                      Bool_t isPbPb = kTRUE,
                      Bool_t isHijing = kFALSE,
                      Bool_t isPID = kFALSE,
                      Double_t whichData = 2015,
                      AliFlowEventCuts::refMultMethod gCentralityEstimator = AliFlowEventCuts::kVZERO,
                      Double_t gPtMin = 0.2,
					  Double_t gPtMax = 3,
                      Int_t gAODfilterBit = 768,
                      Double_t gProbPID = 0.9,
                      UInt_t triggerSelectionString = AliVEvent::kINT7,  //AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral,
                      Double_t gDCAvtxXY = 2.4,
                      Double_t gDCAvtxZ = 3.2,
                      Double_t gCharge = 0,
                      Bool_t doQA = kTRUE,
                      Bool_t usePtWeight = kTRUE,
                      Bool_t useEtaWeight = kFALSE,
                      Bool_t usePhiWeight = kFALSE,
                      Int_t whichCentralityBins = 2, // 1 = 0-30%, 2 = 30-60%, 3 = 60-90% Running 0-90% at one time could use too much memory and crash
                      TString weightDirectory = "") {
    
    // Boolean to use/not use weights for the Q vector
	Bool_t WEIGHTS[] = {usePhiWeight,usePtWeight,useEtaWeight}; //Phi, v'(pt), v'(eta)
  
    //Define centralities
    static const Int_t nCentralitiesMax = 31;  // 24
    Double_t gCentrality[nCentralitiesMax];
    
    if(isHijing) {
        gCentrality[0] = 0;
        gCentrality[1] = 50;
    }
    else {
		if (whichCentralityBins == 1) {
			for(Int_t i = 0; i < 31; i++) 
				gCentrality[i] = 0 + i;
		} else if (whichCentralityBins == 2) {
			for(Int_t i = 0; i < 31; i++) 
				gCentrality[i] = 30 + i;
		} else if (whichCentralityBins == 3) {
			for(Int_t i = 0; i < 31; i++) 
				gCentrality[i] = 60 + i;
		} else {
			cout<<"WARNING: whichCentralityBins is a not valid value"<<endl;
			return NULL;
		}
	}
    
    Int_t nCentralities = 2;
    if((isPbPb)&&(whichData==2011)&&(!isHijing)) nCentralities = nCentralitiesMax;
    else if((isPbPb)&&(whichData!=2011)&&(!isHijing)) nCentralities = 31;
    
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
        cutsRP[iCentralityBin] = createFlowRPCutObject(gCentrality[iCentralityBin],gCentrality[iCentralityBin+1],gPtMin,gPtMax,whichData,gAODfilterBit,gDCAvtxXY,gDCAvtxZ,gCharge,doQA);
        
        //Create the POI cut object
        cutsPOI[iCentralityBin] = createFlowPOICutObject(gCentrality[iCentralityBin],gCentrality[iCentralityBin+1],gPtMin,gPtMax,particleSpecies,isPID,whichData,isPbPb,gAODfilterBit,gProbPID,gDCAvtxXY,gDCAvtxZ,gCharge,doQA);
        
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
    
    //========================WEIGHTS=========================//
    Bool_t useWeights  = WEIGHTS[0] || WEIGHTS[1] || WEIGHTS[2];
	if (useWeights) cout<<"Weights are used"<<endl;
	else cout<<"Weights are not used"<<endl;
	// Open external input files
	//===========================================================================
	//weights: 
	const Int_t nWeightsList = 10; // 0-10,10-20,20-30,30-40,40-50,50-60,60-70,70-80,80-90,90-100
	TString weightsListCentBins[nWeightsList + 1] = {"0","10","20","30","40","50","60","70","80","90","100"};
	TFile *weightsFile = NULL;
	TList *weightsList[nWeightsList];
	
	if (useWeights) {
		//open the file with the weights:
		if (weightDirectory != "") weightsFile = TFile::Open(weightDirectory.Data(),"READ");
		else {cout << "ERROR: PtWeightsHist not found!" << endl; return NULL;}
		if(weightsFile) {
			TList* listTemp = weightsFile->GetListOfKeys();
			for (Int_t i = 0; i < listTemp->GetEntries(); i++) {
				weightsList[i] = NULL;
				TString listName = listTemp->At(i)->GetName();
				weightsFile->GetObject(listName.Data(),weightsList[i]);
				if (!weightsList[i]) {
					cout<<" WARNING: the weight lists in <weights.root> were not available."<<endl;
					return NULL;
				}
			}
		} else {
			cout<<" WARNING: the file <weights.root> with weights from the previous run was not available."<<endl;
			return NULL;
		} 
	}

    //========================FLOWPACKAGE TASKS=========================//
    AliAnalysisDataContainer *cinput1[nCentralitiesMax];
    AliAnalysisDataContainer *coutputFE[nCentralitiesMax];
    AliAnalysisDataContainer *coutputFEQA[nCentralitiesMax];
    AliAnalysisTaskFlowEvent *taskFE[nCentralitiesMax];
    
    AliAnalysisDataContainer *flowEvent[nCentralitiesMax][nHarmonics];
    AliAnalysisTaskFilterFE *tskFilter[nCentralitiesMax][nHarmonics];
    
    AliAnalysisDataContainer *coutputQC[nCentralitiesMax][nHarmonics];
    AliAnalysisTaskQCumulants *taskQC[nCentralitiesMax][nHarmonics];
    AliAnalysisDataContainer *cinputWeights[nCentralitiesMax][nHarmonics];
    
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
        
        mgr = AliAnalysisManager::GetAnalysisManager();
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
            
            
            
            //taskQC[iCentralityBin][nHarmonic-2] = new AliAnalysisTaskQCumulants(Form("TaskQCumulants_%s",outputSlotName[iCentralityBin][nHarmonic-2].Data()),kFALSE);
            taskQC[iCentralityBin][nHarmonic-2] = new AliAnalysisTaskQCumulants(Form("TaskQCumulants_%s",outputSlotName[iCentralityBin][nHarmonic-2].Data()),useWeights);
            taskQC[iCentralityBin][nHarmonic-2]->SetUsePhiWeights(WEIGHTS[0]); 
			taskQC[iCentralityBin][nHarmonic-2]->SetUsePtWeights(WEIGHTS[1]);
			taskQC[iCentralityBin][nHarmonic-2]->SetUseEtaWeights(WEIGHTS[2]); 
            taskQC[iCentralityBin][nHarmonic-2]->SetHarmonic(nHarmonic);
            taskQC[iCentralityBin][nHarmonic-2]->SelectCollisionCandidates(triggerSelectionString);
            taskQC[iCentralityBin][nHarmonic-2]->SetApplyCorrectionForNUA(kTRUE);
            
            //=============Weight container==============
            if (useWeights) {    
				cinputWeights[iCentralityBin][nHarmonic-2] = mgr->CreateContainer(Form("Weights_%s",outputSlotName[iCentralityBin][nHarmonic-2].Data()),
											   TList::Class(),AliAnalysisManager::kInputContainer); 
			}
			//===========================================
            TString outputQC = fileName;
            outputQC += ":outputQCanalysis";
            
            coutputQC[iCentralityBin][nHarmonic-2] = mgr->CreateContainer(Form("%s",outputSlotName[iCentralityBin][nHarmonic-2].Data()),TList::Class(),AliAnalysisManager::kOutputContainer,outputQC);
            
            mgr->AddTask(taskQC[iCentralityBin][nHarmonic-2]);
            
            
            mgr->ConnectInput(taskQC[iCentralityBin][nHarmonic-2],0,flowEvent[iCentralityBin][nHarmonic-2]);
            mgr->ConnectInput(taskQC[iCentralityBin][nHarmonic-2],0,coutputFE[iCentralityBin]);
            mgr->ConnectOutput(taskQC[iCentralityBin][nHarmonic-2],1,coutputQC[iCentralityBin][nHarmonic-2]);
            if (useWeights) {
				mgr->ConnectInput(taskQC[iCentralityBin][nHarmonic-2],1,cinputWeights[iCentralityBin][nHarmonic-2]);
				Int_t ithWeightList = gCentrality[iCentralityBin] / 10;
				cinputWeights[iCentralityBin][nHarmonic-2]->SetData(weightsList[ithWeightList]);
				/*if (ithWeightList < 7.5 ) {
					cinputWeights[iCentralityBin][nHarmonic-2]->SetData(weightsList[ithWeightList]);
				} else if (ithWeightList > 7.5) { // no weights available for centrality > 80%. Use weights for 70-80% instead
					cinputWeights[iCentralityBin][nHarmonic-2]->SetData(weightsList[nWeightsList-1]);
				} */
			}
       	  
		} //loop over the harmonics
        
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
    AliFlowEventCuts *cutsEvent = new AliFlowEventCuts(Form("eventcutsCentrality%fTo%f",gCentralityMin,gCentralityMax));
    if(isPbPb) {
      if(whichData==2015)
	cutsEvent->SetCentralityPercentileRange(gCentralityMin,gCentralityMax,kTRUE);
      else
	cutsEvent->SetCentralityPercentileRange(gCentralityMin,gCentralityMax,kFALSE); 
    cutsEvent->SetCentralityPercentileMethod(gCentralityEstimator);
    if(whichData==2010)
	cutsEvent->SetLHC10h(kTRUE);
      else if (whichData==2011)
	cutsEvent->SetLHC11h(kTRUE);
      cutsEvent->SetCutTPCmultiplicityOutliersAOD(kTRUE);
    }
    else
      cutsEvent->SetCheckPileup(kTRUE);
    
    AliFlowTrackCuts* RefMultCuts = new AliFlowTrackCuts("RefMultCuts");
    RefMultCuts->SetParamType(AliFlowTrackCuts::kAODFilterBit);
    RefMultCuts->SetAODfilterBit(768);
    RefMultCuts->SetMinimalTPCdedx(-999999999);
    RefMultCuts->SetMaxDCAToVertexXY(1000.);
    RefMultCuts->SetMaxDCAToVertexZ(1000.);
    RefMultCuts->SetMinNClustersTPC(70.);
    RefMultCuts->SetMinChi2PerClusterTPC(0.1);
    RefMultCuts->SetMaxChi2PerClusterTPC(4.);
    RefMultCuts->SetPtRange(0.2,20.2);
    RefMultCuts->SetEtaRange(-0.8,0.8);
    RefMultCuts->SetAcceptKinkDaughters(kFALSE);
    cutsEvent->SetRefMultCuts(RefMultCuts);
    cutsEvent->SetRefMultMethod(AliFlowEventCuts::kTPConly);
    
    // vertex-z cut
    cutsEvent->SetPrimaryVertexZrange(gVertexZmin,gVertexZmax);
    // enable the qa plots
    cutsEvent->SetQA(doQA);
    // explicit multiplicity outlier cut
    cutsEvent->SetCutTPCmultiplicityOutliersAOD(kTRUE);
    
    cutsEvent->SetCutSPDTRKVtxZ(kTRUE);
    

    //return the object
    return cutsEvent;
}

//_________________________________________________________//
AliFlowTrackCuts *createFlowRPCutObject(Double_t gCentralityMin = -1,
                                        Double_t gCentralityMax = -1,
                                        Double_t gPtMin = 0.2,
										Double_t gPtMax = 3,
                                        Double_t whichData = 2011,
                                        Int_t gAODfilterBit = 768,
                                        Double_t gDCAvtxXY = -1.,
                                        Double_t gDCAvtxZ = -1.,
                                        Double_t gCharge = 0.,
                                        Bool_t doQA = kFALSE) {
    //Part of the code that creates the RP cut objects
    Double_t gEtaMin = -0.8, gEtaMax = 0.8;
    Int_t gMinTPCdedx = 10;
    
    //Create the event cut objects
    AliFlowTrackCuts *cutsRP = new AliFlowTrackCuts(Form("rpCutsCentrality%fTo%f",gCentralityMin,gCentralityMax));
    
    if (gCharge!=0) cutsRP->SetCharge(gCharge);
		
		cutsRP->SetParamType(AliFlowTrackCuts::kAODFilterBit);
        cutsRP->SetMinNClustersTPC(70);
        cutsRP->SetMinNClustersITS(2);
        cutsRP->SetPtRange(gPtMin,gPtMax);
        cutsRP->SetEtaRange(gEtaMin,gEtaMax);
        cutsRP->SetMinimalTPCdedx(gMinTPCdedx);
        cutsRP->SetAODfilterBit(gAODfilterBit);
        cutsRP->SetQA(doQA);
        cutsRP->SetMinChi2PerClusterTPC(0.1);
		cutsRP->SetMaxChi2PerClusterTPC(4.0);
		cutsRP->SetCutChi2PerClusterITS(36);
		cutsRP->SetAcceptKinkDaughters(kFALSE);
		cutsRP->SetMaxFracSharedTPCCluster(0.4);
		cutsRP->SetMaxFracSharedITSCluster(1);
		
		if(gDCAvtxXY > 0)
			cutsRP->SetMaxDCAToVertexXY(gDCAvtxXY);
		if(gDCAvtxZ > 0)
			cutsRP->SetMaxDCAToVertexZ(gDCAvtxZ);
      
    //return the object
    return cutsRP;
}

//_________________________________________________________//
AliFlowTrackCuts *createFlowPOICutObject(Double_t gCentralityMin = -1,
                                         Double_t gCentralityMax = -1,
                                         Double_t gPtMin = 0.2,
										 Double_t gPtMax = 3,
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
    Int_t bCutTPCbound = 0;
    Bool_t bMimicGlobalCuts=kFALSE;
    Bool_t bRequireTOFSignal=kFALSE;
    Bool_t bPtDepDCAxyCut=kFALSE;
    
    AliFlowTrackCuts *cutsPOI = new AliFlowTrackCuts(Form("poiCutsCentrality%fTo%f",gCentralityMin,gCentralityMax));
    
    cutsPOI->SetParamType(AliFlowTrackCuts::kAODFilterBit);
    
    if (gCharge!=0) cutsPOI->SetCharge(gCharge);
    
    cutsPOI->SetMinNClustersTPC(70);
    cutsPOI->SetMinNClustersITS(2);
    cutsPOI->SetPtRange(gPtMin,gPtMax);
    cutsPOI->SetEtaRange(gEtaMin,gEtaMax);
    cutsPOI->SetMinimalTPCdedx(gMinTPCdedx);
    cutsPOI->SetAODfilterBit(gAODfilterBit);
    cutsPOI->SetMinChi2PerClusterTPC(0.1);
	cutsPOI->SetMaxChi2PerClusterTPC(4.0);
	cutsPOI->SetCutChi2PerClusterITS(36);
	cutsPOI->SetAcceptKinkDaughters(kFALSE);
	cutsPOI->SetMaxFracSharedTPCCluster(0.4);
	cutsPOI->SetMaxFracSharedITSCluster(1);
	
	if(bMimicGlobalCuts) {
      cutsPOI->SetMinNClustersTPC(50);
      cutsPOI->SetCutCrossedTPCRows(70,0.8);
      cutsPOI->SetRequireITSRefit(kTRUE);
      cutsPOI->SetMaxDCAToVertexXYPtDepAOD(kTRUE);
      cutsPOI->SetCutGoldenChi2(kTRUE);
      cutsPOI->SetCutChi2PerClusterITS(36.);
      cutsPOI->SetCutITSClusterGlobal(kTRUE);
    }
    cutsPOI->SetRequireTOFSignal(bRequireTOFSignal);
    
	if(bCutTPCbound==1) cutsPOI->SetCutTPCSecbound(kTRUE,gPtMin); // new cut for LHC15o
    if(bCutTPCbound==2) cutsPOI->SetCutTPCSecboundVar(kTRUE); // new cut for LHC15o
	if(bPtDepDCAxyCut) cutsPOI->SetMaxDCAToVertexXYPtDepAOD(kTRUE);
	
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
            cutsPOI->SetPriors((gCentralityMin + gCentralityMax)*0.5); // set my favourite priors for Bayesian PID (requested if Bayesian PID is used)
    }
    
        //return the object
    return cutsPOI;
}
