void AddTaskPIDFlowSP(TString configFile,
		      TString particleSpecies = "Pion",
		      TString qVector = "Qa",
		      Int_t uptoWhichHarmonics = 5,
		      Double_t gEtaGap = 0.0,
		      Bool_t isVZERO = kFALSE,
		      Bool_t isPbPb = kTRUE,
		      UInt_t triggerSelectionString = "AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral",
		      Bool_t doQA = kTRUE) {
  //Macro to be used for studies of vn for identified particles
  //The macro uses as an input a configuration macro that 
  //creates the AliFlowEventCuts and AliFlowTrackCuts objects
  //Author:Panos.Christakoglou@nikhef.nl
  static const Int_t nCentralitiesMax = 100;
  
  //======================================================================//
  // FLOW cut objects config
  gROOT->LoadMacro(configFile);

  //Define centralities
  Int_t nCentralities = 0;
  Double_t gCentrality[nCentralitiesMax];
  defineCentralities(nCentralities,gCentrality);

  const int nHarmonics = uptoWhichHarmonics - 1;
  AliFlowEventCuts* cutsEvent[nCentralitiesMax];
  AliFlowTrackCuts* cutsRP[nCentralitiesMax];
  AliFlowTrackCuts* cutsPOI[nCentralitiesMax];
  TString suffixName[nCentralitiesMax];
  TString outputSlotName[nCentralitiesMax][nHarmonics];
  for(Int_t iCentralityBin = 0; iCentralityBin < nCentralities - 1; iCentralityBin++) {
    //Create the event cut object
    cutsEvent[iCentralityBin] = createFlowEventCutObject(gCentrality[iCentralityBin],gCentrality[iCentralityBin+1],isPbPb,doQA);
    
    //Create the RP cut object
    cutsRP[iCentralityBin] = createFlowRPCutObject(gCentrality[iCentralityBin],gCentrality[iCentralityBin+1],isVZERO,doQA);        

    //Create the POI cut object
    cutsPOI[iCentralityBin] = createFlowPOICutObject(gCentrality[iCentralityBin],gCentrality[iCentralityBin+1],particleSpecies,qVector,gEtaGap,isVZERO,isPbPb,doQA); 
    
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
