void AddTaskFlowHigherOrdersAllPID(Int_t triggerSelectionString=AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral,
                                   Float_t centrMin=0.,
                                   Float_t centrMax=100.,
                                   Float_t etamin=-0.8,
                                   Float_t etamax=0.8,
                                   Float_t EtaGap=0,
                                   TString fileNameBase="output",
                                   TString uniqueStr="",
                                   Int_t AODfilterBitRP = 768,
                                   Int_t AODfilterBitPOI = 768,
                                   Int_t charge=0,
                                   Bool_t doQA=kTRUE,
                                   Bool_t doPIDQA=kFALSE,
                                   Bool_t isPID = kFALSE,
                                   Bool_t is2011 = kTRUE,
                                   AliPID::EParticleType particleType=AliPID::kPion,
                                   AliFlowTrackCuts::PIDsource sourcePID=AliFlowTrackCuts::kTOFbayesian) {
  // Define the range for eta subevents (for SP method)
  Double_t minA = -0.8;//
  Double_t maxA = 0.8;//
  Double_t minB = -0.8;//
  Double_t maxB = 0.8;//

  if(EtaGap!=0 && EtaGap<=1){
      minA = -0.8;
      maxA = -EtaGap/2;
      minB = EtaGap/2;
      maxB = 0.8;
  }
   
  // AFTERBURNER
  Bool_t useAfterBurner=kFALSE;
  Double_t v1=0.0;
  Double_t v2=0.0;
  Double_t v3=0.0;
  Double_t v4=0.0;
  Int_t numberOfTrackClones=0; //non-flow

  // Define a range of the detector to exclude
  Bool_t ExcludeRegion = kFALSE;
  Double_t excludeEtaMin = -0.;
  Double_t excludeEtaMax = 0.;
  Double_t excludePhiMin = 0.;
  Double_t excludePhiMax = 0.;

  // use physics selection class
  Bool_t  UsePhysicsSelection = kTRUE;

  // QA
//  Bool_t runQAtask=kFALSE;
//  Bool_t FillQAntuple=kFALSE;
//  Bool_t DoQAcorrelations=kFALSE;

  // RUN SETTINGS
  // Flow analysis method can be:(set to kTRUE or kFALSE)
  Bool_t SP       = kTRUE;  // scalar product method (similar to eventplane method)
  Bool_t QC       = kTRUE;  // cumulants using Q vectors
  
  Bool_t METHODS[] = {SP,QC};

  // Boolean to use/not use weights for the Q vector
  Bool_t WEIGHTS[] = {kFALSE,kFALSE,kFALSE}; //Phi, v'(pt), v'(eta)

  // SETTING THE CUTS

  //---------Data selection----------
  //kMC, kGlobal, kESD_TPConly, kESD_SPDtracklet
  AliFlowTrackCuts::trackParameterType rptype = AliFlowTrackCuts::kTPCstandalone;
  AliFlowTrackCuts::trackParameterType poitype = AliFlowTrackCuts::kTPCstandalone;

  //---------Parameter mixing--------
  //kPure - no mixing, kTrackWithMCkine, kTrackWithMCPID, kTrackWithMCpt
  AliFlowTrackCuts::trackParameterMix rpmix = AliFlowTrackCuts::kPure;
  AliFlowTrackCuts::trackParameterMix poimix = AliFlowTrackCuts::kPure;


  const char* rptypestr = AliFlowTrackCuts::GetParamTypeName(rptype);
  const char* poitypestr = AliFlowTrackCuts::GetParamTypeName(poitype);

  //===========================================================================
  // EVENTS CUTS:
  AliFlowEventCuts* cutsEvent = new AliFlowEventCuts("event cuts");
  cutsEvent->SetUsedDataset(is2011);
  cutsEvent->SetCentralityPercentileRange(centrMin,centrMax);
  cutsEvent->SetCentralityPercentileMethod(AliFlowEventCuts::kV0);
//  cutsEvent->SetRefMultMethod(AliFlowEventCuts::kVZERO);
  //cutsEvent->SetCentralityPercentileMethod(AliFlowEventCuts::kSPD1tracklets);
  //cutsEvent->SetNContributorsRange(2);
  cutsEvent->SetPrimaryVertexZrange(-10.,10.);
  cutsEvent->SetQA(doQA);
  cutsEvent->SetCutTPCmultiplicityOutliers();


  // RP TRACK CUTS:
//  AliFlowTrackCuts* cutsRP2 = AliFlowTrackCuts::GetStandardVZEROOnlyTrackCuts();
  AliFlowTrackCuts* cutsRP = new AliFlowTrackCuts("TPConlyRP");
  cutsRP->SetParamType(rptype);
  cutsRP->SetParamMix(rpmix);
  cutsRP->SetPtRange(0.2,5.);
  cutsRP->SetEtaRange(etamin,etamax);
  cutsRP->SetMinNClustersTPC(70);
//  cutsRP->SetMinChi2PerClusterTPC(0.1);//
 // cutsRP->SetMaxChi2PerClusterTPC(4.0);//
  cutsRP->SetMaxDCAToVertexXY(3.0);
  cutsRP->SetMaxDCAToVertexZ(3.0);
  cutsRP->SetAcceptKinkDaughters(kFALSE);
  cutsRP->SetMinimalTPCdedx(10.);
  cutsRP->SetAODfilterBit(AODfilterBitRP);    
  cutsRP->SetQA(doQA);

  // POI TRACK CUTS:
  AliFlowTrackCuts* cutsPOI = new AliFlowTrackCuts("TPConlyPOI");
  cutsPOI->GetBayesianResponse()->ForceOldDedx(); // for 2010 data to use old TPC PID Response instead of the official one
  cutsPOI->SetParamType(poitype);
  cutsPOI->SetParamMix(poimix);
  cutsPOI->SetPtRange(0.2,5.);//
  cutsPOI->SetEtaRange(etamin,etamax);
  //cutsPOI->SetRequireCharge(kTRUE);
  //cutsPOI->SetPID(PdgRP);
  cutsPOI->SetMinNClustersTPC(70);
 // cutsPOI->SetMinChi2PerClusterTPC(0.1); //
 // cutsPOI->SetMaxChi2PerClusterTPC(4.0); //
//  cutsPOI->SetRequireITSRefit(kTRUE);
//  cutsPOI->SetRequireTPCRefit(kTRUE);
//  cutsPOI->SetMinNClustersITS(2);
  //cutsPOI->SetMaxChi2PerClusterITS(1.e+09);
  cutsPOI->SetMaxDCAToVertexXY(3.0);
  cutsPOI->SetMaxDCAToVertexZ(3.0);
  //cutsPOI->SetDCAToVertex2D(kTRUE);
  //cutsPOI->SetMaxNsigmaToVertex(1.e+10);
  //cutsPOI->SetRequireSigmaToVertex(kFALSE);
  cutsPOI->SetAcceptKinkDaughters(kFALSE);
  if(isPID) cutsPOI->SetPID(particleType, sourcePID);//particleType, sourcePID
  if (charge!=0) cutsPOI->SetCharge(charge);
  //cutsPOI->SetAllowTOFmismatch(kFALSE);
  cutsPOI->SetRequireStrictTOFTPCagreement(kTRUE);
  //iexample: francesco's tunig TPC Bethe Bloch for data:
  //cutsPOI->GetESDpid().GetTPCResponse().SetBetheBlochParameters(4.36414e-02,1.75977e+01,1.14385e-08,2.27907e+00,3.36699e+00);
  //cutsPOI->GetESDpid().GetTPCResponse().SetMip(49);
  cutsPOI->SetMinimalTPCdedx(10.);
  cutsPOI->SetAODfilterBit(AODfilterBitPOI);
 // cutsPOI->SetAODfilterBit(768);  
  cutsPOI->SetQA(doQA);
  cutsPOI->SetPriors((centrMin+centrMax)*0.5); // set priors and PID as a function of the centrality

//  Int_t sourcePID = 2;
//  Int_t particleType = 2; //2 or 3 or 4 for kPion, kKaon and kProton 
 
  TString outputSlotName[] = {"","","",""};

  for(int harmonic=2;harmonic<6;harmonic++){  //for v3,v4 and v5
    outputSlotName[harmonic-2]+=uniqueStr;
    outputSlotName[harmonic-2]+=Form("%i",harmonic);
    outputSlotName[harmonic-2]+=cutsRP->GetName();
    outputSlotName[harmonic-2]+="_";
    outputSlotName[harmonic-2]+=cutsPOI->GetName();
    outputSlotName[harmonic-2]+=Form("_%.0f-",centrMin);
    outputSlotName[harmonic-2]+=Form("%.0f_",centrMax);
    if(isPID){
      outputSlotName[harmonic-2]+=AliFlowTrackCuts::PIDsourceName(sourcePID);//sourcePID
      outputSlotName[harmonic-2]+="_";
      outputSlotName[harmonic-2]+=AliPID::ParticleName(particleType);//particleType
    }
    else{
      outputSlotName[harmonic-2]+="AllCharged";
    }
    if (charge<0) outputSlotName[harmonic-2]+="-";
    if (charge>0) outputSlotName[harmonic-2]+="+";
  }
/*
  TString QASlotName;

//QA

    QASlotName+="qa";
    QASlotName+=Form("%i",harmonic);
    QASlotName+="_";
    QASlotName+=Form("_%.0f-",centrMin);
    QASlotName+=Form("%.0f_",centrMax);
    if(isPID){
      QASlotName+=AliFlowTrackCuts::PIDsourceName(sourcePID);//sourcePID
      QASlotName+="_";
      QASlotName+=AliPID::ParticleName(particleType);//particleType
    }
    else{
      QASlotName+="AllCharged";
    }
    if (charge<0) QASlotName+="-";
    if (charge>0) QASlotName+="+";
*/

  TString fileName(fileNameBase);
  fileName.Append(".root");

  Bool_t useWeights  = WEIGHTS[0] || WEIGHTS[1] || WEIGHTS[2];
  if (useWeights) cout<<"Weights are used"<<endl;
  else cout<<"Weights are not used"<<endl;
  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskFlowEvent", "No analysis manager to connect to.");
    return NULL;
  }
  
  // Check the analysis type using the event handlers connected to the analysis
  // manager. The availability of MC handler can also be checked here.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskFlowEvent", "This task requires an input event handler");
    return NULL;
  }  

  // Open external input files
  //===========================================================================
  //weights: 
  TFile *weightsFile = NULL;
  TList *weightsList = NULL;

  if(useWeights) {
    //open the file with the weights:
    weightsFile = TFile::Open("weights.root","READ");
    if(weightsFile) {
      //access the list which holds the histos with weigths:
      weightsList = (TList*)weightsFile->Get("weights");
    }
    else {
      cout<<" WARNING: the file <weights.root> with weights from the previous run was not available."<<endl;
      break;
    } 
  }
   
  // Create the flow event task, add it to the manager.
  //===========================================================================
  AliAnalysisTaskFlowEvent *taskFE[4];
  for(int i=0;i<4;i++){
  if(useAfterBurner){ 
      taskFE[i] = new AliAnalysisTaskFlowEvent(Form("TaskFlowEvent_%s",outputSlotName[i].Data()),"",doQA,1);
      taskFE[i]->SetFlow(v1,v2,v3,v4); 
      taskFE[i]->SetNonFlowNumberOfTrackClones(numberOfTrackClones);
      taskFE[i]->SetAfterburnerOn();
      taskFE[i]->SelectCollisionCandidates(triggerSelectionString);
  } 
  else{
      taskFE[i] = new AliAnalysisTaskFlowEvent(Form("TaskFlowEvent_%s",outputSlotName[i].Data()),"",doQA); 
      taskFE[i]->SelectCollisionCandidates(triggerSelectionString);
  }
  if (ExcludeRegion) {
      taskFE[i]->DefineDeadZone(excludeEtaMin, excludeEtaMax, excludePhiMin, excludePhiMax); 
  }
  taskFE[i]->SetSubeventEtaRange(minA, maxA, minB, maxB);
  if (UsePhysicsSelection) {
    taskFE[i]->SelectCollisionCandidates(AliVEvent::kMB);
    cout<<"Using Physics Selection"<<endl;
  }
  mgr->AddTask(taskFE[i]);
 
  // Pass cuts for RPs and POIs to the task:
  taskFE[i]->SetCutsEvent(cutsEvent);
  taskFE[i]->SetCutsRP(cutsRP);
  taskFE[i]->SetCutsPOI(cutsPOI);
  if (cutsRP->GetParamType()==AliFlowTrackCuts::kVZERO)
  { 
    //TODO: since this is set in a static object all analyses in an analysis train
    //will be affected.
    taskFE[i]->SetHistWeightvsPhiMin(0.);
    taskFE[i]->SetHistWeightvsPhiMax(200.);
  }
  }
  // Create the analysis tasks, add them to the manager.
  //===========================================================================
  AliAnalysisTaskScalarProduct *taskSP[4];
  AliAnalysisTaskQCumulants *taskQC[4];
  AliAnalysisTaskPIDqa *taskPIDQA[4];

  for(int i=0;i<4;i++){
  if (SP){
    taskSP[i] = new AliAnalysisTaskScalarProduct(Form("TaskScalarProduct_%s",outputSlotName[i].Data()),WEIGHTS[0]);
    taskSP[i]->SetHarmonic(i+2);
    taskSP[i]->SelectCollisionCandidates(triggerSelectionString);
    taskSP[i]->SetRelDiffMsub(1.0);
    taskSP[i]->SetApplyCorrectionForNUA(kTRUE);
    mgr->AddTask(taskSP[i]);
  }
  if (QC){
    taskQC[i] = new AliAnalysisTaskQCumulants(Form("TaskQCumulants_%s",outputSlotName[i].Data()),useWeights);
    taskQC[i]->SelectCollisionCandidates(triggerSelectionString);
    taskQC[i]->SetUsePhiWeights(WEIGHTS[0]); 
    taskQC[i]->SetUsePtWeights(WEIGHTS[1]);
    taskQC[i]->SetUseEtaWeights(WEIGHTS[2]); 
    taskQC[i]->SetCalculateCumulantsVsM(kFALSE);
    taskQC[i]->SetnBinsMult(10000);
    taskQC[i]->SetMinMult(0.);
    taskQC[i]->SetMaxMult(10000.);
    taskQC[i]->SetHarmonic(i+2);
    taskQC[i]->SetApplyCorrectionForNUA(kFALSE);
    taskQC[i]->SetFillMultipleControlHistograms(kFALSE);     
    mgr->AddTask(taskQC[i]);
  }
  }
//  if(doPIDQA){
//     gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
//     gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDqa.C");    
//     taskPIDQA = new AliAnalysisTaskPIDqa(Form("TaskPIDQA_%s",QASlotName.Data()));
//     taskPIDQA->SelectCollisionCandidates(triggerSelectionString);
//     AddTaskPIDResponse(kFALSE);
 //    AddTaskPIDqa();     
//     mgr->AddTask(taskPIDQA);

     
//  }

  // Create the output container for the data produced by the task
  // Connect to the input and output containers
  //===========================================================================
  AliAnalysisDataContainer *cinput1[4];
  AliAnalysisDataContainer *coutputFE[4];
  AliAnalysisDataContainer* coutputFEQA[4];
  AliAnalysisDataContainer* coutput_QA[4];
  for(int i=0; i<4; i++) {
//    TString output_QA = fileName;
//    output_QA = ":FlowEventQA"
    cinput1[i] = mgr->GetCommonInputContainer();
    
    coutputFE[i] = mgr->CreateContainer(Form("FlowEventSimple_%s",outputSlotName[i].Data()),AliFlowEventSimple::Class()/*AliVEvent::Class()*/,AliAnalysisManager::kExchangeContainer);
//    coutput_QA[i] = mgr->CreateContainer(Form("FlowEventSimpleQA_%s",outputSlotName[i].Data()),TList::Class(),AliAnalysisManager::kOutputContainer,output_QA);
    mgr->ConnectInput(taskFE[i],0,cinput1[i]); 
    mgr->ConnectOutput(taskFE[i],1,coutputFE[i]);
//    mgr->ConnectOutput(taskFE[i],2,coutput_QA[i]);
    
    if (taskFE[i]->GetQAOn()) {
      TString outputQA = fileName;
      outputQA += ":QA";
      coutputFEQA[i] = mgr->CreateContainer(Form("QA_%s",outputSlotName[i].Data()), TList/*AliVEvent*/::Class(),AliAnalysisManager::kOutputContainer,outputQA);
  //  mgr->ConnectInput(taskPIDQA[i],0,cinput1[i]);  
    mgr->ConnectOutput(taskFE[i],2,coutputFEQA[i]);
    }
  }
  // Create the output containers for the data produced by the analysis tasks
  // Connect to the input and output containers
  //===========================================================================
  AliAnalysisDataContainer *cinputWeights[4];
  AliAnalysisDataContainer *coutputSP[4]; 
  AliAnalysisDataContainer *coutputQC[4]; 
  AliAnalysisDataContainer *coutputQA[4];
  for(int i=0;i<4;i++) {
    if (useWeights) {    
      cinputWeights[i] = mgr->CreateContainer(Form("Weights_%s",outputSlotName[i].Data()),TList::Class(),AliAnalysisManager::kInputContainer); 
    }
    
    if(SP) {
      TString outputSP = fileName;
      outputSP += ":outputSPanalysis";
      outputSP+= rptypestr;
      coutputSP[i] = mgr->CreateContainer(Form("SP_%s",outputSlotName[i].Data()), 
					  TList::Class(),AliAnalysisManager::kOutputContainer,outputSP); 
      mgr->ConnectInput(taskSP[i],0,coutputFE[i]); 
      mgr->ConnectOutput(taskSP[i],1,coutputSP[i]); 
      if (WEIGHTS[0]) {
	mgr->ConnectInput(taskSP[i],1,cinputWeights[i]);
	cinputWeights[i]->SetData(weightsList);
      }
    }
    
    if(QC) {
      TString outputQC = fileName;
      outputQC += ":outputQCanalysis";
      outputQC+= rptypestr;
      
      coutputQC[i] = mgr->CreateContainer(Form("QC_%s",outputSlotName[i].Data()), 
					  TList::Class(),AliAnalysisManager::kOutputContainer,outputQC); 
      mgr->ConnectInput(taskQC[i],0,coutputFE[i]); 
      mgr->ConnectOutput(taskQC[i],1,coutputQC[i]);
      if (useWeights) {
	mgr->ConnectInput(taskQC[i],1,cinputWeights[i]);
	cinputWeights[i]->SetData(weightsList);
      }
    }
/*     if(doPIDQA){
      TString outputPIDQA = fileName;
      outputPIDQA +=":outputQAanalysis";
   //  // cInputFEQA[i] = mgr->CreateContainer(i);
      coutputQA[i] = mgr->CreateContainer(Form("PIDQA_%s",outputSlotName[i].Data()),TList::Class(),AliAnalysisManager::kOutputContainer,outputPIDQA);
      mgr->ConnectInput(taskPIDQA[i],0,mgr->GetCommonInputContainer());
    //  mgr->ConnectInput(taskPIDQA[i],0,cInputFEQA[i]);
      mgr->ConnectOutput(taskPIDQA[i],1,coutputQA[i]);
   }*/
  }
 

}


