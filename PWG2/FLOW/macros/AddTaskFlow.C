/////////////////////////////////////////////////////////////////////////////////////////////
//
// AddTask* macro for flow analysis
// Creates a Flow Event task and adds it to the analysis manager.
// Sets the cuts using the correction framework (CORRFW) classes.
// Also creates Flow Analysis tasks and connects them to the output of the flow event task.
//
/////////////////////////////////////////////////////////////////////////////////////////////


// SETTING THE CUTS

// event selection
const Int_t multminESD = 10;  //used for CORRFW cuts 
const Int_t multmaxESD = 1000000; //used for CORRFW cuts 

const Int_t multmin = 10;     //used for AliFlowEventSimple (to set the centrality)
const Int_t multmax = 1000000;     //used for AliFlowEventSimple (to set the centrality)

// For RP selection
const Double_t ptmin1 = 0.0;
const Double_t ptmax1 = 10.0;
const Double_t ymin1  = -1.;
const Double_t ymax1  = 1.;
const Int_t mintrackrefsTPC1 = 2;
const Int_t mintrackrefsITS1 = 3;
const Int_t charge1 = 1;
Bool_t UsePIDforRP = kFALSE;
const Int_t PDG1 = 211;
const Int_t minclustersTPC1 = 50;
const Int_t maxnsigmatovertex1 = 3;

// For for POI selection
const Double_t ptmin2 = 0.0;
const Double_t ptmax2 = 10.0;
const Double_t ymin2  = -1.;
const Double_t ymax2  = 1.;
const Int_t mintrackrefsTPC2 = 2;
const Int_t mintrackrefsITS2 = 3;
const Int_t charge2 = 1;
Bool_t UsePIDforPOI = kFALSE;
const Int_t PDG2 = 321;
const Int_t minclustersTPC2 = 50;
const Int_t maxnsigmatovertex2 = 3;


AliAnalysisTaskFlowEvent* AddTaskFlow(TString type, Bool_t* METHODS, Bool_t QA, Bool_t* WEIGHTS)
{
  //boleans for the methods
  Bool_t SP       = METHODS[0];
  Bool_t LYZ1SUM  = METHODS[1];
  Bool_t LYZ1PROD = METHODS[2];
  Bool_t LYZ2SUM  = METHODS[3];
  Bool_t LYZ2PROD = METHODS[4];
  Bool_t LYZEP    = METHODS[5];
  Bool_t GFC      = METHODS[6];
  Bool_t QC       = METHODS[7];
  Bool_t FQD      = METHODS[8];
  Bool_t MCEP     = METHODS[9];   
 
  //for using weights
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
  // manager. The availability of MC handler cann also be checked here.
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
    
  if (LYZ2SUM){  
    // read the output file from LYZ1SUM 
    TString inputFileNameLYZ2SUM = "outputLYZ1SUManalysis" ;
    inputFileNameLYZ2SUM += type;
    inputFileNameLYZ2SUM += ".root";
    cout<<"The input file is "<<inputFileNameLYZ2SUM.Data()<<endl;
    TFile* fInputFileLYZ2SUM = new TFile(inputFileNameLYZ2SUM.Data(),"READ");
    if(!fInputFileLYZ2SUM || fInputFileLYZ2SUM->IsZombie()) { 
      cerr << " ERROR: To run LYZ2SUM you need the output file from LYZ1SUM. This file is not there! Please run LYZ1SUM first." << endl ; 
      break;
    }
    else {
      TList* fInputListLYZ2SUM = (TList*)fInputFileLYZ2SUM->Get("cobjLYZ1SUM");
      if (!fInputListLYZ2SUM) {cout<<"list is NULL pointer!"<<endl;}
    }
    cout<<"LYZ2SUM input file/list read..."<<endl;
  }
if (LYZ2PROD){  
    // read the output file from LYZ1PROD 
    TString inputFileNameLYZ2PROD = "outputLYZ1PRODanalysis" ;
    inputFileNameLYZ2PROD += type;
    inputFileNameLYZ2PROD += ".root";
    cout<<"The input file is "<<inputFileNameLYZ2PROD.Data()<<endl;
    TFile* fInputFileLYZ2PROD = new TFile(inputFileNameLYZ2PROD.Data(),"READ");
    if(!fInputFileLYZ2PROD || fInputFileLYZ2PROD->IsZombie()) { 
      cerr << " ERROR: To run LYZ2PROD you need the output file from LYZ1PROD. This file is not there! Please run LYZ1PROD first." << endl ; 
      break;
    }
    else {
      TList* fInputListLYZ2PROD = (TList*)fInputFileLYZ2PROD->Get("cobjLYZ1PROD");
      if (!fInputListLYZ2PROD) {cout<<"list is NULL pointer!"<<endl;}
    }
    cout<<"LYZ2PROD input file/list read..."<<endl;
  }
  
  if (LYZEP) {
    // read the output file from LYZ2SUM
    TString inputFileNameLYZEP = "outputLYZ2SUManalysis" ;
    inputFileNameLYZEP += type;
    inputFileNameLYZEP += ".root";
    cout<<"The input file is "<<inputFileNameLYZEP.Data()<<endl;
    TFile* fInputFileLYZEP = new TFile(inputFileNameLYZEP.Data(),"READ");
    if(!fInputFileLYZEP || fInputFileLYZEP->IsZombie()) { 
      cerr << " ERROR: To run LYZEP you need the output file from LYZ2SUM. This file is not there! Please run LYZ2SUM first." << endl ; 
      break;
    }
    else {
      TList* fInputListLYZEP = (TList*)fInputFileLYZEP->Get("cobjLYZ2SUM");
      if (!fInputListLYZEP) {cout<<"list is NULL pointer!"<<endl;}
    }
    cout<<"LYZEP input file/list read..."<<endl;
  }



  // Create the task, add it to the manager.
  //===========================================================================
  AliAnalysisTaskFlowEvent *taskFE = NULL;
  if (QA) { 
    taskFE = new AliAnalysisTaskFlowEvent("TaskFlowEvent",kTRUE); 
    taskFE->SetAnalysisType(type);
    taskFE->SetMinMult(multmin);
    taskFE->SetMaxMult(multmax);
    mgr->AddTask(taskFE);
  }
  else { 
    taskFE = new AliAnalysisTaskFlowEvent("TaskFlowEvent",kFALSE); 
    taskFE->SetAnalysisType(type);
    taskFE->SetMinMult(multmin);
    taskFE->SetMaxMult(multmax);
    mgr->AddTask(taskFE);
  }
 
  // Create cuts using the correction framework (CORRFW)
  //===========================================================================
  if (QA){
    //Set TList for the QA histograms
    TList* qaRP  = new TList(); 
    TList* qaPOI = new TList();
  }

  //############# event cuts on multiplicity
  AliCFEventGenCuts* mcEventCuts = new AliCFEventGenCuts("mcEventCuts","MC-level event cuts");
  mcEventCuts->SetNTracksCut(multminESD,multmaxESD); 
  if (QA) { 
    mcEventCuts->SetQAOn(qaRP);
  }
  AliCFEventRecCuts* recEventCuts = new AliCFEventRecCuts("recEventCuts","rec-level event cuts");
  recEventCuts->SetNTracksCut(multminESD,multmaxESD); 
  if (QA) { 
    recEventCuts->SetQAOn(qaRP);
  }

  //############# cuts on MC
  AliCFTrackKineCuts* mcKineCuts1 = new AliCFTrackKineCuts("mcKineCuts1","MC-level kinematic cuts");
  mcKineCuts1->SetPtRange(ptmin1,ptmax1);
  mcKineCuts1->SetRapidityRange(ymin1,ymax1);
  //mcKineCuts1->SetChargeMC(charge1);
  if (QA) { 
    mcKineCuts1->SetQAOn(qaRP);
  }

  AliCFTrackKineCuts* mcKineCuts2 = new AliCFTrackKineCuts("mcKineCuts2","MC-level kinematic cuts");
  mcKineCuts2->SetPtRange(ptmin2,ptmax2);
  mcKineCuts2->SetRapidityRange(ymin2,ymax2);
  //mcKineCuts2->SetChargeMC(charge2);
  if (QA) { 
    mcKineCuts2->SetQAOn(qaPOI);
  }
  
  AliCFParticleGenCuts* mcGenCuts1 = new AliCFParticleGenCuts("mcGenCuts1","MC particle generation cuts for RP");
  mcGenCuts1->SetRequireIsPrimary();
  if (UsePIDforRP) {mcGenCuts1->SetRequirePdgCode(PDG1);}
  if (QA) { 
    mcGenCuts1->SetQAOn(qaRP);
  }
  
  AliCFParticleGenCuts* mcGenCuts2 = new AliCFParticleGenCuts("mcGenCuts2","MC particle generation cuts for POI");
  mcGenCuts2->SetRequireIsPrimary();
  if (UsePIDforPOI) {mcGenCuts2->SetRequirePdgCode(PDG2);}
  if (QA) { 
    mcGenCuts2->SetQAOn(qaPOI);
  }

  //############# Acceptance Cuts  
  AliCFAcceptanceCuts *mcAccCuts1 = new AliCFAcceptanceCuts("mcAccCuts1","MC acceptance cuts");
  mcAccCuts1->SetMinNHitITS(mintrackrefsITS1);
  mcAccCuts1->SetMinNHitTPC(mintrackrefsTPC1);
  if (QA) { 
    mcAccCuts1->SetQAOn(qaRP);
  }
  
  AliCFAcceptanceCuts *mcAccCuts2 = new AliCFAcceptanceCuts("mcAccCuts2","MC acceptance cuts");
  mcAccCuts2->SetMinNHitITS(mintrackrefsITS2);
  mcAccCuts2->SetMinNHitTPC(mintrackrefsTPC2);
  if (QA) { 
    mcAccCuts2->SetQAOn(qaPOI);
  }
  //############# Rec-Level kinematic cuts
  AliCFTrackKineCuts *recKineCuts1 = new AliCFTrackKineCuts("recKineCuts1","rec-level kine cuts");
  recKineCuts1->SetPtRange(ptmin1,ptmax1);
  recKineCuts1->SetRapidityRange(ymin1,ymax1);
  //recKineCuts1->SetChargeRec(charge1);
  if (QA) { 
    recKineCuts1->SetQAOn(qaRP);
  }
  
  AliCFTrackKineCuts *recKineCuts2 = new AliCFTrackKineCuts("recKineCuts2","rec-level kine cuts");
  recKineCuts2->SetPtRange(ptmin2,ptmax2);
  recKineCuts2->SetRapidityRange(ymin2,ymax2);
  //recKineCuts2->SetChargeRec(charge2);
  if (QA) { 
    recKineCuts2->SetQAOn(qaPOI);
  }
  
  AliCFTrackQualityCuts *recQualityCuts1 = new AliCFTrackQualityCuts("recQualityCuts1","rec-level quality cuts");
  recQualityCuts1->SetMinNClusterTPC(minclustersTPC1);
  recQualityCuts1->SetStatus(AliESDtrack::kITSrefit);
  if (QA) { 
    recQualityCuts1->SetQAOn(qaRP);
  }
  AliCFTrackQualityCuts *recQualityCuts2 = new AliCFTrackQualityCuts("recQualityCuts2","rec-level quality cuts");
  recQualityCuts2->SetMinNClusterTPC(minclustersTPC2);
  recQualityCuts2->SetStatus(AliESDtrack::kITSrefit);
  if (QA) { 
    recQualityCuts2->SetQAOn(qaPOI);
  }
  
  AliCFTrackIsPrimaryCuts *recIsPrimaryCuts1 = new AliCFTrackIsPrimaryCuts("recIsPrimaryCuts1","rec-level isPrimary cuts");
  recIsPrimaryCuts1->SetMaxNSigmaToVertex(maxnsigmatovertex1);
  if (QA) { 
    recIsPrimaryCuts1->SetQAOn(qaRP);
  }
  
  AliCFTrackIsPrimaryCuts *recIsPrimaryCuts2 = new AliCFTrackIsPrimaryCuts("recIsPrimaryCuts2","rec-level isPrimary cuts");
  recIsPrimaryCuts2->SetMaxNSigmaToVertex(maxnsigmatovertex2);
  if (QA) { 
    recIsPrimaryCuts2->SetQAOn(qaPOI);
  }
  
  int n_species = AliPID::kSPECIES ;
  Double_t* prior = new Double_t[n_species];
  
  prior[0] = 0.0244519 ;
  prior[1] = 0.0143988 ;
  prior[2] = 0.805747  ;
  prior[3] = 0.0928785 ;
  prior[4] = 0.0625243 ;
  
  AliCFTrackCutPid* cutPID1 = NULL;
  if(UsePIDforRP) {
    AliCFTrackCutPid* cutPID1 = new AliCFTrackCutPid("cutPID1","ESD_PID for RP") ;
    cutPID1->SetPriors(prior);
    cutPID1->SetProbabilityCut(0.0);
    cutPID1->SetDetectors("TPC TOF");
    switch(TMath::Abs(PDG1)) {
    case 11   : cutPID1->SetParticleType(AliPID::kElectron, kTRUE); break;
    case 13   : cutPID1->SetParticleType(AliPID::kMuon    , kTRUE); break;
    case 211  : cutPID1->SetParticleType(AliPID::kPion    , kTRUE); break;
    case 321  : cutPID1->SetParticleType(AliPID::kKaon    , kTRUE); break;
    case 2212 : cutPID1->SetParticleType(AliPID::kProton  , kTRUE); break;
    default   : printf("UNDEFINED PID\n"); break;
    }
    if (QA) { 
      cutPID1->SetQAOn(qaRP); 
    }
  }
  
  AliCFTrackCutPid* cutPID2 = NULL;
  if (UsePIDforPOI) {
    AliCFTrackCutPid* cutPID2 = new AliCFTrackCutPid("cutPID2","ESD_PID for POI") ;
    cutPID2->SetPriors(prior);
    cutPID2->SetProbabilityCut(0.0);
    cutPID2->SetDetectors("TPC TOF");
    switch(TMath::Abs(PDG2)) {
    case 11   : cutPID2->SetParticleType(AliPID::kElectron, kTRUE); break;
    case 13   : cutPID2->SetParticleType(AliPID::kMuon    , kTRUE); break;
    case 211  : cutPID2->SetParticleType(AliPID::kPion    , kTRUE); break;
    case 321  : cutPID2->SetParticleType(AliPID::kKaon    , kTRUE); break;
    case 2212 : cutPID2->SetParticleType(AliPID::kProton  , kTRUE); break;
    default   : printf("UNDEFINED PID\n"); break;
    }
    if (QA) { 
      cutPID2->SetQAOn(qaPOI);
    }
  }
  
  printf("CREATE EVENT CUTS\n");
  TObjArray* mcEventList = new TObjArray(0);
  mcEventList->AddLast(mcEventCuts);
    
  TObjArray* recEventList = new TObjArray(0);
  recEventList->AddLast(recEventCuts);

  printf("CREATE MC KINE CUTS\n");
  TObjArray* mcList1 = new TObjArray(0);
  mcList1->AddLast(mcKineCuts1);
  mcList1->AddLast(mcGenCuts1);
  
  TObjArray* mcList2 = new TObjArray(0);
  mcList2->AddLast(mcKineCuts2);
  mcList2->AddLast(mcGenCuts2);
  
  printf("CREATE ACCEPTANCE CUTS\n");
  TObjArray* accList1 = new TObjArray(0) ;
  accList1->AddLast(mcAccCuts1);
  
  TObjArray* accList2 = new TObjArray(0) ;
  accList2->AddLast(mcAccCuts2);
  
  printf("CREATE RECONSTRUCTION CUTS\n");
  TObjArray* recList1 = new TObjArray(0) ;
  recList1->AddLast(recKineCuts1);
  recList1->AddLast(recQualityCuts1);
  recList1->AddLast(recIsPrimaryCuts1);
  
  TObjArray* recList2 = new TObjArray(0) ;
  recList2->AddLast(recKineCuts2);
  recList2->AddLast(recQualityCuts2);
  recList2->AddLast(recIsPrimaryCuts2);
  
  printf("CREATE PID CUTS\n");
  TObjArray* fPIDCutList1 = new TObjArray(0) ;
  if(UsePIDforRP) {fPIDCutList1->AddLast(cutPID1);}
  
  TObjArray* fPIDCutList2 = new TObjArray(0) ;
  if (UsePIDforPOI)  {fPIDCutList2->AddLast(cutPID2);}
  
  printf("CREATE INTERFACE AND CUTS\n");
  AliCFManager* cfmgr1 = new AliCFManager();
  cfmgr1->SetNStepEvent(3);
  cfmgr1->SetEventCutsList(AliCFManager::kEvtGenCuts,mcEventList); 
  cfmgr1->SetEventCutsList(AliCFManager::kEvtRecCuts,recEventList); 
  cfmgr1->SetNStepParticle(4); 
  cfmgr1->SetParticleCutsList(AliCFManager::kPartGenCuts,mcList1);
  cfmgr1->SetParticleCutsList(AliCFManager::kPartAccCuts,accList1);
  cfmgr1->SetParticleCutsList(AliCFManager::kPartRecCuts,recList1);
  cfmgr1->SetParticleCutsList(AliCFManager::kPartSelCuts,fPIDCutList1);
  
  AliCFManager* cfmgr2 = new AliCFManager();
  cfmgr1->SetNStepEvent(3);
  cfmgr1->SetEventCutsList(AliCFManager::kEvtGenCuts,mcEventList); 
  cfmgr1->SetEventCutsList(AliCFManager::kEvtRecCuts,recEventList); 
  cfmgr2->SetNStepParticle(4); 
  cfmgr2->SetParticleCutsList(AliCFManager::kPartGenCuts,mcList2);
  cfmgr2->SetParticleCutsList(AliCFManager::kPartAccCuts,accList2);
  cfmgr2->SetParticleCutsList(AliCFManager::kPartRecCuts,recList2);
  cfmgr2->SetParticleCutsList(AliCFManager::kPartSelCuts,fPIDCutList2);
  
  if (QA) {
    taskFE->SetQAList1(qaRP);
    taskFE->SetQAList2(qaPOI);
  }
  taskFE->SetCFManager1(cfmgr1);
  taskFE->SetCFManager2(cfmgr2);



  // Create the analysis tasks, add them to the manager.
  //===========================================================================
  if (SP){
    AliAnalysisTaskScalarProduct *taskSP = new AliAnalysisTaskScalarProduct("TaskScalarProduct");
    mgr->AddTask(taskSP);
  }
  if (LYZ1SUM){
    AliAnalysisTaskLeeYangZeros *taskLYZ1SUM = new AliAnalysisTaskLeeYangZeros("TaskLeeYangZerosSUM",kTRUE);
    taskLYZ1SUM->SetFirstRunLYZ(kTRUE);
    taskLYZ1SUM->SetUseSumLYZ(kTRUE);
    mgr->AddTask(taskLYZ1SUM);
  }
  if (LYZ1PROD){
    AliAnalysisTaskLeeYangZeros *taskLYZ1PROD = new AliAnalysisTaskLeeYangZeros("TaskLeeYangZerosPROD",kTRUE);
    taskLYZ1PROD->SetFirstRunLYZ(kTRUE);
    taskLYZ1PROD->SetUseSumLYZ(kFALSE);
    mgr->AddTask(taskLYZ1PROD);
  }
  if (LYZ2SUM){
    AliAnalysisTaskLeeYangZeros *taskLYZ2SUM = new AliAnalysisTaskLeeYangZeros("TaskLeeYangZerosSUM",kFALSE);
    taskLYZ2SUM->SetFirstRunLYZ(kFALSE);
    taskLYZ2SUM->SetUseSumLYZ(kTRUE);
    mgr->AddTask(taskLYZ2SUM);
  }
  if (LYZ2PROD){
    AliAnalysisTaskLeeYangZeros *taskLYZ2PROD = new AliAnalysisTaskLeeYangZeros("TaskLeeYangZerosPROD",kFALSE);
    taskLYZ2PROD->SetFirstRunLYZ(kFALSE);
    taskLYZ2PROD->SetUseSumLYZ(kFALSE);
    mgr->AddTask(taskLYZ2PROD);
  }
  if (LYZEP){
    AliAnalysisTaskLYZEventPlane *taskLYZEP = new AliAnalysisTaskLYZEventPlane("TaskLYZEventPlane");
    mgr->AddTask(taskLYZEP);
  }
  if (GFC){
    AliAnalysisTaskCumulants *taskGFC = new AliAnalysisTaskCumulants("TaskCumulants",useWeights);
    taskGFC->SetUsePhiWeights(WEIGHTS[0]); 
    taskGFC->SetUsePtWeights(WEIGHTS[1]);
    taskGFC->SetUseEtaWeights(WEIGHTS[2]); 
    mgr->AddTask(taskGFC);
  }
  if (QC){
    AliAnalysisTaskQCumulants *taskQC = new AliAnalysisTaskQCumulants("TaskQCumulants",useWeights);
    taskQC->SetUsePhiWeights(WEIGHTS[0]); 
    taskQC->SetUsePtWeights(WEIGHTS[1]);
    taskQC->SetUseEtaWeights(WEIGHTS[2]); 
    mgr->AddTask(taskQC);
  }
  if (FQD){
    AliAnalysisTaskFittingQDistribution *taskFQD = new AliAnalysisTaskFittingQDistribution("TaskFittingQDistribution",kFALSE);
    taskFQD->SetUsePhiWeights(WEIGHTS[0]); 
    mgr->AddTask(taskFQD);
  }
  if (MCEP){
    AliAnalysisTaskMCEventPlane *taskMCEP = new AliAnalysisTaskMCEventPlane("TaskMCEventPlane");
    mgr->AddTask(taskMCEP);
  }

  // Create the output container for the data produced by the task
  // Connect to the input and output containers
  //===========================================================================
  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutputFE = mgr->CreateContainer("cobjFlowEventSimple",  AliFlowEventSimple::Class(),AliAnalysisManager::kExchangeContainer);
  mgr->ConnectInput(taskFE,0,cinput1); 
  mgr->ConnectOutput(taskFE,0,coutputFE);

  if (QA) { 
    TString qaNameRPFE = "QAforRP_FE_";
    qaNameRPFE += type;
    qaNameRPFE += ".root";
    AliAnalysisDataContainer *coutputQA1FE = 
      mgr->CreateContainer("QARPFE", TList::Class(),AliAnalysisManager::kOutputContainer,qaNameRPFE);
    
    TString qaNamePOIFE = "QAforPOI_FE_";
    qaNamePOIFE += type;
    qaNamePOIFE += ".root";
    AliAnalysisDataContainer *coutputQA2FE = 
      mgr->CreateContainer("QAPOIFE", TList::Class(),AliAnalysisManager::kOutputContainer,qaNamePOIFE);
    
    mgr->ConnectOutput(taskFE,1,coutputQA1FE);
    mgr->ConnectOutput(taskFE,2,coutputQA2FE); 
  }

  // Create the output containers for the data produced by the analysis tasks
  // Connect to the input and output containers
  //===========================================================================
  if (useWeights) {    
    AliAnalysisDataContainer *cinputWeights = mgr->CreateContainer("cobjWeights",TList::Class(),AliAnalysisManager::kInputContainer); 
  }

  if(SP) {
    TString outputSP = "outputSPanalysis";
    outputSP+= type;
    outputSP+= ".root";
    AliAnalysisDataContainer *coutputSP = mgr->CreateContainer("cobjSP", TList::Class(),AliAnalysisManager::kOutputContainer,outputSP);
    mgr->ConnectInput(taskSP,0,coutputFE); 
    mgr->ConnectOutput(taskSP,0,coutputSP);
  }
  if(LYZ1SUM) {
    TString outputLYZ1SUM = "outputLYZ1SUManalysis";
    outputLYZ1SUM+= type;
    outputLYZ1SUM+= ".root";
    AliAnalysisDataContainer *coutputLYZ1SUM = mgr->CreateContainer("cobjLYZ1SUM", TList::Class(),AliAnalysisManager::kOutputContainer,outputLYZ1SUM);
    mgr->ConnectInput(taskLYZ1SUM,0,coutputFE); 
    mgr->ConnectOutput(taskLYZ1SUM,0,coutputLYZ1SUM);
  }
  if(LYZ1PROD) {
    TString outputLYZ1PROD = "outputLYZ1PRODanalysis";
    outputLYZ1PROD+= type;
    outputLYZ1PROD+= ".root";
    AliAnalysisDataContainer *coutputLYZ1PROD = mgr->CreateContainer("cobjLYZ1PROD", TList::Class(),AliAnalysisManager::kOutputContainer,outputLYZ1PROD);
    mgr->ConnectInput(taskLYZ1PROD,0,coutputFE); 
    mgr->ConnectOutput(taskLYZ1PROD,0,coutputLYZ1PROD);
  }
  if(LYZ2SUM) {
    AliAnalysisDataContainer *cinputLYZ2SUM = mgr->CreateContainer("cobjLYZ2SUMin",TList::Class(),AliAnalysisManager::kInputContainer);
    TString outputLYZ2SUM = "outputLYZ2SUManalysis";
    outputLYZ2SUM+= type;
    outputLYZ2SUM+= ".root";
    AliAnalysisDataContainer *coutputLYZ2SUM = mgr->CreateContainer("cobjLYZ2SUM", TList::Class(),AliAnalysisManager::kOutputContainer,outputLYZ2SUM);
    mgr->ConnectInput(taskLYZ2SUM,0,coutputFE); 
    mgr->ConnectInput(taskLYZ2SUM,1,cinputLYZ2SUM);
    mgr->ConnectOutput(taskLYZ2SUM,0,coutputLYZ2SUM);
    cinputLYZ2SUM->SetData(fInputListLYZ2SUM);
  }
  if(LYZ2PROD) {
    AliAnalysisDataContainer *cinputLYZ2PROD = mgr->CreateContainer("cobjLYZ2PRODin",TList::Class(),AliAnalysisManager::kInputContainer);
    TString outputLYZ2PROD = "outputLYZ2PRODanalysis";
    outputLYZ2PROD+= type;
    outputLYZ2PROD+= ".root";
    AliAnalysisDataContainer *coutputLYZ2PROD = mgr->CreateContainer("cobjLYZ2PROD", TList::Class(),AliAnalysisManager::kOutputContainer,outputLYZ2PROD);
    mgr->ConnectInput(taskLYZ2PROD,0,coutputFE); 
    mgr->ConnectInput(taskLYZ2PROD,1,cinputLYZ2PROD);
    mgr->ConnectOutput(taskLYZ2PROD,0,coutputLYZ2PROD);
    cinputLYZ2PROD->SetData(fInputListLYZ2PROD);
  }
  if(LYZEP) {
    AliAnalysisDataContainer *cinputLYZEP = mgr->CreateContainer("cobjLYZEPin",TList::Class(),AliAnalysisManager::kInputContainer);
    TString outputLYZEP = "outputLYZEPanalysis";
    outputLYZEP+= type;
    outputLYZEP+= ".root";
    AliAnalysisDataContainer *coutputLYZEP = mgr->CreateContainer("cobjLYZEP", TList::Class(),AliAnalysisManager::kOutputContainer,outputLYZEP);
    mgr->ConnectInput(taskLYZEP,0,coutputFE); 
    mgr->ConnectInput(taskLYZEP,1,cinputLYZEP);
    mgr->ConnectOutput(taskLYZEP,0,coutputLYZEP);
    cinputLYZEP->SetData(fInputListLYZEP);
  }
  if(GFC) {
    TString outputGFC = "outputGFCanalysis";
    outputGFC+= type;
    outputGFC+= ".root";
    AliAnalysisDataContainer *coutputGFC = mgr->CreateContainer("cobjGFC", TList::Class(),AliAnalysisManager::kOutputContainer,outputGFC);
    mgr->ConnectInput(taskGFC,0,coutputFE); 
    mgr->ConnectOutput(taskGFC,0,coutputGFC);
    if (useWeights) {
      mgr->ConnectInput(taskGFC,1,cinputWeights);
      cinputWeights->SetData(weightsList);
    } 
  }
  if(QC) {
    TString outputQC = "outputQCanalysis";
    outputQC+= type;
    outputQC+= ".root";
    AliAnalysisDataContainer *coutputQC = mgr->CreateContainer("cobjQC", TList::Class(),AliAnalysisManager::kOutputContainer,outputQC);
    mgr->ConnectInput(taskQC,0,coutputFE); 
    mgr->ConnectOutput(taskQC,0,coutputQC);
    if (useWeights) {
      mgr->ConnectInput(taskQC,1,cinputWeights);
      cinputWeights->SetData(weightsList);
    } 
  }
  if(FQD) {
    TString outputFQD = "outputFQDanalysis";
    outputFQD+= type;
    outputFQD+= ".root";
    AliAnalysisDataContainer *coutputFQD = mgr->CreateContainer("cobjFQD", TList::Class(),AliAnalysisManager::kOutputContainer,outputFQD);
    mgr->ConnectInput(taskFQD,0,coutputFE); 
    mgr->ConnectOutput(taskFQD,0,coutputFQD);
    if(useWeights) {
      mgr->ConnectInput(taskFQD,1,cinputWeights);
      cinputWeights->SetData(weightsList);
    } 
  }
  if(MCEP) {
    TString outputMCEP = "outputMCEPanalysis";
    outputMCEP+= type;
    outputMCEP+= ".root";
    AliAnalysisDataContainer *coutputMCEP = mgr->CreateContainer("cobjMCEP", TList::Class(),AliAnalysisManager::kOutputContainer,outputMCEP);
    mgr->ConnectInput(taskMCEP,0,coutputFE); 
    mgr->ConnectOutput(taskMCEP,0,coutputMCEP);
  }
  

  // Return analysis task
  //===========================================================================
  return taskFE;
  


}





