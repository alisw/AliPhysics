/////////////////////////////////////////////////////////////////////////////////////////////
//
// AddTask* macro for flow analysis
// Creates a Flow Event task and adds it to the analysis manager.
// Sets the cuts using the correction framework (CORRFW) classes.
// Also creates Flow Analysis tasks and connects them to the output of the flow event task.
//
/////////////////////////////////////////////////////////////////////////////////////////////

// Define the range for eta subevents (for SP method)
//(FMD 1.7 - 5.0)
//Double_t minA = -5.0;
//Double_t maxA = -1.7;
//Double_t minB = 1.7;
//Double_t maxB = 5.0;
//(Tracklets 0.9 - 2.0)
Double_t minA = -2.0;
Double_t maxA = -0.9;
Double_t minB = 0.9;
Double_t maxB = 2.0;

// use physics selection class
Bool_t  UsePhysicsSelection = kTRUE;

// SETTING THE CUTS

//----------Event selection----------
Bool_t UseMultCutforESD = kTRUE;
//Bool_t UseMultCutforESD = kFALSE;
const Int_t multminESD = 1;  //used for CORRFW cuts 
const Int_t multmaxESD = 1000000; //used for CORRFW cuts 

Bool_t requireVtxCuts = kTRUE;
const Double_t vertexXmin = -1.e99; 
const Double_t vertexXmax = 1.e99;
const Double_t vertexYmin = -1.e99;
const Double_t vertexYmax = 1.e99;
const Double_t vertexZmin = -1.e99; 
const Double_t vertexZmax = 1.e99; 
const Int_t vertexNContributorsmin = 1;
const Int_t vertexNContributorsmax = 10000;

//Bool_t UseMultCut = kFALSE;
Bool_t UseMultCut = kTRUE;
const Int_t multmin = 1;     //used for AliFlowEventSimple (to set the centrality)
const Int_t multmax = 10000;     //used for AliFlowEventSimple (to set the centrality)
//const Int_t multmin = 10;     //used for AliFlowEventSimple (to set the centrality)
//const Int_t multmax = 1000000;     //used for AliFlowEventSimple (to set the centrality)


//----------For RP selection----------
// Use Global tracks ("Global"), or SPD tracklets ("Tracklet") 
// or FMD hits ("FMD") for the RP selection
const TString rptype = "Global";
//const TString rptype = "Tracklet";
//const TString rptype = "FMD";

//KINEMATICS (on generated and reconstructed tracks)
Bool_t UseKineforRP =  kTRUE;
const Double_t ptminRP = 0.0;
const Double_t ptmaxRP = 10.0;
const Double_t etaminRP  = -0.9;
const Double_t etamaxRP  = 0.9;
const Int_t    chargeRP = 1;  //not used
const Bool_t   isChargedRP = kTRUE;

//PID (on generated and reconstructed tracks)
Bool_t UsePIDforRP = kFALSE;
const Int_t PdgRP = 211;

//TRACK QUALITY (on reconstructed tracks only)
//see /CORRFW/AliCFTrackQualityCuts class
Bool_t UseTrackQualityforRP =  kTRUE;
const Int_t    minClustersTpcRP = 80;           //default = -1; 
const Double_t maxChi2PerClusterTpcRP = 4.0;    //default = 1.e+09;
const UShort_t minDedxClusterTpcRP = 0;
const Int_t    minClustersItsRP = 2;            //panos
const Double_t maxChi2PerClusterItsRP = 1.e+09; 
const Int_t    minClustersTrdRP = -1;
const Int_t    minTrackletTrdRP = -1;
const Int_t    minTrackletTrdPidRP = -1;
const Double_t maxChi2PerClusterTrdRP = 1.e+09;
const ULong_t  statusRP = AliESDtrack::kTPCrefit;   //AliESDtrack::kTPCrefit &  AliESDtrack::kITSrefit 

//PRIMARY (on reconstructed tracks only)
//see /CORRFW/AliCFTrackIsPrimaryCuts class
Bool_t UsePrimariesforRP = kTRUE;
const Bool_t   spdVertexRP = kFALSE;
const Bool_t   tpcVertexRP = kFALSE;
const Float_t  minDcaToVertexXyRP = 0.;
const Float_t  minDcaToVertexZRP = 0.;
const Float_t  maxDcaToVertexXyRP = 2.4;         //default = 1.e+10;  //2.4;
const Float_t  maxDcaToVertexZRP = 3.2;          //default = 1.e+10;  //3.2;
const Bool_t   dcaToVertex2dRP = kFALSE;         //default = kFALSE;
const Bool_t   absDcaToVertexRP = kTRUE;         //default = kTRUE;
const Double_t minNSigmaToVertexRP = 0.;
const Double_t maxNSigmaToVertexRP = 1.e+10; //3.; //1.e+10
const Double_t maxSigmaDcaXySP = 1.e+10;
const Double_t maxSigmaDcaZSP = 1.e+10;
const Bool_t   requireSigmaToVertexSP = kFALSE;
const Bool_t   acceptKinkDaughtersSP = kFALSE;  //default = kTRUE;

//ACCEPTANCE (on generated tracks only : AliMCParticle)
//see /CORRFW/AliCFAcceptanceCuts class
Bool_t UseAcceptanceforRP =  kFALSE; 
const Int_t  minTrackrefsItsRP = 0;//3;
const Int_t  minTrackrefsTpcRP = 0;//2;
const Int_t  minTrackrefsTrdRP = 0; 
const Int_t  minTrackrefsTofRP = 0; 
const Int_t  minTrackrefsMuonRP = 0; 
//default for all is 0

//----------For POI selection----------
//KINEMATICS (on generated and reconstructed tracks)
Bool_t UseKineforPOI = kTRUE;
const Double_t ptminPOI = 0.0;
const Double_t ptmaxPOI = 10.0;
const Double_t etaminPOI  = -0.5;
const Double_t etamaxPOI  = 0.5;
const Int_t    chargePOI = 1;  //not used
const Bool_t   isChargedPOI = kTRUE;

//PID (on generated and reconstructed tracks)
Bool_t UsePIDforPOI = kFALSE;
const Int_t PdgPOI = 321;

//TRACK QUALITY (on reconstructed tracks only)
//see /CORRFW/AliCFTrackQualityCuts class
Bool_t UseTrackQualityforPOI = kTRUE;
const Int_t    minClustersTpcPOI = 80;
const Double_t maxChi2PerClusterTpcPOI = 4.0;    
const UShort_t minDedxClusterTpcPOI = 0;
const Int_t    minClustersItsPOI = 2;
const Double_t maxChi2PerClusterItsPOI = 1.e+09;
const Int_t    minClustersTrdPOI = -1;
const Int_t    minTrackletTrdPOI = -1;
const Int_t    minTrackletTrdPidPOI = -1;
const Double_t maxChi2PerClusterTrdPOI = 1.e+09;
const ULong_t  statusPOI = AliESDtrack::kTPCrefit;   

//PRIMARY (on reconstructed tracks only)
//see /CORRFW/AliCFTrackIsPrimaryCuts class
Bool_t UsePrimariesforPOI = kTRUE;
const Bool_t   spdVertexPOI = kFALSE;
const Bool_t   tpcVertexPOI = kFALSE;
const Float_t  minDcaToVertexXyPOI = 0.;
const Float_t  minDcaToVertexZPOI = 0.;
const Float_t  maxDcaToVertexXyPOI = 2.4;
const Float_t  maxDcaToVertexZPOI = 3.2;
const Bool_t   dcaToVertex2dPOI =  kFALSE;
const Bool_t   absDcaToVertexPOI = kTRUE;
const Double_t minNSigmaToVertexPOI = 0.;
const Double_t maxNSigmaToVertexPOI = 1.e+10;  
const Double_t maxSigmaDcaXyPOI = 1.e+10;
const Double_t maxSigmaDcaZPOI = 1.e+10;
const Bool_t   requireSigmaToVertexPOI = kFALSE;
const Bool_t   acceptKinkDaughtersPOI = kFALSE;

//ACCEPTANCE (on generated tracks only : AliMCParticle)
//see /CORRFW/AliCFAcceptanceCuts class
Bool_t UseAcceptanceforPOI = kFALSE;
const Int_t minTrackrefsItsPOI = 3;
const Int_t minTrackrefsTpcPOI = 2;
const Int_t minTrackrefsTrdPOI = 0; 
const Int_t minTrackrefsTofPOI = 0; 
const Int_t minTrackrefsMuonPOI = 0; 


//----------For Adding Flow to the Event----------
const Bool_t AddToEvent = kFALSE;
Double_t ellipticFlow = 0.05;


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
  Bool_t MH       = METHODS[10];
  Bool_t NL       = METHODS[11];  
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
    
  //LYZ2
  if (LYZ2SUM || LYZ2PROD) {
    //read the outputfile of the first run
    TString outputFileName = "AnalysisResults1.root";
    TString pwd(gSystem->pwd());
    pwd+="/";
    pwd+=outputFileName.Data();
    TFile *outputFile = NULL;
    if(gSystem->AccessPathName(pwd.Data(),kFileExists)) {
      cout<<"WARNING: You do not have an output file:"<<endl;
      cout<<"         "<<pwd.Data()<<endl;
      exit(0);
    } else {
      outputFile = TFile::Open(pwd.Data(),"READ");
    }
    
    if (LYZ2SUM){  
      // read the output directory from LYZ1SUM 
      TString inputFileNameLYZ2SUM = "outputLYZ1SUManalysis" ;
      inputFileNameLYZ2SUM += type;
      cout<<"The input directory is "<<inputFileNameLYZ2SUM.Data()<<endl;
      TFile* fInputFileLYZ2SUM = (TFile*)outputFile->FindObjectAny(inputFileNameLYZ2SUM.Data());
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
      // read the output directory from LYZ1PROD 
      TString inputFileNameLYZ2PROD = "outputLYZ1PRODanalysis" ;
      inputFileNameLYZ2PROD += type;
      cout<<"The input directory is "<<inputFileNameLYZ2PROD.Data()<<endl;
      TFile* fInputFileLYZ2PROD = (TFile*)outputFile->FindObjectAny(inputFileNameLYZ2PROD.Data());
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
  }


  if (LYZEP) {
    //read the outputfile of the second run
    TString outputFileName = "AnalysisResults2.root";
    TString pwd(gSystem->pwd());
    pwd+="/";
    pwd+=outputFileName.Data();
    TFile *outputFile = NULL;
    if(gSystem->AccessPathName(pwd.Data(),kFileExists)) {
      cout<<"WARNING: You do not have an output file:"<<endl;
      cout<<"         "<<pwd.Data()<<endl;
      exit(0);
    } else {
      outputFile = TFile::Open(pwd.Data(),"READ");
    }

    // read the output file from LYZ2SUM
    TString inputFileNameLYZEP = "outputLYZ2SUManalysis" ;
    inputFileNameLYZEP += type;
    cout<<"The input file is "<<inputFileNameLYZEP.Data()<<endl;
    TFile* fInputFileLYZEP = (TFile*)outputFile->FindObjectAny(inputFileNameLYZEP.Data());
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
  
  // Create the FMD task and add it to the manager
  //===========================================================================
  AliFMDAnalysisTaskSE *taskfmd = NULL;
  if (rptype == "FMD") {
    taskfmd = new AliFMDAnalysisTaskSE("TaskFMD");
    mgr->AddTask(taskfmd);
  
    AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
    pars->Init();
    pars->SetProcessPrimary(kTRUE);
    pars->SetProcessHits(kFALSE);
  }


  // Create the task, add it to the manager.
  //===========================================================================
  AliAnalysisTaskFlowEvent *taskFE = NULL;
  if (QA) { 
    if(AddToEvent) { 
      taskFE = new AliAnalysisTaskFlowEvent("TaskFlowEvent",rptype,kTRUE,1);
      taskFE->SetEllipticFlowValue(ellipticFlow); }    //TEST
    else {taskFE = new AliAnalysisTaskFlowEvent("TaskFlowEvent",rptype,kTRUE); }
    taskFE->SetAnalysisType(type);
    taskFE->SetRPType(rptype); //only for ESD
    if (UseMultCut) {
      taskFE->SetMinMult(multmin);
      taskFE->SetMaxMult(multmax);
    }
    taskFE->SetSubeventEtaRange(minA, maxA, minB, maxB);
    if (UsePhysicsSelection) {
      taskFE->SelectCollisionCandidates();
      cout<<"Using Physics Selection"<<endl;
    }
    mgr->AddTask(taskFE);
  }
  else { 
    if(AddToEvent) { 
      taskFE = new AliAnalysisTaskFlowEvent("TaskFlowEvent",rptype,kFALSE,1);
      taskFE->SetEllipticFlowValue(ellipticFlow); }    //TEST
    else {taskFE = new AliAnalysisTaskFlowEvent("TaskFlowEvent",rptype,kFALSE); }
    taskFE->SetAnalysisType(type);
    if (UseMultCut) {
      taskFE->SetMinMult(multmin);
      taskFE->SetMaxMult(multmax);
    }
    taskFE->SetSubeventEtaRange(minA, maxA, minB, maxB);
    if (UsePhysicsSelection) {
      taskFE->SelectCollisionCandidates();
      cout<<"Using Physics Selection"<<endl;
    }
    mgr->AddTask(taskFE);
  }
 
  // Create cuts using the correction framework (CORRFW)
  //===========================================================================
  if (QA){
    //Set TList for the QA histograms
    TList* qaRP  = new TList(); 
    TList* qaPOI = new TList();
  }

  //----------Event cuts----------
  AliCFEventGenCuts* mcEventCuts = new AliCFEventGenCuts("mcEventCuts","MC-level event cuts");
  mcEventCuts->SetNTracksCut(multminESD,multmaxESD); 
  mcEventCuts->SetRequireVtxCuts(requireVtxCuts);
  mcEventCuts->SetVertexXCut(vertexXmin, vertexXmax);
  mcEventCuts->SetVertexYCut(vertexYmin, vertexYmax);
  mcEventCuts->SetVertexZCut(vertexZmin, vertexZmax);
  if (QA) { 
    mcEventCuts->SetQAOn(qaRP);
  }
  AliCFEventRecCuts* recEventCuts = new AliCFEventRecCuts("recEventCuts","rec-level event cuts");
  recEventCuts->SetNTracksCut(multminESD,multmaxESD); 
  recEventCuts->SetRequireVtxCuts(requireVtxCuts);
  recEventCuts->SetVertexXCut(vertexXmin, vertexXmax);
  recEventCuts->SetVertexYCut(vertexYmin, vertexYmax);
  recEventCuts->SetVertexZCut(vertexZmin, vertexZmax);
  recEventCuts->SetVertexNContributors(vertexNContributorsmin,vertexNContributorsmax);
  if (QA) { 
    recEventCuts->SetQAOn(qaRP);
  }
  
  //----------Cuts for RP----------
  //KINEMATICS (MC and reconstructed)
  AliCFTrackKineCuts* mcKineCutsRP = new AliCFTrackKineCuts("mcKineCutsRP","MC-level kinematic cuts");
  mcKineCutsRP->SetPtRange(ptminRP,ptmaxRP);
  mcKineCutsRP->SetEtaRange(etaminRP,etamaxRP);
  //mcKineCutsRP->SetChargeMC(chargeRP);
  mcKineCutsRP->SetRequireIsCharged(isChargedRP);
  if (QA) { 
    mcKineCutsRP->SetQAOn(qaRP);
  }

  AliCFTrackKineCuts *recKineCutsRP = new AliCFTrackKineCuts("recKineCutsRP","rec-level kine cuts");
  recKineCutsRP->SetPtRange(ptminRP,ptmaxRP);
  recKineCutsRP->SetEtaRange(etaminRP,etamaxRP);
  //recKineCutsRP->SetChargeRec(chargeRP);
  recKineCutsRP->SetRequireIsCharged(isChargedRP);
  if (QA) { 
    recKineCutsRP->SetQAOn(qaRP);
  }

  //PID (MC and reconstructed)
  AliCFParticleGenCuts* mcGenCutsRP = new AliCFParticleGenCuts("mcGenCutsRP","MC particle generation cuts for RP");
  mcGenCutsRP->SetRequireIsPrimary();
  if (UsePIDforRP) {mcGenCutsRP->SetRequirePdgCode(PdgRP);}
  if (QA) { 
    mcGenCutsRP->SetQAOn(qaRP);
  }

  int n_species = AliPID::kSPECIES ;
  Double_t* prior = new Double_t[n_species];
  
  prior[0] = 0.0244519 ;
  prior[1] = 0.0143988 ;
  prior[2] = 0.805747  ;
  prior[3] = 0.0928785 ;
  prior[4] = 0.0625243 ;
  
  AliCFTrackCutPid* cutPidRP = NULL;
  if(UsePIDforRP) {
    cutPidRP = new AliCFTrackCutPid("cutPidRP","ESD_PID for RP") ;
    cutPidRP->SetPriors(prior);
    cutPidRP->SetProbabilityCut(0.0);
    cutPidRP->SetDetectors("TPC TOF");
    switch(TMath::Abs(PDG1)) {
    case 11   : cutPidRP->SetParticleType(AliPID::kElectron, kTRUE); break;
    case 13   : cutPidRP->SetParticleType(AliPID::kMuon    , kTRUE); break;
    case 211  : cutPidRP->SetParticleType(AliPID::kPion    , kTRUE); break;
    case 321  : cutPidRP->SetParticleType(AliPID::kKaon    , kTRUE); break;
    case 2212 : cutPidRP->SetParticleType(AliPID::kProton  , kTRUE); break;
    default   : printf("UNDEFINED PID\n"); break;
    }
    if (QA) { 
      cutPidRP->SetQAOn(qaRP); 
    }
  }
  
  //TRACK QUALITY
  AliCFTrackQualityCuts *recQualityCutsRP = new AliCFTrackQualityCuts("recQualityCutsRP","rec-level quality cuts");
  recQualityCutsRP->SetMinNClusterTPC(minClustersTpcRP);
  //recQualityCutsRP->SetMinFoundClusterTPC(minFoundClustersTpcRP); //only for internal TPC QA
  recQualityCutsRP->SetMaxChi2PerClusterTPC(maxChi2PerClusterTpcRP);
  recQualityCutsRP->SetMinNdEdxClusterTPC(minDedxClusterTpcRP);     //to reject secondaries

  recQualityCutsRP->SetMinNClusterITS(minClustersItsRP);
  recQualityCutsRP->SetMaxChi2PerClusterITS(maxChi2PerClusterItsRP);

  recQualityCutsRP->SetMinNClusterTRD(minClustersTrdRP);
  recQualityCutsRP->SetMinNTrackletTRD(minTrackletTrdRP);
  recQualityCutsRP->SetMinNTrackletTRDpid(minTrackletTrdPidRP);
  recQualityCutsRP->SetMaxChi2PerTrackletTRD(maxChi2PerClusterTrdRP);
  recQualityCutsRP->SetStatus(statusRP);  
  if (QA) { 
    recQualityCutsRP->SetQAOn(qaRP);
  }

  /* 
  //How to set this?
  void SetMaxCovDiagonalElements(Float_t c1=1.e+09, Float_t c2=1.e+09, Float_t c3=1.e+09, Float_t c4=1.e+09, Float_t c5=1.e+09)
    {fCovariance11Max=c1;fCovariance22Max=c2;fCovariance33Max=c3;fCovariance44Max=c4;fCovariance55Max=c5;}
  */

  //PRIMARIES
  AliCFTrackIsPrimaryCuts *recIsPrimaryCutsRP = new AliCFTrackIsPrimaryCuts("recIsPrimaryCutsRP","rec-level isPrimary cuts");
  recIsPrimaryCutsRP->UseSPDvertex(spdVertexRP);
  recIsPrimaryCutsRP->UseTPCvertex(tpcVertexRP);
  recIsPrimaryCutsRP->SetMinDCAToVertexXY(minDcaToVertexXyRP); 
  recIsPrimaryCutsRP->SetMinDCAToVertexZ(minDcaToVertexZRP);
  recIsPrimaryCutsRP->SetMaxDCAToVertexXY(maxDcaToVertexXyRP);
  recIsPrimaryCutsRP->SetMaxDCAToVertexZ(maxDcaToVertexZRP); 
  recIsPrimaryCutsRP->SetDCAToVertex2D(dcaToVertex2dRP);
  recIsPrimaryCutsRP->SetAbsDCAToVertex(absDcaToVertexRP);
  recIsPrimaryCutsRP->SetMinNSigmaToVertex(minNSigmaToVertexRP); 
  recIsPrimaryCutsRP->SetMaxNSigmaToVertex(maxNSigmaToVertexRP); 
  recIsPrimaryCutsRP->SetMaxSigmaDCAxy(maxSigmaDcaXySP);
  recIsPrimaryCutsRP->SetMaxSigmaDCAz(maxSigmaDcaZSP);
  recIsPrimaryCutsRP->SetRequireSigmaToVertex(requireSigmaToVertexSP);
  recIsPrimaryCutsRP->SetAcceptKinkDaughters(acceptKinkDaughtersSP);
  if (QA) { 
    recIsPrimaryCutsRP->SetQAOn(qaRP);
  }
  
  //ACCEPTANCE
  AliCFAcceptanceCuts *mcAccCutsRP = new AliCFAcceptanceCuts("mcAccCutsRP","MC acceptance cuts");
  mcAccCutsRP->SetMinNHitITS(minTrackrefsItsRP);
  mcAccCutsRP->SetMinNHitTPC(minTrackrefsTpcRP);
  mcAccCutsRP->SetMinNHitTRD(minTrackrefsTrdRP); 
  mcAccCutsRP->SetMinNHitTOF(minTrackrefsTofRP);
  mcAccCutsRP->SetMinNHitMUON(minTrackrefsMuonRP);
  if (QA) { 
    mcAccCutsRP->SetQAOn(qaRP);
  }

  
  //----------Cuts for POI----------
  //KINEMATICS (MC and reconstructed)
  AliCFTrackKineCuts* mcKineCutsPOI = new AliCFTrackKineCuts("mcKineCutsPOI","MC-level kinematic cuts");
  mcKineCutsPOI->SetPtRange(ptminPOI,ptmaxPOI);
  mcKineCutsPOI->SetEtaRange(etaminPOI,etamaxPOI);
  //mcKineCutsPOI->SetChargeMC(chargePOI);
  mcKineCutsPOI->SetRequireIsCharged(isChargedPOI);
  if (QA) { 
    mcKineCutsPOI->SetQAOn(qaPOI);
  }
  
  AliCFTrackKineCuts *recKineCutsPOI = new AliCFTrackKineCuts("recKineCutsPOI","rec-level kine cuts");
  recKineCutsPOI->SetPtRange(ptminPOI,ptmaxPOI);
  recKineCutsPOI->SetEtaRange(etaminPOI,etamaxPOI);
  //recKineCutsPOI->SetChargeRec(chargePOI);
  recKineCutsPOI->SetRequireIsCharged(isChargedPOI);
  if (QA) { 
    recKineCutsPOI->SetQAOn(qaPOI);
  }
  
  //PID (MC and reconstructed)
  AliCFParticleGenCuts* mcGenCutsPOI = new AliCFParticleGenCuts("mcGenCutsPOI","MC particle generation cuts for POI");
  mcGenCutsPOI->SetRequireIsPrimary();
  if (UsePIDforPOI) {mcGenCutsPOI->SetRequirePdgCode(PdgPOI);}
  if (QA) { 
    mcGenCutsPOI->SetQAOn(qaPOI);
  }

  AliCFTrackCutPid* cutPidPOI = NULL;
  if (UsePIDforPOI) {
    cutPidPOI = new AliCFTrackCutPid("cutPidPOI","ESD_PID for POI") ;
    cutPidPOI->SetPriors(prior);
    cutPidPOI->SetProbabilityCut(0.0);
    cutPidPOI->SetDetectors("TPC TOF");
    switch(TMath::Abs(PDG2)) {
    case 11   : cutPidPOI->SetParticleType(AliPID::kElectron, kTRUE); break;
    case 13   : cutPidPOI->SetParticleType(AliPID::kMuon    , kTRUE); break;
    case 211  : cutPidPOI->SetParticleType(AliPID::kPion    , kTRUE); break;
    case 321  : cutPidPOI->SetParticleType(AliPID::kKaon    , kTRUE); break;
    case 2212 : cutPidPOI->SetParticleType(AliPID::kProton  , kTRUE); break;
    default   : printf("UNDEFINED PID\n"); break;
    }
    if (QA) { 
      cutPidPOI->SetQAOn(qaPOI);
    }
  }

  //TRACK QUALITY
  AliCFTrackQualityCuts *recQualityCutsPOI = new AliCFTrackQualityCuts("recQualityCutsPOI","rec-level quality cuts");
  recQualityCutsPOI->SetMinNClusterTPC(minClustersTpcPOI);
  //recQualityCutsPOI->SetMinFoundClusterTPC(minFoundClustersTpcPOI); //only for internal TPC QA
  recQualityCutsPOI->SetMaxChi2PerClusterTPC(maxChi2PerClusterTpcPOI);
  recQualityCutsPOI->SetMinNdEdxClusterTPC(minDedxClusterTpcPOI);     //to reject secondaries

  recQualityCutsPOI->SetMinNClusterITS(minClustersItsPOI);
  recQualityCutsPOI->SetMaxChi2PerClusterITS(maxChi2PerClusterItsPOI);

  recQualityCutsPOI->SetMinNClusterTRD(minClustersTrdPOI);
  recQualityCutsPOI->SetMinNTrackletTRD(minTrackletTrdPOI);
  recQualityCutsPOI->SetMinNTrackletTRDpid(minTrackletTrdPidPOI);
  recQualityCutsPOI->SetMaxChi2PerTrackletTRD(maxChi2PerClusterTrdPOI);
  recQualityCutsPOI->SetStatus(statusPOI); 
  if (QA) { 
    recQualityCutsPOI->SetQAOn(qaPOI);
  }

  //PRIMARIES
  AliCFTrackIsPrimaryCuts *recIsPrimaryCutsPOI = new AliCFTrackIsPrimaryCuts("recIsPrimaryCutsPOI","rec-level isPrimary cuts");
  recIsPrimaryCutsPOI->UseSPDvertex(spdVertexPOI);
  recIsPrimaryCutsPOI->UseTPCvertex(tpcVertexPOI);
  recIsPrimaryCutsPOI->SetMinDCAToVertexXY(minDcaToVertexXyPOI); 
  recIsPrimaryCutsPOI->SetMinDCAToVertexZ(minDcaToVertexZPOI);
  recIsPrimaryCutsPOI->SetMaxDCAToVertexXY(maxDcaToVertexXyPOI);
  recIsPrimaryCutsPOI->SetMaxDCAToVertexZ(maxDcaToVertexZPOI); 
  recIsPrimaryCutsPOI->SetDCAToVertex2D(dcaToVertex2dPOI);
  recIsPrimaryCutsPOI->SetAbsDCAToVertex(absDcaToVertexPOI);
  recIsPrimaryCutsPOI->SetMinNSigmaToVertex(minNSigmaToVertexPOI); 
  recIsPrimaryCutsPOI->SetMaxNSigmaToVertex(maxNSigmaToVertexPOI); 
  recIsPrimaryCutsPOI->SetMaxSigmaDCAxy(maxSigmaDcaXyPOI);
  recIsPrimaryCutsPOI->SetMaxSigmaDCAz(maxSigmaDcaZPOI);
  recIsPrimaryCutsPOI->SetRequireSigmaToVertex(requireSigmaToVertexPOI);
  recIsPrimaryCutsPOI->SetAcceptKinkDaughters(acceptKinkDaughtersPOI);
  if (QA) { 
    recIsPrimaryCutsPOI->SetQAOn(qaPOI);
  }

  //ACCEPTANCE
  AliCFAcceptanceCuts *mcAccCutsPOI = new AliCFAcceptanceCuts("mcAccCutsPOI","MC acceptance cuts");
  mcAccCutsPOI->SetMinNHitITS(minTrackrefsItsPOI);
  mcAccCutsPOI->SetMinNHitTPC(minTrackrefsTpcPOI);
  mcAccCutsPOI->SetMinNHitTRD(minTrackrefsTrdPOI); 
  mcAccCutsPOI->SetMinNHitTOF(minTrackrefsTofPOI);
  mcAccCutsPOI->SetMinNHitMUON(minTrackrefsMuonPOI);
  if (QA) { 
    mcAccCutsPOI->SetQAOn(qaPOI);
  }

     
  
  //----------Create Cut Lists----------
  printf("CREATE EVENT CUTS\n");
  TObjArray* mcEventList = new TObjArray(0);  
  if (UseMultCutforESD) mcEventList->AddLast(mcEventCuts);//cut on mult and vertex
    
  TObjArray* recEventList = new TObjArray(0);
  if (UseMultCutforESD) recEventList->AddLast(recEventCuts);//cut on mult and vertex

  printf("CREATE MC KINE CUTS\n");
  TObjArray* mcListRP = new TObjArray(0);
  if (UseKineforRP) mcListRP->AddLast(mcKineCutsRP); //cut on pt/eta/phi
  mcListRP->AddLast(mcGenCutsRP); //cut on primary and if (UsePIDforRP) MC PID
  
  TObjArray* mcListPOI = new TObjArray(0);
  if (UseKineforPOI) mcListPOI->AddLast(mcKineCutsPOI); //cut on pt/eta/phi
  mcListPOI->AddLast(mcGenCutsPOI); //cut on primary and if (UsePIDforPOI) MC PID
  
  printf("CREATE MC ACCEPTANCE CUTS\n");
  TObjArray* accListRP = new TObjArray(0) ;
  if (UseAcceptanceforRP) accListRP->AddLast(mcAccCutsRP); //cut on number of track references
  
  TObjArray* accListPOI = new TObjArray(0) ;
  if (UseAcceptanceforPOI) accListPOI->AddLast(mcAccCutsPOI); //cut on number of track references
  
  printf("CREATE ESD RECONSTRUCTION CUTS\n");
  TObjArray* recListRP = new TObjArray(0) ;
  if (UseTrackQualityforRP) recListRP->AddLast(recQualityCutsRP);   //track quality
  if (UsePrimariesforRP)    recListRP->AddLast(recIsPrimaryCutsRP); //cut if it is a primary
  if (UseKineforRP)         recListRP->AddLast(recKineCutsRP);      //cut on pt/eta/phi  

  TObjArray* recListPOI = new TObjArray(0) ;
  if (UseTrackQualityforPOI) recListPOI->AddLast(recQualityCutsPOI);   //track quality
  if (UsePrimariesforPOI)    recListPOI->AddLast(recIsPrimaryCutsPOI); //cut if it is a primary
  if (UseKineforPOI)         recListPOI->AddLast(recKineCutsPOI);      //cut on pt/eta/phi

  printf("CREATE ESD PID CUTS\n");
  TObjArray* fPIDCutListRP = new TObjArray(0) ;
  if(UsePIDforRP) {fPIDCutListRP->AddLast(cutPidRP);} //cut on ESD PID
  
  TObjArray* fPIDCutListPOI = new TObjArray(0) ;
  if (UsePIDforPOI)  {fPIDCutListPOI->AddLast(cutPidPOI);} //cut on ESD PID
  

  //----------Add Cut Lists to the CF Manager----------
  printf("CREATE INTERFACE AND CUTS\n");
  AliCFManager* cfmgrRP = new AliCFManager();
  cfmgrRP->SetNStepEvent(3);
  cfmgrRP->SetEventCutsList(AliCFManager::kEvtGenCuts,mcEventList); 
  cfmgrRP->SetEventCutsList(AliCFManager::kEvtRecCuts,recEventList); 
  cfmgrRP->SetNStepParticle(4); 
  cfmgrRP->SetParticleCutsList(AliCFManager::kPartGenCuts,mcListRP);
  cfmgrRP->SetParticleCutsList(AliCFManager::kPartAccCuts,accListRP);
  cfmgrRP->SetParticleCutsList(AliCFManager::kPartRecCuts,recListRP);
  cfmgrRP->SetParticleCutsList(AliCFManager::kPartSelCuts,fPIDCutListRP);
  
  AliCFManager* cfmgrPOI = new AliCFManager();
  cfmgrPOI->SetNStepEvent(3);
  cfmgrPOI->SetEventCutsList(AliCFManager::kEvtGenCuts,mcEventList); 
  cfmgrPOI->SetEventCutsList(AliCFManager::kEvtRecCuts,recEventList); 
  cfmgrPOI->SetNStepParticle(4); 
  cfmgrPOI->SetParticleCutsList(AliCFManager::kPartGenCuts,mcListPOI);
  cfmgrPOI->SetParticleCutsList(AliCFManager::kPartAccCuts,accListPOI);
  cfmgrPOI->SetParticleCutsList(AliCFManager::kPartRecCuts,recListPOI);
  cfmgrPOI->SetParticleCutsList(AliCFManager::kPartSelCuts,fPIDCutListPOI);
  
  if (QA) {
    taskFE->SetQAList1(qaRP);
    taskFE->SetQAList2(qaPOI);
  }
  taskFE->SetCFManager1(cfmgrRP);
  taskFE->SetCFManager2(cfmgrPOI);



  // Create the analysis tasks, add them to the manager.
  //===========================================================================
  if (SP){
    AliAnalysisTaskScalarProduct *taskSP = new AliAnalysisTaskScalarProduct("TaskScalarProduct",WEIGHTS[0]);
    taskSP->SetRelDiffMsub(0.1);
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
    taskQC->SetnBinsMult(10000);
    taskQC->SetMinMult(0.);
    taskQC->SetMaxMult(10000.);
    mgr->AddTask(taskQC);
  }
  if (FQD){
    AliAnalysisTaskFittingQDistribution *taskFQD = new AliAnalysisTaskFittingQDistribution("TaskFittingQDistribution",kFALSE);
    taskFQD->SetUsePhiWeights(WEIGHTS[0]); 
    taskFQD->SetqMin(0.);
    taskFQD->SetqMax(1000.);
    taskFQD->SetqNbins(10000);
    mgr->AddTask(taskFQD);
  }
  if (MCEP){
    AliAnalysisTaskMCEventPlane *taskMCEP = new AliAnalysisTaskMCEventPlane("TaskMCEventPlane");
    mgr->AddTask(taskMCEP);
  }
  if (MH){
    AliAnalysisTaskMixedHarmonics *taskMH = new AliAnalysisTaskMixedHarmonics("TaskMixedHarmonics",useWeights);
    taskMH->SetCorrelatorInteger(1);
    taskMH->SetNoOfMultipicityBins(10);
    taskMH->SetMultipicityBinWidth(2);
    taskMH->SetMinMultiplicity(3);
    taskMH->SetCorrectForDetectorEffects(kTRUE);
    //taskMH->SetUsePhiWeights(WEIGHTS[0]); 
    //taskMH->SetUsePtWeights(WEIGHTS[1]);
    //taskMH->SetUseEtaWeights(WEIGHTS[2]); 
    mgr->AddTask(taskMH);
  }  
  if (NL){
    AliAnalysisTaskNestedLoops *taskNL = new AliAnalysisTaskNestedLoops("TaskNestedLoops",useWeights);
    //taskNL->SetUsePhiWeights(WEIGHTS[0]); 
    //taskNL->SetUsePtWeights(WEIGHTS[1]);
    //taskNL->SetUseEtaWeights(WEIGHTS[2]); 
    mgr->AddTask(taskNL);
  }
  
  // Create the output container for the data produced by the task
  // Connect to the input and output containers
  //===========================================================================
  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
  
  if (rptype == "FMD") {
    AliAnalysisDataContainer *coutputFMD = 
      mgr->CreateContainer("BackgroundCorrected", TList::Class(), AliAnalysisManager::kExchangeContainer);                        
    //input and output taskFMD     
    mgr->ConnectInput(taskfmd, 0, cinput1);
    mgr->ConnectOutput(taskfmd, 1, coutputFMD);
    //input into taskFE
    mgr->ConnectInput(taskFE,1,coutputFMD);
  }

  AliAnalysisDataContainer *coutputFE = 
    mgr->CreateContainer("cobjFlowEventSimple",  AliFlowEventSimple::Class(),AliAnalysisManager::kExchangeContainer);
  mgr->ConnectInput(taskFE,0,cinput1); 
  mgr->ConnectOutput(taskFE,1,coutputFE);

  if (QA) { 
    TString qaNameRPFE = AliAnalysisManager::GetCommonFileName();
    qaNameRPFE += ":QAforRP_FE_";
    qaNameRPFE += type;

    AliAnalysisDataContainer *coutputQA1FE =
      mgr->CreateContainer("QARPFE", TList::Class(),AliAnalysisManager::kOutputContainer,qaNameRPFE); 
    
    TString qaNamePOIFE = AliAnalysisManager::GetCommonFileName();
    qaNamePOIFE += ":QAforPOI_FE_";
    qaNamePOIFE += type;
        
    AliAnalysisDataContainer *coutputQA2FE =
      mgr->CreateContainer("QAPOIFE", TList::Class(),AliAnalysisManager::kOutputContainer,qaNamePOIFE); 

    mgr->ConnectOutput(taskFE,2,coutputQA1FE); 
    mgr->ConnectOutput(taskFE,3,coutputQA2FE); 
  }

  // Create the output containers for the data produced by the analysis tasks
  // Connect to the input and output containers
  //===========================================================================
  if (useWeights) {    
    AliAnalysisDataContainer *cinputWeights = mgr->CreateContainer("cobjWeights",TList::Class(),AliAnalysisManager::kInputContainer); 
  }

  if(SP) {
    TString outputSP = AliAnalysisManager::GetCommonFileName();
    outputSP += ":outputSPanalysis";
    outputSP+= type;
    
    AliAnalysisDataContainer *coutputSP = mgr->CreateContainer("cobjSP", TList::Class(),AliAnalysisManager::kOutputContainer,outputSP); 
    mgr->ConnectInput(taskSP,0,coutputFE); 
    mgr->ConnectOutput(taskSP,1,coutputSP); 
    if (WEIGHTS[0]) {
      mgr->ConnectInput(taskSP,1,cinputWeights);
      cinputWeights->SetData(weightsList);
    } 
  }
  if(LYZ1SUM) {
    TString outputLYZ1SUM = AliAnalysisManager::GetCommonFileName();
    outputLYZ1SUM += ":outputLYZ1SUManalysis";
    outputLYZ1SUM+= type;
    
    AliAnalysisDataContainer *coutputLYZ1SUM = mgr->CreateContainer("cobjLYZ1SUM", TList::Class(),AliAnalysisManager::kOutputContainer,outputLYZ1SUM); 
    mgr->ConnectInput(taskLYZ1SUM,0,coutputFE); 
    mgr->ConnectOutput(taskLYZ1SUM,1,coutputLYZ1SUM); 
  }
  if(LYZ1PROD) {
    TString outputLYZ1PROD = AliAnalysisManager::GetCommonFileName();
    outputLYZ1PROD += ":outputLYZ1PRODanalysis";
    outputLYZ1PROD+= type;
    
    AliAnalysisDataContainer *coutputLYZ1PROD = mgr->CreateContainer("cobjLYZ1PROD", TList::Class(),AliAnalysisManager::kOutputContainer,outputLYZ1PROD); 
    mgr->ConnectInput(taskLYZ1PROD,0,coutputFE); 
    mgr->ConnectOutput(taskLYZ1PROD,1,coutputLYZ1PROD);
  }
  if(LYZ2SUM) {
    AliAnalysisDataContainer *cinputLYZ2SUM = mgr->CreateContainer("cobjLYZ2SUMin",TList::Class(),AliAnalysisManager::kInputContainer);
    TString outputLYZ2SUM = AliAnalysisManager::GetCommonFileName();
    outputLYZ2SUM += ":outputLYZ2SUManalysis";
    outputLYZ2SUM+= type;
    
    AliAnalysisDataContainer *coutputLYZ2SUM = mgr->CreateContainer("cobjLYZ2SUM", TList::Class(),AliAnalysisManager::kOutputContainer,outputLYZ2SUM); 
    mgr->ConnectInput(taskLYZ2SUM,0,coutputFE); 
    mgr->ConnectInput(taskLYZ2SUM,1,cinputLYZ2SUM);
    mgr->ConnectOutput(taskLYZ2SUM,1,coutputLYZ2SUM); 
    cinputLYZ2SUM->SetData(fInputListLYZ2SUM);
  }
  if(LYZ2PROD) {
    AliAnalysisDataContainer *cinputLYZ2PROD = mgr->CreateContainer("cobjLYZ2PRODin",TList::Class(),AliAnalysisManager::kInputContainer);
    TString outputLYZ2PROD = AliAnalysisManager::GetCommonFileName();
    outputLYZ2PROD += ":outputLYZ2PRODanalysis";
    outputLYZ2PROD+= type;
    
    AliAnalysisDataContainer *coutputLYZ2PROD = mgr->CreateContainer("cobjLYZ2PROD", TList::Class(),AliAnalysisManager::kOutputContainer,outputLYZ2PROD); 
    mgr->ConnectInput(taskLYZ2PROD,0,coutputFE); 
    mgr->ConnectInput(taskLYZ2PROD,1,cinputLYZ2PROD);
    mgr->ConnectOutput(taskLYZ2PROD,1,coutputLYZ2PROD); 
    cinputLYZ2PROD->SetData(fInputListLYZ2PROD);
  }
  if(LYZEP) {
    AliAnalysisDataContainer *cinputLYZEP = mgr->CreateContainer("cobjLYZEPin",TList::Class(),AliAnalysisManager::kInputContainer);
    TString outputLYZEP = AliAnalysisManager::GetCommonFileName();
    outputLYZEP += ":outputLYZEPanalysis";
    outputLYZEP+= type;
    
    AliAnalysisDataContainer *coutputLYZEP = mgr->CreateContainer("cobjLYZEP", TList::Class(),AliAnalysisManager::kOutputContainer,outputLYZEP); 
    mgr->ConnectInput(taskLYZEP,0,coutputFE); 
    mgr->ConnectInput(taskLYZEP,1,cinputLYZEP);
    mgr->ConnectOutput(taskLYZEP,1,coutputLYZEP); 
    cinputLYZEP->SetData(fInputListLYZEP);
  }
  if(GFC) {
    TString outputGFC = AliAnalysisManager::GetCommonFileName();
    outputGFC += ":outputGFCanalysis";
    outputGFC+= type;
    
    AliAnalysisDataContainer *coutputGFC = mgr->CreateContainer("cobjGFC", TList::Class(),AliAnalysisManager::kOutputContainer,outputGFC); 
    mgr->ConnectInput(taskGFC,0,coutputFE); 
    mgr->ConnectOutput(taskGFC,1,coutputGFC);
    if (useWeights) {
      mgr->ConnectInput(taskGFC,1,cinputWeights);
      cinputWeights->SetData(weightsList);
    } 
  }
  if(QC) {
    TString outputQC = AliAnalysisManager::GetCommonFileName();
    outputQC += ":outputQCanalysis";
    outputQC+= type;

    AliAnalysisDataContainer *coutputQC = mgr->CreateContainer("cobjQC", TList::Class(),AliAnalysisManager::kOutputContainer,outputQC); 
    mgr->ConnectInput(taskQC,0,coutputFE); 
    mgr->ConnectOutput(taskQC,1,coutputQC);
    if (useWeights) {
      mgr->ConnectInput(taskQC,1,cinputWeights);
      cinputWeights->SetData(weightsList);
    }    
  }
  if(FQD) {
    TString outputFQD = AliAnalysisManager::GetCommonFileName();
    outputFQD += ":outputFQDanalysis";
    outputFQD+= type;
    
    AliAnalysisDataContainer *coutputFQD = mgr->CreateContainer("cobjFQD", TList::Class(),AliAnalysisManager::kOutputContainer,outputFQD); 
    mgr->ConnectInput(taskFQD,0,coutputFE); 
    mgr->ConnectOutput(taskFQD,1,coutputFQD);
    if(useWeights) {
      mgr->ConnectInput(taskFQD,1,cinputWeights);
      cinputWeights->SetData(weightsList);
    } 
  }
  if(MCEP) {
    TString outputMCEP = AliAnalysisManager::GetCommonFileName();
    outputMCEP += ":outputMCEPanalysis";
    outputMCEP+= type;
    
    AliAnalysisDataContainer *coutputMCEP = mgr->CreateContainer("cobjMCEP", TList::Class(),AliAnalysisManager::kOutputContainer,outputMCEP); 
    mgr->ConnectInput(taskMCEP,0,coutputFE); 
    mgr->ConnectOutput(taskMCEP,1,coutputMCEP); 
  }
  if(MH) {
    TString outputMH = AliAnalysisManager::GetCommonFileName();
    outputMH += ":outputMHanalysis";
    outputMH += type;
        
    AliAnalysisDataContainer *coutputMH = mgr->CreateContainer("cobjMH", TList::Class(),AliAnalysisManager::kOutputContainer,outputMH); 
    mgr->ConnectInput(taskMH,0,coutputFE); 
    mgr->ConnectOutput(taskMH,1,coutputMH); 
    //if (useWeights) {
    //  mgr->ConnectInput(taskMH,1,cinputWeights);
    //  cinputWeights->SetData(weightsList);
    //} 
  }
  if(NL) {
    TString outputNL = AliAnalysisManager::GetCommonFileName();
    outputNL += ":outputNLanalysis";
    outputNL += type;
        
    AliAnalysisDataContainer *coutputNL = mgr->CreateContainer("cobjNL", TList::Class(),AliAnalysisManager::kOutputContainer,outputNL); 
    mgr->ConnectInput(taskNL,0,coutputFE); 
    mgr->ConnectOutput(taskNL,1,coutputNL); 
    //if (useWeights) {
    //  mgr->ConnectInput(taskNL,1,cinputWeights);
    //  cinputWeights->SetData(weightsList);
    //} 
  }

  // Return analysis task
  //===========================================================================
  return taskFE;
  


}





