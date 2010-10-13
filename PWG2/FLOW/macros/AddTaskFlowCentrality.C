/////////////////////////////////////////////////////////////////////////////////////////////
//
// AddTask* macro for flow analysis
// Creates a Flow Event task and adds it to the analysis manager.
// Sets the cuts using the correction framework (CORRFW) classes.
// Also creates Flow Analysis tasks and connects them to the output of the flow event task.
//
/////////////////////////////////////////////////////////////////////////////////////////////

// Define the range for eta subevents (for SP method)
//-----(FMD 1.7 - 5.0)-----
//Double_t minA = -5.0;
//Double_t maxA = -1.7;
//Double_t minB = 1.7;
//Double_t maxB = 5.0;
//-----(Tracklets 0.9 - 2.0)-----
//Double_t minA = -2.0;
//Double_t maxA = -0.9;
//Double_t minB = 0.9;
//Double_t maxB = 2.0;
//-----(Global 0.5 - 0.9)-----
Double_t minA = -0.9;
Double_t maxA = -0.5;
Double_t minB = 0.5;
Double_t maxB = 0.9;

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
const Double_t vertexZmin = -10.; 
const Double_t vertexZmax = 10.; 
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
Bool_t UseKineforRP =  kFALSE;
const Double_t ptminRP = 0.0;
const Double_t ptmaxRP = 10.0;
const Double_t etaminRP  = -0.9;
const Double_t etamaxRP  = 0.9;
Bool_t bSetChargeForRP = kTRUE; // by setting kTRUE you will use only RPs with charge "chargeRP" specified in the line bellow 
const Int_t  chargeRP = 1;  
const Bool_t isChargedRP = kTRUE; // take only charged particle in the analysis

//PID (on generated and reconstructed tracks)
Bool_t bUsePIDforRP = kFALSE; // by setting kTRUE you will use only RPs with PDG code "PdgRP" specified in the line bellow
const Int_t PdgRP = 211; //  pion = 211, kaon = 321, proton = 2212

//TRACK QUALITY (on reconstructed tracks only)
//see /CORRFW/AliCFTrackQualityCuts class
Bool_t UseTrackQualityforRP =  kFALSE;
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
Bool_t UsePrimariesforRP = kFALSE;
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
const Double_t maxSigmaDcaXyRP = 1.e+10;
const Double_t maxSigmaDcaZRP = 1.e+10;
const Bool_t   requireSigmaToVertexRP = kFALSE;
const Bool_t   acceptKinkDaughtersRP = kFALSE;  //default = kTRUE;

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
const Double_t etaminPOI  = -0.9;
const Double_t etamaxPOI  = 0.9;
Bool_t bSetChargeForPOI = kTRUE; // by setting kTRUE you will use only RPs with charge "chargePOI" specified in the line bellow 
const Int_t    chargePOI = -1; 
const Bool_t   isChargedPOI = kTRUE;

//PID (on generated and reconstructed tracks)
Bool_t bUsePIDforPOI = kFALSE; // by setting kTRUE you will use only POIs with PDG code "PdgPOI" specified in the line bellow
const Int_t PdgPOI = 321; //  pion = 211, kaon = 321, proton = 2212

//TRACK QUALITY (on reconstructed tracks only)
//see /CORRFW/AliCFTrackQualityCuts class
Bool_t UseTrackQualityforPOI = kFALSE;
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

void AddTaskFlowCentrality( TString type,
                            Bool_t* METHODS,
                            Bool_t QA,
                            Bool_t* WEIGHTS,
                            Int_t refMultMin=0,
                            Int_t refMultMax=1e10,
                            TString fileName="AnalysisResults.root" )
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
  if (rptype == "FMD") {
    AliFMDAnalysisTaskSE *taskfmd = NULL;
    if (rptype == "FMD") {
      taskfmd = new AliFMDAnalysisTaskSE("TaskFMD");
      mgr->AddTask(taskfmd);
  
      AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
      pars->Init();
      pars->SetProcessPrimary(kTRUE); //for MC only
      pars->SetProcessHits(kFALSE);

      //pars->SetRealData(kTRUE); //for real data
      //pars->SetProcessPrimary(kFALSE); //for real data

    }
  }
  

  // Create the task, add it to the manager.
  //===========================================================================
  AliAnalysisTaskFlowEvent *taskFE = NULL;
  if (QA) { 
    if(useAfterBurner)
    { 
      taskFE = new AliAnalysisTaskFlowEvent("TaskFlowEvent",rptype,kTRUE,1);
      taskFE->SetFlow(v1,v2,v3,v4); 
      taskFE->SetNonFlowNumberOfTrackClones(numberOfTrackClones);
      taskFE->SetAfterburnerOn();
    }
    else {taskFE = new AliAnalysisTaskFlowEvent("TaskFlowEvent",rptype,kTRUE); }
    taskFE->SetAnalysisType(type);
    taskFE->SetRPType(rptype); //only for ESD
    if (UseMultCut) {
      taskFE->SetMinMult(multmin);
      taskFE->SetMaxMult(multmax);
    }
    if (ExcludeRegion) {
      taskFE->DefineDeadZone(excludeEtaMin, excludeEtaMax, excludePhiMin, excludePhiMax); 
    }

    taskFE->SetSubeventEtaRange(minA, maxA, minB, maxB);
    if (UsePhysicsSelection) {
      taskFE->SelectCollisionCandidates();
      cout<<"Using Physics Selection"<<endl;
    }
    mgr->AddTask(taskFE);
  }
  else { 
    if(useAfterBurner)
    { 
      taskFE = new AliAnalysisTaskFlowEvent("TaskFlowEvent",rptype,kFALSE,1);
      taskFE->SetFlow(v1,v2,v3,v4); 
      taskFE->SetNonFlowNumberOfTrackClones(numberOfTrackClones);
      taskFE->SetAfterburnerOn();
    }
    else {taskFE = new AliAnalysisTaskFlowEvent("TaskFlowEvent",rptype,kFALSE); }
    taskFE->SetAnalysisType(type);
    taskFE->SetRPType(rptype); //only for ESD
    if (UseMultCut) {
      taskFE->SetMinMult(multmin);
      taskFE->SetMaxMult(multmax);
    }
    if (ExcludeRegion) {
      taskFE->DefineDeadZone(excludeEtaMin, excludeEtaMax, excludePhiMin, excludePhiMax); 
    }
    taskFE->SetSubeventEtaRange(minA, maxA, minB, maxB);
    if (UsePhysicsSelection) {
      taskFE->SelectCollisionCandidates();
      cout<<"Using Physics Selection"<<endl;
    }
    mgr->AddTask(taskFE);
  }
 
  //----------Add Cut Lists to the CF Manager----------
  printf("CREATE INTERFACE AND CUTS\n");
  
  // EVENTS CUTS:
  AliFlowEventCuts* eventCuts = new AliFlowEventCuts();
  eventCuts->SetRefMultRange(refMultMin,refMultMax);
  taskFE->SetCutsEvent(eventCuts);
  
  // RP CUTS:
  AliFlowTrackCuts* cutsRP = new AliFlowTrackCuts();
  cutsRP->SetPtRange(ptminRP,ptmaxRP);
  cutsRP->SetEtaRange(etaminRP,etamaxRP);
  cutsRP->SetRequireCharge(kTRUE);
  if(bSetChargeForRP){cutsRP->SetCharge(chargeRP);} 
  if(bUsePIDforRP){cutsRP->SetPID(PdgRP);}
      
  if(UseTrackQualityforRP)
  {
   cutsRP->SetMinNClustersTPC(minClustersTpcRP);
   cutsRP->SetMaxChi2PerClusterTPC(maxChi2PerClusterTpcRP);
   cutsRP->SetMinNClustersITS(minClustersItsRP);
   cutsRP->SetMaxChi2PerClusterITS(maxChi2PerClusterItsRP);
   // when compared to original Naomi's implementation we still might need setters for:
   /*
   const UShort_t minDedxClusterTpcRP = 0;
   const Int_t    minClustersTrdRP = -1;
   const Int_t    minTrackletTrdRP = -1;
   const Int_t    minTrackletTrdPidRP = -1;
   const Double_t maxChi2PerClusterTrdRP = 1.e+09;
   const ULong_t  statusRP = AliESDtrack::kTPCrefit;   //AliESDtrack::kTPCrefit &  AliESDtrack::kITSrefit 
   */
  } // end of if(UseTrackQualityforRP)
 
  if(UsePrimariesforRP)
  {
   cutsRP->SetMaxDCAToVertexXY(maxDcaToVertexXyRP);
   cutsRP->SetMaxDCAToVertexZ(maxDcaToVertexZRP);
   cutsRP->SetDCAToVertex2D(dcaToVertex2dRP);
   cutsRP->SetMaxNsigmaToVertex(maxNSigmaToVertexRP);
   cutsRP->SetRequireSigmaToVertex(requireSigmaToVertexRP);
   cutsRP->SetAcceptKinkDaughters(acceptKinkDaughtersRP);
   // when compared to original Naomi's implementation we still might need setters for:
   /*
   const Bool_t   spdVertexRP = kFALSE;
   const Bool_t   tpcVertexRP = kFALSE;
   const Float_t  minDcaToVertexXyRP = 0.;
   const Float_t  minDcaToVertexZRP = 0.;
   const Bool_t   absDcaToVertexRP = kTRUE;  //default = kTRUE;
   const Double_t minNSigmaToVertexRP = 0.;
   const Double_t maxSigmaDcaXyRP = 1.e+10;
   const Double_t maxSigmaDcaZRP = 1.e+10;  
   */
  } // end of if(UsePrimariesforRP)
  
  if(UseAcceptanceforRP)
  {
   // when compared to original Naomi's implementation we still might need setters for:
   /*
   const Int_t minTrackrefsItsRP = 3;
   const Int_t minTrackrefsTpcRP = 2;
   const Int_t minTrackrefsTrdRP = 0; 
   const Int_t minTrackrefsTofRP = 0; 
   const Int_t minTrackrefsMuonRP = 0; 
   */
  }
    
  // POI CUTS:
  AliFlowTrackCuts* cutsPOI = new AliFlowTrackCuts();
  cutsPOI->SetPtRange(ptminPOI,ptmaxPOI);
  cutsPOI->SetEtaRange(etaminPOI,etamaxPOI);
  cutsPOI->SetRequireCharge(kTRUE);
  if(bSetChargeForPOI){cutsPOI->SetCharge(chargePOI);} 
  if(bUsePIDforPOI){cutsPOI->SetPID(PdgPOI);}
      
  if(UseTrackQualityforPOI)
  {
   cutsPOI->SetMinNClustersTPC(minClustersTpcPOI);
   cutsPOI->SetMaxChi2PerClusterTPC(maxChi2PerClusterTpcPOI);
   cutsPOI->SetMinNClustersITS(minClustersItsPOI);
   cutsPOI->SetMaxChi2PerClusterITS(maxChi2PerClusterItsPOI);
   // when compared to original Naomi's implementation we still might need setters for:
   /*
   const UShort_t minDedxClusterTpcPOI = 0;
   const Int_t    minClustersTrdPOI = -1;
   const Int_t    minTrackletTrdPOI = -1;
   const Int_t    minTrackletTrdPidPOI = -1;
   const Double_t maxChi2PerClusterTrdPOI = 1.e+09;
   const ULong_t  statusPOI = AliESDtrack::kTPCrefit;   //AliESDtrack::kTPCrefit &  AliESDtrack::kITSrefit 
   */
  } // end of if(UseTrackQualityforPOI)
 
  if(UsePrimariesforPOI)
  {
   cutsPOI->SetMaxDCAToVertexXY(maxDcaToVertexXyPOI);
   cutsPOI->SetMaxDCAToVertexZ(maxDcaToVertexZPOI);
   cutsPOI->SetDCAToVertex2D(dcaToVertex2dPOI);
   cutsPOI->SetMaxNsigmaToVertex(maxNSigmaToVertexPOI);
   cutsPOI->SetRequireSigmaToVertex(requireSigmaToVertexPOI);
   cutsPOI->SetAcceptKinkDaughters(acceptKinkDaughtersPOI);
   // when compared to original Naomi's implementation we still might need setters for:
   /*
   const Bool_t   spdVertexPOI = kFALSE;
   const Bool_t   tpcVertexPOI = kFALSE;
   const Float_t  minDcaToVertexXyPOI = 0.;
   const Float_t  minDcaToVertexZPOI = 0.;
   const Bool_t   absDcaToVertexPOI = kTRUE;  //default = kTRUE;
   const Double_t minNSigmaToVertexPOI = 0.;
   const Double_t maxSigmaDcaXyPOI = 1.e+10;
   const Double_t maxSigmaDcaZPOI = 1.e+10;  
   */
  } // end of if(UsePrimariesforPOI)

  if(UseAcceptanceforRP)
  {
   // when compared to original Naomi's implementation we still might need setters for:
   /*
   const Int_t minTrackrefsItsRP = 3;
   const Int_t minTrackrefsTpcRP = 2;
   const Int_t minTrackrefsTrdRP = 0; 
   const Int_t minTrackrefsTofRP = 0; 
   const Int_t minTrackrefsMuonRP = 0; 
   */
  }
  
  if (QA) {
    taskFE->SetQAList1(new TList());
    taskFE->SetQAList2(new TList());
  }
  
  // Pass cuts for RPs and POIs to the task:
  taskFE->SetCutsRP(cutsRP);
  taskFE->SetCutsPOI(cutsPOI);

  // Create the analysis tasks, add them to the manager.
  //===========================================================================
  if (SP){
    AliAnalysisTaskScalarProduct *taskSP = new AliAnalysisTaskScalarProduct("TaskScalarProduct",WEIGHTS[0]);
    taskSP->SetRelDiffMsub(1.0);
    taskSP->SetApplyCorrectionForNUA(kFALSE);
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
    taskMH->SetHarmonic(1); // n in cos[n(phi1+phi2-2phi3)] and cos[n(psi1+psi2-2phi3)]
    taskMH->SetNoOfMultipicityBins(10);
    taskMH->SetMultipicityBinWidth(2);
    taskMH->SetMinMultiplicity(3);
    taskMH->SetCorrectForDetectorEffects(kTRUE);
    taskMH->SetEvaluateDifferential3pCorrelator(kFALSE); // evaluate <<cos[n(psi1+psi2-2phi3)]>> (Remark: two nested loops)    
    taskMH->SetOppositeChargesPOI(kFALSE); // POIs psi1 and psi2 in cos[n(psi1+psi2-2phi3)] will have opposite charges  
    mgr->AddTask(taskMH);
  }  
  if (NL){
    AliAnalysisTaskNestedLoops *taskNL = new AliAnalysisTaskNestedLoops("TaskNestedLoops",useWeights);
    taskNL->SetHarmonic(1); // n in cos[n(phi1+phi2-2phi3)] and cos[n(psi1+psi2-2phi3)]
    taskNL->SetEvaluateNestedLoopsForRAD(kTRUE); // RAD = Relative Angle Distribution
    taskNL->SetEvaluateNestedLoopsForMH(kTRUE); // evalaute <<cos[n(phi1+phi2-2phi3)]>> (Remark: three nested loops)   
    taskNL->SetEvaluateDifferential3pCorrelator(kFALSE); // evaluate <<cos[n(psi1+psi2-2phi3)]>>  (Remark: three nested loops)   
    taskNL->SetOppositeChargesPOI(kFALSE); // POIs psi1 and psi2 in cos[n(psi1+psi2-2phi3)] will have opposite charges  
    mgr->AddTask(taskNL);
  }
  
  // Create the output container for the data produced by the task
  // Connect to the input and output containers
  //===========================================================================
  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
  
  if (rptype == "FMD") {
    AliAnalysisDataContainer *coutputFMD = 
      mgr->CreateContainer(Form("BackgroundCorrected_%s",fileName.Data()), TList::Class(), AliAnalysisManager::kExchangeContainer);                        
    //input and output taskFMD     
    mgr->ConnectInput(taskfmd, 0, cinput1);
    mgr->ConnectOutput(taskfmd, 1, coutputFMD);
    //input into taskFE
    mgr->ConnectInput(taskFE,1,coutputFMD);
  }
  
  AliAnalysisDataContainer *coutputFE = 
    mgr->CreateContainer(Form("cobjFlowEventSimple_%s",fileName.Data()),  AliFlowEventSimple::Class(),AliAnalysisManager::kExchangeContainer);
  mgr->ConnectInput(taskFE,0,cinput1); 
  mgr->ConnectOutput(taskFE,1,coutputFE);

  if (QA) { 
    TString qaNameRPFE = fileName;
    qaNameRPFE += ":QAforRP_FE_";
    qaNameRPFE += type;

    AliAnalysisDataContainer *coutputQA1FE =
      mgr->CreateContainer(Form("QARPFE_%s",fileName.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,qaNameRPFE); 
    
    TString qaNamePOIFE = fileName;
    qaNamePOIFE += ":QAforPOI_FE_";
    qaNamePOIFE += type;
        
    AliAnalysisDataContainer *coutputQA2FE =
      mgr->CreateContainer(Form("QAPOIFE_%s",fileName.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,qaNamePOIFE); 

    mgr->ConnectOutput(taskFE,2,coutputQA1FE); 
    mgr->ConnectOutput(taskFE,3,coutputQA2FE); 
  }

  // Create the output containers for the data produced by the analysis tasks
  // Connect to the input and output containers
  //===========================================================================
  if (useWeights) {    
    AliAnalysisDataContainer *cinputWeights = mgr->CreateContainer(Form("cobjWeights_%s",fileName.Data()),TList::Class(),AliAnalysisManager::kInputContainer); 
  }

  if(SP) {
    TString outputSP = fileName;
    outputSP += ":outputSPanalysis";
    outputSP+= type;
    
    AliAnalysisDataContainer *coutputSP = mgr->CreateContainer(Form("cobjSP_%s",fileName.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputSP); 
    mgr->ConnectInput(taskSP,0,coutputFE); 
    mgr->ConnectOutput(taskSP,1,coutputSP); 
    if (WEIGHTS[0]) {
      mgr->ConnectInput(taskSP,1,cinputWeights);
      cinputWeights->SetData(weightsList);
    } 
  }
  if(LYZ1SUM) {
    TString outputLYZ1SUM = fileName;
    outputLYZ1SUM += ":outputLYZ1SUManalysis";
    outputLYZ1SUM+= type;
    
    AliAnalysisDataContainer *coutputLYZ1SUM = mgr->CreateContainer(Form("cobjLYZ1SUM_%s",fileName.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputLYZ1SUM); 
    mgr->ConnectInput(taskLYZ1SUM,0,coutputFE); 
    mgr->ConnectOutput(taskLYZ1SUM,1,coutputLYZ1SUM); 
  }
  if(LYZ1PROD) {
    TString outputLYZ1PROD = fileName;
    outputLYZ1PROD += ":outputLYZ1PRODanalysis";
    outputLYZ1PROD+= type;
    
    AliAnalysisDataContainer *coutputLYZ1PROD = mgr->CreateContainer(Form("cobjLYZ1PROD_%s",fileName.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputLYZ1PROD); 
    mgr->ConnectInput(taskLYZ1PROD,0,coutputFE); 
    mgr->ConnectOutput(taskLYZ1PROD,1,coutputLYZ1PROD);
  }
  if(LYZ2SUM) {
    AliAnalysisDataContainer *cinputLYZ2SUM = mgr->CreateContainer(Form("cobjLYZ2SUMin_%s",fileName.Data()),TList::Class(),AliAnalysisManager::kInputContainer);
    TString outputLYZ2SUM = fileName;
    outputLYZ2SUM += ":outputLYZ2SUManalysis";
    outputLYZ2SUM+= type;
    
    AliAnalysisDataContainer *coutputLYZ2SUM = mgr->CreateContainer(Form("cobjLYZ2SUM_%s",fileName.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputLYZ2SUM); 
    mgr->ConnectInput(taskLYZ2SUM,0,coutputFE); 
    mgr->ConnectInput(taskLYZ2SUM,1,cinputLYZ2SUM);
    mgr->ConnectOutput(taskLYZ2SUM,1,coutputLYZ2SUM); 
    cinputLYZ2SUM->SetData(fInputListLYZ2SUM);
  }
  if(LYZ2PROD) {
    AliAnalysisDataContainer *cinputLYZ2PROD = mgr->CreateContainer(Form("cobjLYZ2PRODin_%s",fileName.Data()),TList::Class(),AliAnalysisManager::kInputContainer);
    TString outputLYZ2PROD = fileName;
    outputLYZ2PROD += ":outputLYZ2PRODanalysis";
    outputLYZ2PROD+= type;
    
    AliAnalysisDataContainer *coutputLYZ2PROD = mgr->CreateContainer(Form("cobjLYZ2PROD_%s",fileName.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputLYZ2PROD); 
    mgr->ConnectInput(taskLYZ2PROD,0,coutputFE); 
    mgr->ConnectInput(taskLYZ2PROD,1,cinputLYZ2PROD);
    mgr->ConnectOutput(taskLYZ2PROD,1,coutputLYZ2PROD); 
    cinputLYZ2PROD->SetData(fInputListLYZ2PROD);
  }
  if(LYZEP) {
    AliAnalysisDataContainer *cinputLYZEP = mgr->CreateContainer(Form("cobjLYZEPin_%s",fileName.Data()),TList::Class(),AliAnalysisManager::kInputContainer);
    TString outputLYZEP = fileName;
    outputLYZEP += ":outputLYZEPanalysis";
    outputLYZEP+= type;
    
    AliAnalysisDataContainer *coutputLYZEP = mgr->CreateContainer(Form("cobjLYZEP_%s",fileName.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputLYZEP); 
    mgr->ConnectInput(taskLYZEP,0,coutputFE); 
    mgr->ConnectInput(taskLYZEP,1,cinputLYZEP);
    mgr->ConnectOutput(taskLYZEP,1,coutputLYZEP); 
    cinputLYZEP->SetData(fInputListLYZEP);
  }
  if(GFC) {
    TString outputGFC = fileName;
    outputGFC += ":outputGFCanalysis";
    outputGFC+= type;
    
    AliAnalysisDataContainer *coutputGFC = mgr->CreateContainer(Form("cobjGFC_%s",fileName.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputGFC); 
    mgr->ConnectInput(taskGFC,0,coutputFE); 
    mgr->ConnectOutput(taskGFC,1,coutputGFC);
    if (useWeights) {
      mgr->ConnectInput(taskGFC,1,cinputWeights);
      cinputWeights->SetData(weightsList);
    } 
  }
  if(QC) {
    TString outputQC = fileName;
    outputQC += ":outputQCanalysis";
    outputQC+= type;

    AliAnalysisDataContainer *coutputQC = mgr->CreateContainer(Form("cobjQC_%s",fileName.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputQC); 
    mgr->ConnectInput(taskQC,0,coutputFE); 
    mgr->ConnectOutput(taskQC,1,coutputQC);
    if (useWeights) {
      mgr->ConnectInput(taskQC,1,cinputWeights);
      cinputWeights->SetData(weightsList);
    }    
  }
  if(FQD) {
    TString outputFQD = fileName;
    outputFQD += ":outputFQDanalysis";
    outputFQD+= type;
    
    AliAnalysisDataContainer *coutputFQD = mgr->CreateContainer(Form("cobjFQD_%s",fileName.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputFQD); 
    mgr->ConnectInput(taskFQD,0,coutputFE); 
    mgr->ConnectOutput(taskFQD,1,coutputFQD);
    if(useWeights) {
      mgr->ConnectInput(taskFQD,1,cinputWeights);
      cinputWeights->SetData(weightsList);
    } 
  }
  if(MCEP) {
    TString outputMCEP = fileName;
    outputMCEP += ":outputMCEPanalysis";
    outputMCEP+= type;
    
    AliAnalysisDataContainer *coutputMCEP = mgr->CreateContainer(Form("cobjMCEP_%s",fileName.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputMCEP); 
    mgr->ConnectInput(taskMCEP,0,coutputFE); 
    mgr->ConnectOutput(taskMCEP,1,coutputMCEP); 
  }
  if(MH) {
    TString outputMH = fileName;
    outputMH += ":outputMHanalysis";
    outputMH += type;
        
    AliAnalysisDataContainer *coutputMH = mgr->CreateContainer(Form("cobjMH_%s",fileName.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputMH); 
    mgr->ConnectInput(taskMH,0,coutputFE); 
    mgr->ConnectOutput(taskMH,1,coutputMH); 
    //if (useWeights) {
    //  mgr->ConnectInput(taskMH,1,cinputWeights);
    //  cinputWeights->SetData(weightsList);
    //} 
  }
  if(NL) {
    TString outputNL = fileName;
    outputNL += ":outputNLanalysis";
    outputNL += type;
        
    AliAnalysisDataContainer *coutputNL = mgr->CreateContainer(Form("cobjNL_%s",fileName.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputNL); 
    mgr->ConnectInput(taskNL,0,coutputFE); 
    mgr->ConnectOutput(taskNL,1,coutputNL); 
    //if (useWeights) {
    //  mgr->ConnectInput(taskNL,1,cinputWeights);
    //  cinputWeights->SetData(weightsList);
    //} 
  }

}





