/************************************************
 Multiplisity Flatuation analysis task

 Auther: Satyajit Jena
 Email:  sjena@cern.ch
 Mon Oct 25 12:47:38 CEST 2010

*************************************************/

AliEbyEMFAnalysisTask *AddTaskMF()
{
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      Error("AddTask", "ERROR: No Analysis Manager");
      return NULL;
   }

   if (!mgr->GetInputEventHandler()) {
     Error("AddTaskMF", "ERROR: No input event handler");
     return NULL;
   }

     
   // AliEbyEEventBase *base =  GetEbyEAnalysisBaseObject();
   AliEbyEMultiplicityFluctuationAnalysis *analysis = GetAnalysisMFObject();
   AliESDtrackCuts *trackCuts = GetTrackCutsObject();
   //  analysis->SetBaseAnalysis(base);   
   analysis->SetAnalysisCutObject(trackCuts);
   
   AliEbyEMFAnalysisTask *taskMF  = new AliEbyEMFAnalysisTask("AliEbyEMFAnalysisTask");
   taskMF->SetAnalysisObject(analysis);	
   mgr->AddTask(taskMF);
   
   AliAnalysisDataContainer *cout 
     = mgr->CreateContainer("MF", TList::Class(), 
			    AliAnalysisManager::kOutputContainer,
			    "AnalysisOutput.mf.root");
   
   mgr->ConnectInput(taskMF, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(taskMF, 1, cout);
   return taskMF;

}

//_______________________________________________________________________________//
AliEbyEFluctuationTask *AddTaskF()
{
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      Error("AddTaskF", "ERROR: No Analysis Manager");
      return NULL;
   }

   if (!mgr->GetInputEventHandler()) {
     Error("AddTaskF", "ERROR: No input event handler");
     return NULL;
   }

     
   //  AliEbyEEventBase *base =  GetEbyEAnalysisBaseObject();
   AliEbyEFluctuationAnalysis *analysis = GetAnalysisFluctuationObject();
   AliESDtrackCuts *trackCuts = GetTrackCutsObject();
   //  analysis->SetBaseAnalysis(base);   
   analysis->SetAnalysisCutObject(trackCuts);
   
   AliEbyEFluctuationTask *taskF  = new AliEbyEFluctuationTask("AliEbyEFluctuationTask");
   taskF->SetAnalysisObject(analysis);	
   mgr->AddTask(taskF);
   
   AliAnalysisDataContainer *cout 
     = mgr->CreateContainer("Fluc", TList::Class(), 
			    AliAnalysisManager::kOutputContainer,
			    "AnalysisOutput.root");
   
   mgr->ConnectInput(taskF, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(taskF, 1, cout);
   return taskF;

}

//_______________________________________________________________________________//


AliEbyEMFAnalysisTaskT *AddTaskMFT()
{
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      Error("AddTask", "ERROR: No Analysis Manager");
      return NULL;
   }

   if (!mgr->GetInputEventHandler()) {
     Error("AddTaskMFT", "ERROR: No input event handler");
     return NULL;
   }

     
   //  AliEbyEEventBase *base =  GetEbyEAnalysisBaseObject();
   AliEbyEMultiplicityFluctuationAnalysis *analysis = GetAnalysisMFObject();
   AliESDtrackCuts *trackCuts = GetTrackCutsObject();

   // analysis->SetBaseAnalysis(base);   
   analysis->SetTreeMode(kTRUE);
   analysis->SetAnalysisCutObject(trackCuts);
   
   AliEbyEMFAnalysisTaskT *taskMF  = new AliEbyEMFAnalysisTaskT("AliEbyEMFAnalysisTaskT");
   taskMF->SetAnalysisObject(analysis);	
   mgr->AddTask(taskMF);
   
   AliAnalysisDataContainer *cout 
     = mgr->CreateContainer("MFT", TList::Class(), 
			    AliAnalysisManager::kOutputContainer,
			    "AnalysisOutput.mf.root");

   AliAnalysisDataContainer *cout1 
     = mgr->CreateContainer("EbyET", TTree::Class(), 
			    AliAnalysisManager::kOutputContainer,
			    "IntermediateMf.root");
   
   mgr->ConnectInput(taskMF, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(taskMF, 1, cout);
   mgr->ConnectOutput(taskMF, 2, cout1);
   return taskMF;

}

//_______________________________________________________________________________//

AliEbyECFAnalysisTask *AddTaskCF()
{
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      Error("AddTaskCF", "ERROR: No Analysis Manager");
      return NULL;
   }

   if (!mgr->GetInputEventHandler()) {
     Error("AddTaskCF", "ERROR: No input event handler");
     return NULL;
   }

   
   AliEbyEChargeFluctuationAnalysis *analysis = GetAnalysisCFObject();
   AliEbyECFAnalysisTask *taskCF = new AliEbyECFAnalysisTask("AliEbyECFAnalysisTask");
  
   taskCF->SetAnalysisObject(analysis);	
   
   mgr->AddTask(taskCF);
   
   AliAnalysisDataContainer *cout 
     = mgr->CreateContainer("CF", TList::Class(), 
			    AliAnalysisManager::kOutputContainer,
			    "AnalysisOutput.cf.root");
   
   mgr->ConnectInput(taskCF, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(taskCF, 1, cout);
   return taskCF;

}


//_______________________________________________________________________________//
AliEbyEEventBase *GetEbyEAnalysisBaseObject(const char* analysisLevel = "ESD",
					    const char* esdAnalysisType = "TPC",
					    const char* CentralityType = "HardType",
					    const char* centEstimator = "V0M") 

{

  printf(" Configuring Base Class \n");
  /*---------------------------------------------------
  | analysisLevel : ESD, AOD, MC
  | esdAnalysisType: TPC, ITC, TPCnITS, FORWARD, GLOBAL
  | CentralityType: 
  | HardType: Hard Coded [ from vzero multiplicity ]
  | Acc:
  | Flat:
  | All :
  | Centrality Estimator: V0M, ZDC, ..... decided later
  | Centrality Bin: 10(10%), 20(5%) , 50(2%)
  ------------------------------------------------------*/
  
  AliEbyEEventBase *base = new AliEbyEEventBase();
  // base->InitQA();
  base->SetAnalysisLevel(analysisLevel);
  
  base->SetCentralityBin(50);
  if(CentralityType == "HardType") {
    
    base->SetCentralityType(AliEbyEEventBase::kHardFlat);
    base->SetCentralityEstimator("V0M");
    
  }

  if (CentralityType == "Central") {
      base->SetCentralityType(AliEbyEEventBase::kFlat);
     base->SetCentralityEstimator("V0M");
  }

  if(analysisLevel == "ESD") {
    switch(esdAnalysisType) {
    case "TPC":
      base->SetAnalysisMode(AliEbyEEventBase::kTPC);
      base->SetPhaseSpace(0.5, 0.3, 1.5, 3.0, 70, 4);
      // eta, ptmin, ptmax, dcaxyz, TPC clusters, kisquare/clustertpc
      break;
      
    case "ITS":
      base->SetAnalysisMode(AliEbyEEventBase::kITS);
      base->SetPhaseSpace(9, -0.9, 0.9, 6, 0.45, 1.05);
      break;
      
    case "TPCnITS":
      base->SetAnalysisMode(AliEbyEEventBase::kTPCnITS);
      base->SetPhaseSpace(9, -0.9, 0.9, 6, 0.45, 1.05);
      break;
      
    case "Global":
      base->SetAnalysisMode(AliEbyEEventBase::kGlobal);
      base->SetPhaseSpace(20, -1.0, 1.0, 48, 0.3, 1.5);
      break;

    case "Forward":
      base->SetAnalysisMode(AliEbyEEventBase::kForward);
      base->SetPhaseSpace(20, -1.0, 1.0, 48, 0.3, 1.5);
      break;
      
    default:
      break;
    }
   base->SetAcceptedVertexDiamond(10.,10.,25.);
     
  }
  if(analysisLevel == "MC") 
    base->SetPhaseSpace(10, -0.5, 0.5, 16, 0.5, 0.9);
  if(analysisLevel == "AOD")
    base->SetPhaseSpace(10, -0.5, 0.5, 16, 0.5, 0.9);
  printf("Event Base Created ---------------- OK\n");
  return base;
}


//_______________________________________________________________________________//
AliEbyEChargeFluctuationAnalysis *GetAnalysisCFObject() {
  AliEbyEEventBase *base =  GetEbyEAnalysisBaseObject();
  AliEbyEChargeFluctuationAnalysis *analysis = new AliEbyEChargeFluctuationAnalysis();
  

  analysis->SetBaseAnalysis(base);
  printf("TaskCF: Charge Fluctuation Base Created ---------------- OK\n");
  return analysis;
}

//_________________________________________________________________________________//
AliEbyEMultiplicityFluctuationAnalysis *GetAnalysisMFObject() {
   AliEbyEEventBase *base =  GetEbyEAnalysisBaseObject();
   AliEbyEMultiplicityFluctuationAnalysis *analysis = new AliEbyEMultiplicityFluctuationAnalysis();
   analysis->SetBaseAnalysis(base);
   printf("TaskMF: Multiplicity Fluctuation Base Created ---------------- OK\n");
   return analysis;
}

//_________________________________________________________________________________//
AliEbyEFluctuationAnalysis *GetAnalysisFluctuationObject() {
   AliEbyEEventBase *base =  GetEbyEAnalysisBaseObject();
   AliEbyEFluctuationAnalysis *analysis = new AliEbyEFluctuationAnalysis();
   analysis->SetBaseAnalysis(base);
   printf("TaskF: Fluctuation Base Created ---------------- OK\n");
   return analysis;
}

//_________________________________________________________________________________//

AliESDtrackCuts *GetTrackCutsObject() {

  AliESDtrackCuts *cuts = new AliESDtrackCuts("bfTrackCuts","bfTrackCuts");
  // cuts->SetMinNClustersTPC(80);
  // cuts->SetMinNClustersITS(2);
  // cuts->SetMaxChi2PerClusterTPC(4.0);
  
  // cuts->SetRequireTPCRefit();
  // cuts->SetRequireITSRefit();
  // cuts->SetAcceptKinkDaughters(kFALSE);
  
  // 

  cuts->SetMaxDCAToVertexXY(0.3);
  cuts->SetMaxDCAToVertexZ(3);
  cuts->SetPtRange(0.0,1.5);
  cuts->SetEtaRange(-0.8,0.8);
  cuts->SaveHistograms("EbyECuts");
  
  return cuts;
}
