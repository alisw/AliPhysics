Bool_t   kSimulation    = kFALSE;
TString  fDataType      = "AOD"; 
Int_t    kYears         = 2011;
TString  kCollisions    = "PbPb";
TString  fTrigger       = "EMCGA";
Bool_t   kEventTriggerAtTaskSE = kFALSE;
Float_t  fMinCen        = -1;
Float_t  fMaxCen        = -1;
TString  fAnaMesonType  = "Pi0";
Bool_t   kAnaPhotonCorr = kTRUE;
Bool_t   kAnaMesonCorr  = kFALSE; 
Bool_t   kTimeCut      = kFALSE;
Bool_t   kDistBC       = kTRUE;
Bool_t   kRecalClus    = kTRUE;
Bool_t   kRecalClusE   = kTRUE;
Bool_t   kRecalClusPos = kTRUE;
Bool_t   kRecalClusSSA = kTRUE;
Bool_t   kNonLin       = kTRUE;
Bool_t   kTM           = kFALSE;
Float_t  fDPhiCut      = 0.03;
Float_t  fDEtaCut      = 0.025;
Bool_t   kExotic       = kTRUE;
Float_t  fExoticFraction = 0.95;///for pp:0.97, for PbPb:0.95
Bool_t   kFidul          = kFALSE;
Bool_t   kReClusterier   = kFALSE;
TString  fName           = "V2";
Float_t  fMinCell        = 0.1;
Float_t  fMinSeed        = 0.3;

AliAnalysisTaskEMCALCaloTrackCorr *AddTaskEmcalCaloTrackCorr(
  const TString  data          = "AOD",
  const TString  coll          = "pp",
  const Bool_t   simulation    = kFALSE,
  const TString  trigger       = "MB", 
  const Bool_t   triggerSE     = kFALSE,
  const Float_t  minCen        = -1,
  const Float_t  maxCen        = -1,
  const Bool_t   anaPhotonCorr = kTRUE,
  const Bool_t   anaMesonCorr  = kFALSE,
  const TString  anaMesonType  = "Pi0",
  const Bool_t   timecut       = kFALSE,
  const Bool_t   tm            = kFALSE,
  const Float_t  dphicut       = 0.03,
  const Float_t  detacut       = 0.025,
  const Bool_t   exotic       = kTRUE,
  const Float_t  exoticFraction = 0.95,
  const Bool_t   reClusterizer  = kFALSE,
  const TString  name          = "V2",
  const Float_t  minCell       = 0.1,
  const Float_t  minSeed       = 0.3)
{
  fDataType      = data;
  kCollisions    = coll;
  kSimulation    = simulation;
  fTrigger       = trigger;
  kEventTriggerAtTaskSE = triggerSE;
  fMinCen        = minCen;
  fMaxCen        = maxCen;
  fAnaMesonType  = anaMesonType;
  kAnaPhotonCorr = anaPhotonCorr;
  kAnaMesonCorr  = anaMesonCorr;
  kTimeCut      = timecut;
  kTM           = tm;
  fDPhiCut      = dphicut;
  fDEtaCut      = detacut;
  kExotic       = exotic;
  fExoticFraction = exoticFraction;
  kReClusterier   = reClusterizer;
  fName           = name;
  fMinCell        = minCell;
  fMinSeed        = minSeed;
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) 
  {
    ::Error("AddTask", "No analysis manager to connect to.");
    return NULL;
  }  

  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskHadCorr", "This task requires an input event handler");
    return NULL;
  }
 
  if(fDataType == "ESD"){ 
   gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/macros/CreateTrackCutsPWGJE.C");
   AliESDtrackCuts * esdTrackCuts = CreateTrackCutsPWGJE(10041004);
   esdTrackCuts->SetEtaRange(-0.8,0.8);
  }

  AliAnalysisTaskEMCALCaloTrackCorr *task = new AliAnalysisTaskEMCALCaloTrackCorr("NeutralCorr"); 
  task->SetMC(kSimulation);
  task->SetDataType(fDataType);
  if(fDataType == "ESD")task->SetTrackCuts(esdTrackCuts);
  if(fDataType == "AOD") task->SwitchOnAODHybridTrackSelection();
  task->SetTrackMatchedDPhiCut(fDPhiCut);
  task->SetTrackMatchedDEtaCut(fDEtaCut);
  task->SetLargeCorrTrigger(5, 50);

  if(kAnaPhotonCorr && !kAnaMesonCorr){
   task->SwitchOnAnaPhotonCorr();
   task->SwitchOffAnaMesonCorr();
   task->SwitchOffFillMesonAOD();
  }
  else if(!kAnaPhotonCorr && kAnaMesonCorr){
   task->SwitchOffAnaPhotonCorr();
   task->SwitchOnAnaMesonCorr();
   task->SetAnaMesonType(fAnaMesonType);
   task->SwitchOnFillMesonAOD();
  }
  else {
   task->SwitchOffAnaPhotonCorr();
   task->SwitchOffAnaMesonCorr();
  } 
  
  task->SetCentralityClass("V0M");
  task->SetCentralityBin(fMinCen,fMaxCen); 
  task->SetEventPlaneMethod("V0");
  task->SetEMCALGeometryName("EMCAL_COMPLETEV1");
  mgr->AddTask(task);

  if(kCollisions =="pp")   task->SwitchOnTrackMultBins();
  if(kCollisions =="PbPb") task->SwitchOffTrackMultBins();

  task->SetDeltaPhiCutRange(TMath::Pi()/2., 3*TMath::Pi()/2.);
  task->SetNTriggPtBins(2);
  Float_t fTriggerPtBins[3]={8,15,25};
  task->SetTriggerBins(fTriggerPtBins);
  task->SetNAssocPtBins(5);
  Float_t fAssociatedPtBins[6]={0.5, 1, 2, 4, 6, 15};
  task->SetAssociatedBins(fAssociatedPtBins);   

  ConfigureTrigger(task);

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
 
 if(kAnaPhotonCorr && !kAnaMesonCorr) {
    AliAnalysisDataContainer *coutputpt1 = mgr->CreateContainer(Form("InclusivePhotonCen%.1f_%.1f",fMinCen, fMaxCen), TList::Class(),
                            AliAnalysisManager::kOutputContainer, "AnalysisResults.root");  
   }
   else if (!kAnaPhotonCorr && kAnaMesonCorr){
    AliAnalysisDataContainer *coutputpt1 = mgr->CreateContainer(Form("Inclusive%sCen%.1f_%.1f", fAnaMesonType.Data(),fMinCen, fMaxCen), TList::Class(),
                             AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
   }
   else {
    AliAnalysisDataContainer *coutputpt1 = mgr->CreateContainer(Form("InclusiveNoCorrCen%.1f_%.1f", fMinCen, fMaxCen), TList::Class(),
                             AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
   }
    mgr->ConnectInput(task, 0, cinput);
    mgr->ConnectOutput(task, 1, coutputpt1);
    return task;

}

void ConfigureTrigger(AliAnalysisTaskEMCALCaloTrackCorr *task1)
{
  if(!kEventTriggerAtTaskSE){
   task1->SwitchOffEventTriggerAtSE();
   if(fTrigger=="EMC7"){
    printf("CaloTrackCorr trigger EMC7\n");
    task1->SetEventTriggerMask(AliVEvent::kEMC7);
   }
   else if (fTrigger=="INT7"){
    printf("CaloTrackCorr trigger INT7\n");
    task1->SetEventTriggerMask(AliVEvent::kINT7);
   }
   else if(fTrigger=="EMC1"){
    printf("CaloTrackCorr trigger EMC1\n");
    task1->SetEventTriggerMask(AliVEvent::kEMC1);
   }
   else if(fTrigger=="MB"){
    printf("CaloTrackCorr trigger MB\n");
    task1->SetEventTriggerMask(AliVEvent::kMB);
   }  
   else if(fTrigger=="AnyINT"){
    printf("CaloTrackCorr trigger AnyINT\n");
    task1->SetEventTriggerMask(AliVEvent::kAnyINT);
   }  
   else if(fTrigger=="EMCEGA"){
    printf("CaloTrackCorr trigger EMC Gamma\n");
    task1->SetEventTriggerMask(AliVEvent::kEMCEGA);
   } 
   else if(fTrigger=="EMCEJE"){
    printf("CaloTrackCorr trigger EMC Jet\n");
    task1->SetEventTriggerMask(AliVEvent::kEMCEJE);
   }
   else if(fTrigger=="Central"){
    printf("CaloTrackCorr trigger Central\n");
    task1->SetEventTriggerMask(AliVEvent::kCentral);
   } 
   else if(fTrigger=="SemiCentral"){
    printf("CaloTrackCorr trigger SemiCentral\n");
    task1->SetEventTriggerMask(AliVEvent::kSemiCentral);
   }
   else if(fTrigger=="SemiOrCentral"){
    printf("CaloTrackCorr trigger SemiCentral Or Central\n");
    task->SetEventTriggerMask(AliVEvent::kSemiCentral | AliVEvent::kCentral);
   }
   else{
    task1->SetEventTriggerMask(AliVEvent::kAny);
   }

  }
  else {
   task1->SwitchOnEventTriggerAtSE();
   if(fTrigger=="EMC7"){
    printf("CaloTrackCorr trigger EMC7\n");
    task1->SelectCollisionCandidates(AliVEvent::kEMC7);
   }
   else if (fTrigger=="INT7"){
    printf("CaloTrackCorr trigger INT7\n");
    task1->SelectCollisionCandidates(AliVEvent::kINT7);
   }
   else if(fTrigger=="EMC1"){
    printf("CaloTrackCorr trigger EMC1\n");
    task1->SelectCollisionCandidates(AliVEvent::kEMC1);
   }
   else if(fTrigger=="MB"){
    printf("CaloTrackCorr trigger MB\n");
    task1->SelectCollisionCandidates(AliVEvent::kMB);
   }  
   else if(fTrigger=="AnyINT"){
    printf("CaloTrackCorr trigger AnyINT\n");
    task1->SelectCollisionCandidates(AliVEvent::kAnyINT);
   }  
   else if(fTrigger=="EMCEGA"){
    printf("CaloTrackCorr trigger EMC Gamma\n");
    task1->SelectCollisionCandidates(AliVEvent::kEMCEGA);
   } 
   else if(fTrigger=="EMCEJE"){
    printf("CaloTrackCorr trigger EMC Jet\n");
    task1->SelectCollisionCandidates(AliVEvent::kEMCEJE);
   }
   else if(fTrigger=="Central"){
    printf("CaloTrackCorr trigger Central\n");
    task1->SelectCollisionCandidates(AliVEvent::kCentral);
   } 
   else if(fTrigger=="SemiCentral"){
    printf("CaloTrackCorr trigger SemiCentral\n");
    task1->SelectCollisionCandidates(AliVEvent::kSemiCentral);
   }
   else if(fTrigger=="SemiOrCentral"){
    printf("CaloTrackCorr trigger SemiCentral Or Central\n");
    task->SelectCollisionCandidates(AliVEvent::kSemiCentral | AliVEvent::kCentral);
   }
   else{
    task1->SelectCollisionCandidates(AliVEvent::kAny);
   }
  }

}
