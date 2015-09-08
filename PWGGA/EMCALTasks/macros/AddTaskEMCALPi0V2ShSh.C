//astahlle
AliAnalysisTaskEMCALPi0V2ShSh *AddTaskEMCALPi0V2ShSh(
                                                     TString & arrayName,
                                                     const Bool_t  bFillAOD   = kFALSE,                                                
                                                     const Int_t   bMC        = kFALSE,
                                                     const Bool_t  exotic     = kTRUE,
                                                     const TString name       = "V1", // V1, V2, NxN, V1Unfold
                                                     const TString trigger    = "", 
                                                     const Bool_t  tm         = kTRUE, 
                                                     const Int_t   minEcell   = 50,
                                                     const Int_t   minEseed   = 100,
                                                     const Int_t   maxDeltaT  = 250, // 0
                                                     const Int_t   timeWindow = 1000, // 30
                                                     const Int_t   minEUnf    = 15, 
                                                     const Int_t   minFrac    = 1,
                                                     const Bool_t  bRecalE    = kTRUE,
                                                     const Bool_t  bBad       = kTRUE,
                                                     const Bool_t  bRecalT    = kTRUE,
                                                     const Bool_t  bNonLine   = kFALSE, // kTRUE
                                                     const Int_t   minCen     = -1,
                                                     const Int_t   maxCen     = -1,
                                                     const Float_t clusterEnergyCutEvent = -1,
                                                     const Int_t   nRowDiff   = 1,
                                                     const Int_t   nColDiff   = 1,
                                                     const Bool_t  skipOrReject = kFALSE,
                                                     const Int_t   debug      = 0
                                                    )
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  
  if (!mgr) {
    ::Error("AddTaskEMCALPi0V2ShSh", "No analysis manager to connect to.");
    return NULL;
  }  
  
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskEventplane", "This task requires an input event handler");
    return NULL;
  }		
  // TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  
  Bool_t ismc = kFALSE;
  ismc = (mgr->GetMCtruthEventHandler())?kTRUE:kFALSE; 
  cout << "AddTaskEMCALPi0V2ShSh - MC config is: " << ismc << endl;  
  if (ismc) return 0;
  
  // UInt_t physMB = 0;
  // UInt_t physJet = 0;
  // UInt_t physGam = 0;
  // UInt_t physEMC = 0; 
  
  // TString sGeomName = AliEMCALGeometry::GetDefaultGeometryName();
  
  // TString arrayName   = "";
  // Bool_t  outAOD      = kFALSE;
  // Bool_t  kMC         = kFALSE; /// With real data kMC = kFALSE
  // Bool_t  exo         = kTRUE;  // Remove exotic cells
  // TString clTrigger   = ""; 
  // Bool_t  clTM        = kTRUE;
  // Int_t   minEcell    = 50;     // 50  MeV (10 MeV used in reconstruction)
  // Int_t   minEseed    = 100;    // 100 MeV
  // Int_t   dTime       = 0;      // default, 250 ns
  // Int_t   wTime       = 30;     // default 425 < T < 825 ns, careful if time calibration is on
  // Int_t   unfMinE     = 15;     // Remove cells with less than 15 MeV from cluster after unfolding
  // Int_t   unfFrac     = 1;      // Remove cells with less than 1% of cluster energy after unfolding
  // Bool_t  calibEE     = kTRUE; // It is set automatically, but here we force to use it or not in any case
  // Bool_t  badMap      = kTRUE; // It is set automatically, but here we force to use it or not in any case  
  // Bool_t  calibTT     = kTRUE; // It is set automatically, but here we force to use it or not in any case
  // Bool_t  clnonlin    = kTRUE;  // Apply non linearity (clusterization)
  // AliAnalysisTaskEMCALClusterize* clv1 = AddTaskEMCALClusterize(arrayName,outAOD,kMC,exo,"V1",clTrigger, clTM,
  //                                                               minEcell,minEseed,dTime,wTime,unfMinE,unfFrac,
  //                                                               calibEE,badMap,calibTT,clnonlin);
  AliAnalysisTaskEMCALClusterize* clust = AddTaskEMCALClusterize(arrayName,bFillAOD,bMC,exotic,name,trigger,tm,
                                                                 minEcell,minEseed,maxDeltaT,timeWindow,minEUnf,minFrac,
                                                                 bRecalE,bBad,bRecalT,bNonLine,
                                                                 minCen,maxCen,clusterEnergyCutEvent,nRowDiff,nColDiff,skipOrReject);    
  printf("Name of clusterizer array: %s\n",arrayName.Data());

  AliAnalysisTaskEMCALPi0V2ShSh *task = new AliAnalysisTaskEMCALPi0V2ShSh("EMCALPi0V2ShSh");
  task->SetEMCALClusterListName(arrayName);
  task->SelectCollisionCandidates(AliVEvent::kSemiCentral);
  task->SetDebugLevel(debug);
  
  if (!ismc) {		
    mgr->AddTask(task);
      
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("hist", TList::Class(),AliAnalysisManager::kOutputContainer, Form("%s",AliAnalysisManager::GetCommonFileName()));
      
    mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task, 1, coutput1);      
  }
  return task;
}
