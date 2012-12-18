// QAtrain_duo - assumes having 2 sets of input data files using 2 suffixes:
//    AliESDs_barrel.root AliESDfriends_barrel.root ...  suffix = _barrel
//    AliESDs_outer.root AliESDfriends_outer.root ...    suffix = _outer
// HOW TO RUN (for run_no=158285):
// QA job: should call:
//    ln -fs AliESDs_barrel.root AliESDs.root
//    ln -fs AliESDfriends_barrel.root AliESDfriends.root
//    ... same for all specific reco files
//    aliroot -b -q QAtrain_duo.C\(\"_barrel\"\,158285\) producing:
//      -> EventStat_temp_barrel.root QAresults_barrel.root 400_400_0_350.stat
//
//    ln -fs AliESDs_outer.root AliESDs.root
//    ln -fs AliESDfriends_outer.root AliESDfriends.root
//    ... same for all specific reco files
//    aliroot QAtrain_duo.C\(\"_outer\"\,158285\) producing:
//      -> EventStat_temp_outer.root QAresults_outer.root 400_400_0_380.stat
//
// Each merging job:
// for stages i < 5
//    aliroot -b -q QAtrain_duo.C\(\"_barrel\"\,158285\,\"wn.xml\"\,i\)
//    aliroot -b -q QAtrain_duo.C\(\"_outer\"\,158285\,\"wn.xml\"\,i\)
//      -> same output files as the QA jobs (except .stat)
// for stage 5
//    aliroot -b -q QAtrain_duo.C\(\"_barrel\"\,158285\,\"Stage_5.xml\"\,5\)
//      -> event_stat_barrel.root trending_barrel.root 145230_145230_0_120345.stat
//    aliroot -b -q QAtrain_duo.C\(\"_outer\"\,158285\,\"Stage_5.xml\"\,5\)
//      -> event_stat_outer.root trending_outer.root 145230_145230_0_132897.stat
#include "Riostream.h"
void LoadLibraries();
void AddAnalysisTasks(const char *, const char *); 
void QAmerge(const char *,const char *, Int_t);

Int_t iCollisionType = 0; // 0=pp, 1=PbPb
// Trigger mask.

UInt_t kTriggerInt = AliVEvent::kAnyINT;
UInt_t kTriggerMuonBarrel = AliVEvent::kMUU7 | AliVEvent::kMuonUnlikeLowPt8 | AliVEvent::kMuonUnlikeLowPt0;
UInt_t kTriggerEMC   = AliVEvent::kEMC7 | AliVEvent::kEMC8 | AliVEvent::kEMCEJE | AliVEvent::kEMCEGA;
UInt_t kTriggerHM   = AliVEvent::kHighMult;
// Main trigger mask used:
UInt_t kTriggerMask = kTriggerInt;

Int_t runNumbers[5] = {158626};

Bool_t doCDBconnect   = 1;
Bool_t doEventStat    = 1;
Bool_t doCentrality   = 0;
Bool_t doQAsym        = 1;
Bool_t doVZERO        = 1;   // there is a 2nd file
Bool_t doVZEROPbPb    = 1; 
Bool_t doVertex       = 1;
Bool_t doSPD          = 1;   // needs RP   
Bool_t doTPC          = 1;
Bool_t doHLT          = 1;
Bool_t doSDD          = 1;   // needs RP
Bool_t doSSDdEdx      = 1;

Bool_t doTRD          = 1;
Bool_t doITS          = 1;
Bool_t doITSsaTracks  = 1; 
Bool_t doITSalign     = 1;  
Bool_t doCALO         = 1;
Bool_t doMUONTrig     = 1;
Bool_t doImpParRes    = 1;
Bool_t doMUON         = 1;
Bool_t doTOF          = 1;
Bool_t doHMPID        = 1;
Bool_t doT0           = 1;
Bool_t doZDC          = 1;
Bool_t doPIDResponse  = 1;
Bool_t doPIDqa        = 1; //new
Bool_t doFMD          = 1; // new
Bool_t doPHOS         = 1; // new
Bool_t doPHOSTrig     = 1; // new
Bool_t doEMCAL        = 1;
Bool_t doFBFqa        = 1; // new - not ported yet to revision

               // Debug level
Int_t       debug_level        = 1;        // Debugging
Int_t       run_number = 0;

void QAtrain_duo(const char *suffix="", Int_t run = 0, 
             const char *xmlfile   = "wn.xml",
             Int_t  stage          = 0, /*0 = QA train, 1...n - merging stage*/
             const char *cdb     = "raw://")
{
  run_number = run;
  TString ss(suffix);
  ss.ToLower();
  Bool_t ibarrel = (ss.Contains("barrel"))?kTRUE:kFALSE;

  TString cdbString(cdb);
  if (cdbString.Contains("raw://"))
  {
    TGrid::Connect("alien://");
    if (!gGrid || !gGrid->IsConnected()) {
      ::Error("QAtrain", "No grid connection");
      return;
    }  
  }
  // Set temporary merging directory to current one
  gSystem->Setenv("TMPDIR", gSystem->pwd());
  // Set temporary compilation directory to current one
  gSystem->SetBuildDir(gSystem->pwd(), kTRUE);
  // Load libraries
  LoadLibraries();
  printf("Include path: %s\n", gSystem->GetIncludePath());
  // Create manager
  AliAnalysisManager *mgr  = new AliAnalysisManager("PilotAnalysis", "Production train");
  mgr->SetRunFromPath(run_number);
  // Input handler
  AliESDInputHandlerRP *esdHandler = new AliESDInputHandlerRP();
  if (ibarrel) {
    esdHandler->SetReadFriends(kTRUE);
    esdHandler->SetActiveBranches("ESDfriend");
  }
  mgr->SetInputEventHandler(esdHandler);
  mgr->SetDebugLevel(debug_level);
  mgr->SetFileInfoLog("fileinfo.log"); 
  
  // AnalysisTasks
  AddAnalysisTasks(suffix, cdb);
  if (stage>0) {
    QAmerge(suffix, xmlfile, stage);
    return;
  }   
  // Input chain
  TChain *chain = new TChain("esdTree");
  chain->Add("AliESDs.root");
  TStopwatch timer;
  timer.Start();
  if (mgr->InitAnalysis()) {                                                                                                              
    mgr->PrintStatus(); 
    mgr->SetSkipTerminate(kTRUE);
//    mgr->SetNSysInfo(1);
    mgr->StartAnalysis("local", chain);
  }
  timer.Print();
}

void LoadLibraries()
{
  gSystem->SetIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TRD -I$ALICE_ROOT/PWGPP -I$ALICE_ROOT/PWGPP/TRD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libTENDER");
  gSystem->Load("libPWGPP.so");
  gSystem->Load("libAliHLTTrigger.so");

  if (doEMCAL || doPHOS || doCALO) {
     gSystem->Load("libEMCALUtils");
     gSystem->Load("libPHOSUtils");
     gSystem->Load("libPWGCaloTrackCorrBase");
     gSystem->Load("libPWGGACaloTrackCorrelations");
     gSystem->Load("libPWGGACaloTasks");
     gSystem->Load("libPWGGAPHOSTasks");
     gSystem->Load("libPWGEMCAL");
     gSystem->Load("libPWGGAEMCALTasks");
  }  
  if(doMUON || doMUONTrig) {
     gSystem->Load("libPWGmuon");
     gSystem->Load("libPWGPPMUONlite");
     gSystem->Load("libPWGmuondep");
  }
  if (doFMD) {
     gSystem->Load("libPWGLFforward2");
  }      
}

void AddAnalysisTasks(const char *suffix, const char *cdb_location)
{
  TString ss(suffix);
  ss.ToLower();
  Bool_t ibarrel = (ss.Contains("barrel"))?kTRUE:kFALSE;
  Bool_t iall = ss.IsNull();
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  mgr->SetCommonFileName(Form("QAresults%s.root",suffix));
  // Statistics task
  mgr->AddStatisticsTask(kTriggerMask);
  //
  // CDB connection
  //
  if (doCDBconnect) {
    gROOT->LoadMacro("$ALICE_ROOT/PWGPP/PilotTrain/AddTaskCDBconnect.C");
    AliTaskCDBconnect *taskCDB = AddTaskCDBconnect();
    if (!taskCDB) return;
    AliCDBManager *cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage(cdb_location);
    taskCDB->SetRunNumber(run_number);
  }    
  
  //
  // Event Statistics (Jan Fiete)
  //
  if (doEventStat) {
      gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
      AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(kFALSE /*MC*/);
      AliAnalysisDataContainer *cstatsout = (AliAnalysisDataContainer*)mgr->GetOutputs()->FindObject("cstatsout");
      cstatsout->SetFileName(Form("EventStat_temp%s.root", suffix));
  }
  
  //
  // Centrality (A. Toia)
  //
  if (doCentrality) {
//     if (!iCollisionType) {
//        printf("Disabling centrality task for p-p\n");
//        doCentrality = kFALSE;
//     } else {           
        gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");
        AliCentralitySelectionTask *taskCentrality = AddTaskCentrality();
      AliAnalysisDataContainer *centralityout = (AliAnalysisDataContainer*)mgr->GetOutputs()->FindObject("CentralityStat");
      centralityout->SetFileName(Form("EventStat_temp%s.root", suffix));
//     }   
  }   
  
  // Vertexing (A. Dainese)
  // 
  if (doVertex && (ibarrel || iall)) {
    gROOT->LoadMacro("$ALICE_ROOT/PWGPP/macros/AddTaskVertexESD.C");
    AliAnalysisTaskVertexESD* taskvertexesd =  AddTaskVertexESD(kFALSE, kTriggerMask);
    taskvertexesd->SelectCollisionCandidates(kTriggerMask);
  }  

  // TPC QA (E. Sicking)
  //
  if (doQAsym && (ibarrel || iall)) {
  // offline trigger in AddTask
    gROOT->LoadMacro("$ALICE_ROOT/PWGPP/PilotTrain/AddTaskQAsym.C");
    AliAnalysisTaskSE * taskqasim = AddTaskQAsym(0, kTriggerMask, kTriggerHM, kTriggerEMC, kTriggerMuonBarrel);
  }  
  //
  // VZERO QA  (C. Cheshkov)
  //
  if (doVZERO && (ibarrel || iall)) {
    gROOT->LoadMacro("$ALICE_ROOT/PWGPP/PilotTrain/AddTaskVZEROQA.C");
    AliAnalysisTaskSE * taskv0qa = AddTaskVZEROQA(0);
//  taskv0qa->SelectCollisionCandidates();
    gROOT->LoadMacro("$ALICE_ROOT/PWGPP/PilotTrain/AddTaskVZEROQATrig.C");
    AliAnaVZEROQA *taskv0qatrig = AddTaskVZEROQATrig(0);
    taskv0qatrig->SelectCollisionCandidates(AliVEvent::kINT8);

  }
  if (doVZEROPbPb && iCollisionType==1) {
    gROOT->LoadMacro("$ALICE_ROOT/PWGPP/VZERO/AddTaskVZEROPbPb.C");
    AliAnaVZEROPbPb* taskV0PbPb = (AliAnaVZEROPbPb*)AddTaskVZEROPbPb(0);
//    taskV0PbPb->SetClassesNames("CTRUE-,C0HWU-,CPBI2WU_B1-,CPBI2_B1-,CPBI1WU-,CPBI1-,CVHNWU-,CVHN-,CVLNWU-,CVLN-");
//    taskV0PbPb->SetClassesNames("CTRUE-,C0HWU-,CPBI2WU,CPBI2,CPBI1WU-,CPBI1-,CVHNWU,CVHN,CVLNWU,CVLN");
//    taskV0PbPb->SetClassesNames("CTRUE-,C0HWU-,CPBI2WU-,CPBI2-,CPBI2WU_B1-,CPBI2_B1-,CPBI1WU-,CPBI1-,CVHNWU-,CVHN-,CVHN_R2-,CVHNWU_R2-,CVLNWU-,CVLN-,CVLN_B2-,CVLNWU_B2-");
//    taskV0PbPb->SetClassesNames("CTRUE-,C0HWU-,CPBI2WU-,CPBI2-,CPBI2WU_B1-,CPBI2_B1-,CPBI1WU-,CPBI1-,CVHNWU-,CVHN-,CVHN_R2-,CVHNWU_R2-,CVLNWU-,CVLN-,CVLN_R1-,CVLN_B2-,CVLNWU_R1-,CVLNWU_B2-");
    taskV0PbPb->SetClassesNames("CTRUE-,C0HWU-,CPBI2WU-,CPBI2-,CPBI2WU_B1-,CPBI2_B1-,CPBI1WU-,CPBI1-,CVHNWU-,CVHN-,CVHN_R2-,CVHNWU_R2-,CVLNWU-,CVLN-,CVLN_R1-,CVLN_B2-,CVLNWU_R1-,CVLNWU_B2-,CSEMI_R1-,CSEMIWU_R1-,CCENT_R2-,CCENTWU_R2-");
  }
  //
  // TPC (Jacek Otwinowski & Michael Knichel)
  //
  //
  // Optionally MC information can be used by setting the 1st argument to true
  // Optionally friends information can be switched off by setting the 2st argument 
  // to false
  // Optionally highMult axis can be used by setting the 3st argument to true (for PbPb)
  if (doTPC && (ibarrel || iall)) {
    gROOT->LoadMacro("$ALICE_ROOT/PWGPP/TPC/macros/AddTaskPerformanceTPCdEdxQA.C");
    AliPerformanceTask *tpcQA = 0;
    if (iCollisionType) {
       // High multiplicity Pb-Pb
       tpcQA = AddTaskPerformanceTPCdEdxQA(kFALSE, kTRUE, kTRUE);
    } else {
      // Low multiplicity (pp)
       tpcQA = AddTaskPerformanceTPCdEdxQA(kFALSE, kTRUE, kFALSE);
    }
    tpcQA->SelectCollisionCandidates(kTriggerMask);
  }  

  // HLT (Alberica Toia)
  if (doHLT && (ibarrel || iall)) {
    gROOT->LoadMacro("$ALICE_ROOT/PWGPP/TPC/macros/AddTaskPerformanceTPCdEdxQA.C");
    AliPerformanceTask *hltQA = AddTaskPerformanceTPCdEdxQA(kFALSE, kTRUE, kFALSE,0,kTRUE);
    hltQA->SelectCollisionCandidates(kTriggerMask);
  }  
  //
  // SPD (A. Mastroserio)
  //
  if (doSPD && (ibarrel || iall)) {
    gROOT->LoadMacro("$ALICE_ROOT/PWGPP/PilotTrain/AddTaskSPDQA.C");
    AliAnalysisTaskSPD* taskspdqa = (AliAnalysisTaskSPD*)AddTaskSPDQA();
    // Request from Annalisa
    if (iCollisionType) taskspdqa->SetHeavyIonMode();
    taskspdqa->SelectCollisionCandidates(kTriggerMask);
    taskspdqa->SetOCDBInfo(run_number, cdb_location);
  }  
  //
  // SDD (F. Prino)
  //
  if (doSDD && (ibarrel || iall)) {
    gROOT->LoadMacro("$ALICE_ROOT/PWGPP/PilotTrain/AddSDDPoints.C");
    AliAnalysisTaskSE* tasksdd = AddSDDPoints();
    tasksdd->SelectCollisionCandidates(kTriggerMask);
  }
  //
  // SSD dEdx (Marek Chojnacki)
  //
  if (doSSDdEdx && (ibarrel || iall)) {
    gROOT->LoadMacro("$ALICE_ROOT/PWGPP/PilotTrain/AddTaskdEdxSSDQA.C");
    AliAnalysisTaskSE* taskssddedx = AddTaskdEdxSSDQA();
    taskssddedx->SelectCollisionCandidates(kTriggerMask);
  }

  //
  // ITS
  //
  if (doITS && (ibarrel || iall)) {
  // hardcoded non-zero trigger mask
      gROOT->LoadMacro("$ALICE_ROOT/PWGPP/macros/AddTaskPerformanceITS.C");
      AliAnalysisTaskITSTrackingCheck *itsQA = 0;
      AliAnalysisTaskITSTrackingCheck *itsQACent0010 = 0;
      AliAnalysisTaskITSTrackingCheck *itsQACent3050 = 0;
      AliAnalysisTaskITSTrackingCheck *itsQACent6080 = 0;
      if(iCollisionType==0) {
        itsQA = AddTaskPerformanceITS(kFALSE);
      } else {
        itsQA = AddTaskPerformanceITS(kFALSE);
        itsQACent0010 = AddTaskPerformanceITS(kFALSE,kFALSE,kFALSE,3500,10000);
        itsQACent3050 = AddTaskPerformanceITS(kFALSE,kFALSE,kFALSE,590,1570);
        itsQACent6080 = AddTaskPerformanceITS(kFALSE,kFALSE,kFALSE,70,310);
      }
  }
  //
  // ITS saTracks, align (F.Prino)
  //
  if (doITSsaTracks && (ibarrel || iall)) {
  // offline trigger in AddTask
     gROOT->LoadMacro("$ALICE_ROOT/PWGPP/macros/AddTaskITSsaTracks.C");
     AliAnalysisTaskITSsaTracks *itssaTracks = AddTaskITSsaTracks(kFALSE,kFALSE);
     itssaTracks->SelectCollisionCandidates(kTriggerMask);
  }   
  if (doITSalign && (ibarrel || iall)) {
  // no offline trigger selection
     gROOT->LoadMacro("$ALICE_ROOT/PWGPP/macros/AddTaskITSAlign.C");
     AliAnalysisTaskITSAlignQA *itsAlign = AddTaskITSAlign(0,2012);
  }   
  //
  // TRD (Alex Bercuci, M. Fasel) 
  //
  if(doTRD && (ibarrel || iall)) {
  // no offline trigger selection
      gROOT->LoadMacro("$ALICE_ROOT/PWGPP/macros/AddTrainPerformanceTRD.C");
      // steer individual TRD tasks
      Bool_t 
      doCheckESD(kTRUE),  // AliTRDcheckESD
      doCheckDET(kTRUE),  // AliTRDcheckDET
      doEffic(kTRUE),     // AliTRDefficiency
      doResolution(kTRUE),// AliTRDresolution
      doCheckPID(kTRUE),  // AliTRDcheckPID
      doV0Monitor(kFALSE);// AliTRDv0Monitor
      AddTrainPerformanceTRD(Translate(doCheckESD, doCheckDET, doEffic, doResolution, doCheckPID, doV0Monitor));
  }

  //
  // ZDC (Chiara Oppedisano) 
  //
  if(doZDC && (ibarrel || iall)) {
  // hardcoded kMB trigger mask
     gROOT->LoadMacro("$ALICE_ROOT/PWGPP/ZDC/AddTaskZDCQA.C");
     AliAnalysisTaskSE *taskZDC = AddTaskZDCQA();
     taskZDC->SelectCollisionCandidates(kTriggerMask);
  }   
  //
  // Calorimetry (Gustavo Conesa)
  //

  if(doCALO) {
  // - In case of MC analysis, it skips the triggered events wagon (given that
  //the normal wagon is always suffixed as "default")
  //- No need to pass the type of data
  //- No need to specify the year. In principle the year is only needed when
  //setting the size of some histograms (change in number of SM along the years),
  //I do not know how to access automatically the run number when the histograms
  // are created in the UserCreateOutput. By default I set them to the maximum
  //expected size, but I still let the possibility to check the year.

  //So the way to initialize the task in trunk is now

  //AddTaskCalorimeterQA(Bool_t kSimulation = kFALSE,
                                   //  const char *suffix="default",
                                   // TString outputFile = "", 
                                   //  Int_t year = 2012, 
                                   //  Bool_t kPrintSettings = kFALSE)

  //So, in principle only the first 2 need to be specified.

  
      gROOT->LoadMacro("$ALICE_ROOT/PWGGA/CaloTrackCorrelations/macros/QA/AddTaskCalorimeterQA.C");
      AliAnalysisTaskCaloTrackCorrelation *taskCaloQA = AddTaskCalorimeterQA(kFALSE, "default");
      taskCaloQA->SetDebugLevel(0);
      // offline mask set in AddTask to kMB
      taskCaloQA->SelectCollisionCandidates(kTriggerMask);
      // Add a new calo task with EMC1 trigger only
      taskCaloQA = AddTaskCalorimeterQA(kFALSE, "trigEMC");
      taskCaloQA->SetDebugLevel(0);
      taskCaloQA->SelectCollisionCandidates(kTriggerEMC);
  }

  //
  // Muon Trigger
  //
  
  if(doMUONTrig) {
  // no offline trigger selection
      gROOT->LoadMacro("$ALICE_ROOT/PWGPP/macros/AddTaskMTRchamberEfficiency.C");
      AliAnalysisTaskTrigChEff *taskMuonTrig = AddTaskMTRchamberEfficiency();
      taskMuonTrig->GetMuonTrackCuts()->SetCustomParamFromRun(160000,2);       
  }

  //
  // Impact parameter resolution (xianbao.yuan@pd.infn.it, andrea.dainese@pd.infn.it)
  //
  if (doImpParRes && (ibarrel || iall)) {
    gROOT->LoadMacro("$ALICE_ROOT/PWGPP/macros/AddTaskImpParRes.C");
    AliAnalysisTaskSE* taskimpparres=0;
    if(iCollisionType==0) {
       taskimpparres= AddTaskImpParRes();
    } else {
       taskimpparres= AddTaskImpParRes(kFALSE,-1,kFALSE,kFALSE);
    }
    taskimpparres->SelectCollisionCandidates(kTriggerMask);
  }  
  //
  // MUON QA (Philippe Pillot)
  //
  if (doMUON) {
  // trigger analysis internal
    gROOT->LoadMacro("$ALICE_ROOT/PWGPP/PilotTrain/AddTaskMuonQA.C");
    AliAnalysisTaskMuonQA* taskmuonqa= AddTaskMuonQA();
    taskmuonqa->GetTrackCuts()->SetCustomParamFromRun(160000,2);
  }  
  //
  // TOF (Francesca Bellini)
  //
  if (doTOF && (ibarrel || iall)) {
    gROOT->LoadMacro("$ALICE_ROOT/PWGPP/TOF/AddTaskTOFQA.C");
    AliAnalysisTaskTOFqa *tofQA = AddTaskTOFQA();
    tofQA->SelectCollisionCandidates(kTriggerMask);
  } 
   //
  // PIDResponse(JENS)
  //
  if (doPIDResponse && (ibarrel || iall)) {
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C"); 
    AliAnalysisTaskPIDResponse *PIDResponse = AddTaskPIDResponse();
    PIDResponse->SelectCollisionCandidates(kTriggerMask);
  }  

  //
  // PIDqa(JENS)
  //
  if (doPIDqa && (ibarrel || iall)) {
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDqa.C");
    AliAnalysisTaskPIDqa *PIDQA = AddTaskPIDqa();
    PIDQA->SelectCollisionCandidates(kTriggerMask);
  }  
 
  //
  // HMPID QA (Giacomo Volpe)
  //
  if (doHMPID && (ibarrel || iall)) {
    gROOT->LoadMacro("$ALICE_ROOT/PWGPP/HMPID/AddTaskHmpidQA.C");
    AliAnalysisTaskSE* taskhmpidqa= AddTaskHmpidQA(kFALSE);
      // offline mask set in AddTask to kMB
    taskhmpidqa->SelectCollisionCandidates(kTriggerMask);
  }      
  
  
  // T0 QA (Alla Mayevskaya)
  if (doT0 && (ibarrel || iall)) {
  // no offline trigger selection
    gROOT->LoadMacro("$ALICE_ROOT/PWGPP/T0/AddTaskT0QA.C");
    AliT0AnalysisTaskQA* taskt0qa= AddTaskT0QA();
    taskt0qa->SelectCollisionCandidates(kTriggerMask);
  }      
  // FMD QA (Christian Holm Christiansen)
  if (doFMD && (ibarrel || iall)) {
    gROOT->LoadMacro("$ALICE_ROOT/PWGLF/FORWARD/analysis2/AddTaskForwardQA.C");
    // Parameters: usemc, usecentrality
    AliForwardQATask *forwardQA = (AliAnalysisTaskSE *)AddTaskForwardQA(kFALSE, kFALSE);
    forwardQA->SetDebug(0);
    // No offline trigger config. needed (see #84077)
  }
   //     
  // PHOS QA (Boris Polishchuk)
  //
  if (doPHOS) {
    gROOT->LoadMacro("$ALICE_ROOT/PWGGA/PHOSTasks/CaloCellQA/macros/AddTaskCaloCellsQA.C");
    AliAnalysisTaskCaloCellsQA *taskPHOSCellQA1 = AddTaskCaloCellsQA(4, 1, NULL,"PHOSCellsQA_AnyInt"); 
    taskPHOSCellQA1->SelectCollisionCandidates(kTriggerMask);
    taskPHOSCellQA1->GetCaloCellsQA()->SetClusterEnergyCuts(0.3,0.3,1.0);
    AliAnalysisTaskCaloCellsQA *taskPHOSCellQA2 = AddTaskCaloCellsQA(4, 1, NULL,"PHOSCellsQA_PHI7"); 
    taskPHOSCellQA2->SelectCollisionCandidates(AliVEvent::kPHI7);
    taskPHOSCellQA2->GetCaloCellsQA()->SetClusterEnergyCuts(0.3,0.3,1.0);
    // Pi0 QA fo PbPb
    if (iCollisionType) {
      gROOT->LoadMacro("$ALICE_ROOT/PWGGA/PHOSTasks/PHOS_PbPbQA/macros/AddTaskPHOSPbPb.C");
      AliAnalysisTaskPHOSPbPbQA* phosPbPb = AddTaskPHOSPbPbQA(0);
    }
  }    
  if (doPHOSTrig) {
     gROOT->LoadMacro("$ALICE_ROOT/PWGGA/PHOSTasks/PHOS_TriggerQA/macros/AddTaskPHOSTriggerQA.C");
     AliAnalysisTaskPHOSTriggerQA *taskPHOSTrig = AddTaskPHOSTriggerQA(NULL);
  }
  //
  // EMCAL QA (Gustavo Conesa)
  //
  if (doEMCAL) {
     gROOT->LoadMacro("$ALICE_ROOT/PWGGA/EMCALTasks/macros/AddTaskEMCALTriggerQA.C");
     AliAnalysisTaskEMCALTriggerQA *emctrig = AddTaskEMCALTriggerQA();
//     emctrig->GetRecoUtils()->SwitchOffBadChannelsRemoval();
  }   
  //     
  // FLOW and BF QA (C.Perez && A.Rodriguez)
  //
  if (doFBFqa && (ibarrel || iall)) {
    gROOT->LoadMacro("$ALICE_ROOT/PWGPP/macros/AddTaskFBFqa.C");
    AliAnalysisTaskSE *qaFBFMB = (AliAnalysisTaskSE*) AddTaskFBFqa("qaFBFmb",kFALSE);
    qaFBFMB->SelectCollisionCandidates(AliVEvent::kMB);
    AliAnalysisTaskSE *qaFBFSC = (AliAnalysisTaskSE*) AddTaskFBFqa("qaFBFsc",kFALSE);
    qaFBFSC->SelectCollisionCandidates(AliVEvent::kSemiCentral);
    AliAnalysisTaskSE *qaFBFCE = (AliAnalysisTaskSE*) AddTaskFBFqa("qaFBFce",kFALSE);
    qaFBFCE->SelectCollisionCandidates(AliVEvent::kCentral);
  }
}

void QAmerge(const char *suffix, const char *dir, Int_t stage)
{
// Merging method
  TStopwatch timer;
  timer.Start();
  TString outputDir = dir;
  TString outputFiles = Form("QAresults%s.root,EventStat_temp%s.root,RecoQAresults%s.root",suffix,suffix,suffix);
  TString mergeExcludes = "";
  TObjArray *list = outputFiles.Tokenize(",");
  TIter *iter = new TIter(list);
  TObjString *str;
  TString outputFile;
  Bool_t merged = kTRUE;
  while((str=(TObjString*)iter->Next())) {
    outputFile = str->GetString();
    // Skip already merged outputs
    if (!gSystem->AccessPathName(outputFile)) {
       printf("Output file <%s> found. Not merging again.",outputFile.Data());
       continue;
    }
    if (mergeExcludes.Contains(outputFile.Data())) continue;
    merged = AliAnalysisAlien::MergeOutput(outputFile, outputDir, 10, stage);
    if (!merged && !outputFile.Contains("RecoQAresults")) {
       printf("ERROR: Cannot merge %s\n", outputFile.Data());
       continue;
    }
  }
  // read the analysis manager from file
  TString infolog = "fileinfo.log";
  AliAnalysisAlien::MergeInfo(infolog, dir); 
  if (!outputDir.Contains("Stage")) {
    ofstream out;
    out.open("outputs_valid", ios::out);
    out.close();    
    return;
  }
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  mgr->SetRunFromPath(mgr->GetRunFromAlienPath(dir));
  mgr->SetSkipTerminate(kFALSE);
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  AliLog::SetGlobalLogLevel(AliLog::kError);
  TTree *tree = NULL;
  gROOT->cd();
  mgr->StartAnalysis("gridterminate", tree);
  if (strlen(suffix)) {
     if (gSystem->Exec(Form("mv trending.root trending%s.root", suffix)))
        ::Error("QAmerge", "File trending.root was not produced");
     if (gSystem->Exec(Form("mv event_stat.root event_stat%s.root", suffix)))
        ::Error("QAmerge", "File trending.root was not produced");
  }   
  ofstream out;
  out.open("outputs_valid", ios::out);
  out.close();
  timer.Print();
}
