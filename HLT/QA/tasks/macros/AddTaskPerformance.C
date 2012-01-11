///////////////////////////////////////////////////////////////////////////////
// Macro to setup AliPerformanceTask for either 
// HLT or offline TPC performance QA to run on PWGPP QA train. 
//
// Input: ESDs, ESDfriends (optional), Kinematics (optional), TrackRefs (optional)
// ESD and MC input handlers must be attached to AliAnalysisManager
// to run default configuration. 
//
// By default 7 performance components are added to 
// the task: 
// 1. AliPerformanceRes (TPC track resolution w.r.t MC at DCA)
// 2. AliPerformanceResTPCInner (TPC track resolution w.r.t MC at inner TPC wall)
// 3. AliPerformanceResTPCOuter (TPC track resolution w.r.t MC at outer TPC wall)
// 4. AliPerformanceEff (TPC track reconstruction efficiency, MC primaries)
// 5. AliPerformanceDEdxTPCInner (TPC dEdx response - track parameters at TPC inner wall)
// 6. AliPerformanceDCA (TPC impact parameters resolution at DCA)
// 7. AliPerformanceTPC (TPC cluster and track and event information)
// 8. AliPerformanceMatch (TPC and ITS/TRD matching and TPC eff w.r.t ITS)
//
// Usage on the analysis train (default configuration):
// gSystem->Load("libANALYSIS");
// gSystem->Load("libANALYSISalice");
// gSystem->Load("libTPCcalib.so");
// gSystem->Load("libTENDER.so");
// gSystem->Load("libPWGPP.so");
//
// gROOT->LoadMacro("$ALICE_ROOT/HLT/QA/tasks/macros/AddTaskPerformanceTPC.C");
// AliPerformanceTask *tpcQA = AddTaskPerformanceTPC("kTRUE","kTRUE", "kTRUE"¸"triggerClass"); 
// 
// Output:
// HLTTPC.Performance.root file with HLT TPC performance components 
// or 
// TPC.Performance.root file with TPC performance components 
//
// Each of the components contains THnSparse generic histograms which 
// have to be analysed (post-analysis) by using Analyse() function. 
// Each component contains such function.
//
//01-09.10 -  J.Otwinowski@gsi.de, jochen@thaeder.de, hege.erdal@ift.uib.no
///////////////////////////////////////////////////////////////////////////////

//__________________________________________________________________________________________________
AliPerformanceTask* AddTaskPerformance(Bool_t bUseMCInfo=kTRUE,
				       Bool_t bUseESDfriend=kTRUE,
				       Bool_t bUseHLT = kFALSE,
				       const char *triggerClass=0) {
  //					const char *triggerClass="CINT1B-ABCE-NOPF-ALL") {
  //
  // Add AliPerformanceTask with HLT/TPC performance components
  //
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr) { 
    Error("AddTaskPerformance","AliAnalysisManager not set!");
    return NULL;
  }
  
  TString type = mgr->GetInputEventHandler()->GetDataType();
  if (!type.Contains("ESD")) {
    Error("AddTaskPerformance", "ESD input handler needed!");
    return NULL;
  }
  
  AliMCEventHandler *mcH = (AliMCEventHandler*)mgr->GetMCtruthEventHandler();
  if (!mcH && bUseMCInfo) {
    Error("AddTaskPerformance", "MC input handler needed!");
    return NULL;
  }

  //============= Add HLT ESD ================================================================
  if(bUseHLT){
    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*>(mgr->GetInputEventHandler());
    esdH->SetReadHLT();
  }
  //==========================================================================================

  //
  // Create task -----------------------------------------------------------------------------------
  //
  if(bUseHLT)
    AliPerformanceTask *task = new AliPerformanceTask("Performance","HLT TPC Performance");
  else
    AliPerformanceTask *task = new AliPerformanceTask("Performance","TPC Performance");

  if (!task) {
    Error("AddTaskPerformance", "Performance task cannot be created!");
    return NULL;
  }

  task->SelectCollisionCandidates();
  task->SetUseMCInfo(bUseMCInfo);
  task->SetUseESDfriend(bUseESDfriend);
  task->SetUseTerminate(kFALSE);
  task->SetUseHLT(bUseHLT);
  //
  // Add task to analysis manager ------------------------------------------------------------------
  //
  
  mgr->AddTask(task);

  //================================================================================================
  // -- Cuts
  //================================================================================================

  //
  // Create TPC-ESD track reconstruction cuts ------------------------------------------------------
  // MB tracks
  //

  AliRecInfoCuts *pRecInfoCutsTPC = new AliRecInfoCuts(); 
  if(pRecInfoCutsTPC) {
    //    pRecInfoCutsTPC->SetMaxDCAToVertexXY(3.0);
    //    pRecInfoCutsTPC->SetMaxDCAToVertexZ(30.0);
    pRecInfoCutsTPC->SetMaxDCAToVertexXY(150.0);
    pRecInfoCutsTPC->SetMaxDCAToVertexZ(150.0);
    pRecInfoCutsTPC->SetRequireSigmaToVertex(kFALSE);
    pRecInfoCutsTPC->SetRequireTPCRefit(kFALSE);

    pRecInfoCutsTPC->SetAcceptKinkDaughters(kTRUE);
    pRecInfoCutsTPC->SetMinNClustersTPC(50);
    pRecInfoCutsTPC->SetMaxChi2PerClusterTPC(1000000.);
    pRecInfoCutsTPC->SetDCAToVertex2D(kFALSE);

    pRecInfoCutsTPC->SetHistogramsOn(kFALSE); 
  } 
  else {
    Error("AddTaskPerformance", "AliRecInfoCutsTPC cannot be created!");
    return NULL;
  }

  //
  // Create TPC-ESD track reconstruction cuts ------------------------------------------------------
  // standard cuts

  AliRecInfoCuts *pRecInfoCuts = new AliRecInfoCuts(); 
  if(pRecInfoCuts) {
    pRecInfoCuts->SetMaxDCAToVertexXY(3.0);
    pRecInfoCuts->SetMaxDCAToVertexZ(3.0);
    pRecInfoCuts->SetMinNClustersTPC(50);
    pRecInfoCuts->SetMinNClustersITS(2);
    pRecInfoCuts->SetHistogramsOn(kFALSE); 
    pRecInfoCuts->SetTPCITSMatchingRadius(70); 
    pRecInfoCuts->SetTPCTRDMatchingRadius(260); 
  } 
  else {
    Error("AddTaskPerformance", "AliRecInfoCuts cannot be created!");
    return NULL;
  }

  //
  // Create TPC-MC track reconstruction cuts -------------------------------------------------------
  //
  AliMCInfoCuts  *pMCInfoCuts = new AliMCInfoCuts();
  if(pMCInfoCuts) {
    pMCInfoCuts->SetMinTrackLength(70);
  } 
  else {
    Error("AddTaskPerformance", "AliMCInfoCuts cannot be created!");
    return NULL;
  }

  //================================================================================================
  // -- Object
  //================================================================================================

  //
  // Create performance objects for TPC and set cuts 
  //
  enum { kTPC = 0, kTPCITS, kConstrained, kTPCInner, kTPCOuter, kTPCSec };

  //
  // Resolution ------------------------------------------------------------------------------------
  //

  AliPerformanceRes *pCompRes0 = new AliPerformanceRes("AliPerformanceRes",
						       "AliPerformanceRes",kTPC,kFALSE); 
  if(!pCompRes0) {
    Error("AddTaskPerformance", "Cannot create AliPerformanceRes");
  }
  pCompRes0->SetAliRecInfoCuts(pRecInfoCuts);
  pCompRes0->SetAliMCInfoCuts(pMCInfoCuts);
  pCompRes0->SetUseTrackVertex(kTRUE);

  //
  // Efficiency ------------------------------------------------------------------------------------
  //

  AliPerformanceEff *pCompEff0 = new AliPerformanceEff("AliPerformanceEff",
						       "AliPerformanceEff",kTPC,kFALSE); 
  if(!pCompEff0) {
    Error("AddTaskPerformance", "Cannot create AliPerformanceEff");
  }
  pCompEff0->SetAliRecInfoCuts(pRecInfoCuts);
  pCompEff0->SetAliMCInfoCuts(pMCInfoCuts);
  pCompEff0->SetUseTrackVertex(kTRUE);

  //
  // dEdx ------------------------------------------------------------------------------------------
  //

  AliPerformanceDEdx *pCompDEdx3 = new AliPerformanceDEdx("AliPerformanceDEdxTPCInner",
							  "AliPerformanceDEdxTPCInner",kTPCInner,kFALSE); 
  if(!pCompDEdx3) {
    Error("AddTaskPerformance", "Cannot create AliPerformanceDEdxTPCInner");
  }
  pCompDEdx3->SetAliRecInfoCuts(pRecInfoCuts);
  pCompDEdx3->SetAliMCInfoCuts(pMCInfoCuts);
  pCompDEdx3->SetUseTrackVertex(kTRUE);

  //
  // DCA -------------------------------------------------------------------------------------------
  //

  AliPerformanceDCA *pCompDCA0 = new AliPerformanceDCA("AliPerformanceDCA",
						       "AliPerformanceDCA",kTPC,kFALSE); 
  if(!pCompDCA0) {
    Error("AddTaskPerformanceTPC", "Cannot create AliPerformanceDCA");
  }
  pCompDCA0->SetAliRecInfoCuts(pRecInfoCuts);
  pCompDCA0->SetAliMCInfoCuts(pMCInfoCuts);
  pCompDCA0->SetUseTrackVertex(kTRUE);

  //
  // TPC performance -------------------------------------------------------------------------------
  //

  AliPerformanceTPC *pCompTPC0 = new AliPerformanceTPC("AliPerformanceTPC",
						       "AliPerformanceTPC",kTPC,kFALSE); 
  if(!pCompTPC0) {
    Error("AddTaskPerformance", "Cannot create AliPerformanceTPC");
  }
  pCompTPC0->SetAliRecInfoCuts(pRecInfoCutsTPC);
  pCompTPC0->SetAliMCInfoCuts(pMCInfoCuts);
  pCompTPC0->SetUseTrackVertex(kTRUE);
  pCompTPC0->SetUseHLT(bUseHLT);

  //
  // Add components to the performance task --------------------------------------------------------
  //

  //01.09.10  At the moment only use pCompTPC0

  //task->AddPerformanceObject( pCompDEdx3 );
  //task->AddPerformanceObject( pCompDCA0 );
  task->AddPerformanceObject( pCompTPC0 );

  /*
    if(bUseMCInfo) { 
    task->AddPerformanceObject( pCompRes0 );
    task->AddPerformanceObject( pCompEff0 );
    }
    
    if(!bUseMCInfo &&  triggerClass) {
    pCompDEdx3->SetTriggerClass(triggerClass);
    pCompDCA0->SetTriggerClass(triggerClass);
    pCompTPC0->SetTriggerClass(triggerClass);
    }
  */
  
  //================================================================================================
  //              data containers
  //================================================================================================
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  
  if(bUseHLT){
    AliAnalysisDataContainer *coutput1 =    
      mgr->CreateContainer("HLTQA", TList::Class(),
			   AliAnalysisManager::kOutputContainer, Form("HLTTPC.%s.root", task->GetName()));
  }
  else{
    AliAnalysisDataContainer *coutput1 = 
      mgr->CreateContainer("TPC", TList::Class(), 
			   AliAnalysisManager::kOutputContainer, Form("TPC.%s.root", task->GetName()));
  }

  // -- connect containers
  mgr->ConnectInput  (task,  0, cinput );
  mgr->ConnectOutput (task,  1, coutput1);

  return task;  
}
