///////////////////////////////////////////////////////////////////////////////
// Macro to setup AliPerformanceTask for 
// TPC performance QA to run on PWG1 QA train. 
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
// 7. AliPerformanceTPC (TPC cluster and track information)
//
// Usage on the analysis train (default configuration):
// gSystem->Load("libANALYSIS");
// gSystem->Load("libANALYSISalice");
// gSystem->Load("libTPCcalib.so");
// gSystem->Load("libPWG1.so");
//
// gROOT->LoadMacro("$ALICE_ROOT/PWG1/macros/AddTaskPerformanceTPC.C");
// AliPerformanceTask *tpcQA = AddTaskPerformanceTPC("kTRUE","kTRUE"); 
// 
// Output:
// TPC.Performance.root file with TPC performance components is created.
//
// Each of the components contains THnSparse generic histograms which 
// have to be analysed (post-analysis) by using Analyse() function. 
// Each component contains such function.
//
//13.10.2009 -  J.Otwinowski@gsi.de
///////////////////////////////////////////////////////////////////////////////

//____________________________________________
AliPerformanceTask* AddTaskPerformanceTPC(Bool_t bUseMCInfo=kTRUE, Bool_t bUseESDfriend=kTRUE)
{
  //
  // Add AliPerformanceTask with TPC performance components
  //
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr) { 
    Error("AddTaskPerformanceTPC","AliAnalysisManager not set!");
    return NULL;
  }
  
  TString type = mgr->GetInputEventHandler()->GetDataType();
  if (!type.Contains("ESD")) {
    Error("AddTaskPerformanceTPC", "ESD input handler needed!");
    return NULL;
  }
  
  AliMCEventHandler *mcH = (AliMCEventHandler*)mgr->GetMCtruthEventHandler();
  if (!mcH && bUseMCInfo) {
    Error("AddTaskPerformanceTPC", "MC input handler needed!");
    return NULL;
  }

  //
  // Create task
  //
  AliPerformanceTask *task = new AliPerformanceTask("Performance","TPC Performance");
  if (!task) {
    Error("AddTaskPerformanceTPC", "TPC performance task cannot be created!");
    return NULL;
  }
  task->SetUseMCInfo(bUseMCInfo);
  task->SetUseESDfriend(bUseESDfriend);

  //
  // Add task to analysis manager
  //
  mgr->AddTask(task);

  //
  // Create TPC-ESD track reconstruction cuts
  //
  AliRecInfoCuts *pRecInfoCuts = new AliRecInfoCuts(); 
  if(pRecInfoCuts) {
    pRecInfoCuts->SetMaxDCAToVertexXY(3.0);
    pRecInfoCuts->SetMaxDCAToVertexZ(3.0);
    pRecInfoCuts->SetMinNClustersTPC(50);
    pRecInfoCuts->SetMinNClustersITS(2);
    pRecInfoCuts->SetHistogramsOn(kFALSE); 
  } 
  else {
    Error("AddTaskPerformanceTPC", "AliRecInfoCuts cannot be created!");
    return NULL;
  }
  //
  // Create TPC-MC track reconstruction cuts
  //
  AliMCInfoCuts  *pMCInfoCuts = new AliMCInfoCuts();
  if(pMCInfoCuts) {
    pMCInfoCuts->SetMinTrackLength(70);
  } 
  else {
    Error("AddTaskPerformanceTPC", "AliMCInfoCuts cannot be created!");
    return NULL;
  }

  //
  // Create performance objects for TPC and set cuts 
  //
  enum { kTPC = 0, kTPCITS, kConstrained, kTPCInner, kTPCOuter, kTPCSec };

  //
  // Resolution
  //
  AliPerformanceRes *pCompRes0 = new AliPerformanceRes("AliPerformanceRes","AliPerformanceRes",kTPC,kFALSE); 
  if(!pCompRes0) {
    Error("AddTaskPerformanceTPC", "Cannot create AliPerformanceRes");
  }
  pCompRes0->SetAliRecInfoCuts(pRecInfoCuts);
  pCompRes0->SetAliMCInfoCuts(pMCInfoCuts);

  AliPerformanceRes *pCompRes3 = new AliPerformanceRes("AliPerformanceResTPCInner","AliPerformanceResTPCInner",kTPCInner,kFALSE); 
  if(!pCompRes3) {
    Error("AddTaskPerformanceTPC", "Cannot create AliPerformanceResTPCInner");
  }
  pCompRes3->SetAliRecInfoCuts(pRecInfoCuts);
  pCompRes3->SetAliMCInfoCuts(pMCInfoCuts);

  AliPerformanceRes *pCompRes4 = new AliPerformanceRes("AliPerformanceResTPCOuter","AliPerformanceResTPCOuter",kTPCOuter,kFALSE); 
  if(!pCompRes4) {
    Error("AddTaskPerformanceTPC", "Cannot create AliPerformanceResTPCOuter");
  }
  pCompRes4->SetAliRecInfoCuts(pRecInfoCuts);
  pCompRes4->SetAliMCInfoCuts(pMCInfoCuts);

  //
  // Efficiency
  //
  AliPerformanceEff *pCompEff0 = new AliPerformanceEff("AliPerformanceEff","AliPerformanceEff",kTPC,kFALSE); 
  if(!pCompEff0) {
    Error("AddTaskPerformanceTPC", "Cannot create AliPerformanceEff");
  }
  pCompEff0->SetAliRecInfoCuts(pRecInfoCuts);
  pCompEff0->SetAliMCInfoCuts(pMCInfoCuts);
  //
  // dEdx
  //
  AliPerformanceDEdx *pCompDEdx3 = new AliPerformanceDEdx("AliPerformanceDEdxTPCInner","AliPerformanceDEdxTPCInner",kTPCInner,kFALSE); 
  if(!pCompDEdx3) {
    Error("AddTaskPerformanceTPC", "Cannot create AliPerformanceDEdxTPCInner");
  }
  pCompDEdx3->SetAliRecInfoCuts(pRecInfoCuts);
  pCompDEdx3->SetAliMCInfoCuts(pMCInfoCuts);
  //
  // DCA
  //
  AliPerformanceDCA *pCompDCA0 = new AliPerformanceDCA("AliPerformanceDCA","AliPerformanceDCA",kTPC,kFALSE); 
  if(!pCompDCA0) {
    Error("AddTaskPerformanceTPC", "Cannot create AliPerformanceDCA");
  }
  pCompDCA0->SetAliRecInfoCuts(pRecInfoCuts);
  pCompDCA0->SetAliMCInfoCuts(pMCInfoCuts);
  //
  // TPC performance
  //
  AliPerformanceTPC *pCompTPC0 = new AliPerformanceTPC("AliPerformanceTPC","AliPerformanceTPC",kTPC,kFALSE); 
  if(!pCompTPC0) {
    Error("AddTaskPerformanceTPC", "Cannot create AliPerformanceTPC");
  }
  pCompTPC0->SetAliRecInfoCuts(pRecInfoCuts);
  pCompTPC0->SetAliMCInfoCuts(pMCInfoCuts);
 
  //
  // add components to the performance task
  //
  task->AddPerformanceObject( pCompRes0 );
  task->AddPerformanceObject( pCompRes3 );
  task->AddPerformanceObject( pCompRes4 );
  task->AddPerformanceObject( pCompEff0 );
  task->AddPerformanceObject( pCompDEdx3 );
  task->AddPerformanceObject( pCompDCA0 );
  task->AddPerformanceObject( pCompTPC0 );

  //
  // Create containers for input
  //
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

  //
  // Create containers for output
  //
  AliAnalysisDataContainer *coutput_tpc = mgr->CreateContainer("TPC", TList::Class(), AliAnalysisManager::kOutputContainer, Form("TPC.%s.root", task->GetName()));
  mgr->ConnectOutput(task, 0, coutput_tpc);

return task;  
}
