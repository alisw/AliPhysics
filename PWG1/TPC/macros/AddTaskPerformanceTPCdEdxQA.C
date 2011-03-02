///////////////////////////////////////////////////////////////////////////////
// Macro to setup AliPerformanceTask for 
// TPC performance QA to run on PWG1 QA train. 
//
// Input: ESDs, ESDfriends (optional), Kinematics (optional), TrackRefs (optional)
// ESD and MC input handlers must be attached to AliAnalysisManager
// to run default configuration. 
//
// By default 1 performance components are added to 
// the task: 
// 1. AliPerformanceTPC (TPC cluster and track and event information)
// 2. AliPerformancedEdx (TPC dEdx information)
//
// Usage on the analysis train (default configuration):
// gSystem->Load("libANALYSIS");
// gSystem->Load("libANALYSISalice");
// gSystem->Load("libTPCcalib.so");
// gSystem->Load("libTENDER.so");
// gSystem->Load("libPWG1.so");
//
// gROOT->LoadMacro("$ALICE_ROOT/PWG1/TPC/macros/AddTaskPerformanceTPCdEdxQA.C");
// AliPerformanceTask *tpcQA = AddTaskPerformanceTPCdEdxQA("kFALSE","kTRUE","kFALSE","triggerClass"); 
// 
// Output:
// TPC.Performance.root file with TPC performance components is created.
//
// Each of the components contains THnSparse generic histograms which 
// have to be analysed (post-analysis) by using Analyse() function. 
// Each component contains such function.
//
//30.09.2010 -  J.Otwinowski@gsi.de
///////////////////////////////////////////////////////////////////////////////

//____________________________________________
AliPerformanceTask* AddTaskPerformanceTPCdEdxQA(Bool_t bUseMCInfo=kFALSE, Bool_t bUseESDfriend=kTRUE, Bool_t highMult = kFALSE, const char *triggerClass=0)
{
  //
  // Add AliPerformanceTask with TPC performance components
  //
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr) { 
    Error("AddTaskPerformanceTPCdEdxQA","AliAnalysisManager not set!");
    return NULL;
  }
  
  TString type = mgr->GetInputEventHandler()->GetDataType();
  if (!type.Contains("ESD")) {
    Error("AddTaskPerformanceTPCdEdxQA", "ESD input handler needed!");
    return NULL;
  }
  
  AliMCEventHandler *mcH = (AliMCEventHandler*)mgr->GetMCtruthEventHandler();
  if (!mcH && bUseMCInfo) {
    Error("AddTaskPerformanceTPCdEdxQA", "MC input handler needed!");
    return NULL;
  }

  //
  // Create task
  //
  AliPerformanceTask *task = new AliPerformanceTask("PerformanceQA","TPC Performance");
  if (!task) {
    Error("AddTaskPerformanceTPCdEdxQA", "TPC performance task cannot be created!");
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
  // MB tracks
  //
  AliRecInfoCuts *pRecInfoCutsTPC = new AliRecInfoCuts(); 
  if(pRecInfoCutsTPC) {
    pRecInfoCutsTPC->SetMaxDCAToVertexXY(3.0);
    pRecInfoCutsTPC->SetMaxDCAToVertexZ(3.0);
    pRecInfoCutsTPC->SetRequireSigmaToVertex(kFALSE);
    pRecInfoCutsTPC->SetRequireTPCRefit(kFALSE);
    pRecInfoCutsTPC->SetAcceptKinkDaughters(kTRUE);
    pRecInfoCutsTPC->SetMinNClustersTPC(70);
    pRecInfoCutsTPC->SetMaxChi2PerClusterTPC(1000000.);
    pRecInfoCutsTPC->SetDCAToVertex2D(kFALSE);

    pRecInfoCutsTPC->SetHistogramsOn(kFALSE); 
  } 
  else {
    Error("AddTaskPerformanceTPC", "AliRecInfoCutsTPC cannot be created!");
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
  // TPC performance
  //
  AliPerformanceTPC *pCompTPC0 = new AliPerformanceTPC("AliPerformanceTPC","AliPerformanceTPC",kTPC,kFALSE,-1,highMult); 
  if(!pCompTPC0) {
    Error("AddTaskPerformanceTPC", "Cannot create AliPerformanceTPC");
  }
  pCompTPC0->SetAliRecInfoCuts(pRecInfoCutsTPC);
  pCompTPC0->SetAliMCInfoCuts(pMCInfoCuts);
  pCompTPC0->SetUseTrackVertex(kFALSE);
  pCompTPC0->SetUseHLT(kFALSE);
  //
  // dEdx
  //
  AliPerformanceDEdx *pCompDEdx3 = new AliPerformanceDEdx("AliPerformanceDEdxTPCInner","AliPerformanceDEdxTPCInner",kTPCInner,kFALSE); 
  if(!pCompDEdx3) {
    Error("AddTaskPerformanceTPC", "Cannot create AliPerformanceDEdxTPCInner");
  }
  pCompDEdx3->SetAliRecInfoCuts(pRecInfoCutsTPC);
  pCompDEdx3->SetAliMCInfoCuts(pMCInfoCuts);
  pCompDEdx3->SetUseTrackVertex(kFALSE);


  //
  // Add components to the performance task
  //
  if(!bUseMCInfo) { 
    pCompTPC0->SetTriggerClass(triggerClass);
    pCompDEdx3->SetTriggerClass(triggerClass);
  }
  task->AddPerformanceObject( pCompTPC0 );
  task->AddPerformanceObject( pCompDEdx3 );

  //
  // Create containers for input
  //
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

  //
  // Create containers for output
  //
  AliAnalysisDataContainer *coutput_tpc = mgr->CreateContainer("TPCQA", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:TPC_%s", mgr->GetCommonFileName(), task->GetName()));

  AliAnalysisDataContainer *coutput2_tpc = mgr->CreateContainer("TPCQASummary", TTree::Class(), AliAnalysisManager::kParamContainer, "trending.root:SummaryTPCQA"); 

  mgr->ConnectOutput(task, 1, coutput_tpc);
  mgr->ConnectOutput(task, 0, coutput2_tpc);

return task;  
}
