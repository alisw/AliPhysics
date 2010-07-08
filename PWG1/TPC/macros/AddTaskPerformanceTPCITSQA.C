///////////////////////////////////////////////////////////////////////////////
// Macro to setup AliPerformanceTask for 
// TPC and TPC+ITS performance QA to run on PWG1 QA train. 
//
// Input: ESDs, ESDfriends (optional), Kinematics (optional), TrackRefs (optional)
// ESD and MC input handlers must be attached to AliAnalysisManager
// to run default configuration. 
//
// By default 1 performance components are added to 
// the task: 
// 1. AliPerformanceTPC (TPC cluster and track and event information)
//
// Usage on the analysis train (default configuration):
// gSystem->Load("libANALYSIS");
// gSystem->Load("libANALYSISalice");
// gSystem->Load("libTPCcalib.so");
// gSystem->Load("libTENDER.so");
// gSystem->Load("libPWG1.so");
//
// gROOT->LoadMacro("$ALICE_ROOT/PWG1/macros/AddTaskPerformanceTPCITSQA.C");
// AliPerformanceTask *tpcQA = AddTaskPerformanceTPC("kFALSE","kTRUE","triggerClass"); 
// 
// Output:
// TPCITS.Performance.root file with TPC and TPCITS performance components is created.
//
// Each of the components contains THnSparse generic histograms which 
// have to be analysed (post-analysis) by using Analyse() function. 
// Each component contains such function.
//
//14.12.2009 -  J.Otwinowski@gsi.de
///////////////////////////////////////////////////////////////////////////////

//____________________________________________
AliPerformanceTask* AddTaskPerformanceTPCITSQA(Bool_t bUseMCInfo=kFALSE, Bool_t bUseESDfriend=kTRUE, const char *triggerClass=0)
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
  AliPerformanceTask *task = new AliPerformanceTask("PerformanceQA","TPC Performance");
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
  // Create TPC-ITS ESD track reconstruction cuts
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
    pRecInfoCutsTPC->SetMinNClustersITS(2);

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
  // Create performance objects for TPC and TPCITS set cuts 
  //
  enum { kTPC = 0, kTPCITS};

  //
  // TPC performance
  //
  AliPerformanceTPC *pCompTPC0 = new AliPerformanceTPC("AliPerformanceTPC","AliPerformanceTPC",kTPC,kFALSE); 
  if(!pCompTPC0) {
    Error("AddTaskPerformanceTPC", "Cannot create AliPerformanceTPC");
  }
  pCompTPC0->SetAliRecInfoCuts(pRecInfoCutsTPC);
  pCompTPC0->SetAliMCInfoCuts(pMCInfoCuts);

  //
  // TPC+ITS performance
  //
  AliPerformanceTPC *pCompTPC1 = new AliPerformanceTPC("AliPerformanceTPCITS","AliPerformanceTPCITS",kTPCITS,kFALSE); 
  if(!pCompTPC1) {
    Error("AddTaskPerformanceTPC", "Cannot create AliPerformanceTPCITS");
  }
  pCompTPC1->SetAliRecInfoCuts(pRecInfoCutsTPC);
  pCompTPC1->SetAliMCInfoCuts(pMCInfoCuts);

  //
  // Add components to the performance task
  //
  if(!bUseMCInfo) { 
    pCompTPC0->SetTriggerClass(triggerClass);
    pCompTPC1->SetTriggerClass(triggerClass);
  }
  task->AddPerformanceObject( pCompTPC0 );
  task->AddPerformanceObject( pCompTPC1 );

  //
  // Create containers for input
  //
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

  //
  // Create containers for output
  //
  AliAnalysisDataContainer *coutput_tpc = mgr->CreateContainer("TPCITSQA", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:TPCITS_%s", mgr->GetCommonFileName(), task->GetName()));
  mgr->ConnectOutput(task, 1, coutput_tpc);

return task;  
}
