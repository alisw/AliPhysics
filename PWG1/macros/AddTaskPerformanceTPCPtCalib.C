///////////////////////////////////////////////////////////////////////////////
// Macro to setup AliPerformanceTask for 
// TPC performance QA to run on PWG1 QA train. 
//
// By default 1 performance component is added to 
// the task: 
// 1. AliPerformancePtCalib
// or AliPerformancePtCalibMC if bUseMCinfo = kTRUE (use MC info)

// Usage on the analysis train (default configuration):
// gSystem->Load("libANALYSIS");
// gSystem->Load("libANALYSISalice");
// gSystem->Load("libTPCcalib.so");
// gSystem->Load("libPWG1.so");
//
// gROOT->LoadMacro("$ALICE_ROOT/PWG1/macros/AddTaskPerformanceTPCPtCalib.C");
// AliPerformanceTask *tpcQA = AddTaskPerformanceTPCPtCalib("kTRUE","kTRUE","CINT1B-ABCE-NOPF-ALL"); 
// 
// Output:
// TPCPtCalib.Performance.root file with TPC performance components is created.
//
// Each of the components contains THnSparse generic histograms which 
// have to be analysed (post-analysis) by using Analyse() function. 
// Each component contains such function.
//
// June 2010 -  Simone Schuchmann sschuchm@ikf.uni-frankfurt.de
///////////////////////////////////////////////////////////////////////////////

//____________________________________________
AliPerformanceTask* AddTaskPerformanceTPCPtCalib(Bool_t bUseMCInfo=kFALSE, Bool_t bUseESDfriend=kFALSE, const char *triggerClass=0)
{

//
  // Create physics trigger selection class
  //
    AliPhysicsSelection *physTrigSel =  new AliPhysicsSelection();
  
  //
  // Add AliPerformanceTask with TPC performance components
  //
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr) { 
     Error("AddTaskPerformanceTPCPtCalib","AliAnalysisManager not set!");
    return NULL;
  }
  
  TString type = mgr->GetInputEventHandler()->GetDataType();
  if (!type.Contains("ESD")) {
     Error("AddTaskPerformanceTPCPtCalib", "ESD input handler needed!");
    return NULL;
  }
  
  AliMCEventHandler *mcH = (AliMCEventHandler*)mgr->GetMCtruthEventHandler();
  if (!mcH && bUseMCInfo) {
    Error("AddTaskPerformanceTPCPtCalib", "MC input handler needed!");
    return NULL;
  }

  //
  // Create task
  //
  AliPerformanceTask *task = new AliPerformanceTask("TPCPerformanceInvPt","TPC Performance PtCalib");
  if (!task) {
    Error("AddTaskPerformanceTPCPtCalib", "TPC performance task cannot be created!");
    return NULL;
  }
  task->SetUseMCInfo(bUseMCInfo);
  task->SetUseESDfriend(kFALSE);
  task->SelectCollisionCandidates();

  //
  // Add physics selection task to analysis manager
  //
  //    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  //    AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection();
  //    mgr->AddTask(physSelTask);

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
    Error("AddTaskPerformanceTPCPtCalib", "AliRecInfoCuts cannot be created!");
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
    Error("AddTaskPerformanceTPCPtCalib", "AliMCInfoCuts cannot be created!");
    return NULL;
  }

  //
  // Create performance objects for TPC and set cuts 
  //
  enum { kTPC = 0, kTPCITS, kConstrained, kTPCInner, kTPCOuter, kTPCSec };



  AliPerformancePtCalib *ptCalib =  NULL;
  AliPerformancePtCalibMC *ptCalibMC = NULL;

  if(bUseMCInfo){
     ptCalibMC = new AliPerformancePtCalibMC("AliPerformancePtCalibMC","AliPerformancePtCalibMC");
     if(!ptCalibMC) {
	Error("AddTaskPerformanceTPCPtCalib", "Cannot create AliPerformancePtCalibMC");
     }
     // physTrigSel->SetAnalyzeMC();
     // ptCalibMC->SetPhysicsTriggerSelection(physTrigSel); 
     ptCalibMC->SetAliRecInfoCuts(pRecInfoCuts);
     ptCalibMC->SetReadTPCTracks(kTRUE);  
     AliESDtrackCuts* esdtrackCuts = new AliESDtrackCuts("AliESDtrackCutsPtMC");
     esdtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2009(kTRUE);
     ptCalibMC->SetAliESDtrackCuts(esdtrackCuts);
     //ptCalibMC->SetEtaRange(0.9);
  }
  else{

     ptCalib =  new AliPerformancePtCalib("AliPerformancePtCalib","AliPerformancePtCalib");
     if(!ptCalib) {
	Error("AddTaskPerformanceTPCPtCalib", "Cannot create AliPerformancePtCalib");
     }
     ptCalib->SetAliRecInfoCuts(pRecInfoCuts);
     ptCalib->SetReadTPCTracks(kTRUE);  
     ptCalib->SetEtaRange(0.8);
     // ptCalib->SetAliMCInfoCuts(pMCInfoCut);

     //if(triggerClass) ptCalib->SetTriggerClass(triggerClass);
     //ptCalib->SetPhysicsTriggerSelection(physTrigSel);
     //ptCalib->SetTrigger(AliTriggerAnalysis::kMB1);
     AliESDtrackCuts* esdtrackCuts = new AliESDtrackCuts("AliESDtrackCutsPt");
     esdtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2009(kTRUE);
     ptCalib->SetAliESDtrackCuts(esdtrackCuts);
      
  }
  
  // add components to the performance task
  
  if(bUseMCInfo) task->AddPerformanceObject(ptCalibMC);
  else task->AddPerformanceObject(ptCalib);
  
  // Create containers for input
  //
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  
  //
  // Create containers for output
  //
  AliAnalysisDataContainer *coutput_tpcptcalib = mgr->CreateContainer("TPCPtCalib", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s.Performance.root", task->GetName()));
  mgr->ConnectOutput(task, 1, coutput_tpcptcalib);

return task;  
}
