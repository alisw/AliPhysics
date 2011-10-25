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
// 0. AliPerformanceTPC (TPC cluster and track and event information)
// 1. AliPerformanceMatch (TPC and ITS/TRD matching and TPC eff w.r.t ITS)
// 2. AliPerformanceMatch (TPC and ITS/TRD matching and TPC eff w.r.t TPC)
// 3. AliPerformancedEdx (TPC dEdx information)
// 4. AliPerformanceRes (TPC track resolution w.r.t MC at DCA)
// 5. AliPerformanceEff (TPC track reconstruction efficiency, MC primaries)
//
// Usage on the analysis train (default configuration):
// gSystem->Load("libANALYSIS");
// gSystem->Load("libANALYSISalice");
// gSystem->Load("libTPCcalib.so");
// gSystem->Load("libTENDER.so");
// gSystem->Load("libPWG1.so");
//
// gROOT->LoadMacro("$ALICE_ROOT/PWG1/TPC/macros/AddTaskPerformanceTPCdEdxQA.C");
// AliPerformanceTask *tpcQA = AddTaskPerformanceTPCdEdxQA("kFALSE","kTRUE","kFALSE","triggerClass",kFALSE); 
// 
// Output:
// TPC.Performance.root file with TPC performance components is created.
//
// Each of the components contains THnSparse generic histograms which 
// have to be analysed (post-analysis) by using Analyse() function. 
// Each component contains such function.
//
//30.09.2010 -  J.Otwinowski@gsi.de
//22.09.2011 -  jochen@thaeder.de - Updated
///////////////////////////////////////////////////////////////////////////////

//____________________________________________
AliPerformanceTask* AddTaskPerformanceTPCdEdxQA(Bool_t bUseMCInfo=kFALSE, Bool_t bUseESDfriend=kTRUE, 
						Bool_t highMult = kFALSE, const char *triggerClass=0, 
						Bool_t bUseHLT = kFALSE)
{
  Char_t *taskName[] = {"TPC", "HLT"};
  Int_t idx = 0;
  if (bUseHLT) idx = 1;
  
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
  // Add HLT Event
  //
  if (bUseHLT) {
    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*>(mgr->GetInputEventHandler());
    esdH->SetReadHLT();
  }

  //
  // Create task
  //
  AliPerformanceTask *task = new AliPerformanceTask("PerformanceQA",Form("%s Performance",taskName[idx]));
  if (!task) {
    Error("AddTaskPerformanceTPCdEdxQA", Form("%s performance task cannot be created!",taskName[idx]));
    return NULL;
  }
  task->SetUseMCInfo(bUseMCInfo);
  task->SetUseESDfriend(bUseESDfriend);
  //  task->SetUseTerminate(kFALSE);
  task->SetUseHLT(bUseHLT);

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
    pRecInfoCutsTPC->SetAcceptKinkDaughters(kFALSE);
    pRecInfoCutsTPC->SetMinNClustersTPC(70);
    pRecInfoCutsTPC->SetMaxChi2PerClusterTPC(4.);
    pRecInfoCutsTPC->SetDCAToVertex2D(kTRUE);

    pRecInfoCutsTPC->SetHistogramsOn(kFALSE); 
  } 
  else {
    Error("AddTaskPerformanceTPCdEdxQA", "AliRecInfoCutsTPC cannot be created!");
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
    Error("AddTaskPerformanceTPCdEdxQA", "AliMCInfoCuts cannot be created!");
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
    Error("AddTaskPerformanceTPCdEdxQA", "Cannot create AliPerformanceTPC");
  }
  pCompTPC0->SetAliRecInfoCuts(pRecInfoCutsTPC);
  pCompTPC0->SetAliMCInfoCuts(pMCInfoCuts);
  //  pCompTPC0->SetUseTrackVertex(kFALSE);
  pCompTPC0->SetUseTrackVertex(kTRUE);
  pCompTPC0->SetUseHLT(bUseHLT);
  pCompTPC0->SetUseTOFBunchCrossing(kTRUE);
  
  //
  // TPC ITS match
  //
  AliPerformanceMatch *pCompMatch1 = new AliPerformanceMatch("AliPerformanceMatchTPCITS","AliPerformanceMatchTPCITS",0,kFALSE); 
  if(!pCompMatch1) {
    Error("AddTaskPerformanceTPCdEdxQA", "Cannot create AliPerformanceMatchTPCITS");
  }
  pCompMatch1->SetAliRecInfoCuts(pRecInfoCutsTPC);
  pCompMatch1->SetAliMCInfoCuts(pMCInfoCuts);
  pCompMatch1->SetUseTOFBunchCrossing(kTRUE);


  //
  // ITS TPC match
  //
  AliPerformanceMatch *pCompMatch2 = new AliPerformanceMatch("AliPerformanceMatchITSTPC","AliPerformanceMatchITSTPC",1,kFALSE); 
  if(!pCompMatch2) {
    Error("AddTaskPerformanceTPCdEdxQA", "Cannot create AliPerformanceMatchITSTPC");  }
  pCompMatch2->SetAliRecInfoCuts(pRecInfoCutsTPC);
  pCompMatch2->SetAliMCInfoCuts(pMCInfoCuts);
  pCompMatch2->SetUseTOFBunchCrossing(kTRUE);

  //
  // dEdx
  //
  AliPerformanceDEdx *pCompDEdx3 = new AliPerformanceDEdx("AliPerformanceDEdxTPCInner","AliPerformanceDEdxTPCInner",kTPCInner,kFALSE); 
  if(!pCompDEdx3) {
    Error("AddTaskPerformanceTPCdEdxQA", "Cannot create AliPerformanceDEdxTPCInner");
  }
  pCompDEdx3->SetAliRecInfoCuts(pRecInfoCutsTPC);
  pCompDEdx3->SetAliMCInfoCuts(pMCInfoCuts);
  //pCompDEdx3->SetUseTrackVertex(kFALSE);
  pCompDEdx3->SetUseTrackVertex(kTRUE);

  //
  // Resolution ------------------------------------------------------------------------------------
  //

  AliPerformanceRes *pCompRes4 = new AliPerformanceRes("AliPerformanceRes",
						       "AliPerformanceRes",kTPC,kFALSE); 
  if(!pCompRes4) {
    Error("AddTaskPerformanceTPC", "Cannot create AliPerformanceRes");
  }


  pCompRes4->SetAliRecInfoCuts(pRecInfoCutsTPC);
  pCompRes4->SetAliMCInfoCuts(pMCInfoCuts);
  pCompRes4->SetUseTrackVertex(kTRUE);

  //
  // Efficiency ------------------------------------------------------------------------------------
  //

  AliPerformanceEff *pCompEff5 = new AliPerformanceEff("AliPerformanceEff",
						       "AliPerformanceEff",kTPC,kFALSE); 
  if(!pCompEff5) {
    Error("AddTaskPerformanceTPC", "Cannot create AliPerformanceEff");
  }

  pCompEff5->SetAliRecInfoCuts(pRecInfoCutsTPC);
  pCompEff5->SetAliMCInfoCuts(pMCInfoCuts);
  pCompEff5->SetUseTrackVertex(kTRUE);



  //
  // Add components to the performance task
  //
  if(!bUseMCInfo) { 
    pCompTPC0->SetTriggerClass(triggerClass);
    pCompMatch1->SetTriggerClass(triggerClass);
    pCompMatch2->SetTriggerClass(triggerClass);
    pCompDEdx3->SetTriggerClass(triggerClass);
  }
  task->AddPerformanceObject( pCompTPC0 );
  task->AddPerformanceObject( pCompMatch1 );
  task->AddPerformanceObject( pCompMatch2 );
  task->AddPerformanceObject( pCompDEdx3 );
  if(bUseMCInfo)   {
      task->AddPerformanceObject( pCompRes4 );
      task->AddPerformanceObject( pCompEff5 );
  }

  //
  // Create containers for input
  //
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

  //
  // Create containers for output
  //
   
  AliAnalysisDataContainer *coutput = mgr->CreateContainer(Form("%sQA", taskName[idx]), TList::Class(), 
							   AliAnalysisManager::kOutputContainer, 
							   Form("%s:%s_%s", mgr->GetCommonFileName(), taskName[idx],task->GetName()));


  AliAnalysisDataContainer *coutput2  = mgr->CreateContainer(Form("%sQASummary", taskName[idx]), TTree::Class(), 
							     AliAnalysisManager::kParamContainer, 
							     Form("trending.root:Summary%sQA", taskName[idx])); 
  
  mgr->ConnectOutput(task, 1, coutput);
  mgr->ConnectOutput(task, 0, coutput2);

return task;  
}
