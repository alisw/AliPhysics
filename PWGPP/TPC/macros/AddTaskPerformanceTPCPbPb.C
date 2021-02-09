///////////////////////////////////////////////////////////////////////////////
// Macro to setup AliPerformanceTask for 
// TPC performance QA to run on PWGPP QA train. 
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
// gSystem->Load("libTPCcalib");
// gSystem->Load("libTender");
// gSystem->Load("libPWGPP");
//
// gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/TPC/macros/AddTaskPerformanceTPCPbPb.C");
// AliPerformanceTask *tpcQA = AddTaskPerformanceTPCPbPb("kFALSE","kTRUE","triggerClass"); 
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
AliPerformanceTask* AddTaskPerformanceTPCPbPb(Bool_t bUseMCInfo=kFALSE, Bool_t bUseESDfriend=kTRUE, const char *triggerClass=0)
{
  //
  // Add AliPerformanceTask with TPC performance components
  //
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr) { 
    Error("AddTaskPerformanceTPCPbPb","AliAnalysisManager not set!");
    return NULL;
  }
  
  TString type = mgr->GetInputEventHandler()->GetDataType();
  if (!type.Contains("ESD")) {
    Error("AddTaskPerformanceTPCPbPb", "ESD input handler needed!");
    return NULL;
  }
  
  AliMCEventHandler *mcH = (AliMCEventHandler*)mgr->GetMCtruthEventHandler();
  if (!mcH && bUseMCInfo) {
    Error("AddTaskPerformanceTPCPbPb", "MC input handler needed!");
    return NULL;
  }

  //
  // Create task
  //
  AliPerformanceTask *task = new AliPerformanceTask("PerformanceQA","TPC Performance");
  if (!task) {
    Error("AddTaskPerformanceTPCPbPb", "TPC performance task cannot be created!");
    return NULL;
  }
  task->SetUseMCInfo(bUseMCInfo);
  task->SetUseESDfriend(bUseESDfriend);

  task->SetUseTerminate(kFALSE);

  task->SetUseCentrality(1); // 0=off,1=Vzero,2=SPD
  Int_t cent0 = 0;
  Int_t cent1 = 20;
  Int_t cent2 = 40;

  //
  // Add task to analysis manager
  //
  mgr->AddTask(task);

  //
  // Create TPC-ESD track reconstruction cuts
  // MB tracks
  //
  AliRecInfoCuts *pRecInfoCutsTPC = new AliRecInfoCuts("pRecInfoCutsTPC"); 
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
  AliMCInfoCuts  *pMCInfoCuts = new AliMCInfoCuts("pMCInfoCuts");
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
  // TPC performance  ------------------------------------------------------------------------------------
  //
  
  AliPerformanceTPC *pCompTPC0 = NULL;
  AliPerformanceTPC *pCompTPCc0 = NULL;
  AliPerformanceTPC *pCompTPCc1 = NULL;
  AliPerformanceTPC *pCompTPCc2 = NULL;

  if( !task->GetUseCentrality() ) {
    
    pCompTPC0 = new AliPerformanceTPC("AliPerformanceTPC","AliPerformanceTPC",kTPC,kFALSE); 
    if(!pCompTPC0) {
      Error("AddTaskPerformanceTPC", "Cannot create AliPerformanceTPC");
    }
    pCompTPC0->SetAliRecInfoCuts(pRecInfoCutsTPC);
    pCompTPC0->SetAliMCInfoCuts(pMCInfoCuts);
    pCompTPC0->SetUseTrackVertex(kFALSE);
    pCompTPC0->SetUseHLT(kFALSE);
  }
  else {
    pCompTPCc0 = new AliPerformanceTPC("AliPerformanceTPC",Form("AliPerformanceTPC_cent_%d",cent0),kTPC,kFALSE); 
    if(!pCompTPCc0) {
      Error(Form("AddTaskPerformanceTPC_cent_%d",cent0), "Cannot create AliPerformanceTPC");
    }
    pCompTPCc0->SetAliRecInfoCuts(pRecInfoCutsTPC);
    pCompTPCc0->SetAliMCInfoCuts(pMCInfoCuts);
    pCompTPCc0->SetUseTrackVertex(kFALSE);
    pCompTPCc0->SetUseHLT(kFALSE);
    pCompTPCc0->SetUseCentralityBin(cent0);

    pCompTPCc1 = new AliPerformanceTPC("AliPerformanceTPC",Form("AliPerformanceTPC_cent_%d",cent1),kTPC,kFALSE); 
    if(!pCompTPCc1) {
      Error(Form("AddTaskPerformanceTPC_cent_%d",cent1), "Cannot create AliPerformanceTPC");
    }
    pCompTPCc1->SetAliRecInfoCuts(pRecInfoCutsTPC);
    pCompTPCc1->SetAliMCInfoCuts(pMCInfoCuts);
    pCompTPCc1->SetUseTrackVertex(kFALSE);
    pCompTPCc1->SetUseHLT(kFALSE);
    pCompTPCc1->SetUseCentralityBin(cent1);

    pCompTPCc2 = new AliPerformanceTPC("AliPerformanceTPC",Form("AliPerformanceTPC_cent_%d",cent2),kTPC,kFALSE); 
    if(!pCompTPCc2) {
      Error(Form("AddTaskPerformanceTPC_cent_%d",cent2), "Cannot create AliPerformanceTPC");
    }
    pCompTPCc2->SetAliRecInfoCuts(pRecInfoCutsTPC);
    pCompTPCc2->SetAliMCInfoCuts(pMCInfoCuts);
    pCompTPCc2->SetUseTrackVertex(kFALSE);
    pCompTPCc2->SetUseHLT(kFALSE);
    pCompTPCc2->SetUseCentralityBin(cent2);
  }
  
  //
  // dEdx ------------------------------------------------------------------------------------
  //

  AliPerformanceDEdx *pCompDEdx0 = NULL;
  AliPerformanceDEdx *pCompDEdxc0 = NULL;
  AliPerformanceDEdx *pCompDEdxc1 = NULL;
  AliPerformanceDEdx *pCompDEdxc2 = NULL;

  if( !task->GetUseCentrality() ) {
    pCompDEdx0 = new AliPerformanceDEdx("AliPerformanceDEdxTPCInner","AliPerformanceDEdxTPCInner",kTPCInner,kFALSE); 
    if(!pCompDEdx0) {
      Error("AddTaskPerformanceTPC", "Cannot create AliPerformanceDEdxTPCInner");
    }
    pCompDEdx0->SetAliRecInfoCuts(pRecInfoCutsTPC);
    pCompDEdx0->SetAliMCInfoCuts(pMCInfoCuts);
    pCompDEdx0->SetUseTrackVertex(kFALSE);
  }
  else {
    pCompDEdxc0 = new AliPerformanceDEdx("AliPerformanceDEdxTPCInner",
					 Form("AliPerformanceDEdxTPCInner_cent_%d",cent0),kTPCInner,kFALSE); 
    if(!pCompDEdxc0) {
      Error(Form("AddTaskPerformanceTPC_cent_%d",cent0), "Cannot create AliPerformanceDEdxTPCInner");
    }
    pCompDEdxc0->SetAliRecInfoCuts(pRecInfoCutsTPC);
    pCompDEdxc0->SetAliMCInfoCuts(pMCInfoCuts);
    pCompDEdxc0->SetUseTrackVertex(kFALSE);   
    pCompDEdxc0->SetUseCentralityBin(cent0);

    pCompDEdxc1 = new AliPerformanceDEdx("AliPerformanceDEdxTPCInner",
					 Form("AliPerformanceDEdxTPCInner_cent_%d",cent1),kTPCInner,kFALSE); 
    if(!pCompDEdxc1) {
      Error(Form("AddTaskPerformanceTPC_cent_%d",cent1), "Cannot create AliPerformanceDEdxTPCInner");
    }
    pCompDEdxc1->SetAliRecInfoCuts(pRecInfoCutsTPC);
    pCompDEdxc1->SetAliMCInfoCuts(pMCInfoCuts);
    pCompDEdxc1->SetUseTrackVertex(kFALSE);   
    pCompDEdxc1->SetUseCentralityBin(cent1);

    pCompDEdxc2 = new AliPerformanceDEdx("AliPerformanceDEdxTPCInner",
					 Form("AliPerformanceDEdxTPCInner_cent_%d",cent2),kTPCInner,kFALSE); 
    if(!pCompDEdxc2) {
      Error(Form("AddTaskPerformanceTPC_cent_%d",cent2), "Cannot create AliPerformanceDEdxTPCInner");
    }
    pCompDEdxc2->SetAliRecInfoCuts(pRecInfoCutsTPC);
    pCompDEdxc2->SetAliMCInfoCuts(pMCInfoCuts);
    pCompDEdxc2->SetUseTrackVertex(kFALSE);
    pCompDEdxc2->SetUseCentralityBin(cent2);  
  }

  //
  // Efficiency ------------------------------------------------------------------------------------
  //

  AliPerformanceEff *pCompEff0  = NULL;
  AliPerformanceEff *pCompEffc0 = NULL;
  AliPerformanceEff *pCompEffc1 = NULL;
  AliPerformanceEff *pCompEffc2 = NULL;

  if( !task->GetUseCentrality() ) {
    pCompEff0 = new AliPerformanceEff("AliPerformanceEff", "AliPerformanceEff",kTPC,kFALSE); 
    if(!pCompEff0) {
      Error("AddTaskPerformanceTPC", "Cannot create AliPerformanceEff");
    }
    pCompEff0->SetAliRecInfoCuts(pRecInfoCutsTPC);
    pCompEff0->SetAliMCInfoCuts(pMCInfoCuts);
    pCompEff0->SetUseTrackVertex(kFALSE);
    pCompEff0->InitHighMult();
  }
  else {
    pCompEffc0 = new AliPerformanceEff("AliPerformanceEff", Form("AliPerformanceEff_cent_%d",cent0),kTPC,kFALSE); 
    if(!pCompEffc0) {
      Error(Form("AddTaskPerformanceTPC_cent_%d",cent0), "Cannot create AliPerformanceEff");
    }
    pCompEffc0->SetAliRecInfoCuts(pRecInfoCutsTPC);
    pCompEffc0->SetAliMCInfoCuts(pMCInfoCuts);
    pCompEffc0->SetUseTrackVertex(kFALSE);
    pCompEffc0->InitHighMult();
    pCompEffc0->SetUseCentralityBin(cent0);  

    pCompEffc1 = new AliPerformanceEff("AliPerformanceEff", Form("AliPerformanceEff_cent_%d",cent1),kTPC,kFALSE); 
    if(!pCompEffc1) {
      Error(Form("AddTaskPerformanceTPC_cent_%d",cent1), "Cannot create AliPerformanceEff");
    }
    pCompEffc1->SetAliRecInfoCuts(pRecInfoCutsTPC);
    pCompEffc1->SetAliMCInfoCuts(pMCInfoCuts);
    pCompEffc1->SetUseTrackVertex(kFALSE);
    pCompEffc1->InitHighMult();
    pCompEffc1->SetUseCentralityBin(cent1);  

    pCompEffc2 = new AliPerformanceEff("AliPerformanceEff", Form("AliPerformanceEff_cent_%d",cent2),kTPC,kFALSE); 
    if(!pCompEffc2) {
      Error(Form("AddTaskPerformanceTPC_cent_%d",cent2), "Cannot create AliPerformanceEff");
    }
    pCompEffc2->SetAliRecInfoCuts(pRecInfoCutsTPC);
    pCompEffc2->SetAliMCInfoCuts(pMCInfoCuts);
    pCompEffc2->SetUseTrackVertex(kFALSE);
    pCompEffc2->InitHighMult();
    pCompEffc2->SetUseCentralityBin(cent2);  
  }

  //
  // Resolution ------------------------------------------------------------------------------------
  //

  AliPerformanceRes *pCompRes0 = NULL;
  AliPerformanceRes *pCompResc0 = NULL;
  AliPerformanceRes *pCompResc1 = NULL;
  AliPerformanceRes *pCompResc2 = NULL;

  if( !task->GetUseCentrality() ) {
    pCompRes0 = new AliPerformanceRes("AliPerformanceRes", "AliPerformanceRes",kTPC,kFALSE); 
    if(!pCompRes0) {
      Error("AddTaskPerformanceTPC", "Cannot create AliPerformanceRes");
    }
    pCompRes0->SetAliRecInfoCuts(pRecInfoCutsTPC);
    pCompRes0->SetAliMCInfoCuts(pMCInfoCuts);
    pCompRes0->SetUseTrackVertex(kFALSE);
    pCompRes0->InitHighMult();
  }
  else {
    pCompResc0 = new AliPerformanceRes("AliPerformanceRes", Form("AliPerformanceRes_cent_%d",cent0), kTPC, kFALSE); 
    if(!pCompResc0) {
      Error(Form("AddTaskPerformanceTPC_cent_%d",cent0), "Cannot create AliPerformanceRes");
    }
    pCompResc0->SetAliRecInfoCuts(pRecInfoCutsTPC);
    pCompResc0->SetAliMCInfoCuts(pMCInfoCuts);
    pCompResc0->SetUseTrackVertex(kFALSE);
    pCompResc0->InitHighMult();
    pCompResc0->SetUseCentralityBin(cent0); 

    pCompResc1 = new AliPerformanceRes("AliPerformanceRes", Form("AliPerformanceRes_cent_%d",cent1), kTPC, kFALSE); 
    if(!pCompResc1) {
      Error(Form("AddTaskPerformanceTPC_cent_%d",cent1), "Cannot create AliPerformanceRes");
    }
    pCompResc1->SetAliRecInfoCuts(pRecInfoCutsTPC);
    pCompResc1->SetAliMCInfoCuts(pMCInfoCuts);
    pCompResc1->SetUseTrackVertex(kFALSE);
    pCompResc1->InitHighMult();
    pCompResc1->SetUseCentralityBin(cent1); 


    pCompResc2 = new AliPerformanceRes("AliPerformanceRes", Form("AliPerformanceRes_cent_%d",cent2), kTPC, kFALSE); 
    if(!pCompResc2) {
      Error(Form("AddTaskPerformanceTPC_cent_%d",cent2), "Cannot create AliPerformanceRes");
    }
    pCompResc2->SetAliRecInfoCuts(pRecInfoCutsTPC);
    pCompResc2->SetAliMCInfoCuts(pMCInfoCuts);
    pCompResc2->SetUseTrackVertex(kFALSE);
    pCompResc2->InitHighMult();
    pCompResc2->SetUseCentralityBin(cent2); 
  }

  // ----------------------------------------------------------------------------------------------------------------------------
  
  //
  // Add components to the performance task  ------------------------------------------------------------------------------------
  //
  if( !task->GetUseCentrality() ) {
    task->AddPerformanceObject( pCompTPC0 );
    task->AddPerformanceObject( pCompDEdx0 );
    task->AddPerformanceObject( pCompEff0 );
    task->AddPerformanceObject( pCompRes0 );

    if(!bUseMCInfo) { 
      pCompTPC0->SetTriggerClass(triggerClass);
      pCompDEdx0->SetTriggerClass(triggerClass);
      pCompEff0->SetTriggerClass(triggerClass);
      pCompRes0->SetTriggerClass(triggerClass);
    }
  }
  else {
    task->AddPerformanceObject( pCompTPCc0 );
    task->AddPerformanceObject( pCompTPCc1 );
    task->AddPerformanceObject( pCompTPCc2 );
    task->AddPerformanceObject( pCompDEdxc0 );
    task->AddPerformanceObject( pCompDEdxc1 );
    task->AddPerformanceObject( pCompDEdxc2 );
    task->AddPerformanceObject( pCompEffc0 );
    task->AddPerformanceObject( pCompEffc1 );
    task->AddPerformanceObject( pCompEffc2 );
    task->AddPerformanceObject( pCompResc0 );
    task->AddPerformanceObject( pCompResc1 );
    task->AddPerformanceObject( pCompResc2 );

    if(!bUseMCInfo) { 
      pCompTPCc0->SetTriggerClass(triggerClass);
      pCompDEdxc0->SetTriggerClass(triggerClass);
      pCompTPCc1->SetTriggerClass(triggerClass);
      pCompDEdxc1->SetTriggerClass(triggerClass);
      pCompTPCc2->SetTriggerClass(triggerClass);
      pCompDEdxc2->SetTriggerClass(triggerClass);
      pCompEffc0->SetTriggerClass(triggerClass);
      pCompEffc1->SetTriggerClass(triggerClass);
      pCompEffc2->SetTriggerClass(triggerClass);
      pCompResc0->SetTriggerClass(triggerClass);
      pCompResc1->SetTriggerClass(triggerClass);
      pCompResc2->SetTriggerClass(triggerClass);
    }
  }

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
  mgr->ConnectOutput(task, 2, coutput2_tpc);

return task;  
}
