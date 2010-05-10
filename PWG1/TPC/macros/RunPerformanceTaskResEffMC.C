// Macro to run TPC performance task (locally, proof).
//
// By default 2 performance components are added to 
// AliPerformanceRes - resloution at DCA to prim. vertex
// AliPerformanceEff - efficiency for prim. tracks
//

/*
 
  //1. Run locally 
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libTENDER");
  gSystem->Load("libCORRFW");
  gSystem->Load("libPWG1");

  gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C");
  TChain* chain = CreateESDChain("esds.txt",10, 0);
  //TChain* chain = CreateESDChain("/u/jacek/alice/dNdPt/input/LHC10a8/esds_104867_MC_LHC10a8.txt",3, 0);
  chain->Lookup();

  gROOT->LoadMacro("$ALICE_ROOT/PWG1/TPC/macros/RunPerformanceTaskResEffMC.C");
  RunPerformanceTaskResEffMC(chain, kTRUE, kFALSE, kFALSE);

  //2. Make final spectra and store them in the
  // output folder and generate control pictures

  aliroot -b
  TFile f("TPC.MC.Performance.root");

  AliPerformanceRes * compObjRes = (AliPerformanceRes*)TPCMC->FindObject("AliPerformanceRes");
  compObjRes->Analyse();
  compObjRes->GetAnalysisFolder()->ls("*");
  // create control pictures
  compObjRes->PrintHisto(kTRUE,"TPC.MC.PerformanceRes.Pict.ps");
  // store output QA histograms in file
  TFile fout("TPC.MC.PerformanceRes.Histo.root","recreate");
  compObjRes->GetAnalysisFolder()->Write();
  fout.Close();

  //
  f->cd();
  AliPerformanceEff * compObjEff = (AliPerformanceEff*)TPCMC->FindObject("AliPerformanceEff");
  compObjEff->Analyse();
  compObjEff->GetAnalysisFolder()->ls("*");
  // create control pictures
  compObjEff->PrintHisto(kTRUE,"TPC.MC.PerformanceEff.Pict.ps");
  // store output QA histograms in file
  TFile fout("TPC.MC.PerformanceEff.Histo.root","recreate");
  compObjEff->GetAnalysisFolder()->Write();
  fout.Close();

  f.Close();

*/

//_____________________________________________________________________________
void RunPerformanceTaskResEffMC(TChain *chain, Bool_t bUseMCInfo=kTRUE, Bool_t bUseESDfriend=kTRUE,  Bool_t bProof=kTRUE)
{
  if(!chain) 
  {
    AliDebug(AliLog::kError, "ERROR: No input chain available");
    return;
  }
  //
  // Swtich off all AliInfo (too much output!)
  //
  AliLog::SetGlobalLogLevel(AliLog::kError);

  //
  // Create analysis manager
  //
  AliAnalysisManager *mgr = new AliAnalysisManager;
  if(!mgr) { 
    Error("runTPCQA","AliAnalysisManager not set!");
    return;
  }

  //
  // Set ESD input handler
  //
  AliESDInputHandler* esdH = new AliESDInputHandler;
  if(!esdH) { 
    Error("runTPCQA","AliESDInputHandler not created!");
    return;
  }
  if(bUseESDfriend) esdH->SetActiveBranches("ESDfriend");
  mgr->SetInputEventHandler(esdH);

  //
  // Set MC input handler
  //
  if(bUseMCInfo) {
    AliMCEventHandler* mcH = new AliMCEventHandler;
    if(!esdH) { 
      Error("runTPCQA","AliMCEventHandler not created!");
      return;
    }
    mcH->SetReadTR(kTRUE);
    mgr->SetMCtruthEventHandler(mcH);
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
  // MB tracks
  //
  AliRecInfoCuts *pRecInfoCutsTPC = new AliRecInfoCuts(); 
  if(pRecInfoCutsTPC) {
    pRecInfoCutsTPC->SetMaxDCAToVertexXY(3.0);
    pRecInfoCutsTPC->SetMaxDCAToVertexZ(30.0);
    //pRecInfoCutsTPC->SetMaxDCAToVertexZ(3.0);
    pRecInfoCutsTPC->SetRequireSigmaToVertex(kFALSE);
    pRecInfoCutsTPC->SetRequireTPCRefit(kFALSE);
    pRecInfoCutsTPC->SetAcceptKinkDaughters(kTRUE);
    pRecInfoCutsTPC->SetMinNClustersTPC(50);
    pRecInfoCutsTPC->SetMaxChi2PerClusterTPC(1000000.);
    pRecInfoCutsTPC->SetDCAToVertex2D(kFALSE);

    pRecInfoCutsTPC->SetHistogramsOn(kFALSE); 
  } 
  else {
    Error("AddTaskPerformanceTPC", "AliRecInfoCutsTPC cannot be created!");
    return NULL;
  }

  //
  // Create TPC-ESD track reconstruction cuts
  // MATCH tracks
  //
  AliRecInfoCuts *pRecInfoCutsMATCH = new AliRecInfoCuts(); 
  if(pRecInfoCutsMATCH) {
    pRecInfoCutsMATCH->SetMaxDCAToVertexXY(3.0);
    pRecInfoCutsMATCH->SetMaxDCAToVertexZ(3.0);
    pRecInfoCutsMATCH->SetRequireSigmaToVertex(kFALSE);
    pRecInfoCutsMATCH->SetRequireTPCRefit(kFALSE);
    pRecInfoCutsMATCH->SetAcceptKinkDaughters(kTRUE);
    pRecInfoCutsMATCH->SetMinNClustersTPC(50);
    pRecInfoCutsMATCH->SetMaxChi2PerClusterTPC(1000000.);
    pRecInfoCutsMATCH->SetDCAToVertex2D(kFALSE);
    pRecInfoCutsMATCH->SetTPCITSMatchingRadius(70); 
    pRecInfoCutsMATCH->SetTPCTRDMatchingRadius(260); 
    pRecInfoCutsMATCH->SetMinNClustersITS(3);

    pRecInfoCutsMATCH->SetHistogramsOn(kFALSE); 
  } 
  else {
    Error("AddTaskPerformanceTPC", "AliRecInfoCutsTPC cannot be created!");
    return NULL;
  }

  //
  // Create TPC-ESD track reconstruction cuts
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
  pCompTPC0->SetAliRecInfoCuts(pRecInfoCutsTPC);
  pCompTPC0->SetAliMCInfoCuts(pMCInfoCuts);
  //
  // TPC+ITS matching performance
  //
  AliPerformanceMatch *pCompMatch0 = new AliPerformanceMatch("AliPerformanceMatchTPCITS","AliPerformanceMatchTPCITS",0,kFALSE); 
  if(!pCompMatch0) {
    Error("AddTaskPerformanceMatch", "Cannot create AliPerformanceMatchTPCITS");
  }
  pCompMatch0->SetAliRecInfoCuts(pRecInfoCutsMATCH);
  pCompMatch0->SetAliMCInfoCuts(pMCInfoCuts);
  //
  // TPC+TRD matching performance
  //
  AliPerformanceMatch *pCompMatch1 = new AliPerformanceMatch("AliPerformanceMatchTPCTRD","AliPerformanceMatchTPCTRD",1,kFALSE); 
  if(!pCompMatch1) {
    Error("AddTaskPerformanceMatch", "Cannot create AliPerformanceMatchTPCTRD");
  }
  pCompMatch1->SetAliRecInfoCuts(pRecInfoCutsMATCH);
  pCompMatch1->SetAliMCInfoCuts(pMCInfoCuts);
 
  AliPerformanceMatch *pCompMatch2 = new AliPerformanceMatch("AliPerformanceMatchTPCEFF","AliPerformanceMatchTPCEFF",2,kFALSE); 
  if(!pCompMatch2) {
    Error("AddTaskPerformanceMatch", "Cannot create AliPerformanceMatchTPCEFF");
  }
  pCompMatch2->SetAliRecInfoCuts(pRecInfoCutsMATCH);
  pCompMatch2->SetAliMCInfoCuts(pMCInfoCuts);

  //
  // Add components to the performance task
  //
  //task->AddPerformanceObject( pCompDEdx3 );
  //task->AddPerformanceObject( pCompDCA0 );
  //task->AddPerformanceObject( pCompTPC0 );
  //task->AddPerformanceObject( pCompMatch0 );
  //task->AddPerformanceObject( pCompMatch1 );
  //task->AddPerformanceObject( pCompMatch2 );

  //
  if(bUseMCInfo) { 
     task->AddPerformanceObject( pCompRes0 );
     //task->AddPerformanceObject( pCompRes3 );
     //task->AddPerformanceObject( pCompRes4 );
     task->AddPerformanceObject( pCompEff0 );
  }

  //
  if(!bUseMCInfo) {
    pCompDEdx3->SetTriggerClass(triggerClass);
    pCompDCA0->SetTriggerClass(triggerClass);
    pCompTPC0->SetTriggerClass(triggerClass);
    pCompMatch0->SetTriggerClass(triggerClass);
    pCompMatch1->SetTriggerClass(triggerClass);
    pCompMatch2->SetTriggerClass(triggerClass);
  }

  //
  // Create containers for input
  //
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

  //
  // Create containers for output
  //
  AliAnalysisDataContainer *coutput_tpc = mgr->CreateContainer("TPCMC", TList::Class(), AliAnalysisManager::kOutputContainer, Form("TPC.MC.%s.root", task->GetName()));
  mgr->ConnectOutput(task, 1, coutput_tpc);

  // Enable debug printouts
  mgr->SetDebugLevel(0);

  if (!mgr->InitAnalysis())
    return;

  mgr->PrintStatus();

  if(bProof) mgr->StartAnalysis("proof",chain);
  else mgr->StartAnalysis("local",chain);
}

