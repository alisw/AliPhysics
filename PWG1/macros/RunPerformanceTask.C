// Macro to run TPC performance
// By default 6 performance components are added to 
// the output: 
// 1. AliPerformanceRes (TPC track resolution at DCA)
// 2. AliPerformanceResTPCInner (TPC track resolution at inner TPC wall)
// 3. AliPerformanceEff (TPC track reconstruction efficieny)
// 4. AliPerformanceDEdxTPCInner (TPC dEdxresponse - track parameters at TPC inner wall)
// 5. AliPerformanceDCA (TPC impact parameter resolution at DCA)
// 6. AliPerformanceTPC (TPC cluster and track information)

/*
  //1. Run locally e.g.

  gROOT->LoadMacro("$ALICE_ROOT/PWG1/macros/LoadMyLibs.C");

  gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C");
  TChain* chain = CreateESDChain("list_PbPb_EMCAL.txt", 1, 0);
  chain->Lookup();

  // set magnetic field
  TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", 2, 1., 1., 10., AliMagF::k5kG));

  gROOT->LoadMacro("$ALICE_ROOT/PWG1/macros/RunPerformanceTask.C");
  RunPerformanceTask(chain, kTRUE, kFALSE, kFALSE);

  //2. Run on PROOF Lite e.g.

  TProof::Open(""); 

  gROOT->LoadMacro("$ALICE_ROOT/PWG1/macros/ProofEnableAliRoot.C");
  ProofEnableAliRoot("/u/jacek/alice/AliRoot/HEAD/");

  gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C");
  TChain* chain = CreateESDChain("list_PbPb_EMCAL.txt",20, 0);
  chain->Lookup();

  // set magnetic field
  TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", 2, 1., 1., 10., AliMagF::k5kG));

  gROOT->LoadMacro("$ALICE_ROOT/PWG1/macros/RunPerformanceTask.C");
  RunPerformanceTask(chain, kTRUE, kFALSE, kTRUE);


  //3. Run only on static PROOF at GSI e.g.

  TProof::Reset("jacek@lxgrid5.gsi.de");
  TProofMgr * proofmgr = TProof::Mgr("jacek@lxgrid5.gsi.de");
  proofmgr->SetROOTVersion("523-04");
  TProof * proof = proofmgr->CreateSession();
  proof->SetParameter("PROOF_MaxSlavesPerNode", (Long_t)10000);

  gROOT->LoadMacro("$ALICE_ROOT/PWG1/macros/ProofEnableAliRoot.C");
  ProofEnableAliRoot("/u/jacek/alice/AliRoot/HEAD/");

  gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C");
  TChain* chain = CreateESDChain("esd_v4-16-Rev-08-grid.txt", 200, 0);
  chain->Lookup();

  // set magnetic field
  TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", 2, 1., 1., 10., AliMagF::k5kG));

  gROOT->LoadMacro("$ALICE_ROOT/PWG1/macros/RunPerformanceTask.C");
  RunPerformanceTask(chain, kTRUE, kFALSE, kTRUE); 

  //4. Make final spectra and store them in the
  // output folder and generate control pictures e.g.

  TFile f("TPC.Performance.root");
  AliPerformanceEff * compObjEff = (AliPerformanceEff*)coutput->FindObject("AliPerformanceEff");
  compObjEff->Analyse();
  compObjEff->GetAnalysisFolder()->ls("*");
  // create pictures
  compObjEff->PrintHisto(kTRUE,"PerformanceEffQA.ps");
  // store output QA histograms in file
  TFile fout("PerformanceEffQAHisto.root","recreate");
  compObjEff->GetAnalysisFolder()->Write();
  fout.Close();
  f.Close();

*/

//_____________________________________________________________________________
void RunPerformanceTask(TChain *chain, Bool_t bUseMCInfo=kTRUE, Bool_t bUseFriend=kTRUE,  Bool_t bProof=kTRUE)
{
  if(!chain) 
  {
    AliDebug(AliLog::kError, "ERROR: No input chain available");
    return;
  }

  // set Proof
  if(bProof) { 

    /*
    // ONLY FOR GSI
    TProof::Reset("jacek@lxgrid5.gsi.de");
    TProofMgr * proofmgr = TProof::Mgr("jacek@lxgrid5.gsi.de");
    proofmgr->SetROOTVersion("523-04");
    TProof * proof = proofmgr->CreateSession();
    proof->SetParameter("PROOF_MaxSlavesPerNode", (Long_t)10000);
    */

    //cout << "*** START PROOF Lite SESSION ***" << endl;
    TProof::Open(""); 
    gROOT->LoadMacro("ProofEnableAliRoot.C");
    ProofEnableAliRoot("/u/jacek/alice/AliRoot/HEAD/");

  }

  //
  // Create global cuts objects 
  //

  // Create ESD track reconstruction cuts
  AliRecInfoCuts *pRecInfoCuts = new AliRecInfoCuts(); 
  if(pRecInfoCuts) {
    pRecInfoCuts->SetMaxDCAToVertexXY(3.0);
    pRecInfoCuts->SetMaxDCAToVertexZ(3.0);
    pRecInfoCuts->SetMinNClustersTPC(50);
    pRecInfoCuts->SetMinNClustersITS(2);
    pRecInfoCuts->SetHistogramsOn(kFALSE); 
  } else {
    AliDebug(AliLog::kError, "ERROR: Cannot create AliRecInfoCuts object");
  }

  // Create MC track reconstruction cuts
  AliMCInfoCuts  *pMCInfoCuts = new AliMCInfoCuts();
  if(pMCInfoCuts) {
    pMCInfoCuts->SetMinTrackLength(50);
  } else {
    AliDebug(AliLog::kError, "ERROR: Cannot AliMCInfoCuts object");
  }

  //
  // Create performance objects and set cuts 
  //
  enum { kTPC = 0, kTPCITS, kConstrained, kTPCInner, kTPCOuter, kTPCSec };

  //
  // Resolution
  //
  AliPerformanceRes *pCompRes0 = new AliPerformanceRes("AliPerformanceRes","AliPerformanceRes",kTPC,kFALSE); 
  if(!pCompRes0) {
    AliDebug(AliLog::kError, "ERROR: Cannot create AliPerformanceRes object");
  }
  pCompRes0->SetAliRecInfoCuts(pRecInfoCuts);
  pCompRes0->SetAliMCInfoCuts(pMCInfoCuts);

  AliPerformanceRes *pCompRes3 = new AliPerformanceRes("AliPerformanceResTPCInner","AliPerformanceResTPCInner",kTPCInner,kFALSE); 
  if(!pCompRes3) {
    AliDebug(AliLog::kError, "ERROR: Cannot create AliPerformanceResInnerTPC object");
  }
  pCompRes3->SetAliRecInfoCuts(pRecInfoCuts);
  pCompRes3->SetAliMCInfoCuts(pMCInfoCuts);

  AliPerformanceRes *pCompRes4 = new AliPerformanceRes("AliPerformanceResTPCOuter","AliPerformanceResTPCOuter",kTPCOuter,kFALSE); 
  if(!pCompRes4) {
    AliDebug(AliLog::kError, "ERROR: Cannot create AliPerformanceResOuterTPC object");
  }
  pCompRes4->SetAliRecInfoCuts(pRecInfoCuts);
  pCompRes4->SetAliMCInfoCuts(pMCInfoCuts);

  //
  // Efficiency
  //
  AliPerformanceEff *pCompEff0 = new AliPerformanceEff("AliPerformanceEff","AliPerformanceEff",kTPC,kFALSE); 
  if(!pCompEff0) {
    AliDebug(AliLog::kError, "ERROR: Cannot create AliPerformanceEff object");
  }
  pCompEff0->SetAliRecInfoCuts(pRecInfoCuts);
  pCompEff0->SetAliMCInfoCuts(pMCInfoCuts);
  //
  // dEdx
  //
  AliPerformanceDEdx *pCompDEdx3 = new AliPerformanceDEdx("AliPerformanceDEdxTPCInner","AliPerformanceDEdxTPCInner",kTPCInner,kFALSE); 
  if(!pCompDEdx3) {
    AliDebug(AliLog::kError, "ERROR: Cannot create AliPerformanceDEdxTPCInner object");
  }
  pCompDEdx3->SetAliRecInfoCuts(pRecInfoCuts);
  pCompDEdx3->SetAliMCInfoCuts(pMCInfoCuts);
  //
  // DCA
  //
  AliPerformanceDCA *pCompDCA0 = new AliPerformanceDCA("AliPerformanceDCA","AliPerformanceDCA",kTPC,kFALSE); 
  if(!pCompDCA0) {
    AliDebug(AliLog::kError, "ERROR: Cannot create AliPerformanceDCA object");
  }
  pCompDCA0->SetAliRecInfoCuts(pRecInfoCuts);
  pCompDCA0->SetAliMCInfoCuts(pMCInfoCuts);
  //
  // TPC performance
  //
  AliPerformanceTPC *pCompTPC0 = new AliPerformanceTPC("AliPerformanceTPC","AliPerformanceTPC",kTPC,kFALSE); 
  if(!pCompTPC0) {
    AliDebug(AliLog::kError, "ERROR: Cannot create AliPerformanceTPC object");
  }
  pCompTPC0->SetAliRecInfoCuts(pRecInfoCuts);
  pCompTPC0->SetAliMCInfoCuts(pMCInfoCuts);
 
  // Swtich off all AliInfo (too much output!!!)
  AliLog::SetGlobalLogLevel(AliLog::kError);

  // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager;

  // Create task
  AliPerformanceTask *task = new AliPerformanceTask("Performance","TPC Performance");
  if (bUseMCInfo) task->SetUseMCInfo(kTRUE);
  if (bUseFriend) task->SetUseESDfriend(kTRUE);
  task->AddPerformanceObject( pCompRes0 );
  task->AddPerformanceObject( pCompRes3 );
  if(bUseFriend) task->AddPerformanceObject( pCompRes4 );
  task->AddPerformanceObject( pCompEff0 );
  task->AddPerformanceObject( pCompDEdx3 );
  task->AddPerformanceObject( pCompDCA0 );
  task->AddPerformanceObject( pCompTPC0 );

  // Add task
  mgr->AddTask(task);

  // Add ESD handler
  AliESDInputHandler* esdH = new AliESDInputHandler;
  mgr->SetInputEventHandler(esdH);

  if(bUseMCInfo) {
  // Enable MC event handler
    AliMCEventHandler* handler = new AliMCEventHandler;
    handler->SetReadTR(kTRUE);
    mgr->SetMCtruthEventHandler(handler);
  }
 
  // Create containers for input
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  mgr->ConnectInput(task, 0, cinput);

  // Create containers for output
  AliAnalysisDataContainer *coutput = mgr->CreateContainer("coutput", TList::Class(), AliAnalysisManager::kOutputContainer, Form("TPC.%s.root", task->GetName()));
  mgr->ConnectOutput(task, 0, coutput);

  // Enable debug printouts
  mgr->SetDebugLevel(0);

  if (!mgr->InitAnalysis())
    return;

  mgr->PrintStatus();

  if(bProof) mgr->StartAnalysis("proof",chain);
  else mgr->StartAnalysis("local",chain);
}

