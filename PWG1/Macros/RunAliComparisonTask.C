void RunAliComparisonTask(TChain  *chain = 0, Bool_t aProof = kTRUE, Bool_t aDebug = kFALSE)
{
  //
  // Set mag field map (needed to propagate track to the DCA)
  //
  //Int_t magField = 2;  // 0 - 0.2 T, 1 = 0.4 T, 2  - 0.5 T
  //magFMap = new AliMagFMaps("Maps","Maps", 2, 1., 10., magField);
  //AliTracker::SetFieldMap(magFMap,kFALSE);

  //AliMagWrapCheb* field = 0x0;
  //field = new AliMagWrapCheb("Maps","Maps", 2, 1., 10., AliMagWrapCheb::k5kG,kTRUE,"$(ALICE_ROOT)/data/maps/mfchebKGI_sym.root");
  //Bool_t uniform=kFALSE;
  //AliTracker::SetFieldMap(field,uniform);  // tracking with the real map

  TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", 2, 1., 1., 10., AliMagF::k5kG));

  //
  // Create global cuts objects 
  //

  // Create ESD track reconstruction cuts
  AliRecInfoCuts *pRecInfoCuts = new AliRecInfoCuts(); 
  if(pRecInfoCuts) {
    pRecInfoCuts->SetPtRange(0.20,200.0);
    //pRecInfoCuts->SetEtaRange(-0.9,0.9);
    pRecInfoCuts->SetMaxDCAToVertexXY(3.0);
    pRecInfoCuts->SetMaxDCAToVertexZ(3.0);
    pRecInfoCuts->SetMinNClustersTPC(50);
    pRecInfoCuts->SetMinNClustersITS(2);
    pRecInfoCuts->SetMinTPCsignalN(50);

	pRecInfoCuts->SetHistogramsOn(kFALSE); 
  } else {
    AliDebug(AliLog::kError, "ERROR: Cannot create AliRecInfoCuts object");
  }

  // Create MC track reconstruction cuts
  AliMCInfoCuts  *pMCInfoCuts = new AliMCInfoCuts();
  if(pMCInfoCuts) {
    pMCInfoCuts->SetMinRowsWithDigits(50);
    pMCInfoCuts->SetMaxR(0.025); // from diamond xy size (pp@10TeV) 
    pMCInfoCuts->SetMaxVz(10.);  // from diamond z size  (pp@10TeV)
    pMCInfoCuts->SetRangeTPCSignal(0.5,1.4); 
  } else {
    AliDebug(AliLog::kError, "ERROR: Cannot AliMCInfoCuts object");
  }

  //
  // Create comparison objects and set cuts 
  //
  const Int_t kTPC = 0;
  const Int_t kTPCITS = 1;
  const Int_t kConstrained = 2;

  // Resolution
  AliComparisonRes *pCompRes0 = new AliComparisonRes("AliComparisonResTPC","AliComparisonResTPC",kTPC,kFALSE); 
  if(!pCompRes0) {
    AliDebug(AliLog::kError, "ERROR: Cannot create AliComparisonRes0 object");
  }
  pCompRes0->SetAliRecInfoCuts(pRecInfoCuts);
  pCompRes0->SetAliMCInfoCuts(pMCInfoCuts);

  AliComparisonRes *pCompRes1 = new AliComparisonRes("AliComparisonResTPCITS","AliComparisonResTPCITS",kTPCITS,kFALSE); 
  if(!pCompRes1) {
    AliDebug(AliLog::kError, "ERROR: Cannot create AliComparisonRes1 object");
  }
  pCompRes1->SetAliRecInfoCuts(pRecInfoCuts);
  pCompRes1->SetAliMCInfoCuts(pMCInfoCuts);

  AliComparisonRes *pCompRes2 = new AliComparisonRes("AliComparisonResConstrained","AliComparisonResConstrained",kConstrained,kFALSE); 
  if(!pCompRes2) {
    AliDebug(AliLog::kError, "ERROR: Cannot create AliComparisonRes2 object");
  }
  pCompRes2->SetAliRecInfoCuts(pRecInfoCuts);
  pCompRes2->SetAliMCInfoCuts(pMCInfoCuts);

  // Efficiency
  AliComparisonEff *pCompEff0 =  new AliComparisonEff("AliComparisonEffTPC","AliComparisonEffTPC",kTPC,kFALSE);
  if(!pCompEff0) {
    AliDebug(AliLog::kError, "ERROR: Cannot create AliComparisonEff object");
  }
  pCompEff0->SetAliRecInfoCuts(pRecInfoCuts);
  pCompEff0->SetAliMCInfoCuts(pMCInfoCuts);

  AliComparisonEff *pCompEff1 =  new AliComparisonEff("AliComparisonEffTPCITS","AliComparisonEffTPCITS",kTPCITS,kFALSE);
  if(!pCompEff1) {
    AliDebug(AliLog::kError, "ERROR: Cannot create AliComparisonEff object");
  }
  pCompEff1->SetAliRecInfoCuts(pRecInfoCuts);
  pCompEff1->SetAliMCInfoCuts(pMCInfoCuts);

  AliComparisonEff *pCompEff2 =  new AliComparisonEff("AliComparisonEffConstrained","AliComparisonEffConstrained",kConstrained,kFALSE);
  if(!pCompEff2) {
    AliDebug(AliLog::kError, "ERROR: Cannot create AliComparisonEff object");
  }
  pCompEff2->SetAliRecInfoCuts(pRecInfoCuts);
  pCompEff2->SetAliMCInfoCuts(pMCInfoCuts);

  // dE/dx
  AliComparisonDEdx *pCompDEdx0 = new AliComparisonDEdx("AliComparisonDEdxTPC","AliComparisonDEdxTPC",kTPC,kFALSE);
  if(!pCompDEdx0) {
    AliDebug(AliLog::kError, "ERROR: Cannot create AliComparisonDEdx object");
  }
  pCompDEdx0->SetAliRecInfoCuts(pRecInfoCuts);
  pCompDEdx0->SetAliMCInfoCuts(pMCInfoCuts);

  // DCA
  AliComparisonDCA *pCompDCA0 = new AliComparisonDCA("AliComparisonDCATPC","AliComparisonDCATPC",kTPC,kFALSE);
  if(!pCompDCA0) {
    AliDebug(AliLog::kError, "ERROR: Cannot create AliComparisonDCATPC object");
  }
  pCompDCA0->SetAliRecInfoCuts(pRecInfoCuts);
  pCompDCA0->SetAliMCInfoCuts(pMCInfoCuts);

  AliComparisonDCA *pCompDCA1 = new AliComparisonDCA("AliComparisonDCATPCITS","AliComparisonDCAITS",kTPCITS,kFALSE);
  if(!pCompDCA1) {
    AliDebug(AliLog::kError, "ERROR: Cannot create AliComparisonDCAITS");
  }
  pCompDCA1->SetAliRecInfoCuts(pRecInfoCuts);
  pCompDCA1->SetAliMCInfoCuts(pMCInfoCuts);

  // Create the analysis manager
  mgr = new AliAnalysisManager("testAnalysis");

  // Create, add task
  task = new AliComparisonTask;

  task->AddComparisonObject( pCompRes0 );
  task->AddComparisonObject( pCompRes1 );
  task->AddComparisonObject( pCompRes2 );
  task->AddComparisonObject( pCompEff0 );
  task->AddComparisonObject( pCompEff1 );
  task->AddComparisonObject( pCompEff2 );
  task->AddComparisonObject( pCompDEdx0 );
  task->AddComparisonObject( pCompDCA0 );
  task->AddComparisonObject( pCompDCA1 );
  
  mgr->AddTask(task);

  // Attach input
  cInput  = mgr->CreateContainer("cInput", TChain::Class(), AliAnalysisManager::kInputContainer);
  mgr->ConnectInput(task, 0, cInput);

  // Attach output
  cOutput = mgr->CreateContainer("cOutput", TList::Class(), AliAnalysisManager::kOutputContainer,"Output.root");
  mgr->ConnectOutput(task, 0, cOutput);

  // Enable debug printouts
  if (aDebug)
    mgr->SetDebugLevel(2);

  // Run analysis
  mgr->InitAnalysis();
  mgr->PrintStatus();

  if(chain) {
    mgr->StartAnalysis((aProof) ? "proof" : "local", chain);
  } else {
    AliDebug(AliLog::kError, "ERROR: No chain available");
  }
}
