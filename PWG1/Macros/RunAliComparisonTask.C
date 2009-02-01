void RunAliComparisonTask(TChain  *chain = 0, Bool_t aProof = kTRUE, Bool_t aDebug = kFALSE)
{
  //
  // Set mag field map (needed to propagate track to the DCA)
  //
  TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", 2, 1., 1., 10., AliMagF::k5kG));

  //
  // Create global cuts objects 
  //

  // Create ESD track reconstruction cuts
  AliRecInfoCuts *pRecInfoCuts = new AliRecInfoCuts(); 
  if(pRecInfoCuts) {
    pRecInfoCuts->SetPtRange(0.15,200.0);
    pRecInfoCuts->SetMaxAbsTanTheta(1.0);
    pRecInfoCuts->SetMinNClustersTPC(10);
    pRecInfoCuts->SetMinTPCsignalN(50);

	pRecInfoCuts->SetHistogramsOn(kFALSE); 
  } else {
    AliDebug(AliLog::kError, "ERROR: Cannot create AliRecInfoCuts object");
  }

  // Create MC track reconstruction cuts
  AliMCInfoCuts  *pMCInfoCuts = new AliMCInfoCuts();
  if(pMCInfoCuts) {
    pMCInfoCuts->SetMinRowsWithDigits(50);
    pMCInfoCuts->SetMaxR(0.001);  
    pMCInfoCuts->SetMaxVz(0.001); 
    pMCInfoCuts->SetRangeTPCSignal(0.5,1.4); 
  } else {
    AliDebug(AliLog::kError, "ERROR: Cannot AliMCInfoCuts object");
  }

  //
  // Create comparison objects and set cuts 
  //

  // Resolution
  AliComparisonRes *pCompRes = new AliComparisonRes(); 
  if(!pCompRes) {
    AliDebug(AliLog::kError, "ERROR: Cannot create AliComparisonRes object");
  }
  pCompRes->SetAliRecInfoCuts(pRecInfoCuts);
  pCompRes->SetAliMCInfoCuts(pMCInfoCuts);

  // Efficiency
  AliComparisonEff *pCompEff =  new AliComparisonEff();
  if(!pCompEff) {
    AliDebug(AliLog::kError, "ERROR: Cannot create AliComparisonEff object");
  }
  pCompEff->SetAliRecInfoCuts(pRecInfoCuts);
  pCompEff->SetAliMCInfoCuts(pMCInfoCuts);

  // dE/dx
  AliComparisonDEdx *pCompDEdx = new AliComparisonDEdx();
  if(!pCompDEdx) {
    AliDebug(AliLog::kError, "ERROR: Cannot create AliComparisonDEdx object");
  }
  pCompDEdx->SetAliRecInfoCuts(pRecInfoCuts);
  pCompDEdx->SetAliMCInfoCuts(pMCInfoCuts);
  pCompDEdx->SetMCPtMin(0.5);
  pCompDEdx->SetMCAbsTanThetaMax(0.5);
  pCompDEdx->SetMCPdgCode(pMCInfoCuts->GetPiP()); // only pi+ particles

  // DCA
  AliComparisonDCA *pCompDCA = new AliComparisonDCA();
  if(!pCompDCA) {
    AliDebug(AliLog::kError, "ERROR: Cannot create AliComparisonDCA object");
  }
  pCompDCA->SetAliRecInfoCuts(pRecInfoCuts);
  pCompDCA->SetAliMCInfoCuts(pMCInfoCuts);

  // Create the analysis manager
  mgr = new AliAnalysisManager("testAnalysis");

  // Create, add task
  task = new AliComparisonTask;

  task->AddComparisonObject( pCompRes );
  task->AddComparisonObject( pCompEff );
  task->AddComparisonObject( pCompDEdx );
  task->AddComparisonObject( pCompDCA );
  
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
