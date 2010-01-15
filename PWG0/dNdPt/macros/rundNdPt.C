void rundNdPt(const char *fileList ="inputList.txt",const char *outFile = "outputFile.root", Int_t NumberOfFiles=10,Int_t fromFile=0, Int_t nEvents=1000, Int_t firstEvent =0, Bool_t bUseMCInfo=kTRUE, Float_t zvWindow=20., Int_t cutMode=9,Float_t etaWindow=0.9, Float_t ptMin=0.15,AlidNdPtHelper::AnalysisMode analysisMode = AlidNdPtHelper::kTPC, AlidNdPtHelper::ParticleMode particleMode = AlidNdPtHelper::kAllPart,  AlidNdPtHelper::OutputObject outputObject=AlidNdPtHelper::kCorrection,const char *corrFile = "corrMatricesFile.root",Bool_t bProof=kFALSE)
{
  // set Proof
  if(bProof) { 
    
    //cout << "*** START PROOF Lite SESSION ***" << endl;
    //TProof::Open(""); // 1. Enter your username here
    
    //TProofMgr * proofmgr = TProof::Mgr("lxialpod2.gsi.de:21001");
    //TProofMgr * proofmgr = TProof::Mgr("lxialpod.gsi.de:21001");
    TProofMgr * proofmgr = TProof::Mgr("lxialpod.gsi.de:21004");
    //TProofMgr * proofmgr = TProof::Mgr("lxial39.gsi.de:21001");
    TProof * proof = proofmgr->CreateSession();
    proof->SetParameter("PROOF_MaxSlavesPerNode", (Long_t)10000);

    // -- Load AliRoot Libraries
    gROOT->LoadMacro("ProofEnableAliRootGSI.C");
    ProofEnableAliRootGSI("/u/jacek/alice/AliRoot/trunk");
  }

  // Swtich off all AliInfo (too much output!!!)
  AliLog::SetGlobalLogLevel(AliLog::kError);

  // Create analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager;

  //
  // Create physics trigger selection class
  //
  AliPhysicsSelection *physTrigSel =  new AliPhysicsSelection();

  //
  // Create event cuts
  //
  AlidNdPtEventCuts *evtCuts = new AlidNdPtEventCuts("AlidNdPtEventCuts","Event cuts");
  evtCuts->SetZvRange(-zvWindow,zvWindow);
  evtCuts->SetMeanXYZv(0.0,0.0,0.0);
  evtCuts->SetSigmaMeanXYZv(1.0,1.0,10.0);
  evtCuts->SetTriggerRequired(kTRUE);
  //evtCuts->SetTriggerRequired(kFALSE);

  // Create geom. acceptance cuts
  AlidNdPtAcceptanceCuts *accCuts = new AlidNdPtAcceptanceCuts("AlidNdPtAcceptanceCuts","Geom. acceptance cuts");
  accCuts->SetEtaRange(-etaWindow,etaWindow);
  accCuts->SetPtRange(ptMin,1.e10);
  accCuts->SetMaxDCAr(3.0);
  accCuts->SetMaxDCAz(30.0);

  // Create standard esd track cuts
  gROOT->LoadMacro("$ALICE_ROOT/PWG0/dNdPt/macros/CreatedNdPtTrackCuts.C");
  AliESDtrackCuts* esdTrackCuts = CreatedNdPtTrackCuts(cutMode);
  if (!esdTrackCuts) {
    printf("ERROR: esdTrackCuts could not be created\n");
    return;
  } else {
    esdTrackCuts->SetHistogramsOn(kTRUE);
  }

  //
  // Create task
  //
  AlidNdPtTask *task = new AlidNdPtTask("AlidNdPtTask");
  if (bUseMCInfo) task->SetUseMCInfo(kTRUE);

  // create cut analysis object
  if(outputObject==AlidNdPtHelper::kCutAnalysis) 
  {
    AlidNdPtCutAnalysis *fdNdPtCutAnalysis = new AlidNdPtCutAnalysis("dNdPtCutAnalysis","dN/dPt Cut Analysis");
    fdNdPtCutAnalysis->SetEventCuts(evtCuts);
    fdNdPtCutAnalysis->SetAcceptanceCuts(accCuts);
    fdNdPtCutAnalysis->SetTrackCuts(esdTrackCuts);
    fdNdPtCutAnalysis->SetAnalysisMode(analysisMode); 
    fdNdPtCutAnalysis->SetParticleMode(particleMode); 
    if(bUseMCInfo) 
    {
       //fdNdPtCutAnalysis->SetTrigger(AliTriggerAnalysis::kMB1); 
       physTrigSel->SetAnalyzeMC();
       fdNdPtCutAnalysis->SetPhysicsTriggerSelection(physTrigSel); 
       fdNdPtCutAnalysis->SetUseMCInfo(kTRUE);
    }
    else { // online trigger
    fdNdPtCutAnalysis->SetPhysicsTriggerSelection(physTrigSel); 
    //fdNdPtCutAnalysis->SetTriggerClass("CBEAMB-ABCE-NOPF-ALL"); 
    //fdNdPtCutAnalysis->SetTriggerClass("CINT1B-ABCE-NOPF-ALL");
    //fdNdPtCutAnalysis->SetTriggerClass("CINT1A-ABCE-NOPF-ALL"); 
    //fdNdPtCutAnalysis->SetTriggerClass("CINT1C-ABCE-NOPF-ALL"); 
    //fdNdPtCutAnalysis->SetTriggerClass("CINT1-E-NOPF-ALL"); 
    }

    task->AddAnalysisObject( fdNdPtCutAnalysis );
  }

  // create analysis object
  if(outputObject==AlidNdPtHelper::kAnalysis) 
  {
    AlidNdPtAnalysis *fdNdPtAnalysis = new AlidNdPtAnalysis("dNdPtAnalysis","dN/dPt Analysis");
    fdNdPtAnalysis->SetEventCuts(evtCuts);
    fdNdPtAnalysis->SetAcceptanceCuts(accCuts);
    fdNdPtAnalysis->SetTrackCuts(esdTrackCuts);
    fdNdPtAnalysis->SetAnalysisMode(analysisMode); 
    fdNdPtAnalysis->SetParticleMode(particleMode); 
    if(bUseMCInfo) 
    {
       //fdNdPtAnalysis->SetTrigger(AliTriggerAnalysis::kMB1); 
       physTrigSel->SetAnalyzeMC();
       fdNdPtAnalysis->SetPhysicsTriggerSelection(physTrigSel); 

       fdNdPtAnalysis->SetUseMCInfo(kTRUE);
       //fdNdPtAnalysis->SetHistogramsOn(kTRUE);
       fdNdPtAnalysis->SetHistogramsOn(kFALSE);
    }
    else { // online trigger
    fdNdPtAnalysis->SetPhysicsTriggerSelection(physTrigSel); 
    //fdNdPtAnalysis->SetTriggerClass("CBEAMB-ABCE-NOPF-ALL"); 
    //fdNdPtAnalysis->SetTriggerClass("CINT1B-ABCE-NOPF-ALL");
    //fdNdPtAnalysis->SetTriggerClass("CINT1A-ABCE-NOPF-ALL"); 
    //fdNdPtAnalysis->SetTriggerClass("CINT1C-ABCE-NOPF-ALL"); 
    //fdNdPtAnalysis->SetTriggerClass("CINT1-E-NOPF-ALL"); 
    }

    task->AddAnalysisObject( fdNdPtAnalysis );
  }

  // create correction object
  if(outputObject == AlidNdPtHelper::kCorrection) 
  {
    AlidNdPtCorrection *fdNdPtCorrection = new AlidNdPtCorrection("dNdPtCorrection","dN/dPt Correction", corrFile);

    fdNdPtCorrection->SetEventCuts(evtCuts);
    fdNdPtCorrection->SetAcceptanceCuts(accCuts);
    fdNdPtCorrection->SetTrackCuts(esdTrackCuts);
    fdNdPtCorrection->SetAnalysisMode(analysisMode); 
    fdNdPtCorrection->SetParticleMode(particleMode); 
    if(bUseMCInfo) 
    {
       //fdNdPtCorrection->SetTrigger(AliTriggerAnalysis::kMB1); 
       physTrigSel->SetAnalyzeMC();
       fdNdPtCorrection->SetPhysicsTriggerSelection(physTrigSel); 
       fdNdPtCorrection->SetUseMCInfo(kTRUE);
    }
    else { // online trigger
    fdNdPtCorrection->SetPhysicsTriggerSelection(physTrigSel); 
    //fdNdPtCorrection->SetTriggerClass("CBEAMB-ABCE-NOPF-ALL"); 
    //fdNdPtCorrection->SetTriggerClass("CINT1B-ABCE-NOPF-ALL");
    //fdNdPtCorrection->SetTriggerClass("CINT1A-ABCE-NOPF-ALL"); 
    //fdNdPtCorrection->SetTriggerClass("CINT1C-ABCE-NOPF-ALL"); 
    //fdNdPtCorrection->SetTriggerClass("CINT1-E-NOPF-ALL"); 
    }

    task->AddAnalysisObject( fdNdPtCorrection );
  }

  // Add task
  mgr->AddTask(task);

  // Add ESD handler
  AliESDInputHandler* esdH = new AliESDInputHandler;
  //esdH->SetInactiveBranches("*");
  mgr->SetInputEventHandler(esdH);

  if(bUseMCInfo) {
  // Enable MC event handler
    AliMCEventHandler* handler = new AliMCEventHandler;
    handler->SetReadTR(kFALSE);
    //handler->SetReadTR(kTRUE);
    mgr->SetMCtruthEventHandler(handler);
  }

  // Create input chain
  gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C");
  TChain* chain = CreateESDChain(fileList, NumberOfFiles, fromFile);
  if(!chain) {
    printf("ERROR: chain cannot be created\n");
    return;
  }

  // Create containers for input
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  mgr->ConnectInput(task, 0, cinput);

  // Create containers for output
  AliAnalysisDataContainer *coutput = mgr->CreateContainer("coutput", TList::Class(), AliAnalysisManager::kOutputContainer, outFile);
  mgr->ConnectOutput(task, 0, coutput);

  // Enable debug printouts
  mgr->SetDebugLevel(0);

  if (!mgr->InitAnalysis())
    return;

  mgr->PrintStatus();

  if(bProof) mgr->StartAnalysis("proof",chain, nEvents, firstEvent);
  else mgr->StartAnalysis("local",chain,nEvents, firstEvent);
}

