void Load(const char* taskName, Bool_t debug)
{
  TString compileTaskName;
  compileTaskName.Form("%s.cxx++", taskName);
  if (debug)
    compileTaskName += "g";

  if (gProof) {
    gProof->Load(compileTaskName);
  } else
    gROOT->Macro(compileTaskName);

  // Enable debug printouts
  if (debug)
  {
    AliLog::SetClassDebugLevel(taskName, AliLog::kDebug+2);
  }
  else
    AliLog::SetClassDebugLevel(taskName, AliLog::kWarning);
}

void run(Int_t runWhat, const Char_t* data, Int_t nRuns=20, Int_t offset=0, Bool_t aDebug = kFALSE, Int_t aProof = kFALSE, Bool_t mc = kTRUE, const char* option = "")
{
  // runWhat options: 0 = AlidNdEtaTask
  //                  1 = AlidNdEtaCorrectionTask
  //                  2 = both
  //
  // aProof option: 0 no proof
  //                1 proof with chain
  //                2 proof with dataset
  
  TString taskName;
  if (runWhat == 0 || runWhat == 2)
  {
    Printf("Running AlidNdEtaTask");
  }
  if (runWhat == 1 || runWhat == 2)
  {
    Printf("Running AlidNdEtaCorrectionTask");
    if (!mc)
    {
      Printf("AlidNdEtaCorrectionTask needs MC. Exiting...");
      return;
    }
  }

  if (nRuns < 0)
    nRuns = 1234567890;

  if (aProof)
  {
    TProof::Open("alicecaf");
    //gProof->SetParallel(1);

    // Enable the needed package
    if (1)
    {
      gProof->UploadPackage("$ALICE_ROOT/STEERBase");
      gProof->EnablePackage("$ALICE_ROOT/STEERBase");
      gProof->UploadPackage("$ALICE_ROOT/ESD");
      gProof->EnablePackage("$ALICE_ROOT/ESD");
      gProof->UploadPackage("$ALICE_ROOT/AOD");
      gProof->EnablePackage("$ALICE_ROOT/AOD");
      gProof->UploadPackage("$ALICE_ROOT/ANALYSIS");
      gProof->EnablePackage("$ALICE_ROOT/ANALYSIS");
      gProof->UploadPackage("$ALICE_ROOT/ANALYSISalice");
      gProof->EnablePackage("$ALICE_ROOT/ANALYSISalice");
    }
    else
    {
      gProof->UploadPackage("/afs/cern.ch/alice/caf/sw/ALICE/PARs/v4-16-Release/AF-v4-16");
      gProof->EnablePackage("AF-v4-16");
    }

    gProof->UploadPackage("$ALICE_ROOT/PWG0base");
    gProof->EnablePackage("$ALICE_ROOT/PWG0base");
  }
  else
  {
    gSystem->AddIncludePath("-I${ALICE_ROOT}/include/ -I${ALICE_ROOT}/PWG0/ -I${ALICE_ROOT}/PWG0/dNdEta/"); 
    gSystem->Load("libVMC");
    gSystem->Load("libTree");
    gSystem->Load("libProof");
    gSystem->Load("libSTEERBase");
    gSystem->Load("libESD");
    gSystem->Load("libAOD");
    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice");
    gSystem->Load("libPWG0base");
  }

  // Create the analysis manager
  mgr = new AliAnalysisManager;

  // Add ESD handler
  AliESDInputHandler* esdH = new AliESDInputHandler;
  esdH->SetInactiveBranches("AliESDACORDE FMD ALIESDTZERO ALIESDVZERO ALIESDZDC AliRawDataErrorLogs CaloClusters Cascades EMCALCells EMCALTrigger ESDfriend Kinks Kinks Cascades AliESDTZERO ALIESDACORDE MuonTracks TrdTracks CaloClusters");
  mgr->SetInputEventHandler(esdH);

  AliPWG0Helper::AnalysisMode analysisMode = AliPWG0Helper::kSPD;
  AliPWG0Helper::Trigger      trigger      = AliPWG0Helper::kMB1;

  AliPWG0Helper::PrintConf(analysisMode, trigger);

  AliESDtrackCuts* esdTrackCuts = 0;
  if (analysisMode != AliPWG0Helper::kSPD)
  {
    // selection of esd tracks
    gROOT->ProcessLine(".L ../CreateStandardCuts.C");
    esdTrackCuts = CreateTrackCuts(analysisMode);
    if (!esdTrackCuts)
    {
      printf("ERROR: esdTrackCuts could not be created\n");
      return;
    }
    esdTrackCuts->SetHistogramsOn(kTRUE);
  }

  cInput  = mgr->GetCommonInputContainer();

  // Create, add task
  if (runWhat == 0 || runWhat == 2)
  {
    Load("AlidNdEtaTask", aDebug);
    task = new AlidNdEtaTask(option);

    if (mc)
      task->SetReadMC();

    // syst. error flags
    //task->SetUseMCVertex();
    //task->SetUseMCKine();
    //task->SetOnlyPrimaries();
    //task->SetFillPhi();
    
    task->SetTrigger(trigger);
    task->SetAnalysisMode(analysisMode);
    task->SetTrackCuts(esdTrackCuts);
    task->SetDeltaPhiCut(0.05);

    mgr->AddTask(task);

    // Attach input
    mgr->ConnectInput(task, 0, cInput);

    // Attach output
    cOutput = mgr->CreateContainer("cOutput", TList::Class(), AliAnalysisManager::kOutputContainer);
    mgr->ConnectOutput(task, 0, cOutput);
  }
  if (runWhat == 1 || runWhat == 2)
  {
    Load("AlidNdEtaCorrectionTask", aDebug);
    task2 = new AlidNdEtaCorrectionTask(option);

    // syst. error flags
    //task2->SetFillPhi();
    //task2->SetOnlyPrimaries();

    task2->SetTrigger(trigger);
    task2->SetAnalysisMode(analysisMode);
    task2->SetTrackCuts(esdTrackCuts);
    //task2->SetDeltaPhiCut(0.05);

    mgr->AddTask(task2);

    // Attach input
    mgr->ConnectInput(task2, 0, cInput);

    // Attach output
    cOutput = mgr->CreateContainer("cOutput2", TList::Class(), AliAnalysisManager::kOutputContainer);
    mgr->ConnectOutput(task2, 0, cOutput);
  }

  if (mc) {
    // Enable MC event handler
    AliMCEventHandler* handler = new AliMCEventHandler;
    handler->SetReadTR(kFALSE);
    mgr->SetMCtruthEventHandler(handler);
  }

  // Enable debug printouts
  if (aDebug)
    mgr->SetDebugLevel(2);

  // Run analysis
  mgr->InitAnalysis();
  mgr->PrintStatus();

  if (aProof == 2)
  {
    // process dataset

    mgr->StartAnalysis("proof", data, nRuns, offset);
  }
  else if (aProof == 3)
  {
    gROOT->ProcessLine(".L CreateChainFromDataSet.C");
    ds = gProof->GetDataSet(data)->GetStagedSubset();
    chain = CreateChainFromDataSet(ds);
    mgr->StartAnalysis("local", chain, nRuns, offset);
  }
  else
  {
    // Create chain of input files
    gROOT->LoadMacro("../CreateESDChain.C");

    chain = CreateESDChain(data, nRuns, offset);
    //chain = CreateChain("TE", data, nRuns, offset);

    mgr->StartAnalysis((aProof > 0) ? "proof" : "local", chain);
  }
}

