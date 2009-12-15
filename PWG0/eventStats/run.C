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

void run(const Char_t* data, Int_t nRuns=20, Int_t offset=0, Bool_t aDebug = kFALSE, Int_t aProof = kFALSE, const char* option = "")
{
  // aProof option: 0 no proof
  //                1 proof with chain
  //                2 proof with dataset
  //
  // option is passed to the task(s)

  if (nRuns < 0)
    nRuns = 1234567890;

  if (aProof > 0)
  {
    TProof::Open("alicecaf"); 
    //gProof->SetParallel(2);

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
      gProof->UploadPackage("$ALICE_ROOT/AF-v4-18-12-AN.par");
      gProof->EnablePackage("AF-v4-18-12-AN");
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
    gSystem->Load("libPWG0dep");
  }
  
  // Create the analysis manager
  mgr = new AliAnalysisManager;

  // Add ESD handler
  AliESDInputHandler* esdH = new AliESDInputHandler;
  esdH->SetInactiveBranches("AliRawDataErrorLogs CaloClusters Cascades EMCALCells EMCALTrigger ESDfriend Kinks Kinks Cascades AliESDTZERO MuonTracks TrdTracks CaloClusters");
  mgr->SetInputEventHandler(esdH);

  cInput = mgr->GetCommonInputContainer();
  
  Load("AliEventStatsTask", aDebug);
  TString optStr(option);
  
  // remove SAVE option if set
  Bool_t save = kFALSE;
  if (optStr.Contains("SAVE"))
  {
    optStr = optStr(0,optStr.Index("SAVE")) + optStr(optStr.Index("SAVE")+4, optStr.Length());
    save = kTRUE;
  }
  
  task = new AliEventStatsTask(optStr);
  physicsSelection = new AliPhysicsSelection;
  task->SetPhysicsSelection(physicsSelection);
  
  mgr->AddTask(task);

  // Attach input
  mgr->ConnectInput(task, 0, cInput);

  // Attach output
  cOutput = mgr->CreateContainer("cOutput", TList::Class(), AliAnalysisManager::kOutputContainer);
  mgr->ConnectOutput(task, 1, cOutput);

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

    if (save)
    {
      TString path("maps/");
      path += TString(data).Tokenize("/")->Last()->GetName();
      
      gSystem->mkdir(path, kTRUE);
      gSystem->Rename("event_stats.root", path + "/event_stats.root");
      
      Printf(">>>>> Moved files to %s", path.Data());
    }  
  }
  else if (aProof == 3)
  {
    gROOT->ProcessLine(".L CreateChainFromDataSet.C");
    ds = gProof->GetDataSet(data)->GetStagedSubset();
    file = TFile::Open("dataset.root", "RECREATE");
    ds->Write("dataset");
    file->Close();
    chain = CreateChainFromDataSet(ds);
    mgr->StartAnalysis("local", chain, nRuns, offset);
  }
  else if (aProof == -1)
  {
    gROOT->ProcessLine(".L CreateChainFromDataSet.C");
    TFile::Open("dataset.root");
    ds = (TFileCollection*) gFile->Get("dataset");
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
