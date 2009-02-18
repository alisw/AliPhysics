void run(Char_t* data, Long64_t nRuns = -1, Long64_t offset = 0, Bool_t aDebug = kFALSE, Int_t aProof = 0, Bool_t mc = kTRUE, const char* option = "")
{
  // aProof option: 0 no proof
  //                1 proof with chain
  //                2 proof with dataset

  if (nRuns < 0)
    nRuns = 1234567890;

  if (aProof)
  {
    TProof::Open("lxb6046");

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
      gProof->UploadPackage("$ALICE_ROOT/AF-v4-12");
      gProof->EnablePackage("$ALICE_ROOT/AF-v4-12");
    }

    gProof->UploadPackage("$ALICE_ROOT/PWG0base");
    gProof->EnablePackage("$ALICE_ROOT/PWG0base");
  }
  else
  {
    gSystem->Load("libVMC");
    gSystem->Load("libTree");
    gSystem->Load("libSTEERBase");
    gSystem->Load("libESD");
    gSystem->Load("libAOD");
    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice");
    gSystem->Load("libPWG0base");
  }

  // Create the analysis manager
  mgr = new AliAnalysisManager;

  AliPWG0Helper::AnalysisMode analysisMode = AliPWG0Helper::kTPC;
  AliPWG0Helper::Trigger      trigger      = AliPWG0Helper::kMB1;

  AliPWG0Helper::PrintConf(analysisMode, trigger);

  TString taskName("AliMultiplicityTask.cxx+");
  if (aDebug)
    taskName += "+g";

  // Create, add task
  if (aProof > 0) {
    gProof->Load(taskName);
  } else
    gROOT->Macro(taskName);

  task = new AliMultiplicityTask(option);

  if (analysisMode != AliPWG0Helper::kSPD)
  {
    // selection of esd tracks
    gROOT->ProcessLine(".L ../CreateStandardCuts.C");
    AliESDtrackCuts* esdTrackCuts = CreateTrackCuts(analysisMode);
    if (!esdTrackCuts)
    {
      printf("ERROR: esdTrackCuts could not be created\n");
      return;
    }

    task->SetTrackCuts(esdTrackCuts);
  }

  task->SetAnalysisMode(analysisMode);
  task->SetTrigger(trigger);

  if (mc)
    task->SetReadMC();

  //task->SetUseMCVertex();

  mgr->AddTask(task);

  if (mc) {
    // Enable MC event handler
    AliMCEventHandler* handler = new AliMCEventHandler;
    handler->SetReadTR(kFALSE);
    mgr->SetMCtruthEventHandler(handler);
  }

  // pt study
  TString optionStr(option);
  if (optionStr.Contains("pt-spectrum-func"))
  {
    //TF1* func = new TF1("func", "0.7 + x", 0, 0.3);
    //TF1* func = new TF1("func", "1.3 - x", 0, 0.3);
    //TF1* func = new TF1("func", "1", 0, 0.3);
    //new TCanvas; func->Draw();
    //inputList.Add(func->GetHistogram()->Clone("pt-spectrum"));

    TFile* file = TFile::Open("ptspectrum_fit.root");
    if (!file)
    {
      Printf("Could not open ptspectrum_fit.root");
      return;
    }

    TString subStr(optionStr(optionStr.Index("pt-spectrum-func")+17, 3));
    TString histName(Form("ptspectrum_%s", subStr.Data()));
    Printf("Pt-Spectrum modification. Using %s.", histName.Data());
    TH1* hist = (TH1*) file->Get(histName);
    if (!hist)
    {
      Printf("Could not read histogram.");
      return;
    }

    new TCanvas; hist->Draw();
    task->SetPtSpectrum((TH1*) hist->Clone("pt-spectrum"));
  }

  // Add ESD handler
  AliESDInputHandler* esdH = new AliESDInputHandler;
  mgr->SetInputEventHandler(esdH);

  // Attach input
  cInput  = mgr->GetCommonInputContainer();
  mgr->ConnectInput(task, 0, cInput);

  // Attach output
  cOutput = mgr->CreateContainer("cOutput", TList::Class(), AliAnalysisManager::kOutputContainer);
  //cOutput->SetDataOwned(kTRUE);
  mgr->ConnectOutput(task, 0, cOutput);

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
  else
  {
    // Create chain of input files
    gROOT->LoadMacro("../CreateESDChain.C");
    chain = CreateESDChain(data, nRuns, offset);

    mgr->StartAnalysis((aProof > 0) ? "proof" : "local", chain);
  }
}
