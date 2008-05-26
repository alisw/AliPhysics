void run(Int_t runWhat, const Char_t* data, Int_t nRuns=20, Int_t offset=0, Bool_t aDebug = kFALSE, Int_t aProof = kFALSE, Bool_t mc = kTRUE, const char* option = "")
{
  // runWhat options: 0 = AlidNdEtaTask
  //                  1 = AlidNdEtaCorrectionTask
  //
  // aProof option: 0 no proof
  //                1 proof with chain
  //                2 proof with dataset

  TString taskName;
  if (runWhat == 0)
  {
    taskName = "AlidNdEtaTask";
  }
  else if (runWhat == 1)
  {
    taskName = "AlidNdEtaCorrectionTask";
    if (!mc)
    {
      Printf("%s needs MC. Exiting...", taskName.Data());
      return;
    }
  }
  else
  {
    Printf("Do not know what to run. Exiting...");
    return;
  }

  Printf("Processing task: %s", taskName.Data());

  if (nRuns < 0)
    nRuns = 1234567890;

  if (aProof)
  {
    TProof::Open("lxb6046");
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
    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice");
    gSystem->Load("libPWG0base");
  }

  // Create the analysis manager
  mgr = new AliAnalysisManager;

  TString compileTaskName;
  compileTaskName.Form("%s.cxx+", taskName.Data());
  if (aDebug)
    compileTaskName += "+g";

  // Create, add task
  if (aProof) {
    gProof->Load(compileTaskName);
  } else
    gROOT->Macro(compileTaskName);

  AliPWG0Helper::AnalysisMode analysisMode = AliPWG0Helper::kSPD;
  AliPWG0Helper::Trigger      trigger = AliPWG0Helper::kMB1;

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
  }

  if (runWhat == 0)
  {
    task = new AlidNdEtaTask(option);

    if (mc)
      task->SetReadMC();

    //task->SetUseMCVertex();
    //task->SetUseMCKine();
  }
  else if (runWhat == 1)
  {
    task = new AlidNdEtaCorrectionTask(option);

    //task->SetOnlyPrimaries();
  }

  task->SetTrigger(trigger);
  task->SetAnalysisMode(analysisMode);
  task->SetTrackCuts(esdTrackCuts);

  mgr->AddTask(task);

  if (mc) {
    // Enable MC event handler
    AliMCEventHandler* handler = new AliMCEventHandler;
    handler->SetReadTR(kFALSE);
    mgr->SetMCtruthEventHandler(handler);
  }

  // Add ESD handler
  AliESDInputHandler* esdH = new AliESDInputHandler;
  //esdH->SetInactiveBranches("*");
  mgr->SetInputEventHandler(esdH);

  // Attach input
  cInput  = mgr->CreateContainer("cInput", TChain::Class(), AliAnalysisManager::kInputContainer);
  mgr->ConnectInput(task, 0, cInput);

  // Attach output
  cOutput = mgr->CreateContainer("cOutput", TList::Class(), AliAnalysisManager::kOutputContainer);
  mgr->ConnectOutput(task, 0, cOutput);

  // Enable debug printouts
  if (aDebug)
  {
    mgr->SetDebugLevel(2);
    AliLog::SetClassDebugLevel(taskName, AliLog::kDebug+2);
  }
  else
    AliLog::SetClassDebugLevel(taskName, AliLog::kWarning);

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

void loadlibs()
{
  gSystem->Load("libTree");
  gSystem->Load("libVMC");

  gSystem->Load("libSTEERBase");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libPWG0base");
}

void FinishAnalysisAll(const char* dataInput = "analysis_esd_raw.root", const char* dataOutput = "analysis_esd.root", const char* correctionMapFile = "correction_map.root", const char* correctionMapFolder = "dndeta_correction")
{
  loadlibs();

  AlidNdEtaCorrection* dNdEtaCorrection = new AlidNdEtaCorrection(correctionMapFolder, correctionMapFolder);
  TFile::Open(correctionMapFile);
  dNdEtaCorrection->LoadHistograms();

  TFile* file = TFile::Open(dataInput);

  if (!file)
  {
    cout << "Error. File not found" << endl;
    return;
  }

    // Note: the last parameter does not define which analysis is going to happen, the histograms will be overwritten when loading from the f
  dNdEtaAnalysis* fdNdEtaAnalysis = new dNdEtaAnalysis("dndetaNSD", "dndetaNSD");
  fdNdEtaAnalysis->LoadHistograms("fdNdEtaAnalysisESD");
  fdNdEtaAnalysis->Finish(dNdEtaCorrection, 0.3, AlidNdEtaCorrection::kNSD, "ESD -> NSD");
  //fdNdEtaAnalysis->DrawHistograms(kTRUE);
  TFile* file2 = TFile::Open(dataOutput, "RECREATE");
  fdNdEtaAnalysis->SaveHistograms();

  file->cd();
  dNdEtaAnalysis* fdNdEtaAnalysis = new dNdEtaAnalysis("dndeta", "dndeta");
  fdNdEtaAnalysis->LoadHistograms("fdNdEtaAnalysisESD");
  fdNdEtaAnalysis->Finish(dNdEtaCorrection, 0.3, AlidNdEtaCorrection::kINEL, "ESD -> full inelastic");
  //fdNdEtaAnalysis->DrawHistograms(kTRUE);
  file2->cd();
  fdNdEtaAnalysis->SaveHistograms();

  file->cd();
  fdNdEtaAnalysis = new dNdEtaAnalysis("dndetaTr", "dndetaTr");
  fdNdEtaAnalysis->LoadHistograms("fdNdEtaAnalysisESD");
  fdNdEtaAnalysis->Finish(dNdEtaCorrection, 0.3, AlidNdEtaCorrection::kVertexReco, "ESD -> minimum bias");
  //fdNdEtaAnalysis->DrawHistograms(kTRUE);
  file2->cd();
  fdNdEtaAnalysis->SaveHistograms();

  file->cd();
  fdNdEtaAnalysis = new dNdEtaAnalysis("dndetaTrVtx", "dndetaTrVtx");
  fdNdEtaAnalysis->LoadHistograms("fdNdEtaAnalysisESD");
  fdNdEtaAnalysis->Finish(dNdEtaCorrection, 0.3, AlidNdEtaCorrection::kTrack2Particle, "ESD -> MB with vertex");
  //fdNdEtaAnalysis->DrawHistograms(kTRUE);
  file2->cd();
  fdNdEtaAnalysis->SaveHistograms();

  file->cd();
  fdNdEtaAnalysis = new dNdEtaAnalysis("dndetaTracks", "dndetaTracks");
  fdNdEtaAnalysis->LoadHistograms("fdNdEtaAnalysisESD");
  fdNdEtaAnalysis->Finish(0, 0.3, AlidNdEtaCorrection::kNone, "ESD raw");
  //fdNdEtaAnalysis->DrawHistograms(kTRUE);
  file2->cd();
  fdNdEtaAnalysis->SaveHistograms();

  file2->Write();
  file2->Close();
}

void* FinishAnalysis(const char* analysisFile = "analysis_esd_raw.root", const char* analysisDir = "fdNdEtaAnalysisESD", const char* correctionMapFile = "correction_map.root", const char* correctionMapFolder = "dndeta_correction", Bool_t useUncorrected = kFALSE, Bool_t simple = kFALSE)
{
  loadlibs();

  TFile* file = TFile::Open(analysisFile);

  dNdEtaAnalysis* fdNdEtaAnalysis = new dNdEtaAnalysis(analysisDir, analysisDir);
  fdNdEtaAnalysis->LoadHistograms();

  if (correctionMapFile)
  {
    AlidNdEtaCorrection* dNdEtaCorrection = new AlidNdEtaCorrection(correctionMapFolder, correctionMapFolder);
    TFile::Open(correctionMapFile);
    dNdEtaCorrection->LoadHistograms();

    fdNdEtaAnalysis->Finish(dNdEtaCorrection, 0.3, AlidNdEtaCorrection::kINEL);
    //fdNdEtaAnalysis->Finish(dNdEtaCorrection, 0, AlidNdEtaCorrection::kINEL);
    //fdNdEtaAnalysis->Finish(dNdEtaCorrection, 0, AlidNdEtaCorrection::kTrack2Particle);
  }
  else
    fdNdEtaAnalysis->Finish(0, 0.3, AlidNdEtaCorrection::kNone);

  fdNdEtaAnalysis->DrawHistograms(simple);

  TH1* hist = fdNdEtaAnalysis->GetdNdEtaHistogram(1);
  Int_t binLeft = hist->GetXaxis()->FindBin(-0.5);
  Int_t binRight = hist->GetXaxis()->FindBin(0.5);
  Float_t value1 = hist->Integral(binLeft, binRight);

  hist = fdNdEtaAnalysis->GetdNdEtaHistogram(2);
  Float_t value2 = hist->Integral(binLeft, binRight);

  if (value2 > 0)
    printf("Ratio is %f, values are %f %f\n", value1 / value2, value1, value2);

  return fdNdEtaAnalysis;
}
