void run(Char_t* data, Int_t nRuns=20, Int_t offset=0, Bool_t aDebug = kFALSE, Bool_t aProof = kFALSE, Bool_t mc = kTRUE, const char* option = "")
{
  if (aProof)
  {
    TProof::Open("lxb6046");

    // Enable the needed package
    gProof->UploadPackage("STEERBase");
    gProof->EnablePackage("STEERBase");
    gProof->UploadPackage("ESD");
    gProof->EnablePackage("ESD");
    gProof->UploadPackage("AOD");
    gProof->EnablePackage("AOD");
    gProof->UploadPackage("ANALYSIS");
    gProof->EnablePackage("ANALYSIS");
    gProof->UploadPackage("PWG0base");
    gProof->EnablePackage("PWG0base");
  }
  else
  {
    gSystem->Load("libVMC");
    gSystem->Load("libTree");
    gSystem->Load("libSTEERBase");
    gSystem->Load("libESD");
    gSystem->Load("libANALYSIS");
    gSystem->Load("libPWG0base");
  }

  // Create chain of input files
  gROOT->LoadMacro("../CreateESDChain.C");
  chain = CreateESDChain(data, nRuns, offset);

  // Create the analysis manager
  mgr = new AliAnalysisManager;

  TString taskName("AlidNdEtaTask.cxx+");
  if (aDebug)
    taskName += "+g";

  // Create, add task
  if (aProof) {
    gProof->Load(taskName);
  } else
    gROOT->Macro(taskName);

  task = new AlidNdEtaTask(option);

  AliPWG0Helper::AnalysisMode analysisMode = AliPWG0Helper::kTPC;
  task->SetAnalysisMode(analysisMode);

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

  if (mc)
    task->SetReadMC();

  mgr->AddTask(task);

  if (mc) {
    // Enable MC event handler
    AliMCEventHandler* handler = new AliMCEventHandler;
    handler->SetReadTR(kFALSE);
    mgr->SetMCtruthEventHandler(handler);
  }

  // Add ESD handler
  AliESDInputHandler* esdH = new AliESDInputHandler;
  mgr->SetInputEventHandler(esdH);

  // Attach input
  cInput  = mgr->CreateContainer("cInput", TChain::Class(), AliAnalysisManager::kInputContainer);
  mgr->ConnectInput(task, 0, cInput);

  // Attach output
  cOutput = mgr->CreateContainer("cOutput", TList::Class(), AliAnalysisManager::kOutputContainer);
  mgr->ConnectOutput(task, 0, cOutput);

  // Enable debug printouts
  if (aDebug)
    mgr->SetDebugLevel(2);

  // Run analysis
  mgr->InitAnalysis();
  mgr->PrintStatus();

  mgr->StartAnalysis((aProof) ? "proof" : "local", chain);
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

  dNdEtaAnalysis* fdNdEtaAnalysis = new dNdEtaAnalysis("dndeta", "dndeta");
  fdNdEtaAnalysis->LoadHistograms("fdNdEtaAnalysisESD");
  fdNdEtaAnalysis->Finish(dNdEtaCorrection, 0.3, AlidNdEtaCorrection::kINEL, 1);
  //fdNdEtaAnalysis->DrawHistograms(kTRUE);
  TFile* file2 = TFile::Open(dataOutput, "RECREATE");
  fdNdEtaAnalysis->SaveHistograms();

  file->cd();
  fdNdEtaAnalysis = new dNdEtaAnalysis("dndetaTr", "dndetaTr");
  fdNdEtaAnalysis->LoadHistograms("fdNdEtaAnalysisESD");
  fdNdEtaAnalysis->Finish(dNdEtaCorrection, 0.3, AlidNdEtaCorrection::kVertexReco, 1);
  //fdNdEtaAnalysis->DrawHistograms(kTRUE);
  file2->cd();
  fdNdEtaAnalysis->SaveHistograms();

  file->cd();
  fdNdEtaAnalysis = new dNdEtaAnalysis("dndetaTrVtx", "dndetaTrVtx");
  fdNdEtaAnalysis->LoadHistograms("fdNdEtaAnalysisESD");
  fdNdEtaAnalysis->Finish(dNdEtaCorrection, 0.3, AlidNdEtaCorrection::kTrack2Particle, 1);
  //fdNdEtaAnalysis->DrawHistograms(kTRUE);
  file2->cd();
  fdNdEtaAnalysis->SaveHistograms();

  file->cd();
  fdNdEtaAnalysis = new dNdEtaAnalysis("dndetaTracks", "dndetaTracks");
  fdNdEtaAnalysis->LoadHistograms("fdNdEtaAnalysisESD");
  fdNdEtaAnalysis->Finish(0, 0.3, AlidNdEtaCorrection::kNone, 1);
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
