AliAnalysisTaskSEDStarCharmFraction *AddTaskDStarCharmFraction(TString cutsFileName, Bool_t readMC = kFALSE, TString suffix = "", Double_t peakSigmaCut = 3., Double_t sidebandSigmaCut = 6., Bool_t singleSideband = kFALSE, Double_t impParMin = 0., TString fileOut = "AnalysisResults.root") {
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    AliFatal("No analysis manager to connect to");
    return NULL;
  }

  Bool_t standardCuts = kFALSE;

  TFile *cutsFile;
  if (cutsFileName.EqualTo("")) {
    standardCuts = kTRUE; 
  }
  else {
    cutsFile = TFile::Open(cutsFileName.Data());
    if (!cutsFile || (cutsFile && !cutsFile->IsOpen())) {
      AliFatal("Input cut file not found");
      return NULL;
    }
  }


  AliRDHFCutsDStartoKpipi *cuts;
  if (standardCuts) {
    cuts = new AliRDHFCutsDStartoKpipi();
    cuts->SetStandardCutsPP2010();
  }
  else {
    cuts = (AliRDHFCutsDStartoKpipi*) cutsFile->Get("DStartoKpipiCuts");
  }
  cuts->SetName("DStartoKpipiCuts");

  if (!cuts) {
    AliFatal("Cut object not found");
    return NULL;
  }

  printf("CREATE TASK\n");

  AliAnalysisTaskSEDStarCharmFraction *task = new AliAnalysisTaskSEDStarCharmFraction("AliAnalysisTaskSEDStarSpectra", cuts);
  task->SetReadMC(readMC);
  task->SetSingleSideband(singleSideband);
  task->SetImpParCut(impParMin);

  // Legacy: set constant width
  Int_t nPtBins = cuts->GetNPtBins();
  Double_t *sigmas = new Double_t[nPtBins];
  Double_t *peakCuts = new Double_t[nPtBins];
  Double_t *sidebandCuts = new Double_t[nPtBins];
  Double_t *sidebandWindows = new Double_t[nPtBins];

  for (Int_t i=0;i<nPtBins;i++) {
    sigmas[i] = 0.00060;
    peakCuts[i] = peakSigmaCut*sigmas[i];
    sidebandCuts[i] = sidebandSigmaCut*sigmas[i]; //sideband selection starts at "sidebandSigmaCut" sigma from PDG mass
    sidebandWindows[i] = 9.*sigmas[i]; //sideband selection has width of 9 sigma (so stops at sidebandSigmaCut+9 sigma from PDG mass)

    cout << "peakCuts[" << i << "]=" << peakCuts[i] << endl;
  }

  /*Double_t sigma = 0.0005;
  task->SetPeakCut(3.*sigma);
  task->SetSidebandCut(6.*sigma);
  task->SetSidebandWindow(9.*sigma);*/

  task->SetPeakCut(peakCuts);
  task->SetSidebandCut(sidebandCuts);
  task->SetSidebandWindow(sidebandWindows);

  mgr->AddTask(task);

  AliAnalysisDataContainer *cInput = mgr->GetCommonInputContainer();
  mgr->ConnectInput(task, 0, cInput);

  AliAnalysisDataContainer *cOutputNEvents = mgr->CreateContainer(Form("cNEvents%s", suffix.Data()), TH1D::Class(), AliAnalysisManager::kOutputContainer, fileOut.Data());
  AliAnalysisDataContainer *cOutputListCandidate = mgr->CreateContainer(Form("cListCandidate%s", suffix.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileOut.Data());
  AliAnalysisDataContainer *cOutputListSignal = mgr->CreateContainer(Form("cListSignal%s", suffix.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileOut.Data());
  AliAnalysisDataContainer *cOutputListSignalPrompt = mgr->CreateContainer(Form("cListSignalPrompt%s", suffix.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileOut.Data());
  AliAnalysisDataContainer *cOutputListSignalFromB = mgr->CreateContainer(Form("cListSignalFromB%s", suffix.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileOut.Data());
  AliAnalysisDataContainer *cOutputListBackground = mgr->CreateContainer(Form("cListBackground%s", suffix.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileOut.Data());
  AliAnalysisDataContainer *cOutputCuts = mgr->CreateContainer(Form("cCuts%s", suffix.Data()), AliRDHFCutsDStartoKpipi::Class(), AliAnalysisManager::kOutputContainer, fileOut.Data());
  AliAnalysisDataContainer *cOutputNormCount = mgr->CreateContainer(Form("cNormCount%s", suffix.Data()), AliNormalizationCounter::Class(), AliAnalysisManager::kOutputContainer, fileOut.Data());
  AliAnalysisDataContainer *cOutputTreeCandidate = mgr->CreateContainer(Form("cTreeCandidate%s", suffix.Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, fileOut.Data());

  mgr->ConnectOutput(task, 1, cOutputNEvents);
  mgr->ConnectOutput(task, 2, cOutputListCandidate);
  mgr->ConnectOutput(task, 3, cOutputListSignal);
  mgr->ConnectOutput(task, 4, cOutputListSignalPrompt);
  mgr->ConnectOutput(task, 5, cOutputListSignalFromB);
  mgr->ConnectOutput(task, 6, cOutputListBackground);
  mgr->ConnectOutput(task, 7, cOutputCuts);
  mgr->ConnectOutput(task, 8, cOutputNormCount);
  mgr->ConnectOutput(task, 9, cOutputTreeCandidate);

  return task;
}
