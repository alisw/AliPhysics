AliAnalysisTaskSEDStarCharmFraction *AddTaskDStarCharmFraction(TString cutsFileName, Bool_t readMC = kFALSE, TString suffix = "", Double_t peakSigmaCut = 3., Double_t sidebandSigmaCut = 6., Double_t sidebandSigmaWidth = 9., Bool_t singleSideband = kFALSE, Double_t impParMin = 0., TString fileOut = "AnalysisResults.root") {
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

  /*  1., 2., 3., 4., 5., 6., 7., 8., 10., 12., 16., 20., 24., 36. */
  /*  1., 2., 3., 4., 5., 6., 7., 8.,      12., 16.,      24.      */

  /*
     1- 2: 0.000701202
     2- 3: 0.000630536
     3- 4: 0.000558947
     4- 5: 0.000559443
     5- 6: 0.000562800
     6- 7: 0.000552661
     7- 8: 0.000552661
     8-10: 0.000578663
    10-12: 0.000578663
    12-16: 0.000650723
    16-20: 0.000755267
    20-24: 0.000755267
    24-36: 0.000755267
  */
  //                       1-2       2-3     3-4      4-5      5-6      6-7      7-8      8-10    10-12    12-16    16-20    20-24    24-36
  Double_t sigmas[13] = {0.00070, 0.00063, 0.00056, 0.00056, 0.00056, 0.00055, 0.00055, 0.00058, 0.00058, 0.00065, 0.00076, 0.00076, 0.00076};
  Double_t peakCuts[13];
  Double_t sidebandCuts[13];
  Double_t sidebandWindows[13];

  for (Int_t i=0;i<13;i++) {
    peakCuts[i] = peakSigmaCut*sigmas[i];
    sidebandCuts[i] = sidebandSigmaCut*sigmas[i]; //sideband selection starts at "sidebandSigmaCut" (default: 6) sigma from PDG mass
    sidebandWindows[i] = sidebandSigmaWidth*sigmas[i]; //sideband selection has width of "sidebandSigmaWidth" (default: 9) sigma (so stops at sidebandSigmaCut+sidebandSigmaWidth sigma from PDG mass)

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