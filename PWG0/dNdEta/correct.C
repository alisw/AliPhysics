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
  fdNdEtaAnalysis->Finish(dNdEtaCorrection, 0.21, AlidNdEtaCorrection::kNSD, "ESD -> NSD");
  //fdNdEtaAnalysis->DrawHistograms(kTRUE);
  TFile* file2 = TFile::Open(dataOutput, "RECREATE");
  fdNdEtaAnalysis->SaveHistograms();

  file->cd();
  dNdEtaAnalysis* fdNdEtaAnalysis = new dNdEtaAnalysis("dndeta", "dndeta");
  fdNdEtaAnalysis->LoadHistograms("fdNdEtaAnalysisESD");
  fdNdEtaAnalysis->Finish(dNdEtaCorrection, 0.21, AlidNdEtaCorrection::kINEL, "ESD -> full inelastic");
  //fdNdEtaAnalysis->DrawHistograms(kTRUE);
  file2->cd();
  fdNdEtaAnalysis->SaveHistograms();

  file->cd();
  fdNdEtaAnalysis = new dNdEtaAnalysis("dndetaTr", "dndetaTr");
  fdNdEtaAnalysis->LoadHistograms("fdNdEtaAnalysisESD");
  fdNdEtaAnalysis->Finish(dNdEtaCorrection, 0.21, AlidNdEtaCorrection::kVertexReco, "ESD -> minimum bias");
  //fdNdEtaAnalysis->DrawHistograms(kTRUE);
  file2->cd();
  fdNdEtaAnalysis->SaveHistograms();

  file->cd();
  fdNdEtaAnalysis = new dNdEtaAnalysis("dndetaTrVtx", "dndetaTrVtx");
  fdNdEtaAnalysis->LoadHistograms("fdNdEtaAnalysisESD");
  fdNdEtaAnalysis->Finish(dNdEtaCorrection, 0.21, AlidNdEtaCorrection::kTrack2Particle, "ESD -> MB with vertex");
  //fdNdEtaAnalysis->DrawHistograms(kTRUE);
  file2->cd();
  fdNdEtaAnalysis->SaveHistograms();

  file->cd();
  fdNdEtaAnalysis = new dNdEtaAnalysis("dndetaTracks", "dndetaTracks");
  fdNdEtaAnalysis->LoadHistograms("fdNdEtaAnalysisESD");
  fdNdEtaAnalysis->Finish(0, 0.21, AlidNdEtaCorrection::kNone, "ESD raw with pt cut");
  //fdNdEtaAnalysis->DrawHistograms(kTRUE);
  file2->cd();
  fdNdEtaAnalysis->SaveHistograms();

  file->cd();
  fdNdEtaAnalysis = new dNdEtaAnalysis("dndetaTracksAll", "dndetaTracksAll");
  fdNdEtaAnalysis->LoadHistograms("fdNdEtaAnalysisESD");
  fdNdEtaAnalysis->Finish(0, -1, AlidNdEtaCorrection::kNone, "ESD raw");
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

    fdNdEtaAnalysis->Finish(dNdEtaCorrection, 0.21, AlidNdEtaCorrection::kINEL);
    //fdNdEtaAnalysis->Finish(dNdEtaCorrection, 0, AlidNdEtaCorrection::kINEL);
    //fdNdEtaAnalysis->Finish(dNdEtaCorrection, 0, AlidNdEtaCorrection::kTrack2Particle);
  }
  else
    fdNdEtaAnalysis->Finish(0, 0.21, AlidNdEtaCorrection::kNone);

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

void correct(Bool_t onlyESD = kFALSE)
{
  FinishAnalysisAll();
  gROOT->ProcessLine(".L $ALICE_ROOT/PWG0/dNdEta/drawPlots.C");
  dNdEta(onlyESD);
}


