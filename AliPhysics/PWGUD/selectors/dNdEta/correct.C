void loadlibs()
{
  gSystem->Load("libTree");
  gSystem->Load("libVMC");

  gSystem->Load("libSTEERBase");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libPWG0base");
}

void FinishAnalysisAll(const char* dataInput = "analysis_esd_raw.root", const char* dataOutput = "analysis_esd.root", const char* correctionMapFile = "correction_map2.root", const char* correctionMapFolder = "dndeta_correction_ua5")
{
  loadlibs();

  AlidNdEtaCorrection* dNdEtaCorrection = 0;
  
  if (correctionMapFile)
  {
    dNdEtaCorrection = new AlidNdEtaCorrection(correctionMapFolder, correctionMapFolder);
    if (!TFile::Open(correctionMapFile))
      return;
    dNdEtaCorrection->LoadHistograms();
  }
  
  TFile* file = TFile::Open(dataInput);

  if (!file)
  {
    cout << "Error. File not found" << endl;
    return;
  }
  
  Int_t backgroundEvents = 0;
  
  //backgroundEvents = 1162+434; // Michele for MB1, run 104892, 15.02.10
  //backgroundEvents = 842;    // JF estimate for V0 systematic case 1
  //backgroundEvents = 6;          // Michele for V0AND, run 104892, 15.02.10
  
  //backgroundEvents = 1758+434; // MB1, run 104867-92
  
  //backgroundEvents = 4398+961;   // Michele for MB1, run 104824-52, 16.02.10
  //backgroundEvents = 19;         // Michele for V0AND, run 104824-52, 16.02.10
  
  backgroundEvents = -1;    // use 0 bin from MC! for 2.36 TeV
  //backgroundEvents = 900; // my estimate for 2.36 TeV
  
  //backgroundEvents = 918;   // V0OR for run 114786 w/o bunch intensities w/o proper 0 checking!
  //backgroundEvents = 723; // V0OR for run 114798 w/o bunch intensities w/o proper 0 checking!
  
  Printf("Subtracting %d background events!!!", backgroundEvents);
  gSystem->Sleep(1000);
  
  TH1* combinatoricsCorrection = 0;
  if (1)
  {
    TFile::Open("corrComb.root");
    combinatoricsCorrection = (TH1*) gFile->Get("ratioofratios");
    if (!combinatoricsCorrection)
      combinatoricsCorrection = (TH1*) gFile->Get("correctiondNdEta");
    file->cd();
  }
  
  // Note: the last parameter does not define which analysis is going to happen, the histograms will be overwritten when loading from the file
  dNdEtaAnalysis* fdNdEtaAnalysis = new dNdEtaAnalysis("dndetaNSD", "dndetaNSD");
  fdNdEtaAnalysis->LoadHistograms("fdNdEtaAnalysisESD");
  fdNdEtaAnalysis->Finish(dNdEtaCorrection, 0.151, AlidNdEtaCorrection::kNSD, "ESD -> NSD", backgroundEvents, combinatoricsCorrection);
  //fdNdEtaAnalysis->DrawHistograms(kTRUE);
  TFile* file2 = TFile::Open(dataOutput, "RECREATE");
  fdNdEtaAnalysis->SaveHistograms();
  
  file->cd();
  dNdEtaAnalysis* fdNdEtaAnalysis = new dNdEtaAnalysis("dndeta", "dndeta");
  fdNdEtaAnalysis->LoadHistograms("fdNdEtaAnalysisESD");
  fdNdEtaAnalysis->Finish(dNdEtaCorrection, 0.151, AlidNdEtaCorrection::kINEL, "ESD -> full inelastic", backgroundEvents, combinatoricsCorrection);
  //fdNdEtaAnalysis->DrawHistograms(kTRUE);
  file2->cd();
  fdNdEtaAnalysis->SaveHistograms();

  file->cd();
  fdNdEtaAnalysis = new dNdEtaAnalysis("dndetaTr", "dndetaTr");
  fdNdEtaAnalysis->LoadHistograms("fdNdEtaAnalysisESD");
  fdNdEtaAnalysis->Finish(dNdEtaCorrection, 0.151, AlidNdEtaCorrection::kVertexReco, "ESD -> minimum bias", backgroundEvents, combinatoricsCorrection);
  //fdNdEtaAnalysis->DrawHistograms(kTRUE);
  file2->cd();
  fdNdEtaAnalysis->SaveHistograms();

  file->cd();
  fdNdEtaAnalysis = new dNdEtaAnalysis("dndetaTrVtx", "dndetaTrVtx");
  fdNdEtaAnalysis->LoadHistograms("fdNdEtaAnalysisESD");
  fdNdEtaAnalysis->Finish(dNdEtaCorrection, 0.151, AlidNdEtaCorrection::kTrack2Particle, "ESD -> MB with vertex", backgroundEvents, combinatoricsCorrection);
  //fdNdEtaAnalysis->DrawHistograms(kTRUE);
  file2->cd();
  fdNdEtaAnalysis->SaveHistograms();

  file->cd();
  fdNdEtaAnalysis = new dNdEtaAnalysis("dndetaTracks", "dndetaTracks");
  fdNdEtaAnalysis->LoadHistograms("fdNdEtaAnalysisESD");
  fdNdEtaAnalysis->Finish(0, 0.151, AlidNdEtaCorrection::kNone, "ESD raw with pt cut", backgroundEvents, combinatoricsCorrection);
  //fdNdEtaAnalysis->DrawHistograms(kTRUE);
  file2->cd();
  fdNdEtaAnalysis->SaveHistograms();

  file->cd();
  fdNdEtaAnalysis = new dNdEtaAnalysis("dndetaTracksAll", "dndetaTracksAll");
  fdNdEtaAnalysis->LoadHistograms("fdNdEtaAnalysisESD");
  fdNdEtaAnalysis->Finish(0, -1, AlidNdEtaCorrection::kNone, "ESD raw", backgroundEvents, combinatoricsCorrection);
  //fdNdEtaAnalysis->DrawHistograms(kTRUE);
  file2->cd();
  fdNdEtaAnalysis->SaveHistograms();

  file->cd();
  fdNdEtaAnalysis = new dNdEtaAnalysis("dndetaOnePart", "dndetaOnePart");
  fdNdEtaAnalysis->LoadHistograms("fdNdEtaAnalysisESD");
  fdNdEtaAnalysis->Finish(dNdEtaCorrection, 0.151, AlidNdEtaCorrection::kOnePart, "ESD -> OnePart", backgroundEvents, combinatoricsCorrection);
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

    fdNdEtaAnalysis->Finish(dNdEtaCorrection, 0.151, AlidNdEtaCorrection::kINEL);
    //fdNdEtaAnalysis->Finish(dNdEtaCorrection, 0, AlidNdEtaCorrection::kINEL);
    //fdNdEtaAnalysis->Finish(dNdEtaCorrection, 0, AlidNdEtaCorrection::kTrack2Particle);
  }
  else
    fdNdEtaAnalysis->Finish(0, 0.151, AlidNdEtaCorrection::kNone);

  return fdNdEtaAnalysis;
  
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

void FinishMC()
{
  loadlibs();
  
  result1 = (dNdEtaAnalysis*) FinishAnalysis("analysis_mc.root", "dndeta", 0, 0, 1);
  result2 = (dNdEtaAnalysis*) FinishAnalysis("analysis_mc.root", "dndetaNSD", 0, 0, 1);
  
  file = TFile::Open("out.root", "RECREATE");
  result1->SaveHistograms();
  result2->SaveHistograms();
  file->Close();
}

void correct(Bool_t onlyESD = kFALSE, Bool_t mergedXSections = kTRUE)
{
  gSystem->Unlink("analysis_esd.root");
  if (mergedXSections)
    FinishAnalysisAll("analysis_esd_raw.root", "analysis_esd.root", "correction_map2.root", "dndeta_correction_ua5");
  else
    FinishAnalysisAll("analysis_esd_raw.root", "analysis_esd.root", "correction_map.root", "dndeta_correction");
  
  gROOT->ProcessLine(".L $ALICE_ROOT/PWG0/dNdEta/drawPlots.C");
  dNdEta(onlyESD);
}


