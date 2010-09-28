const char* objectName = "PWG4_LeadingTrackUE/histosLeadingTrackUE";
//const char* objectName = "PWG4_LeadingTrackUE/histosLeadingTrackUE_filterbit32";

void SetRanges(TAxis* axis)
{
  if (strcmp(axis->GetTitle(), "leading p_{T} (GeV/c)") == 0)
    axis->SetRangeUser(0, 10);

  if (strcmp(axis->GetTitle(), "multiplicity") == 0)
    axis->SetRangeUser(0, 10);
}

void SetRanges(TH1* hist)
{
  SetRanges(hist->GetXaxis());
  SetRanges(hist->GetYaxis());
  SetRanges(hist->GetZaxis());
}

void Prepare1DPlot(TH1* hist)
{
  hist->SetLineWidth(2);
  hist->SetStats(kFALSE);

  hist->GetXaxis()->SetLabelOffset(0.02);
  hist->GetXaxis()->SetTitleOffset(1.3);
  hist->GetYaxis()->SetTitleOffset(1.3);

  SetRanges(hist);
}

void DrawRatio(TH1* corr, TH1* mc, const char* name)
{
  mc->SetLineColor(2);

  TCanvas* canvas3 = new TCanvas(name, name, 600, 600);
  canvas3->Range(0, 0, 1, 1);

  TPad* pad1 = new TPad(Form("%s_1", name), "", 0, 0.5, 0.98, 0.98);
  pad1->Draw();

  TPad* pad2 = new TPad(Form("%s_2", name), "", 0, 0.02, 0.98, 0.5);
  pad2->Draw();

  pad1->SetRightMargin(0.01);
  pad2->SetRightMargin(0.01);
  pad1->SetTopMargin(0.05);
  pad1->SetLeftMargin(0.13);
  pad2->SetLeftMargin(0.13);
  pad2->SetBottomMargin(0.22);
  
  // no border between them
  pad1->SetBottomMargin(0);
  pad2->SetTopMargin(0);

  pad1->cd();
  pad1->SetGridx();
  pad1->SetGridy();

  TLegend* legend = new TLegend(0.15, 0.65, 0.55, 0.90);
  legend->SetFillColor(0);
  legend->AddEntry(corr, "Corrected");
  legend->AddEntry(mc, "MC prediction");
  legend->SetTextSize(0.08);

  Prepare1DPlot(corr);

  TH2F* dummy = new TH2F("dummy", "", 100, corr->GetXaxis()->GetBinLowEdge(1), corr->GetXaxis()->GetBinUpEdge(corr->GetNbinsX()), 1000, 0, TMath::Max(corr->GetMaximum(), mc->GetMaximum()) * 1.1);
  dummy->SetStats(kFALSE);
  dummy->SetXTitle(corr->GetXaxis()->GetTitle());
  dummy->SetYTitle(corr->GetYaxis()->GetTitle());
  dummy->GetYaxis()->SetTitleOffset(1);
  Prepare1DPlot(dummy);

  dummy->GetXaxis()->SetLabelSize(0.08);
  dummy->GetYaxis()->SetLabelSize(0.08);
  dummy->GetXaxis()->SetTitleSize(0.08);
  dummy->GetYaxis()->SetTitleSize(0.08);
  dummy->GetYaxis()->SetTitleOffset(0.8);
  dummy->DrawCopy();

  corr->Draw("SAME");
  mc->Draw("SAME");

  legend->Draw();

  pad2->cd();
  //pad2->SetBottomMargin(0.15);
  //pad2->SetGridx();
  //pad2->SetGridy();

  TH1* ratio = (TH1*) mc->Clone("ratio");
  ratio->Divide(corr);

  Float_t minR = TMath::Min(0.91, ratio->GetMinimum() * 0.95);
  Float_t maxR = TMath::Max(1.09, ratio->GetMaximum() * 1.05);
  
  minR = 0.8;
  maxR = 1.2;

  TH1F dummy3("dummy3", ";;Ratio: MC / corr", 100, corr->GetXaxis()->GetBinLowEdge(1), corr->GetXaxis()->GetBinUpEdge(corr->GetNbinsX()));
  dummy3.SetXTitle(corr->GetXaxis()->GetTitle());
  Prepare1DPlot(&dummy3);
  dummy3.SetStats(kFALSE);
  for (Int_t i=1; i<=100; ++i)
  	dummy3.SetBinContent(i, 1);
  dummy3.GetYaxis()->SetRangeUser(minR, maxR);
  dummy3.SetLineWidth(2);
  dummy3.GetXaxis()->SetLabelSize(0.08);
  dummy3.GetYaxis()->SetLabelSize(0.08);
  dummy3.GetXaxis()->SetTitleSize(0.08);
  dummy3.GetYaxis()->SetTitleSize(0.08);
  dummy3.GetYaxis()->SetTitleOffset(0.8);
  dummy3.DrawCopy();

  ratio->Draw("SAME");

  canvas3->Modified();
}

void loadlibs()
{
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libJETAN");
  gSystem->Load("libPWG4JetTasks");
}

void CompareBias(const char* mcFile = "PWG4_JetTasksOutput.root", Int_t region)
{
  loadlibs();

  TFile::Open(mcFile);
  list = (TList*) gFile->Get(objectName);
  AliUEHistograms* h = (AliUEHistograms*) list->FindObject("AliUEHistograms");
  h->SetEtaRange(-0.79, 0.79);
  h->SetPtRange(0.51, 100);
  h->SetCombineMinMax(kTRUE);
  
  biasFromData = (TH1*) h->GetNumberDensitypT()->GetBias(AliUEHist::kCFStepBiasStudy, AliUEHist::kCFStepReconstructed, region, "z")->Clone("biasFromData");
  biasFromData2 = (TH1*) h->GetNumberDensitypT()->GetBias(AliUEHist::kCFStepBiasStudy2, AliUEHist::kCFStepReconstructed, region, "z")->Clone("biasFromData2");
  //biasFromData = (TH1*) h->GetNumberDensitypT()->GetBias(AliUEHist::kCFStepReconstructed, AliUEHist::kCFStepRealLeading, region, "z")->Clone("biasFromData");
  biasFromMC   = (TH1*) h->GetNumberDensitypT()->GetBias(AliUEHist::kCFStepReconstructed, AliUEHist::kCFStepTracked, region, "z")->Clone("biasFromMC");
  //biasFromMC   = (TH1*) h->GetNumberDensitypT()->GetBias(AliUEHist::kCFStepReconstructed, AliUEHist::kCFStepRealLeading, region, "z")->Clone("biasFromMC");
  
/*  biasFromData->Draw();
  biasFromMC->SetLineColor(2);
  biasFromMC->Draw("SAME");
  return;*/
  
  DrawRatio(biasFromData, biasFromMC, "bias: data vs MC");
  DrawRatio(biasFromData, biasFromData2, "bias: data vs data two step");
}
  
void CompareBiasWithData(const char* mcFile = "PWG4_JetTasksOutput.root", const char* dataFile = "esd.root")
{
  loadlibs();

  file1 = TFile::Open(mcFile);
  list = (TList*) file1->Get(objectName);
  AliUEHistograms* h = (AliUEHistograms*) list->FindObject("AliUEHistograms");
  h->SetEtaRange(-0.79, 0.79);
  h->SetPtRange(0.51, 100);
  
  biasFromMC   = (TH1*) h->GetNumberDensitypT()->GetBias(AliUEHist::kCFStepReconstructed, AliUEHist::kCFStepTracked, "z")->Clone("biasFromMC");
  
  file2 = TFile::Open(dataFile);
  list = (TList*) file2->Get(objectName);
  AliUEHistograms* h2 = (AliUEHistograms*) list->FindObject("AliUEHistograms");
  h2->SetEtaRange(-0.79, 0.79);
  h2->SetPtRange(0.51, 100);
  
  biasFromData = (TH1*) h2->GetNumberDensitypT()->GetBias(AliUEHist::kCFStepBiasStudy, AliUEHist::kCFStepReconstructed, "z")->Clone("biasFromData");
  biasFromData2 = (TH1*) h2->GetNumberDensitypT()->GetBias(AliUEHist::kCFStepBiasStudy2, AliUEHist::kCFStepReconstructed, "z")->Clone("biasFromData2");
  
  DrawRatio(biasFromData, biasFromMC, "bias: data vs MC");
  DrawRatio(biasFromData, biasFromData2, "bias: data vs data two step");
}

void Compare(const char* fileName1, const char* fileName2, Int_t step1, Int_t step2)
{
  loadlibs();

  file1 = TFile::Open(fileName1);
  list = (TList*) file1->Get(objectName);
  AliUEHistograms* h = (AliUEHistograms*) list->FindObject("AliUEHistograms");
  h->SetEtaRange(-0.79, 0.79);
  h->SetPtRange(0.51, 100);

  file2 = TFile::Open(fileName2);
  list = (TList*) file2->Get(objectName);
  AliUEHistograms* h2 = (AliUEHistograms*) list->FindObject("AliUEHistograms");
  h2->SetEtaRange(-0.79, 0.79);
  h2->SetPtRange(0.51, 100);

  Int_t region = 0;
    
  TH1* hist1 = h->GetNumberDensitypT()->GetUEHist(step1, region);
  TH1* hist2 = h2->GetNumberDensitypT()->GetUEHist(step2, region);

  DrawRatio(hist1, hist2, "compare");
}  

void CompareRaw(const char* fileName1, const char* fileName2)
{
  loadlibs();

  Compare(fileName1, fileName2, AliUEHist::kCFStepReconstructed, AliUEHist::kCFStepReconstructed);
}

void ProfileMultiplicity(const char* fileName = "PWG4_JetTasksOutput.root")
{
  loadlibs();

  TFile::Open(fileName);
  list = (TList*) gFile->Get(objectName);
  AliUEHistograms* h = (AliUEHistograms*) list->FindObject("AliUEHistograms");

  new TCanvas;
  h->GetCorrelationMultiplicity()->Draw("colz");
  gPad->SetLogz();

  new TCanvas;
  h->GetCorrelationMultiplicity()->ProfileX()->DrawCopy()->Fit("pol1", "", "", 1, 10);
}

void SetupRanges(void* obj)
{
  ((AliUEHistograms*) obj)->SetEtaRange(-0.79, 0.79);
  ((AliUEHistograms*) obj)->SetPtRange(0.51, 100);
  ((AliUEHistograms*) obj)->SetCombineMinMax(kTRUE);
}

void DrawRatios(const char* name, void* correctedVoid, void* comparisonVoid, Int_t compareStep = -1, Int_t compareRegion = 2)
{
  AliUEHist* corrected = (AliUEHist*) correctedVoid;
  AliUEHist* comparison = (AliUEHist*) comparisonVoid;

  Int_t beginStep = AliUEHist::kCFStepAll;
  Int_t endStep = AliUEHist::kCFStepReconstructed;
  
  if (compareStep != -1)
  {
    beginStep = compareStep;
    endStep = compareStep;
  }
  
  Int_t beginRegion = 0;
  Int_t endRegion = 2;
  
  if (compareRegion != -1)
  {
    beginRegion = compareRegion;
    endRegion = compareRegion;
  }

  for (Int_t step=beginStep; step<=endStep; step++)
  {
    for (Int_t region=beginRegion; region <= endRegion; region++)
    {
      TH1* corrHist = corrected->GetUEHist(step, region);
      TH1* mcHist   = comparison->GetUEHist(step, region);
      
      DrawRatio(corrHist, mcHist, Form("%s: step %d %s %s", name, step, corrected->GetStepTitle(step), corrected->GetRegionTitle(region)));
    }
  }
}

void DrawRatios(void* correctedVoid, void* comparisonVoid, Int_t compareStep = -1, Int_t compareRegion = 2)
{
  AliUEHistograms* corrected = (AliUEHistograms*) correctedVoid;
  AliUEHistograms* comparison = (AliUEHistograms*) comparisonVoid;

  DrawRatios("Number density", corrected->GetNumberDensitypT(), comparison->GetNumberDensitypT(), compareStep, compareRegion);
  if (compareStep != -1)
  {
    DrawRatios("Pt sum", corrected->GetSumpT(), comparison->GetSumpT(), compareStep, compareRegion);
    //DrawRatios("Phi correlation", corrected->GetNumberDensityPhi(), comparison->GetNumberDensityPhi(), compareStep, compareRegion);
  }
}

void correctMC(const char* fileNameCorrections, const char* fileNameESD = 0, Int_t compareStep = -1, Int_t compareRegion = 2)
{
  // corrects the reconstructed step in fileNameESD with fileNameCorr
  // if fileNameESD is 0 data from fileNameCorr is taken
  // afterwards the corrected distributions are compared with the MC stored in fileNameESD
  
  loadlibs();
  
  TFile::Open(fileNameCorrections);
  list = (TList*) gFile->Get(objectName);
  AliUEHistograms* corr = (AliUEHistograms*) list->FindObject("AliUEHistograms");
  
  if (fileNameESD)
  {
    file2 = TFile::Open(fileNameESD);
    list = (TList*) file2->Get(objectName);
  }
  AliUEHistograms* testSample = (AliUEHistograms*) list->FindObject("AliUEHistograms");
    
  // copy to esd object
  AliUEHistograms* esd = new AliUEHistograms;
  esd->CopyReconstructedData(testSample);
  
  SetupRanges(corr);
  SetupRanges(testSample);
  SetupRanges(esd);
  
  esd->Correct(corr);
  
  DrawRatios(testSample, esd, compareStep, compareRegion);
}

// function to compare only final step for all regions and distributions

void correctData(const char* fileNameCorrections, const char* fileNameESD, Int_t compareStep = -1, Int_t compareRegion = 2)
{
  // corrects fileNameESD with fileNameCorrections and compares the two
  
  loadlibs();
  
  TFile::Open(fileNameCorrections);
  list = (TList*) gFile->Get(objectName);
  AliUEHistograms* corr = (AliUEHistograms*) list->FindObject("AliUEHistograms");
  
  file2 = TFile::Open(fileNameESD);
  list = (TList*) file2->Get(objectName);
  AliUEHistograms* esd = (AliUEHistograms*) list->FindObject("AliUEHistograms");
  
  SetupRanges(corr);
  SetupRanges(esd);
  
  esd->Correct(corr);
  
  file3 = TFile::Open("corrected.root", "RECREATE");
  esd->Write();
  file3->Close();
  
  DrawRatios(corr, esd, compareStep, compareRegion);
}
