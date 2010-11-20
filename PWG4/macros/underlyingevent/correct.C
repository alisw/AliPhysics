const char* objectName = "PWG4_LeadingTrackUE/histosLeadingTrackUE";
Float_t gLeadingpTMin = 0.51;
Int_t gUEHist = 0;
Bool_t gCache = 0;
void* gFirst = 0;
void* gSecond = 0;
Float_t gForceRange = -1;
Int_t gEnergy = 900;
Int_t gRegion = 0;

//const char* objectName = "PWG4_LeadingTrackUE/histosLeadingTrackUE_filterbit32";

void SetRanges(TAxis* axis)
{
  if (strcmp(axis->GetTitle(), "leading p_{T} (GeV/c)") == 0)
    axis->SetRangeUser(0, (gEnergy == 900) ? 10 : 25);

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

TH1* GetSystematicUncertainty(TH1* corr, TH1* trackHist)
{
  if (!trackHist)
  {
    systError = (TH1*) corr->Clone("systError");
    systError->Reset();
  }
  else  // for dphi evaluation
    systError = new TH1F("systError", "", 100, 0, 50);
  
  Float_t constantUnc = 0;
  
  // particle composition
  constantUnc += 0.8 ** 2;
  
  // tpc efficiency
  if (gEnergy == 900 && gLeadingpTMin < 1.0)
    constantUnc += 1.0 ** 2;
  else if (gEnergy == 900 && gLeadingpTMin < 1.5)
    constantUnc += 0.5 ** 2;
  if (gEnergy == 7000 && gLeadingpTMin < 1.0)
    constantUnc += 1.0 ** 2;
  else if (gEnergy == 7000 && gLeadingpTMin < 1.5)
    constantUnc += 0.6 ** 2;
  
  // track cuts
  if (gEnergy == 900 && gLeadingpTMin < 1.0)
    constantUnc += 2.5 ** 2;
  else if (gEnergy == 900 && gLeadingpTMin < 1.5)
    constantUnc += 2.0 ** 2;
  if (gEnergy == 7000)
    constantUnc += 3.0 ** 2;

  // difference corrected with pythia and phojet
  if (gEnergy == 900 && gLeadingpTMin < 1.0)
    constantUnc += 0.6 ** 2;
  else if (gEnergy == 900 && gLeadingpTMin < 1.5)
    constantUnc += 0.8 ** 2;
  
  if (gEnergy == 7000 && gLeadingpTMin < 1.0)
  {
    if (gUEHist == 0)
      constantUnc += 0.6 ** 2;
    if (gUEHist == 1)
      constantUnc += 0.8 ** 2;
    if (gUEHist == 2)
      constantUnc += 1.0 ** 2;
  }
  else if (gEnergy == 7000 && gLeadingpTMin < 1.5)
    constantUnc += 1.0 ** 2;
    
  for (Int_t bin=1; bin<=systError->GetNbinsX(); bin++)
    systError->SetBinContent(bin, constantUnc);

  // mis id bias
  if (gUEHist == 0 || gUEHist == 2)
    systError->Fill(0.75, 4.0 ** 2);
  if (gUEHist == 1)
    systError->Fill(0.75, 5.0 ** 2);

  if (gEnergy == 900)
  {
    if (gLeadingpTMin < 1.0)
      systError->Fill(1.25, 1.0 ** 2);
    else if (gLeadingpTMin < 1.5)
      systError->Fill(1.25, 2.0 ** 2);
  }
  
  // non-closure in MC
  if (gEnergy == 900)
    for (Int_t bin=1; bin<=systError->GetNbinsX(); bin++)
      systError->Fill(systError->GetXaxis()->GetBinCenter(bin), 1.0 ** 2);
      
  if (gEnergy == 7000)
  {
    if (gUEHist == 0 && gUEHist == 1)
      systError->Fill(0.75, 2.0 ** 2);
    if (gUEHist == 2)
      systError->Fill(0.75, 1.2 ** 2);
  }
  
  // vertex efficiency
  systError->Fill(0.75, 1.0 ** 2);

  // strangeness
  for (Int_t bin=1; bin<=systError->GetNbinsX(); bin++)
  {
    if (gEnergy == 900)
      systError->Fill(systError->GetXaxis()->GetBinCenter(bin), 0.5 ** 2);
    if (gEnergy == 7000 && systError->GetXaxis()->GetBinCenter(bin) < 1.5)
      systError->Fill(systError->GetXaxis()->GetBinCenter(bin), 2.0 ** 2);
    else if (gEnergy == 7000)
      systError->Fill(systError->GetXaxis()->GetBinCenter(bin), 1.0 ** 2);
  }  
    
  for (Int_t bin=1; bin<=systError->GetNbinsX(); bin++)
    systError->SetBinContent(bin, TMath::Sqrt(systError->GetBinContent(bin)));
  
  if (trackHist)
  {
    //new TCanvas; trackHist->Draw();
    //new TCanvas; systError->DrawCopy("");
    
    Float_t uncFlat = 0;
    for (Int_t i=1; i<=trackHist->GetNbinsX(); i++)
      uncFlat += trackHist->GetBinContent(i) * systError->GetBinContent(systError->FindBin(trackHist->GetBinCenter(i)));
    uncFlat /= trackHist->Integral();
    
    systError = (TH1F*) corr->Clone("systError");
    systError->Reset();
  
    for (Int_t i=1; i<=systError->GetNbinsX(); i++)
      systError->SetBinContent(i, uncFlat);
      
    //new TCanvas; systError->DrawCopy("");
  }
  
  systError->SetFillColor(kGray);
  systError->SetFillStyle(1001);
  systError->SetMarkerStyle(0);
  systError->SetLineColor(0);
  
  return systError;
}

void DrawRatio(TH1* corr, TH1* mc, const char* name, TH1* syst = 0)
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

  if (gUEHist != 2)
    TLegend* legend = new TLegend(0.15, 0.65, 0.55, 0.90);
  else
    TLegend* legend = new TLegend(0.55, 0.65, 0.95, 0.90);

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
  
  if (gForceRange > 0)
    dummy->GetYaxis()->SetRangeUser(0, gForceRange);
  
  dummy->DrawCopy();
  
  if (gUEHist != 2)
  {
    const char* regionStr[] = { "Toward Region", "Away Region", "Transverse Region" };
    latex = new TLatex(0.65, 0.1, regionStr[gRegion]);
    latex->SetTextSize(0.075);
    latex->SetNDC();
    latex->Draw();
  }
  
  if (syst)
  {
    systError = (TH1*) syst->Clone("corrSystError");
    for (Int_t bin=1; bin<=systError->GetNbinsX(); bin++)
    {
      systError->SetBinError(bin, corr->GetBinContent(bin) * syst->GetBinContent(bin) / 100);
      systError->SetBinContent(bin, corr->GetBinContent(bin));
    }
    
    systError->Draw("E2 ][ SAME");
  }
  
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
  
  minR = 0.6;
  maxR = 1.4;

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
  
  if (syst)
  {
    // for the ratio add in quadrature
    for (Int_t bin=1; bin<=syst->GetNbinsX(); bin++)
    {
      if (corr->GetBinError(bin) > 0)
        syst->SetBinError(bin, TMath::Sqrt(TMath::Power(syst->GetBinContent(bin) / 100, 2) + TMath::Power(corr->GetBinError(bin) / corr->GetBinContent(bin), 2)));
      else
        syst->SetBinError(bin, 0);
      syst->SetBinContent(bin, 1);
    }
    
    syst->Draw("E2 ][ SAME");
    dummy3.DrawCopy("SAME");
  }

  ratio->Draw("SAME");
  
  ratio->Fit("pol0", "N");

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

void CompareBias(const char* mcFile = "PWG4_JetTasksOutput.root", Int_t region, Int_t ueHist)
{
  loadlibs();

  TFile::Open(mcFile);
  list = (TList*) gFile->Get(objectName);
  AliUEHistograms* h = (AliUEHistograms*) list->FindObject("AliUEHistograms");
  SetupRanges(h);
  
  const char* axis = "z";
  Float_t ptLeadMin = 0;
  Float_t ptLeadMax = -1;
  
  if (ueHist == 2)
  {
    ptLeadMin = 0.51 + 0;
    ptLeadMax = 1.99 + 0;
    axis = "y";
  }
  
  biasFromData = (TH1*) h->GetUEHist(ueHist)->GetBias(AliUEHist::kCFStepBiasStudy, AliUEHist::kCFStepReconstructed, region, axis, ptLeadMin, ptLeadMax)->Clone("biasFromData");
  biasFromData2 = (TH1*) h->GetUEHist(ueHist)->GetBias(AliUEHist::kCFStepBiasStudy2, AliUEHist::kCFStepReconstructed, region, axis, ptLeadMin, ptLeadMax)->Clone("biasFromData2");
  //biasFromData = (TH1*) h->GetUEHist(ueHist)->GetBias(AliUEHist::kCFStepReconstructed, AliUEHist::kCFStepRealLeading, region, axis)->Clone("biasFromData");
  biasFromMC   = (TH1*) h->GetUEHist(ueHist)->GetBias(AliUEHist::kCFStepReconstructed, AliUEHist::kCFStepTracked, region, axis, ptLeadMin, ptLeadMax)->Clone("biasFromMC");
  //biasFromMC   = (TH1*) h->GetUEHist(ueHist)->GetBias(AliUEHist::kCFStepReconstructed, AliUEHist::kCFStepRealLeading, region, axis)->Clone("biasFromMC");
  
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
  SetupRanges(h);
  
  biasFromMC   = (TH1*) h->GetNumberDensitypT()->GetBias(AliUEHist::kCFStepReconstructed, AliUEHist::kCFStepTracked, "z")->Clone("biasFromMC");
  
  file2 = TFile::Open(dataFile);
  list = (TList*) file2->Get(objectName);
  AliUEHistograms* h2 = (AliUEHistograms*) list->FindObject("AliUEHistograms");
  SetupRanges(h2);
  
  biasFromData = (TH1*) h2->GetNumberDensitypT()->GetBias(AliUEHist::kCFStepBiasStudy, AliUEHist::kCFStepReconstructed, "z")->Clone("biasFromData");
  biasFromData2 = (TH1*) h2->GetNumberDensitypT()->GetBias(AliUEHist::kCFStepBiasStudy2, AliUEHist::kCFStepReconstructed, "z")->Clone("biasFromData2");
  
  DrawRatio(biasFromData, biasFromMC, "bias: data vs MC");
  DrawRatio(biasFromData, biasFromData2, "bias: data vs data two step");
}

void Compare(const char* fileName1, const char* fileName2, Int_t id, Int_t step1, Int_t step2, Int_t region, Float_t ptLeadMin = -1, Float_t ptLeadMax = -1)
{
  loadlibs();
  
  if (!gCache || !gFirst)
  {
    file1 = TFile::Open(fileName1);
    list = (TList*) file1->Get(objectName);
    AliUEHistograms* h = (AliUEHistograms*) list->FindObject("AliUEHistograms");
    SetupRanges(h);
  
    file2 = TFile::Open(fileName2);
    list = (TList*) file2->Get(objectName);
    AliUEHistograms* h2 = (AliUEHistograms*) list->FindObject("AliUEHistograms");
    SetupRanges(h2);
  }
  
  if (gCache)
  {
    if (!gFirst)
    {
      gFirst = h;
      gSecond = h2;
    }
    else
    {
      AliUEHistograms* h = (AliUEHistograms*) gFirst;
      AliUEHistograms* h2 = (AliUEHistograms*) gSecond;
    }
  }
    
  TH1* hist1 = h->GetUEHist(id)->GetUEHist(step1, region, ptLeadMin, ptLeadMax);
  TH1* hist2 = h2->GetUEHist(id)->GetUEHist(step2, region, ptLeadMin, ptLeadMax);

  TH1* trackHist = 0;
  if (id == 2)
  {
    trackHist = h->GetUEHist(id)->GetTrackHist(region)->ShowProjection(2, 0);
    // only keep bins under consideration
    for (Int_t bin=1; bin<=trackHist->GetNbinsX(); bin++)
      if (bin < trackHist->FindBin(ptLeadMin) || bin > trackHist->FindBin(ptLeadMax))
        trackHist->SetBinContent(bin, 0);
  }

  // systematic uncertainty
  TH1* syst = GetSystematicUncertainty(hist1, trackHist);
  
  DrawRatio(hist1, hist2, Form("%s_%d_%d_%d_%.2f_%.2f", TString(fileName1).Tokenize(".")->First()->GetName(), id, step1, region, ptLeadMin, ptLeadMax), syst);
}  

void CompareEventHist(const char* fileName1, const char* fileName2, Int_t id, Int_t step, Int_t var)
{
  loadlibs();

  file1 = TFile::Open(fileName1);
  list = (TList*) file1->Get(objectName);
  AliUEHistograms* h = (AliUEHistograms*) list->FindObject("AliUEHistograms");
  SetupRanges(h);

  file2 = TFile::Open(fileName2);
  list = (TList*) file2->Get(objectName);
  AliUEHistograms* h2 = (AliUEHistograms*) list->FindObject("AliUEHistograms");
  SetupRanges(h2);

  TH1* hist1 = h->GetUEHist(id)->GetEventHist()->ShowProjection(var, step);
  TH1* hist2 = h2->GetUEHist(id)->GetEventHist()->ShowProjection(var, step);
 
  DrawRatio(hist1, hist2, "compare");
}

void CompareStep(const char* fileName1, const char* fileName2, Int_t id, Int_t step, Int_t region, Float_t ptLeadMin = -1, Float_t ptLeadMax = -1)
{
  // fileName1 is labelled Corrected in the plot

  loadlibs();

  gUEHist = id;
  gRegion = region;
  Compare(fileName1, fileName2, id, step, step, region, ptLeadMin, ptLeadMax);
}

void DrawStep(const char* fileName, Int_t id, Int_t step, Int_t region, Float_t ptLeadMin = -1, Float_t ptLeadMax = -1)
{
  loadlibs();

  file1 = TFile::Open(fileName);
  list = (TList*) file1->Get(objectName);
  AliUEHistograms* h = (AliUEHistograms*) list->FindObject("AliUEHistograms");
  SetupRanges(h);

  new TCanvas;
  h->GetUEHist(id)->GetUEHist(step, region, ptLeadMin, ptLeadMax)->Draw();
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
  ((AliUEHistograms*) obj)->SetPtRange(gLeadingpTMin, 100);
  ((AliUEHistograms*) obj)->SetCombineMinMax(kTRUE);
}

void DrawRatios(const char* name, void* correctedVoid, void* comparisonVoid, Int_t compareStep = -1, Int_t compareRegion = 2, Float_t ptLeadMin = -1, Float_t ptLeadMax = -1)
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
      Printf("%f %f", ptLeadMin, ptLeadMax);
      TH1* corrHist = corrected->GetUEHist(step, region, ptLeadMin, ptLeadMax);
      TH1* mcHist   = comparison->GetUEHist(step, region, ptLeadMin, ptLeadMax);
      
      DrawRatio(corrHist, mcHist, TString(Form("%s: step %d %s %s", name, step, corrected->GetStepTitle(step), corrected->GetRegionTitle(region))));
    }
  }
}

void DrawRatios(void* correctedVoid, void* comparisonVoid, Int_t compareStep = -1, Int_t compareRegion = 2, Int_t compareUEHist = 0)
{
  AliUEHistograms* corrected = (AliUEHistograms*) correctedVoid;
  AliUEHistograms* comparison = (AliUEHistograms*) comparisonVoid;

  if (1 && compareUEHist == 2)
  {
    for (Float_t ptLeadMin = 0.51; ptLeadMin < 10; ptLeadMin += 1.5)
      DrawRatios(TString(Form("UE %d pT %f", compareUEHist, ptLeadMin)), corrected->GetUEHist(compareUEHist), comparison->GetUEHist(compareUEHist), compareStep, compareRegion, ptLeadMin, ptLeadMin + 0.48);      
    return;
  }

  if (compareUEHist == -1)
  {
    for (Int_t i=0; i<2; i++)
      DrawRatios(TString(Form("UE %d", i)), corrected->GetUEHist(i), comparison->GetUEHist(i), compareStep, compareRegion);
  }
  else
    DrawRatios(TString(Form("UE %d", compareUEHist)), corrected->GetUEHist(compareUEHist), comparison->GetUEHist(compareUEHist), compareStep, compareRegion);
}

void CompareEventsTracks(void* corrVoid, void* mcVoid, Int_t compareStep, Int_t compareRegion, Int_t compareUEHist = 0)
{
  AliUEHist* corr = ((AliUEHistograms*) corrVoid)->GetUEHist(compareUEHist);
  AliUEHist* mc   = ((AliUEHistograms*) mcVoid)->GetUEHist(compareUEHist);

  Float_t ptLeadMin = 0;
  Float_t ptLeadMax = -1;
  Int_t axis = 2;
  
  if (compareUEHist == 2)
  {
    ptLeadMin = 1.01;
    ptLeadMax = 1.49;
    axis = 4;
  }

  TH1* corrHist = corr->GetUEHist(compareStep, compareRegion, ptLeadMin, ptLeadMax);
  TH1* mcHist   = mc  ->GetUEHist(compareStep, compareRegion, ptLeadMin, ptLeadMax);
  
  DrawRatio(corrHist, mcHist, Form("check"));
  
  corr->SetBinLimits(corr->GetTrackHist(compareRegion)->GetGrid(compareStep));
  mc->SetBinLimits(corr->GetTrackHist(compareRegion)->GetGrid(compareStep));

  corrHist =  corr->GetTrackHist(compareRegion)->GetGrid(compareStep)->Project(axis);
  mcHist   =  mc  ->GetTrackHist(compareRegion)->GetGrid(compareStep)->Project(axis);
  DrawRatio(corrHist, mcHist, Form("check2"));
  
  corrHist =  corr->GetEventHist()->GetGrid(compareStep)->Project(0);
  mcHist   =  mc  ->GetEventHist()->GetGrid(compareStep)->Project(0);
  DrawRatio(corrHist, mcHist, Form("check3"));
}

void correctMC(const char* fileNameCorrections, const char* fileNameESD = 0, Int_t compareStep = -1, Int_t compareRegion = 2, Int_t compareUEHist = 0)
{
  // corrects the reconstructed step in fileNameESD with fileNameCorr
  // if fileNameESD is 0 data from fileNameCorr is taken
  // afterwards the corrected distributions are compared with the MC stored in fileNameESD
  
  loadlibs();
  
  TFile::Open(fileNameCorrections);
  list = (TList*) gFile->Get(objectName);
  AliUEHistograms* corr = (AliUEHistograms*) list->FindObject("AliUEHistograms");
  SetupRanges(corr);
  corr->ExtendTrackingEfficiency();
  
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
  esd->GetUEHist(2)->AdditionalDPhiCorrection(0);
  
  DrawRatios(esd, testSample, compareStep, compareRegion, compareUEHist);
  
  //CompareEventsTracks(esd, testSample, compareStep, compareRegion, compareUEHist);
}

// function to compare only final step for all regions and distributions

void correctData(const char* fileNameCorrections, const char* fileNameESD, const char* contEnhancement, Float_t contEncUpTo = 1.0, Int_t compareStep = -1, Int_t compareRegion = 2, Int_t compareUEHist = 0)
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
  
  corr->ExtendTrackingEfficiency();
  
  TFile::Open(contEnhancement);
  contEncHist = (TH1*) gFile->Get("histo");
  contEncHistFullRange = (TH1*) corr->GetUEHist(0)->GetTrackingEfficiency(1)->Clone("contEncHistFullRange");
  
  contEncHistFullRange->Reset();
  for (Int_t i=1; i<=contEncHistFullRange->GetNbinsX(); i++)
  {
    contEncHistFullRange->SetBinContent(i, 1);
    if (i <= contEncHist->GetNbinsX() && contEncHist->GetXaxis()->GetBinCenter(i) < contEncUpTo && contEncHist->GetBinContent(i) > 0)
      contEncHistFullRange->SetBinContent(i, contEncHist->GetBinContent(i));
  }
  corr->SetContaminationEnhancement((TH1F*) contEncHistFullRange);
  
  esd->Correct(corr);
  esd->GetUEHist(2)->AdditionalDPhiCorrection(0);
  
  file3 = TFile::Open("corrected.root", "RECREATE");
  file3->mkdir("PWG4_LeadingTrackUE");
  file3->cd("PWG4_LeadingTrackUE");
  list->Write(0, TObject::kSingleKey);
  file3->Close();
  
  DrawRatios(esd, corr, compareStep, compareRegion, compareUEHist);
}

void ITSTPCEfficiency(const char* fileNameData, Int_t id, Int_t itsTPC = 0)
{
  // its = 0; tpc = 1

  // uncertainty from dN/dpT paper
  Double_t pTBins[] =  {0.0, 0.1, 0.15,  0.2,  0.25,  0.3,   0.35,  0.4,   0.45,  0.5,   0.6,   0.7,   0.8,   0.9,   1.0,   1.5,   2.0,   2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 100.0};
  Float_t effITS[] =   {0.,  0.,  0.995, 0.98, 0.986, 0.996, 1,     1,     1,     1,     1,     1,     1,     1,     1,     1,     1,    };  // the last three are the same because i don't have entries
  Float_t effTPC[] =   {0.,  0,   1.042, 1.026,1.021, 1.018, 1.015, 1.015, 1.012, 1.012, 1.007, 1.0075,1.006, 1.006, 1.004, 1.004, 1.009 }; // the last bins put as if they were the same

  TH1F* effHist = new TH1F("effHist", "effHist", 39, pTBins);
  for (Int_t i=0; i<39; i++)
  {
    Int_t bin = i;
    if (i > 16)
      bin = 16;
    effHist->SetBinContent(i+1, (itsTPC == 0) ? effITS[bin] : effTPC[bin]);
  }

  new TCanvas; effHist->Draw();

  EffectOfModifiedTrackingEfficiency(fileNameData, id, 2, effHist, 1, -1, (itsTPC == 0) ? "ITS" : "TPC");
} 


void EffectOfModifiedTrackingEfficiency(const char* fileNameData, Int_t id, Int_t region, TH1* trackingEff, Int_t axis1, Int_t axis2 = -1, const char* name = "EffectOfModifiedTrackingEfficiency")
{
  // trackingEff should contain the change in tracking efficiency, i.e. between before and after in the eta-pT plane

  loadlibs();
  
  TFile::Open(fileNameData);
  list = (TList*) gFile->Get(objectName);
  AliUEHistograms* corrected = (AliUEHistograms*) list->FindObject("AliUEHistograms");
  SetupRanges(corrected);
  
  AliUEHist* ueHist = corrected->GetUEHist(id);
  
  Float_t ptLeadMin = -1;
  Float_t ptLeadMax = -1;
  if (id == 2)
  {
    // the uncertainty is flat in delta phi, so use this trick to get directly the uncertainty as function of leading pT
    //ptLeadMin = 1.01;
    //ptLeadMax = 1.99;
  }
    
  // histogram before
  TH1* before = ueHist->GetUEHist(AliUEHist::kCFStepAll, region, ptLeadMin, ptLeadMax);

  // copy histogram
  // the CFStepTriggered step is overwritten here and cannot be used for comparison afterwards anymore
  ueHist->CorrectTracks(AliUEHist::kCFStepAll, AliUEHist::kCFStepTriggered, (TH1*) 0, 0, -1);

  // reapply tracking efficiency
  ueHist->CorrectTracks(AliUEHist::kCFStepTriggered, AliUEHist::kCFStepAll, trackingEff, axis1, axis2);

  // histogram after
  TH1* after = ueHist->GetUEHist(AliUEHist::kCFStepAll, region, ptLeadMin, ptLeadMax);
  
  DrawRatio(before, after, name);
  gPad->GetCanvas()->SaveAs(Form("%s.png", name));
}

void EffectOfTrackCuts(const char* fileNameData, Int_t id, const char* systFile)
{
  loadlibs();

  TFile::Open(fileNameData);
  list = (TList*) gFile->Get(objectName);
  AliUEHistograms* corrected = (AliUEHistograms*) list->FindObject("AliUEHistograms");
  effHist = (TH2D*) corrected->GetUEHist(0)->GetTrackingEfficiency()->Clone("effHist");

  file = TFile::Open(systFile);

  Int_t maxSyst = 3;
  const char* systNames[] = { "NClusTPC", "Chi2TPC", "SigmaDCA" };

  for (Int_t i=0; i<maxSyst; i++)
  {
    for (Int_t j=0; j<2; j++)
    {
      TString histName;
      histName.Form("%s_syst_%s", systNames[i], (j == 0) ? "up" : "down");
      systEffect = (TH2*) file->Get(histName);

      // rebin
      effHist->Reset();
      for (Int_t x=1; x <= effHist->GetNbinsX(); x++)
        for (Int_t y=1; y <= effHist->GetNbinsY(); y++)
          effHist->SetBinContent(x, y, 1);

      for (Int_t x=1; x <= systEffect->GetNbinsX(); x++)
        for (Int_t y=1; y <= systEffect->GetNbinsY(); y++)
          if (systEffect->GetBinContent(x, y) != 0)
            effHist->SetBinContent(effHist->GetXaxis()->FindBin(systEffect->GetYaxis()->GetBinCenter(y)), effHist->GetYaxis()->FindBin(systEffect->GetXaxis()->GetBinCenter(x)), systEffect->GetBinContent(x, y));
           
   
      //new TCanvas; systEffect->Draw("COLZ"); new TCanvas; effHist->Draw("COLZ");
 
      EffectOfModifiedTrackingEfficiency(fileNameData, id, (id == 2) ? 0 : 2, effHist, 0, 1, histName);

      //return;
    }
  } 
}

void ModifyComposition(const char* fileNameData, const char* fileNameCorrections, Int_t id, Bool_t verbose = kFALSE)
{
  loadlibs();
  
  TFile::Open(fileNameData);
  list = (TList*) gFile->Get(objectName);
  AliUEHistograms* corrected = (AliUEHistograms*) list->FindObject("AliUEHistograms");
  SetupRanges(corrected);
  
  file2 = TFile::Open(fileNameCorrections);
  list2 = (TList*) file2->Get(objectName);
  AliUEHistograms* corrections = (AliUEHistograms*) list2->FindObject("AliUEHistograms");
  SetupRanges(corrections);
  
  ueHistData        = (AliUEHist*) corrected->GetUEHist(id);
  ueHistCorrections = (AliUEHist*) corrections->GetUEHist(id);
  
  // copy histogram
  // the CFStepTriggered step is overwritten here and cannot be used for comparison afterwards anymore
  ueHistData->CorrectTracks(AliUEHist::kCFStepAll, AliUEHist::kCFStepTriggered, (TH1*) 0, 0);
  
  Int_t maxRegion = 3;
  Float_t ptLeadMin = -1;
  Float_t ptLeadMax = -1;
  if (id == 2)
  {
    maxRegion = 1;
    // the uncertainty is flat in delta phi, so use this trick to get directly the uncertainty as function of leading pT
    //ptLeadMin = 1.01;
    //ptLeadMax = 1.99;
  }
  
  // histogram before
  TH1* before[3];
  for (Int_t region=0; region<maxRegion; region++)
    before[region] = ueHistData->GetUEHist(AliUEHist::kCFStepAll, region, ptLeadMin, ptLeadMax);
  
  //defaultEff = ueHistCorrections->GetTrackingEfficiency();
  defaultEff = ueHistCorrections->GetTrackingCorrection();
  //defaultEffpT = ueHistCorrections->GetTrackingEfficiency(1);
  defaultEffpT = ueHistCorrections->GetTrackingCorrection(1);
  defaultContainer = ueHistCorrections->GetTrackHistEfficiency();
  
  c = new TCanvas;
  defaultEffpT->Draw("");
  
  Float_t largestDeviation[3];
  for (Int_t i=0; i<maxRegion; i++)
    largestDeviation[i] = 0;  
  
  for (Int_t i=0; i<7; i++)
  {
    // case 0: // default
    // case 1: // + 30% kaons
    // case 2: // - 30% kaons
    // case 3: // + 30% protons
    // case 4: // - 30% protons
    // case 5: // + 30% others
    // case 6: // - 30% others
    Int_t correctionIndex = (i+1) / 2 + 1; // bin 1 == protons, bin 2 == kaons, ...
    Double_t scaleFactor = (i % 2 == 1) ? 1.3 : 0.7;
    if (i == 0)
      scaleFactor = 1;
    
    newContainer = (AliCFContainer*) defaultContainer->Clone();
    
    // modify, change all steps
    for (Int_t j=0; j<newContainer->GetNStep(); j++)
    {
      THnSparse* grid = newContainer->GetGrid(j)->GetGrid();
      
      for (Int_t binIdx = 0; binIdx < grid->GetNbins(); binIdx++)
      {
        Int_t bins[5];
        Double_t value = grid->GetBinContent(binIdx, bins);
        Double_t error = grid->GetBinError(binIdx);
        
        if (bins[2] != correctionIndex)
          continue;
    
        value *= scaleFactor;
        error *= scaleFactor;
    
        grid->SetBinContent(bins, value);
        grid->SetBinError(bins, error);      
      }
    }
    
    // put in corrections
    ueHistCorrections->SetTrackHistEfficiency(newContainer);
    
    // ratio
    //modifiedEff = ueHistCorrections->GetTrackingEfficiency();
    modifiedEff = ueHistCorrections->GetTrackingCorrection();
    modifiedEff->Divide(modifiedEff, defaultEff);
    //modifiedEff->Draw("COLZ");
    
    c->cd();
    //modifiedEffpT = ueHistCorrections->GetTrackingEfficiency(1);
    modifiedEffpT = ueHistCorrections->GetTrackingCorrection(1);
    modifiedEffpT->SetLineColor(i+1);
    modifiedEffpT->Draw("SAME");
    
    // apply change in tracking efficiency
    ueHistData->CorrectTracks(AliUEHist::kCFStepTriggered, AliUEHist::kCFStepAll, modifiedEff, 0, 1);
  
    for (Int_t region=0; region<maxRegion; region++)
    {
      // histogram after
      TH1* after = ueHistData->GetUEHist(AliUEHist::kCFStepAll, region, ptLeadMin, ptLeadMax);
      
      if (verbose)
        DrawRatio(before[region], (TH1*) after->Clone(), Form("Region %d Composition %d", region, i));
      
      // ratio is flat, extract deviation
      after->Divide(before[region]);
      after->Fit("pol0", "0");
      Float_t deviation = 100.0 - 100.0 * after->GetFunction("pol0")->GetParameter(0);
      Printf("Deviation for region %d case %d is %.2f %%", region, i, deviation);
      
      if (TMath::Abs(deviation) > largestDeviation[region])
        largestDeviation[region] = TMath::Abs(deviation);
    }
    //return;
  }
  
  for (Int_t i=0; i<maxRegion; i++)
    Printf("Largest deviation in region %d is %f", i, largestDeviation[i]);
}    

void MergeList()
{
  loadlibs();

  ifstream in;
  in.open("list");

  TFileMerger m;

  TString line;
  while (in.good())
  {
    in >> line;

    if (line.Length() == 0)
      continue;

    TString fileName;
    fileName.Form("%s/%s/PWG4_JetTasksOutput.root", "maps", line.Data());
    Printf("%s", fileName.Data());
    
    m.AddFile(fileName);
  }
  
  m.SetFastMethod();
  m.OutputFile("merged.root");
  m.Merge();
}

void MergeList2(const char* listFile, Bool_t onlyPrintEvents = kFALSE)
{
  loadlibs();

  ifstream in;
  in.open(listFile);
  
  AliUEHistograms* final = 0;
  TList* finalList = 0;

  TString line;
  while (in.good())
  {
    in >> line;

    if (line.Length() == 0)
      continue;

    TString fileName;
    fileName.Form("%s/%s/PWG4_JetTasksOutput.root", "maps", line.Data());
    Printf("%s", fileName.Data());
    
    file = TFile::Open(fileName);
    if (!file)
      continue;
      
    list = (TList*) file->Get(objectName);
    
    if (!final)
    {
      final = (AliUEHistograms*) list->FindObject("AliUEHistograms");
      //final->GetEventCount()->Draw(); return;
      Printf("Events: %d", (Int_t) final->GetEventCount()->ProjectionX()->GetBinContent(4));
      finalList = list;
    }
    else
    {
      additional = (AliUEHistograms*) list->FindObject("AliUEHistograms");
      Printf("Events: %d", (Int_t) additional->GetEventCount()->ProjectionX()->GetBinContent(4));
      
      if (!onlyPrintEvents)
      {
        TList list2;
        list2.Add(additional);
        final->Merge(&list2);
      }
      delete additional;
      file->Close();
    }
  }
  
  if (onlyPrintEvents)
    return;
    
  Printf("Total events (at step 0): %d", (Int_t) final->GetEventCount()->ProjectionX()->GetBinContent(4));
  
  file3 = TFile::Open("merged.root", "RECREATE");
  file3->mkdir("PWG4_LeadingTrackUE");
  file3->cd("PWG4_LeadingTrackUE");
  finalList->Write(0, TObject::kSingleKey);
  file3->Close();
}

void PlotAll(const char* correctedFile, const char* mcFile)
{
  gCache = 1;
  
  if (gEnergy == 900)
  {
    Float_t range[] = { 1.5, 2 };
  }
  else
  {
    Float_t range[] = { 3, 10 };
  }
  
  for (Int_t id=0; id<3; id++)
  {
    if (id < 2)
      gForceRange = range[id];
    else
      gForceRange = -1;
      
    if (id < 2)
    {
      for (Int_t region=0; region<3; region++)
      {
        CompareStep(correctedFile, mcFile, id, 0, region);
        gPad->GetCanvas()->SaveAs(Form("%s_%s_%d_%d.png", TString(correctedFile).Tokenize(".")->First()->GetName(), TString(mcFile).Tokenize(".")->First()->GetName(), id, region));
      }
    }
    else
    {
      Float_t leadingPtArr[] = { 0.50, 2.0, 4.0, 6.0, 10.0, 20.0, 50.0 };
      for (Int_t leadingPtID=0; leadingPtID<6; leadingPtID++)
      {
        CompareStep(correctedFile, mcFile, id, 0, 0, leadingPtArr[leadingPtID] + 0.01, leadingPtArr[leadingPtID+1] - 0.01);
        gPad->GetCanvas()->SaveAs(Form("%s_%s_%d_%.2f_%.2f.png", TString(correctedFile).Tokenize(".")->First()->GetName(), TString(mcFile).Tokenize(".")->First()->GetName(), id, leadingPtArr[leadingPtID], leadingPtArr[leadingPtID+1]));
      }
    }
  }
}
 