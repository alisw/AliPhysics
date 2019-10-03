Float_t gpTMin = 0.51;
Float_t gpTMax = 49.99;
Int_t gUEHist = 0;
Bool_t gCache = 0;
void* gFirst = 0;
void* gSecond = 0;
Float_t gForceRange = -1;
Int_t gEnergy = 900;
Int_t gRegion = 0;
const char* gText = "";

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
  if (!hist)
    return;

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
  if (gEnergy == 900 && gpTMin < 1.0)
    constantUnc += 1.0 ** 2;
  else if (gEnergy == 900 && gpTMin < 1.5)
    constantUnc += 0.5 ** 2;
  if (gEnergy == 7000 && gpTMin < 1.0)
    constantUnc += 1.0 ** 2;
  else if (gEnergy == 7000 && gpTMin < 1.5)
    constantUnc += 0.6 ** 2;
  
  // track cuts
  if (gEnergy == 900 && gpTMin < 1.0)
    constantUnc += 2.5 ** 2;
  else if (gEnergy == 900 && gpTMin < 1.5)
    constantUnc += 2.0 ** 2;
  if (gEnergy == 7000)
    constantUnc += 3.0 ** 2;

  // difference corrected with pythia and phojet
  if (gEnergy == 900 && gpTMin < 1.0)
    constantUnc += 0.6 ** 2;
  else if (gEnergy == 900 && gpTMin < 1.5)
    constantUnc += 0.8 ** 2;
  
  if (gEnergy == 7000 && gpTMin < 1.0)
  {
    if (gUEHist == 0)
      constantUnc += 0.6 ** 2;
    if (gUEHist == 1)
      constantUnc += 0.8 ** 2;
    if (gUEHist == 2)
      constantUnc += 1.0 ** 2;
  }
  else if (gEnergy == 7000 && gpTMin < 1.5)
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
    if (gpTMin < 1.0)
      systError->Fill(1.25, 1.0 ** 2);
    else if (gpTMin < 1.5)
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
    if (trackHist->Integral() > 0)
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
  
  if (strlen(gText) > 0)
  {
    latex = new TLatex(0.2, 0.2, gText);
    latex->SetNDC();
    latex->SetTextSize(0.06);
    latex->Draw();
  }

  //legend->Draw();

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

void* GetUEHistogram(const char* fileName, TList** listRef = 0)
{
  file = TFile::Open(fileName);
  if (!file)
    return 0;
  
  list = (TList*) gFile->Get("PWG4_LeadingTrackUE/histosLeadingTrackUE");
  if (!list)
    list = (TList*) gFile->Get("PWG4_PhiCorrelations/histosPhiCorrelations");
    
  if (!list)
    return 0;
    
  if (listRef)
    *listRef = list;
    
  if (list->FindObject("AliUEHistograms"))
    return list->FindObject("AliUEHistograms");
    
  return list->FindObject("AliUEHistogramsSame");
}

void CompareBias(const char* mcFile = "PWG4_JetTasksOutput.root", Int_t region, Int_t ueHist)
{
  loadlibs();

  AliUEHistograms* h = (AliUEHistograms*) GetUEHistogram(mcFile);
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
  
void CompareBiasWithData(const char* mcFile = "PWG4_JetTasksOutput.root", const char* dataFile = "esd.root", Int_t region = 2)
{
  loadlibs();

  AliUEHistograms* h = (AliUEHistograms*) GetUEHistogram(mcFile);
  SetupRanges(h);
  
  biasFromMC   = (TH1*) h->GetNumberDensitypT()->GetBias(AliUEHist::kCFStepReconstructed, AliUEHist::kCFStepTracked, region, "z")->Clone("biasFromMC");
  
  AliUEHistograms* h2 = (AliUEHistograms*) GetUEHistogram(dataFile);
  SetupRanges(h2);
  
  biasFromData = (TH1*) h2->GetNumberDensitypT()->GetBias(AliUEHist::kCFStepBiasStudy, AliUEHist::kCFStepReconstructed, region, "z")->Clone("biasFromData");
  biasFromData2 = (TH1*) h2->GetNumberDensitypT()->GetBias(AliUEHist::kCFStepBiasStudy2, AliUEHist::kCFStepReconstructed, region, "z")->Clone("biasFromData2");
  
  DrawRatio(biasFromData, biasFromMC, "bias: data vs MC");
  DrawRatio(biasFromData, biasFromData2, "bias: data vs data two step");
}

Int_t count = 0;

void Compare(const char* fileName1, const char* fileName2, Int_t id, Int_t step1, Int_t step2, Int_t region, Float_t ptLeadMin = -1, Float_t ptLeadMax = -1)
{
  loadlibs();
  
  if (!gCache || !gFirst)
  {
    AliUEHistograms* h = (AliUEHistograms*) GetUEHistogram(fileName1);
    AliUEHistograms* h2 = (AliUEHistograms*) GetUEHistogram(fileName2);
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
  
  SetupRanges(h);
  SetupRanges(h2);
  
  TH1* hist1 = h->GetUEHist(id)->GetUEHist(step1, region, ptLeadMin, ptLeadMax);
  //TH1* hist1 = h->GetUEHist(id)->GetUEHist(step1, region, ptLeadMin, ptLeadMax);
  TH1* hist2 = h2->GetUEHist(id)->GetUEHist(step2, region, ptLeadMin, ptLeadMax);

  //hist1->Scale(1.0 / hist1->Integral());
  //hist2->Scale(1.0 / hist2->Integral());

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
  
  DrawRatio(hist1, hist2, Form("%d_%s_%d_%d_%d_%.2f_%.2f", count++, TString(gSystem->BaseName(fileName1)).Tokenize(".")->First()->GetName(), id, step1, region, ptLeadMin, ptLeadMax), syst);
}  

void CompareEventHist(const char* fileName1, const char* fileName2, Int_t id, Int_t step, Int_t var)
{
  loadlibs();

  AliUEHistograms* h = (AliUEHistograms*) GetUEHistogram(fileName1);
  SetupRanges(h);

  AliUEHistograms* h2 = (AliUEHistograms*) GetUEHistogram(fileName2);
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

  AliUEHistograms* h = (AliUEHistograms*) GetUEHistogram(fileName);
  SetupRanges(h);

  new TCanvas;
  h->GetUEHist(id)->GetUEHist(step, region, ptLeadMin, ptLeadMax)->Draw();
}

void ProfileMultiplicity(const char* fileName = "PWG4_JetTasksOutput.root")
{
  loadlibs();

  AliUEHistograms* h = (AliUEHistograms*) GetUEHistogram(fileName);

  new TCanvas;
  h->GetCorrelationMultiplicity()->Draw("colz");
  gPad->SetLogz();

  new TCanvas;
  h->GetCorrelationMultiplicity()->ProfileX()->DrawCopy()->Fit("pol1", "", "", 1, 10);
}

void SetupRanges(void* obj)
{
  ((AliUEHistograms*) obj)->SetEtaRange(-0.79, 0.79);
  ((AliUEHistograms*) obj)->SetPtRange(gpTMin, gpTMax);
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
    for (Float_t ptLeadMin = 0.51; ptLeadMin < 10; ptLeadMin += 3)
      DrawRatios(TString(Form("UE %d pT %f", compareUEHist, ptLeadMin)), corrected->GetUEHist(compareUEHist), comparison->GetUEHist(compareUEHist), compareStep, compareRegion, ptLeadMin, ptLeadMin + 1.48);      
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
  
  AliUEHistograms* corr = (AliUEHistograms*) GetUEHistogram(fileNameCorrections);
  SetupRanges(corr);
  corr->ExtendTrackingEfficiency();
  
  AliUEHistograms* testSample = corr;
  if (fileNameESD)
    testSample = (AliUEHistograms*) GetUEHistogram(fileNameESD);
      
  // copy to esd object
  AliUEHistograms* esd = (AliUEHistograms*) corr->Clone();
  esd->Reset();
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
  
  AliUEHistograms* corr = (AliUEHistograms*) GetUEHistogram(fileNameCorrections);
  
  AliUEHistograms* esd = (AliUEHistograms*) GetUEHistogram(fileNameESD);
  
  SetupRanges(corr);
  SetupRanges(esd);
  
  corr->ExtendTrackingEfficiency();
  
  if (contEnhancement)
  {
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
  }
  
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
  
  AliUEHistograms* corrected = (AliUEHistograms*) GetUEHistogram(fileNameData);
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

  AliUEHistograms* corrected = (AliUEHistograms*) GetUEHistogram(fileNameData);
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
  
  AliUEHistograms* corrected = (AliUEHistograms*) GetUEHistogram(fileNameData);
  SetupRanges(corrected);
  
  AliUEHistograms* corrections = (AliUEHistograms*) GetUEHistogram(fileNameCorrections);
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

void MergeList2(const char* listFile, const char* dir, Bool_t onlyPrintEvents = kFALSE, const char* targetDir = "PWG4_LeadingTrackUE")
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
    fileName.Form("%s/%s/PWG4_JetTasksOutput.root", dir, line.Data());
    Printf("%s", fileName.Data());
    
    TList* list = 0;
    AliUEHistograms* obj = (AliUEHistograms*) GetUEHistogram(fileName, &list);
    if (!obj)
      continue;
    
    if (!final)
    {
      final = (AliUEHistograms*) obj;
      //final->GetEventCount()->Draw(); return;
      Printf("Events: %d", (Int_t) final->GetEventCount()->ProjectionX()->GetBinContent(4));
      finalList = list;
    }
    else
    {
      additional = (AliUEHistograms*) obj;
      Printf("Events: %d", (Int_t) additional->GetEventCount()->ProjectionX()->GetBinContent(4));
      
      if (!onlyPrintEvents)
      {
        TList list2;
        list2.Add(additional);
        final->Merge(&list2);
      }
      delete additional;
      gFile->Close();
    }
  }
  
  if (onlyPrintEvents)
    return;
    
  Printf("Total events (at step 0): %d", (Int_t) final->GetEventCount()->ProjectionX()->GetBinContent(4));
  
  file3 = TFile::Open("merged.root", "RECREATE");
  file3->mkdir(targetDir);
  file3->cd(targetDir);
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

void RemoveBaseLine(TH1* hist)
{
  if (!hist)
    return;
    
  //hist->Rebin(2);
  //hist->Scale(0.5);

  hist->Fit("pol0", "0", "", 1.07, 2.07);
  
  if (!hist->GetFunction("pol0"))
    return;
  
  Float_t zyam = hist->GetFunction("pol0")->GetParameter(0);
  
  if (zyam <= 0)
    return;
    
  //return;
  
  for (Int_t i=1; i<=hist->GetNbinsX(); i++)
    hist->SetBinContent(i, hist->GetBinContent(i) - zyam);
}
 
void PlotDeltaPhiDistributions(const char* fileName1, const char* fileName2, Float_t yMax = 0.1)
{
  loadlibs();

  Int_t maxLeadingPt = 5;
  Float_t leadingPtArr[] = { 1.0, 2.0, 6.0, 10.0, 20.0, 40.0 };
  Int_t maxAssocPt = 5;
  Float_t assocPtArr[] =   { 1.0, 2.0, 6.0, 10.0, 20.0, 40.0 };
  
  AliUEHistograms* h = (AliUEHistograms*) GetUEHistogram(fileName1);
  AliUEHistograms* h2 = (AliUEHistograms*) GetUEHistogram(fileName2);
  
  TCanvas* canvas = new TCanvas("DeltaPhi", "DeltaPhi", 1000, 700);
  canvas->Divide(maxAssocPt, maxLeadingPt);
  
  TCanvas* canvas2 = new TCanvas("Centrality", "Centrality", 800, 600);
  h->GetCentralityDistribution()->Draw();
  gPad->SetLogy();
  
  TLegend* legend = new TLegend(0.2, 0.5, 0.95, 0.90);
  TLegend* legend2 = new TLegend(0.5, 0.63, 0.95, 0.90);
  legend2->SetFillColor(0);
  legend2->SetTextSize(0.04);
  
  for (Int_t i=0; i<maxLeadingPt; i++)
    for (Int_t j=0; j<maxAssocPt; j++)
    {
      TString str;
      str.Form("%.1f < p_{T,trig} < %.1f", leadingPtArr[i], leadingPtArr[i+1]);
      
      if (j == 0)
      {
        canvas2->cd();
        h->GetUEHist(2)->GetEventHist()->GetGrid(6)->SetRangeUser(0, leadingPtArr[i] + 0.01, leadingPtArr[i+1] - 0.01);
        centralityHist = h->GetUEHist(2)->GetEventHist()->ShowProjection(1, 6);
        centralityHist->SetLineColor(i+2);
        centralityHist->DrawCopy("SAME");
        h->GetUEHist(2)->GetEventHist()->GetGrid(6)->SetRangeUser(0, 0, -1);
        legend2->AddEntry(centralityHist, str);
      }
    
      canvas->cd(j+1 + i * maxAssocPt);
      gPad->SetLeftMargin(0.15);
      gPad->SetBottomMargin(0.2);
      gPad->SetTopMargin(0.01);
      gPad->SetRightMargin(0.01);
      
      if (i == 0 && j == 3)
        legend->Draw();
      
      if (assocPtArr[j] > leadingPtArr[i])
        continue;
    
      gpTMin = assocPtArr[j] + 0.01;
      gpTMax = assocPtArr[j+1] - 0.01;
      
      SetupRanges(h);
      // for HI file do not set range in eta anymore after it was changed to delta eta axis
      h->SetEtaRange(0, 0);

      SetupRanges(h2);
      
//       SPD clusters 1
//       centrality multiplicity
//       5 4349.5
//       10 3509.5
//       20 2369.5
//       30 1569.5
//       40 989.5
//       50 589.5
//       60 309.5
//       70 149.5
//       80 69.5
//       90 29.5 
      
      // 0-5%
      // 0-5% --> 10, 15
      // 0-10% --> 8, 15
      // 0-20% --> 6, 15
      // 20-40% --> 3, 5
      TH1* hist1 = h->GetUEHist(2)->GetUEHist(6, 0, leadingPtArr[i] + 0.01, leadingPtArr[i+1] - 0.01, 10, 15);
      TH1* hist2 = h->GetUEHist(2)->GetUEHist(6, 0, leadingPtArr[i] + 0.01, leadingPtArr[i+1] - 0.01, 3, 5);
      TH1* hist2b = h->GetUEHist(2)->GetUEHist(6, 0, leadingPtArr[i] + 0.01, leadingPtArr[i+1] - 0.01, 1, 2);
      TH1* hist3 = h2->GetUEHist(2)->GetUEHist(0, 0, leadingPtArr[i] + 0.01, leadingPtArr[i+1] - 0.01);
      
      RemoveBaseLine(hist1);
      RemoveBaseLine(hist2);
      RemoveBaseLine(hist2b);
      RemoveBaseLine(hist3);
      
      if (0)
      {
        hist1->Scale(1.0 / hist1->Integral());
        hist2->Scale(1.0 / hist2->Integral());
        hist3->Scale(1.0 / hist3->Integral());
      }
    
      if (i == 0 && j == 0)
      {
        legend->SetFillColor(0);
        legend->AddEntry(hist1, "Pb+Pb 0-5%");
        if (hist2)
          legend->AddEntry(hist2, "Pb+Pb 20-40%");
        if (hist2b)
          legend->AddEntry(hist2b, "Pb+Pb >40%");
        legend->AddEntry(hist3, "p+p 7 TeV");
        legend->SetTextSize(0.08);
      }
    
      Prepare1DPlot(hist1);
      Prepare1DPlot(hist2);
      Prepare1DPlot(hist2b);
      Prepare1DPlot(hist3);
      
      Double_t yMin = 0.01;
      Double_t yMax2 = yMax;
      
      if (yMax < 0)
      {
        yMin = -0.1; //TMath::Min(hist1->GetMinimum(), hist2->GetMinimum()) * 0.97;
        yMax2 = TMath::Max((0 && hist3) ? hist3->GetMaximum() : 0.0, TMath::Max(hist1->GetMaximum(), (hist2) ? hist2->GetMaximum() : 0.0)) * 1.03;
      }
    
      TH2F* dummy = new TH2F("dummy", "", 100, hist1->GetXaxis()->GetBinLowEdge(1), hist1->GetXaxis()->GetBinUpEdge(hist1->GetNbinsX()), 1000, yMin, yMax2); //TMath::Max(hist1->GetMaximum(), hist2->GetMaximum()) * 1.1);
      dummy->SetStats(kFALSE);
      dummy->SetXTitle(hist1->GetXaxis()->GetTitle());
      dummy->SetYTitle(hist1->GetYaxis()->GetTitle());
      dummy->GetYaxis()->SetTitleOffset(1);
      Prepare1DPlot(dummy);
    
      dummy->GetXaxis()->SetLabelSize(0.06);
      dummy->GetYaxis()->SetLabelSize(0.06);
      dummy->GetXaxis()->SetTitleSize(0.06);
      dummy->GetYaxis()->SetTitleSize(0.06);
      /*
      dummy->GetYaxis()->SetTitleOffset(0.8);
      */
      
      dummy->DrawCopy();
      
      hist1->Draw("SAME");
      if (hist2)
      {
        hist2->SetLineColor(2);
        hist2->Draw("SAME");
      }
      if (hist2b)
      {
        hist2b->SetLineColor(3);
        hist2b->Draw("SAME");
      }
      if (hist3)
      {
        hist3->SetLineColor(4);
        hist3->Draw("SAME");
      }
      
      latex = new TLatex(0.55, 0.8, str);
      latex->SetNDC();
      latex->SetTextSize(0.06);
      latex->Draw();
       
      str.Form("%.1f < p_{T,ass} < %.1f",  assocPtArr[j], assocPtArr[j+1]);
      latex = new TLatex(0.55, 0.88, str);
      latex->SetNDC();
      latex->SetTextSize(0.06);
      latex->Draw();
      
      //return;
    }

  canvas->SaveAs(Form("DeltaPhi_%.2f.png", yMax));
  
  canvas2->cd();
  legend2->Draw();

  //TString name;
  //name.Form("%s_%.2f_%.2f_%.2f_%.2f.png", TString(gSystem->BaseName(fileName1)).Tokenize(".")->First()->GetName(), leadingPtArr[i], leadingPtArr[i+1], assocPtArr[j], assocPtArr[j+1]);
      

}
  
void PlotDeltaEtaDeltaPhi(const char* fileName)
{
  loadlibs();
  
  AliUEHistograms* h = (AliUEHistograms*) GetUEHistogram(fileName);
  h->SetPtRange(3.01, 9.99);
  
  h->GetUEHist(2)->GetUEHist(6, 0, 4.01, 9.99, 10, 15, kTRUE)->Draw("SURF1");
}
