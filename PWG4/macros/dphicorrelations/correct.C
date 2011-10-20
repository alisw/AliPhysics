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
  dummy->SetYTitle(corr->GetYaxis()->GetTitle());
  
  //dummy = new TH2F("dummy", "", 100, corr->GetXaxis()->GetBinLowEdge(1), corr->GetXaxis()->GetBinUpEdge(corr->GetNbinsX()), 1000, 13.5, 20.5);
  //dummy->SetYTitle("1/N_{trig} dN/d#Delta#phi");
  
  dummy->SetStats(kFALSE);
  dummy->SetXTitle(corr->GetXaxis()->GetTitle());
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
  
  if (0 && gUEHist != 2)
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
  //canvas3->SaveAs(Form("%s.eps", canvas3->GetTitle()));
}

void loadlibs()
{
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libJETAN");
  gSystem->Load("libPWG4JetTasks");
}

const char* lastFileName = 0;
void* cacheSameEvent = 0;
void* cacheMixedEvent = 0;

void* GetUEHistogram(const char* fileName, TList** listRef = 0, Bool_t mixed = kFALSE)
{
  if (!lastFileName || strcmp(lastFileName, fileName) != 0)
  {
    lastFileName = fileName;
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
      
    cacheMixedEvent = list->FindObject("AliUEHistogramsMixed");
    cacheSameEvent = list->FindObject("AliUEHistogramsSame");

    if (mixed)
      return cacheMixedEvent;
    
    if (list->FindObject("AliUEHistograms"))
      return list->FindObject("AliUEHistograms");
      
    return cacheSameEvent;
  }
  else
  {
    Printf("GetUEHistogram --> Using cache for %s", fileName);
    
    if (mixed)
      return cacheMixedEvent;
    else
      return cacheSameEvent;
  }
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

void Compare(const char* fileName1, const char* fileName2, Int_t id, Int_t step1, Int_t step2, Int_t region, Float_t ptLeadMin = -1, Float_t ptLeadMax = -1, Int_t centralityBegin = 0, Int_t centralityEnd = -1)
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
  
  TH1* hist1 = h->GetUEHist(id)->GetUEHist(step1, region, ptLeadMin, ptLeadMax, centralityBegin, centralityEnd);
  //TH1* hist1 = h->GetUEHist(id)->GetUEHist(step1, region, ptLeadMin, ptLeadMax);
  TH1* hist2 = h2->GetUEHist(id)->GetUEHist(step2, region, ptLeadMin, ptLeadMax, centralityBegin, centralityEnd);

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
  TH1* syst = 0; //GetSystematicUncertainty(hist1, trackHist);
  
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

void CompareEventHist(const char* fileName1, const char* fileName2, Int_t id, Int_t step1, Int_t step2, Int_t var)
{
  loadlibs();

  AliUEHistograms* h = (AliUEHistograms*) GetUEHistogram(fileName1);
  SetupRanges(h);

  AliUEHistograms* h2 = (AliUEHistograms*) GetUEHistogram(fileName2);
  SetupRanges(h2);
  
//   h->GetUEHist(id)->GetEventHist()->GetGrid(step1)->GetGrid()->GetAxis(1)->SetRange(1, 5);
//   h2->GetUEHist(id)->GetEventHist()->GetGrid(step2)->GetGrid()->GetAxis(1)->SetRange(1, 5);
  
  TH1* hist1 = h->GetUEHist(id)->GetEventHist()->ShowProjection(var, step1);
  TH1* hist2 = h2->GetUEHist(id)->GetEventHist()->ShowProjection(var, step2);
 
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

void CompareStep(const char* fileName1, const char* fileName2, Int_t id, Int_t step1, Int_t step2, Int_t region, Float_t ptLeadMin = -1, Float_t ptLeadMax = -1, Int_t centralityBegin = 0, Int_t centralityEnd = -1)
{
  // fileName1 is labelled Corrected in the plot

  loadlibs();

  gUEHist = id;
  gRegion = region;
  Compare(fileName1, fileName2, id, step1, step2, region, ptLeadMin, ptLeadMax, centralityBegin, centralityEnd);
}

TH1* DrawStep(const char* fileName, Int_t id, Int_t step, Int_t region, Float_t ptLeadMin = -1, Float_t ptLeadMax = -1, Int_t centralityBegin = 0, Int_t centralityEnd = -1, Int_t twoD = 0, Bool_t mixed = kFALSE)
{
  loadlibs();

  AliUEHistograms* h = (AliUEHistograms*) GetUEHistogram(fileName, 0, mixed);
  SetupRanges(h);

  new TCanvas;
  return h->GetUEHist(id)->GetUEHist(step, region, ptLeadMin, ptLeadMax, centralityBegin, centralityEnd, twoD)->DrawCopy();
}

void DrawExample(const char* fileName, const char* fileNamePbPbMix)
{
  gpTMin = 1.01;
  gpTMax = 1.99;
  
 loadlibs();
  
  AliUEHistograms* h = (AliUEHistograms*) GetUEHistogram(fileName);
  hMixed = (AliUEHistograms*) GetUEHistogram(fileNamePbPbMix, 0, kTRUE);
  
  SetupRanges(h);
  SetupRanges(hMixed);

  TH1* hist1 = 0;
  
  GetDistAndFlow(h, hMixed, &hist1,  0, 6, 0,  2, 2.01, 3.99, 1, kTRUE, 0, kFALSE); 
  
  ((TH2*) hist1)->Rebin2D(2, 2);
  hist1->Scale(0.25);
  
  new TCanvas("c", "c", 800, 800);
  gPad->SetLeftMargin(0.15);
  hist1->SetTitle("");
  hist1->GetYaxis()->SetRangeUser(-1.79, 1.79);
  hist1->GetXaxis()->SetTitleOffset(1.5);
  hist1->GetYaxis()->SetTitleOffset(2);
  hist1->SetStats(kFALSE);
  hist1->Draw("SURF1");
  
  latex = new TLatex(0.82, 0.74, "ALICE performance");
  latex->SetTextSize(0.02);
  latex->SetTextAlign(22);
  latex->SetNDC();
  latex->Draw();
  latex = new TLatex(0.82, 0.72, "Pb-Pb 2.76 TeV");
  latex->SetTextSize(0.02);
  latex->SetTextAlign(22);
  latex->SetNDC();
  latex->Draw();
  latex = new TLatex(0.82, 0.70, "28.09.11");
  latex->SetTextSize(0.02);
  latex->SetTextAlign(22);
  latex->SetNDC();
  latex->Draw();

  DrawALICELogo(0.75, 0.75, 0.9, 0.9);
}

void DrawValidation(const char* fileName1, const char* fileName2)
{
  gpTMin = 1.01;
  gpTMax = 3.99;
  
  CompareStep(fileName1, fileName2, 2, 6, 4, 0, 4.01, 19.99);
  CompareStep(fileName1, fileName2, 2, 4, 4, 0, 4.01, 19.99);

  CompareStep(fileName1, fileName2, 2, 4, 2, 0, 4.01, 19.99);
  CompareStep(fileName1, fileName2, 2, 2, 2, 0, 4.01, 19.99);
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
  ((AliUEHistograms*) obj)->SetEtaRange(0, 0);
//   ((AliUEHistograms*) obj)->SetEtaRange(-0.99, 0.99); Printf("WARNING: Setting eta Range!");
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
    if (compareStep == -1 && (step == AliUEHist::kCFStepAnaTopology || step == AliUEHist::kCFStepTriggered))
      continue;
      
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
    for (Float_t ptMin = 2.01; ptMin < 8; ptMin += 2)
    {
      ((AliUEHistograms*) corrected)->SetPtRange(ptMin, ptMin + 1.98);
      ((AliUEHistograms*) comparison)->SetPtRange(ptMin, ptMin + 1.98);
      
      DrawRatios(TString(Form("Dphi %d pT %f", compareUEHist, ptMin)), corrected->GetUEHist(compareUEHist), comparison->GetUEHist(compareUEHist), compareStep, compareRegion, 8.01, 19.99);      
    }
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
  
  list = new TList;
  list->Add(esd);
  
  file3 = TFile::Open("correctedMC.root", "RECREATE");
  file3->mkdir("PWG4_PhiCorrelations");
  file3->cd("PWG4_PhiCorrelations");
  list->Write("histosPhiCorrelations", TObject::kSingleKey);
  file3->Close();
  
  if (1)
    DrawRatios(esd, testSample, compareStep, compareRegion, compareUEHist);
  
  if (1)
  {
    esd->SetPtRange(2.01, 3.99);
    corrected = esd->GetUEHist(2)->GetUEHist(0, 0, 4.01, 7.99, 0, -1, 1);
    testSample->SetPtRange(2.01, 3.99);
    mc = testSample->GetUEHist(2)->GetUEHist(0, 0, 4.01, 7.99, 0, -1, 1);
    new TCanvas; corrected->DrawCopy("SURF1");
    new TCanvas; mc->DrawCopy("SURF1");
    new TCanvas; mc->DrawCopy("SURF1")->Divide(corrected);
  }
  
  //CompareEventsTracks(esd, testSample, compareStep, compareRegion, compareUEHist);
}

// function to compare only final step for all regions and distributions

void correctData(const char* fileNameCorrections, const char* fileNameESD, const char* contEnhancement = 0, Float_t contEncUpTo = 1.0, Int_t compareStep = 0, Int_t compareRegion = 0, Int_t compareUEHist = 2)
{
  // corrects fileNameESD with fileNameCorrections and compares the two
  
  loadlibs();
  
  AliUEHistograms* corr = (AliUEHistograms*) GetUEHistogram(fileNameCorrections);
  
  TList* list = 0;
  AliUEHistograms* esd = (AliUEHistograms*) GetUEHistogram(fileNameESD, &list);
  
  SetupRanges(corr);
  SetupRanges(esd);
  
  corr->SetEtaRange(-0.99, 0.99);
  corr->ExtendTrackingEfficiency(0);
  
//   return;
  
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
  //esd->GetUEHist(2)->AdditionalDPhiCorrection(0);
  
  file3 = TFile::Open("corrected.root", "RECREATE");
  file3->mkdir("PWG4_PhiCorrelations");
  file3->cd("PWG4_PhiCorrelations");
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

/*
TGraph* GetFlow2040()
{
  // from first ALICE flow paper (provided by Raimond)
  // http://www-library.desy.de/spires/find/hep/www?eprint=arXiv:1011.3914
  
  // centrality 20-30% 
  Double_t xCumulant4th2030ALICE[] = {0.050000,0.150000,0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,0.850000,0.950000,
  1.100000,1.300000,1.500000,1.700000,1.900000,2.250000,2.750000,3.250000,3.750000,4.500000,
  5.500000,7.000000,9.000000};
  Double_t yCumulant4th2030ALICE[] = {0.000000,0.000000,0.030926,0.041076,0.052063,0.059429,0.070500,0.084461,0.086745,0.099254,
  0.109691,0.116398,0.130831,0.141959,0.158932,0.169680,0.171387,0.178858,0.171475,0.140358,
  0.000000,0.000000,0.000000};
  Double_t xErrCumulant4th2030ALICE[23] = {0.};
  Double_t yErrCumulant4th2030ALICE[] = {0.000000,0.000000,0.002857,0.003451,0.003567,0.003859,0.004609,0.004976,0.005412,0.006277,
  0.004748,0.005808,0.006896,0.007987,0.008683,0.008080,0.013278,0.018413,0.024873,0.026057,
  0.000000,0.000000,0.000000};
  Int_t nPointsCumulant4th2030ALICE = sizeof(xCumulant4th2030ALICE)/sizeof(Double_t);                                      
  TGraphErrors *Cumulant4th2030ALICE = new TGraphErrors(nPointsCumulant4th2030ALICE,xCumulant4th2030ALICE,yCumulant4th2030ALICE,
                                                        xErrCumulant4th2030ALICE,yErrCumulant4th2030ALICE);
  Cumulant4th2030ALICE->SetMarkerStyle(kFullSquare);
  Cumulant4th2030ALICE->SetMarkerColor(kRed);
  Cumulant4th2030ALICE->SetMarkerSize(1.2);
  Cumulant4th2030ALICE->SetFillStyle(1001);
  Cumulant4th2030ALICE->SetFillColor(kRed-10);
  
  //===================================================================================================================
  // centrality 30-40% 
  Double_t xCumulant4th3040ALICE[] = {0.050000,0.150000,0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,0.850000,0.950000,
  1.100000,1.300000,1.500000,1.700000,1.900000,2.250000,2.750000,3.250000,3.750000,4.500000,
  5.500000,7.000000,9.000000};
  Double_t yCumulant4th3040ALICE[] = {0.000000,0.000000,0.037071,0.048566,0.061083,0.070910,0.078831,0.091396,0.102026,0.109691,
  0.124449,0.139819,0.155561,0.165701,0.173678,0.191149,0.202015,0.204540,0.212560,0.195885,
  0.000000,0.000000,0.000000};
  Double_t xErrCumulant4th3040ALICE[23] = {0.};
  Double_t yErrCumulant4th3040ALICE[] = {0.000000,0.000000,0.002992,0.003364,0.003669,0.003931,0.004698,0.005261,0.005446,0.006151,
  0.004980,0.005741,0.007198,0.008576,0.010868,0.009926,0.015269,0.020691,0.027601,0.031834,
  0.000000,0.000000,0.000000};
  Int_t nPointsCumulant4th3040ALICE = sizeof(xCumulant4th3040ALICE)/sizeof(Double_t);
  TGraphErrors *Cumulant4th3040ALICE = new TGraphErrors(nPointsCumulant4th3040ALICE,xCumulant4th3040ALICE,yCumulant4th3040ALICE,
                                                        xErrCumulant4th3040ALICE,yErrCumulant4th3040ALICE);
  Cumulant4th3040ALICE->SetMarkerStyle(kFullTriangleUp);
  Cumulant4th3040ALICE->SetMarkerColor(kGreen+2);
  Cumulant4th3040ALICE->SetMarkerSize(1.2);
  Cumulant4th3040ALICE->SetFillStyle(1001);
  Cumulant4th3040ALICE->SetFillColor(kGreen+2);
  
  // build average between the two (for class 20-40%)
  Double_t* yAverage = new Double_t[nPointsCumulant4th3040ALICE];
  for (Int_t i=0; i<nPointsCumulant4th3040ALICE; i++)
    yAverage[i] = (yCumulant4th2030ALICE[i] + yCumulant4th3040ALICE[i]) / 2;
    
  // assume flow constant above highest pT; not neccessarily physically sound ;)
  if (1)
  {
    yAverage[20] = yAverage[19];
    xCumulant4th3040ALICE[20] = 100;
    nPointsCumulant4th3040ALICE -= 2;
  }
  
  TGraph *flow2040 = new TGraph(nPointsCumulant4th3040ALICE,xCumulant4th3040ALICE,yAverage);
  
  if (0)
  {
    flow2040->Draw("*A");
    Cumulant4th2030ALICE->Draw("PSAME");
    Cumulant4th3040ALICE->Draw("PSAME");
  }
  
  return flow2040;
}
*/

/*
TGraph* GetFlow1020()
{
  // from first ALICE flow paper (provided by Raimond)
  // http://www-library.desy.de/spires/find/hep/www?eprint=arXiv:1011.3914
  
  //===================================================================================================================
  // centrality 10-20% 
  Double_t xCumulant4th1020ALICE[] = {0.050000,0.150000,0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,0.850000,0.950000,
  1.100000,1.300000,1.500000,1.700000,1.900000,2.250000,2.750000,3.250000,3.750000,4.500000,
  5.500000,7.000000,9.000000};
  Double_t yCumulant4th1020ALICE[] = {0.000000,0.000000,0.024075,0.031505,0.040413,0.044981,0.055358,0.060563,0.063378,0.070030,
  0.082692,0.091611,0.099641,0.107223,0.122376,0.131240,0.137425,0.146050,0.131365,0.124708,
  0.000000,0.000000,0.000000};
  Double_t xErrCumulant4th1020ALICE[23] = {0.};
  Double_t yErrCumulant4th1020ALICE[] = {0.000000,0.000000,0.002413,0.002931,0.003444,0.003950,0.004338,0.004835,0.005059,0.005586,
  0.004521,0.005278,0.005999,0.007072,0.008260,0.007279,0.011897,0.017409,0.023995,0.025701,
  0.000000,0.000000,0.000000};
  Int_t nPointsCumulant4th1020ALICE = sizeof(xCumulant4th1020ALICE)/sizeof(Double_t);                                      
  
  // assume flow constant above highest pT; not neccessarily physically sound ;)
  if (1)
  {
    yCumulant4th1020ALICE[20] = yCumulant4th1020ALICE[19];
    xCumulant4th1020ALICE[20] = 100;
    nPointsCumulant4th1020ALICE -= 2;
  }
  
  TGraphErrors *Cumulant4th1020ALICE = new TGraphErrors(nPointsCumulant4th1020ALICE,xCumulant4th1020ALICE,yCumulant4th1020ALICE,
                                                        xErrCumulant4th1020ALICE,yErrCumulant4th1020ALICE);
  
 Cumulant4th1020ALICE->SetMarkerStyle(kFullCircle);
 Cumulant4th1020ALICE->SetMarkerColor(kBlue);
 Cumulant4th1020ALICE->SetMarkerSize(1.2);
 Cumulant4th1020ALICE->SetFillStyle(1001);
 Cumulant4th1020ALICE->SetFillColor(kBlue-10);
  
  TGraph *flow1020 = new TGraph(nPointsCumulant4th1020ALICE,xCumulant4th1020ALICE,yCumulant4th1020ALICE);
  
  if (0)
  {
    flow1020->Draw("*A");
    Cumulant4th1020ALICE->Draw("PSAME");
  }
  
  return flow1020;
}
*/

/*
TGraph* GetFlow05()
{
  // takes flow measurement from 10-20% and scales by a factor extracted from Fig.3 in the ALICE flow paper
  // factor = integrated flow in 0-5% / integrated flow in 10-20%
  
  graph = GetFlow1020();
  for (Int_t i=0; i<graph->GetN(); i++)
    graph->GetY()[i] *= 0.016 / 0.055;
    
  return graph;
}
*/
TGraphErrors* GetFlow01_Rap02(Int_t n)
{
  // private communication 19.04.11, Raimond / Ante

  if (n == 2)
  {
    //  v2{SP}(pt) for 0-1%, rapidity gap = 0.2:
    const Int_t nPointsSP_0001ALICE_v2_etaGap02 = 18;
    Double_t xSP_0001ALICE_v2_etaGap02[nPointsSP_0001ALICE_v2_etaGap02] = {0.250000,0.350000,0.450000,0.550000,0.650000,
    0.750000,0.850000,0.950000,1.100000,1.300000,1.500000,1.700000,1.900000,2.250000,2.750000,3.250000,3.750000,4.500000};
    Double_t ySP_0001ALICE_v2_etaGap02[nPointsSP_0001ALICE_v2_etaGap02] = {0.009235,0.014105,0.017274,0.018245,0.023190,
    0.024871,0.028216,0.031506,0.034179,0.035337,0.039836,0.045261,0.043393,0.050693,0.056469,0.055247,0.049718,0.052748};
    Double_t xErrSP_0001ALICE_v2_etaGap02[nPointsSP_0001ALICE_v2_etaGap02] = {0.};
    Double_t yErrSP_0001ALICE_v2_etaGap02[nPointsSP_0001ALICE_v2_etaGap02] = {0.000515,0.000504,0.000532,0.000585,0.000641,
    0.000709,0.000788,0.000876,0.000723,0.000891,0.001098,0.001354,0.001671,0.001485,0.002463,0.004038,0.006441,0.008091};
    TGraphErrors *SP_0001ALICE_v2_etaGap02 = new TGraphErrors(nPointsSP_0001ALICE_v2_etaGap02,xSP_0001ALICE_v2_etaGap02,
							      ySP_0001ALICE_v2_etaGap02,xErrSP_0001ALICE_v2_etaGap02,yErrSP_0001ALICE_v2_etaGap02);
    SP_0001ALICE_v2_etaGap02->SetMarkerStyle(kFullCircle);
    SP_0001ALICE_v2_etaGap02->SetMarkerColor(kBlue);  
    SP_0001ALICE_v2_etaGap02->SetFillStyle(1001);
    SP_0001ALICE_v2_etaGap02->SetFillColor(kBlue-10);  
    
    return SP_0001ALICE_v2_etaGap02;
  }
  
  if (n == 3)
  {
    //  v3{SP}(pt) for 0-1%, rapidity gap = 0.2:
    const Int_t nPointsSP_0001ALICE_v3_etaGap02 = 18;
    Double_t xSP_0001ALICE_v3_etaGap02[nPointsSP_0001ALICE_v3_etaGap02] = {0.250000,0.350000,0.450000,0.550000,0.650000,
    0.750000,0.850000,0.950000,1.100000,1.300000,1.500000,1.700000,1.900000,2.250000,2.750000,3.250000,3.750000,4.500000};
    Double_t ySP_0001ALICE_v3_etaGap02[nPointsSP_0001ALICE_v3_etaGap02] = {0.005688,0.007222,0.010305,0.013795,0.016077,
    0.018693,0.022310,0.026991,0.030162,0.035119,0.043097,0.048201,0.058249,0.063273,0.079233,0.083465,0.087807,0.069577};
    Double_t xErrSP_0001ALICE_v3_etaGap02[nPointsSP_0001ALICE_v3_etaGap02] = {0.};
    Double_t yErrSP_0001ALICE_v3_etaGap02[nPointsSP_0001ALICE_v3_etaGap02] = {0.000585,0.000582,0.000614,0.000667,0.000734,
    0.000811,0.000898,0.000989,0.000817,0.001000,0.001234,0.001517,0.001874,0.001669,0.002765,0.004528,0.007202,0.009066};
    TGraphErrors *SP_0001ALICE_v3_etaGap02 = new TGraphErrors(nPointsSP_0001ALICE_v3_etaGap02,xSP_0001ALICE_v3_etaGap02,ySP_0001ALICE_v3_etaGap02,
							      xErrSP_0001ALICE_v3_etaGap02,yErrSP_0001ALICE_v3_etaGap02);
    SP_0001ALICE_v3_etaGap02->SetMarkerStyle(kFullTriangleUp);
    SP_0001ALICE_v3_etaGap02->SetMarkerSize(1.4);  
    SP_0001ALICE_v3_etaGap02->SetMarkerColor(kGreen+2);
    SP_0001ALICE_v3_etaGap02->SetFillStyle(1001);
    SP_0001ALICE_v3_etaGap02->SetFillColor(kGreen-10);     
    
    return SP_0001ALICE_v3_etaGap02;
  }
   
  if (n == 4)
  {
    //  v4{SP}(pt) for 0-1%, rapidity gap = 0.2:
    const Int_t nPointsSP_0001ALICE_v4_etaGap02 = 18;
    Double_t xSP_0001ALICE_v4_etaGap02[nPointsSP_0001ALICE_v4_etaGap02] = {0.250000,0.350000,0.450000,0.550000,0.650000,
    0.750000,0.850000,0.950000,1.100000,1.300000,1.500000,1.700000,1.900000,2.250000,2.750000,3.250000,3.750000,4.500000};
    Double_t ySP_0001ALICE_v4_etaGap02[nPointsSP_0001ALICE_v4_etaGap02] = {0.001189,0.003540,0.004682,0.004210,0.007032,
    0.008627,0.010226,0.013671,0.016214,0.020054,0.023878,0.033939,0.033693,0.040006,0.055612,0.066287,0.074857,0.078751};
    Double_t xErrSP_0001ALICE_v4_etaGap02[nPointsSP_0001ALICE_v4_etaGap02] = {0.};
    Double_t yErrSP_0001ALICE_v4_etaGap02[nPointsSP_0001ALICE_v4_etaGap02] = {0.001035,0.001017,0.001081,0.001187,0.001299,
    0.001432,0.001590,0.001757,0.001443,0.001769,0.002175,0.002674,0.003296,0.002924,0.004844,0.007921,0.012524,0.015771};
    TGraphErrors *SP_0001ALICE_v4_etaGap02 = new TGraphErrors(nPointsSP_0001ALICE_v4_etaGap02,xSP_0001ALICE_v4_etaGap02,ySP_0001ALICE_v4_etaGap02,
							      xErrSP_0001ALICE_v4_etaGap02,yErrSP_0001ALICE_v4_etaGap02);
    SP_0001ALICE_v4_etaGap02->SetMarkerStyle(kFullSquare);
    SP_0001ALICE_v4_etaGap02->SetMarkerColor(kRed);
    SP_0001ALICE_v4_etaGap02->SetFillStyle(1001);
    SP_0001ALICE_v4_etaGap02->SetFillColor(kRed-10);  
    
    return SP_0001ALICE_v4_etaGap02;
  }
  
  return 0;
}

TGraphErrors* GetFlow01_Rap10(Int_t n)
{
  if (n == 2)
  {
    //  v2{SP}(pt) for 0-1%, rapidity gap = 1.0:
    const Int_t nPointsSP_0001ALICE_v2_etaGap10 = 21;
    Double_t xSP_0001ALICE_v2_etaGap10[nPointsSP_0001ALICE_v2_etaGap10] = {0.250000,0.350000,0.450000,0.550000,0.650000,
    0.750000,0.850000,0.950000,1.100000,1.300000,1.500000,1.700000,1.900000,2.100000,2.300000,2.500000,2.700000,2.900000,
    3.250000,3.750000,4.500000};
    Double_t ySP_0001ALICE_v2_etaGap10[nPointsSP_0001ALICE_v2_etaGap10] = {0.009129,0.013461,0.017567,0.018041,0.020384,
    0.023780,0.021647,0.029543,0.028912,0.029464,0.037016,0.044131,0.043135,0.047286,0.051983,0.049311,0.050472,0.046569,
    0.036905,0.054836,0.030527};
    Double_t xErrSP_0001ALICE_v2_etaGap10[nPointsSP_0001ALICE_v2_etaGap10] = {0.};
    Double_t yErrSP_0001ALICE_v2_etaGap10[nPointsSP_0001ALICE_v2_etaGap10] = {0.001179,0.001152,0.001219,0.001339,0.001480,
    0.001644,0.001831,0.002016,0.001662,0.002033,0.002497,0.003056,0.003777,0.004645,0.005713,0.007069,0.008540,0.010447,
    0.009145,0.014749,0.018698};
    TGraphErrors *SP_0001ALICE_v2_etaGap10 = new TGraphErrors(nPointsSP_0001ALICE_v2_etaGap10,xSP_0001ALICE_v2_etaGap10,
							      ySP_0001ALICE_v2_etaGap10,xErrSP_0001ALICE_v2_etaGap10,yErrSP_0001ALICE_v2_etaGap10);
    SP_0001ALICE_v2_etaGap10->SetMarkerStyle(kOpenCircle);
    SP_0001ALICE_v2_etaGap10->SetMarkerColor(kBlue);  
    SP_0001ALICE_v2_etaGap10->SetFillStyle(1001);
    SP_0001ALICE_v2_etaGap10->SetFillColor(kBlue-10);  
    
    return SP_0001ALICE_v2_etaGap10;
  }

  if (n == 3)
  {
    //  v3{SP}(pt) for 0-1%, rapidity gap = 1.0:
    const Int_t nPointsSP_0001ALICE_v3_etaGap10 = 18;
    Double_t xSP_0001ALICE_v3_etaGap10[nPointsSP_0001ALICE_v3_etaGap10] = {0.250000,0.350000,0.450000,0.550000,0.650000,
    0.750000,0.850000,0.950000,1.100000,1.300000,1.500000,1.700000,1.900000,2.250000,2.750000,3.250000,3.750000,4.500000};
    Double_t ySP_0001ALICE_v3_etaGap10[nPointsSP_0001ALICE_v3_etaGap10] = {0.006373,0.008403,0.010848,0.011505,0.016728,
    0.018519,0.020163,0.027119,0.029315,0.036832,0.040974,0.043287,0.054395,0.060676,0.081763,0.074333,0.096016,0.074909};
    Double_t xErrSP_0001ALICE_v3_etaGap10[nPointsSP_0001ALICE_v3_etaGap10] = {0.};
    Double_t yErrSP_0001ALICE_v3_etaGap10[nPointsSP_0001ALICE_v3_etaGap10] = {0.001286,0.001269,0.001346,0.001474,0.001620,
    0.001796,0.001991,0.002187,0.001803,0.002203,0.002697,0.003316,0.004078,0.003640,0.006050,0.009873,0.015824,0.020174};
    TGraphErrors *SP_0001ALICE_v3_etaGap10 = new TGraphErrors(nPointsSP_0001ALICE_v3_etaGap10,xSP_0001ALICE_v3_etaGap10,ySP_0001ALICE_v3_etaGap10,
							      xErrSP_0001ALICE_v3_etaGap10,yErrSP_0001ALICE_v3_etaGap10);
    SP_0001ALICE_v3_etaGap10->SetMarkerStyle(kOpenTriangleUp);
    SP_0001ALICE_v3_etaGap10->SetMarkerSize(1.2);  
    SP_0001ALICE_v3_etaGap10->SetMarkerColor(kGreen+2);
    SP_0001ALICE_v3_etaGap10->SetFillStyle(1001);
    SP_0001ALICE_v3_etaGap10->SetFillColor(kGreen-10);     
    
    return SP_0001ALICE_v3_etaGap10;
  }

  if (n == 4)
  {
    //  v4{SP}(pt) for 0-1%, rapidity gap = 1.0:
    const Int_t nPointsSP_0001ALICE_v4_etaGap10 = 11;
    Double_t xSP_0001ALICE_v4_etaGap10[nPointsSP_0001ALICE_v4_etaGap10] = {0.300000,0.500000,0.700000,0.900000,1.200000,1.600000,2.000000,2.400000,2.800000,3.500000,4.500000};
    Double_t ySP_0001ALICE_v4_etaGap10[nPointsSP_0001ALICE_v4_etaGap10] = {-0.000458,0.006444,0.005490,0.010870,0.018866,0.024370,0.029703,0.052505,0.060334,0.048189,0.128184};
    Double_t xErrSP_0001ALICE_v4_etaGap10[nPointsSP_0001ALICE_v4_etaGap10] = {0.};
    Double_t yErrSP_0001ALICE_v4_etaGap10[nPointsSP_0001ALICE_v4_etaGap10] = {0.001901,0.002012,0.002477,0.003014,0.002852,0.004297,0.006491,0.009846,0.014623,0.017120,0.040568};
    TGraphErrors *SP_0001ALICE_v4_etaGap10 = new TGraphErrors(nPointsSP_0001ALICE_v4_etaGap10,xSP_0001ALICE_v4_etaGap10,ySP_0001ALICE_v4_etaGap10,
							      xErrSP_0001ALICE_v4_etaGap10,yErrSP_0001ALICE_v4_etaGap10);
    SP_0001ALICE_v4_etaGap10->SetMarkerStyle(kOpenSquare);
    SP_0001ALICE_v4_etaGap10->SetMarkerColor(kRed);
    SP_0001ALICE_v4_etaGap10->SetFillStyle(1001);
    SP_0001ALICE_v4_etaGap10->SetFillColor(kRed-10);  
    
    return SP_0001ALICE_v4_etaGap10;
  }
  
  if (n == 5)
  {
    //  v5{SP}(pt) for 0-1%, rapidity gap = 1.0:
    const Int_t nPointsSP_0001ALICE_v5_etaGap10 = 11;
    Double_t xSP_0001ALICE_v5_etaGap10[nPointsSP_0001ALICE_v5_etaGap10] = {0.300000,0.500000,0.700000,0.900000,1.200000,1.600000,2.000000,2.400000,2.800000,3.500000,4.500000};
    Double_t ySP_0001ALICE_v5_etaGap10[nPointsSP_0001ALICE_v5_etaGap10] = {0.007022,0.001344,0.008380,0.004298,-0.001444,0.014114,0.015012,0.041880,0.019820,0.042083,0.015268};
    Double_t xErrSP_0001ALICE_v5_etaGap10[nPointsSP_0001ALICE_v5_etaGap10] = {0.};
    Double_t yErrSP_0001ALICE_v5_etaGap10[nPointsSP_0001ALICE_v5_etaGap10] = {0.002713,0.003167,0.003741,0.004650,0.004525,0.006578,0.009986,0.015185,0.022535,0.026356,0.064773};
    TGraphErrors *SP_0001ALICE_v5_etaGap10 = new TGraphErrors(nPointsSP_0001ALICE_v5_etaGap10,xSP_0001ALICE_v5_etaGap10,ySP_0001ALICE_v5_etaGap10,
							      xErrSP_0001ALICE_v5_etaGap10,yErrSP_0001ALICE_v5_etaGap10);
    return SP_0001ALICE_v5_etaGap10;
  }
}

TGraphErrors* GetFlow02_Rap10(Int_t n)
{
  // private communication 20.04.11, Ante B. / Raimond

  if (n == 2)
  {
     //  v2{SP}(pt) for 00-02%, eta gap = 1.0:
    const Int_t nPointsSP_0002_v2_etaGap10 = 15;
    Double_t xSP_0002_v2_etaGap10[nPointsSP_0002_v2_etaGap10] = {0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,
    0.900000,1.100000,1.350000,1.650000,1.950000,2.250000,2.700000,3.500000,4.500000};
    Double_t ySP_0002_v2_etaGap10[nPointsSP_0002_v2_etaGap10] = {0.010171,0.013190,0.017342,0.020629,0.022617,0.026549,
    0.027423,0.032261,0.037467,0.041001,0.045763,0.049327,0.049688,0.051480,0.038527};
    Double_t xErrSP_0002_v2_etaGap10[nPointsSP_0002_v2_etaGap10] = {0.};
    Double_t yErrSP_0002_v2_etaGap10[nPointsSP_0002_v2_etaGap10] = {0.000600,0.000590,0.000625,0.000683,0.000757,0.000839,
    0.000692,0.000848,0.000888,0.001209,0.001653,0.002252,0.002465,0.003968,0.009391};
    TGraphErrors *SP_0002_v2_etaGap10 = new TGraphErrors(nPointsSP_0002_v2_etaGap10,xSP_0002_v2_etaGap10,ySP_0002_v2_etaGap10,
                                                  xErrSP_0002_v2_etaGap10,yErrSP_0002_v2_etaGap10);
						  
    return SP_0002_v2_etaGap10;
  }
  
  if (n == 3)
  {
    const Int_t nPointsSP_0002_v3_etaGap10 = 15;
    Double_t xSP_0002_v3_etaGap10[nPointsSP_0002_v3_etaGap10] = {0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,
    0.900000,1.100000,1.350000,1.650000,1.950000,2.250000,2.700000,3.500000,4.500000};
    Double_t ySP_0002_v3_etaGap10[nPointsSP_0002_v3_etaGap10] = {0.006592,0.007286,0.012180,0.012242,0.017416,0.018393,
    0.024716,0.030980,0.037703,0.046558,0.051285,0.064613,0.074831,0.077093,0.082442};
    Double_t xErrSP_0002_v3_etaGap10[nPointsSP_0002_v3_etaGap10] = {0.};
    Double_t yErrSP_0002_v3_etaGap10[nPointsSP_0002_v3_etaGap10] = {0.000682,0.000676,0.000713,0.000782,0.000860,0.000953,
    0.000782,0.000957,0.001002,0.001361,0.001862,0.002541,0.002767,0.004466,0.010586};
    TGraphErrors *SP_0002_v3_etaGap10 = new TGraphErrors(nPointsSP_0002_v3_etaGap10,xSP_0002_v3_etaGap10,ySP_0002_v3_etaGap10,
							  xErrSP_0002_v3_etaGap10,yErrSP_0002_v3_etaGap10);    
							  
    return SP_0002_v3_etaGap10;
  }
  
  if (n == 4)
  {
    const Int_t nPointsSP_0002_v4_etaGap10 = 15;
    Double_t xSP_0002_v4_etaGap10[nPointsSP_0002_v4_etaGap10] = {0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,
    0.900000,1.100000,1.350000,1.650000,1.950000,2.250000,2.700000,3.500000,4.500000};
    Double_t ySP_0002_v4_etaGap10[nPointsSP_0002_v4_etaGap10] = {-0.000533,0.001167,0.002081,0.005218,0.006826,0.008440,
    0.013009,0.014812,0.017125,0.030106,0.038279,0.050488,0.067640,0.071637,0.084239};
    Double_t xErrSP_0002_v4_etaGap10[nPointsSP_0002_v4_etaGap10] = {0.};
    Double_t yErrSP_0002_v4_etaGap10[nPointsSP_0002_v4_etaGap10] = {0.001427,0.001398,0.001482,0.001594,0.001758,0.001945,
    0.001593,0.001951,0.002046,0.002787,0.003802,0.005182,0.005663,0.009064,0.021449};
    TGraphErrors *SP_0002_v4_etaGap10 = new TGraphErrors(nPointsSP_0002_v4_etaGap10,xSP_0002_v4_etaGap10,ySP_0002_v4_etaGap10,
                                                      xErrSP_0002_v4_etaGap10,yErrSP_0002_v4_etaGap10);
    return SP_0002_v4_etaGap10;
  }
  
  if (n == 5)
  {
    //  v5{SP}(pt) for 00-02%, eta gap = 0.2:
    const Int_t nPointsSP_0002_v5_etaGap02 = 13;
    Double_t xSP_0002_v5_etaGap02[nPointsSP_0002_v5_etaGap02] = {0.300000,0.500000,0.700000,0.900000,1.100000,1.300000,1.500000,
    1.700000,2.000000,2.550000,3.250000,3.950000,4.650000};
    Double_t ySP_0002_v5_etaGap02[nPointsSP_0002_v5_etaGap02] = {0.000570,0.002922,0.002151,0.005256,0.006287,0.005849,0.009399,
    0.011420,0.012455,0.032134,0.057009,0.020607,0.013551};
    Double_t xErrSP_0002_v5_etaGap02[nPointsSP_0002_v5_etaGap02] = {0.};
    Double_t yErrSP_0002_v5_etaGap02[nPointsSP_0002_v5_etaGap02] = {0.001074,0.001155,0.001433,0.001725,0.002123,0.002608,0.003196,
    0.003930,0.003755,0.004869,0.009719,0.018353,0.031814};
    TGraphErrors *SP_0002_v5_etaGap02 = new TGraphErrors(nPointsSP_0002_v5_etaGap02,xSP_0002_v5_etaGap02,ySP_0002_v5_etaGap02,
							  xErrSP_0002_v5_etaGap02,yErrSP_0002_v5_etaGap02);
    return SP_0002_v5_etaGap02;
  }
}

TGraphErrors* GetFlow02(Int_t n)
{
  // private communication 28.01.11, Ante B. / Raimond

  if (n == 2)
  {
    // v2{2}(pt):
    Double_t xCumulant2nd0002ALICE_v2[] = {0.050000,0.150000,0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,0.850000,0.950000,
    1.050000,1.150000,1.250000,1.350000,1.450000,1.550000,1.650000,1.750000,1.850000,1.950000,
    2.250000,2.750000,3.250000,3.750000,4.250000,4.750000,5.500000,6.500000,7.500000,8.500000,
    9.500000};
    Double_t yCumulant2nd0002ALICE_v2[] = {0.000000,0.000000,0.012173,0.015186,0.018580,0.021114,0.024646,0.027040,0.030269,0.032677,
    0.035332,0.037382,0.039228,0.040614,0.042460,0.044658,0.046246,0.050392,0.051436,0.054669,
    0.057330,0.063439,0.067425,0.060144,0.071260,0.070206,0.000000,0.000000,0.000000,0.000000,
    0.000000};
    Double_t xErrCumulant2nd0002ALICE_v2[31] = {0.};
    Double_t yErrCumulant2nd0002ALICE_v2[] = {0.000000,0.000000,0.000256,0.000259,0.000271,0.000296,0.000322,0.000357,0.000397,0.000438,
    0.000483,0.000529,0.000590,0.000639,0.000713,0.000793,0.000877,0.000976,0.001070,0.001197,
    0.000725,0.001265,0.002069,0.003156,0.004605,0.006543,0.000000,0.000000,0.000000,0.000000,
    0.000000};
    Int_t nPointsCumulant2nd0002ALICE_v2 = sizeof(xCumulant2nd0002ALICE_v2)/sizeof(Double_t);                                      
    TGraphErrors *Cumulant2nd0002ALICE_v2 = new TGraphErrors(nPointsCumulant2nd0002ALICE_v2,xCumulant2nd0002ALICE_v2,yCumulant2nd0002ALICE_v2,
                                                          xErrCumulant2nd0002ALICE_v2,yErrCumulant2nd0002ALICE_v2);
    Cumulant2nd0002ALICE_v2->SetMarkerStyle(kFullCircle);
    Cumulant2nd0002ALICE_v2->SetMarkerColor(kBlue);
    
    return Cumulant2nd0002ALICE_v2;
  }
  
  if (n == 3)
  {
    // v3{2}(pt):
    Double_t xCumulant2nd0002ALICE_v3[] = {0.050000,0.150000,0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,0.850000,0.950000,
    1.050000,1.150000,1.250000,1.350000,1.450000,1.550000,1.650000,1.750000,1.850000,1.950000,
    2.250000,2.750000,3.250000,3.750000,4.250000,4.750000,5.500000,6.500000,7.500000,8.500000,
    9.500000};
    Double_t yCumulant2nd0002ALICE_v3[] = {0.000000,0.000000,0.007696,0.008994,0.010544,0.013269,0.016330,0.019234,0.023465,0.026803,
    0.029906,0.032211,0.035300,0.038158,0.041861,0.046002,0.049382,0.053574,0.055773,0.059420,
    0.069373,0.079922,0.090265,0.103583,0.111358,0.090740,0.000000,0.000000,0.000000,0.000000,
    0.000000};
    Double_t xErrCumulant2nd0002ALICE_v3[31] = {0.};
    Double_t yErrCumulant2nd0002ALICE_v3[] = {0.000000,0.000000,0.000318,0.000317,0.000333,0.000360,0.000392,0.000431,0.000476,0.000523,
    0.000575,0.000637,0.000707,0.000785,0.000878,0.000964,0.001064,0.001175,0.001320,0.001459,
    0.000889,0.001539,0.002530,0.003826,0.005614,0.007892,0.000000,0.000000,0.000000,0.000000,
    0.000000};
    Int_t nPointsCumulant2nd0002ALICE_v3 = sizeof(xCumulant2nd0002ALICE_v3)/sizeof(Double_t);                                      
    TGraphErrors *Cumulant2nd0002ALICE_v3 = new TGraphErrors(nPointsCumulant2nd0002ALICE_v3,xCumulant2nd0002ALICE_v3,yCumulant2nd0002ALICE_v3,
                                                          xErrCumulant2nd0002ALICE_v3,yErrCumulant2nd0002ALICE_v3);
    Cumulant2nd0002ALICE_v3->SetMarkerStyle(kFullTriangleUp);
    Cumulant2nd0002ALICE_v3->SetMarkerSize(1.2);
    Cumulant2nd0002ALICE_v3->SetMarkerColor(kGreen+2); 
    
    return Cumulant2nd0002ALICE_v3;
  }
  
  if (n == 4)
  {
    // v4{2}(pt):
    Double_t xCumulant2nd0002ALICE_v4[] = {0.050000,0.150000,0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,0.850000,0.950000,
    1.050000,1.150000,1.250000,1.350000,1.450000,1.550000,1.650000,1.750000,1.850000,1.950000,
    2.250000,2.750000,3.250000,3.750000,4.250000,4.750000,5.500000,6.500000,7.500000,8.500000,
    9.500000};
    Double_t yCumulant2nd0002ALICE_v4[] = {0.000000,0.000000,0.005710,0.006014,0.004483,0.005453,0.007714,0.006837,0.009721,0.011288,
    0.012531,0.016461,0.016606,0.018587,0.022722,0.025497,0.025832,0.030994,0.030349,0.034730,
    0.045529,0.061153,0.074238,0.079307,0.088885,0.085218,0.000000,0.000000,0.000000,0.000000,
    0.000000};
    Double_t xErrCumulant2nd0002ALICE_v4[31] = {0.};
    Double_t yErrCumulant2nd0002ALICE_v4[] = {0.000000,0.000000,0.000488,0.000493,0.000523,0.000571,0.000609,0.000678,0.000742,0.000805,
    0.000903,0.000985,0.001100,0.001219,0.001352,0.001503,0.001682,0.001847,0.002060,0.002303,
    0.001400,0.002431,0.003974,0.006040,0.008901,0.012343,0.000000,0.000000,0.000000,0.000000,
    0.000000};
    Int_t nPointsCumulant2nd0002ALICE_v4 = sizeof(xCumulant2nd0002ALICE_v4)/sizeof(Double_t);                                      
    TGraphErrors *Cumulant2nd0002ALICE_v4 = new TGraphErrors(nPointsCumulant2nd0002ALICE_v4,xCumulant2nd0002ALICE_v4,yCumulant2nd0002ALICE_v4,
                                                          xErrCumulant2nd0002ALICE_v4,yErrCumulant2nd0002ALICE_v4);
    Cumulant2nd0002ALICE_v4->SetMarkerStyle(kFullSquare);
    Cumulant2nd0002ALICE_v4->SetMarkerColor(kRed);  
    
    return Cumulant2nd0002ALICE_v4;
  }
  
  return 0;
}

/* results up to 5 GeV/c

TGraphErrors* GetFlow05(Int_t n)
{
  // private communication 02.02.11, Ante B. / Raimond

  if (n == 2)
  {
    // v2{2}(pt) for 0-5%:
    Double_t xCumulant2nd0005ALICE_v2[] = {0.050000,0.150000,0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,0.850000,0.950000,
    1.100000,1.300000,1.500000,1.700000,1.900000,2.250000,2.750000,3.250000,3.750000,4.500000};
    Double_t yCumulant2nd0005ALICE_v2[] = {0.000000,0.000000,0.013891,0.017693,0.021693,0.025323,0.029131,0.032443,0.035781,0.038256,
    0.042801,0.047705,0.053229,0.057387,0.062677,0.068815,0.077695,0.082058,0.082511,0.079791};
    Double_t xErrCumulant2nd0005ALICE_v2[20] = {0.};
    Double_t yErrCumulant2nd0005ALICE_v2[] = {0.000000,0.000000,0.000149,0.000150,0.000160,0.000174,0.000191,0.000211,0.000233,0.000257,
    0.000208,0.000254,0.000311,0.000377,0.000464,0.000419,0.000726,0.001180,0.001791,0.002131};
    Int_t nPointsCumulant2nd0005ALICE_v2 = sizeof(xCumulant2nd0005ALICE_v2)/sizeof(Double_t);                                      
    TGraphErrors *Cumulant2nd0005ALICE_v2 = new TGraphErrors(nPointsCumulant2nd0005ALICE_v2,xCumulant2nd0005ALICE_v2,yCumulant2nd0005ALICE_v2,
                                                          xErrCumulant2nd0005ALICE_v2,yErrCumulant2nd0005ALICE_v2);
    Cumulant2nd0005ALICE_v2->SetMarkerStyle(kFullCircle);
    Cumulant2nd0005ALICE_v2->SetMarkerColor(kBlue);    
    
    return Cumulant2nd0005ALICE_v2;
  }
  
  if (n == 3)
  {
    // v3{2}(pt) for 0-5%:
    Double_t xCumulant2nd0005ALICE_v3[] = {0.050000,0.150000,0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,0.850000,0.950000,
    1.100000,1.300000,1.500000,1.700000,1.900000,2.250000,2.750000,3.250000,3.750000,4.500000};
    Double_t yCumulant2nd0005ALICE_v3[] = {0.000000,0.000000,0.007788,0.009472,0.011596,0.014618,0.017540,0.021020,0.024946,0.028004,
    0.032330,0.039491,0.046368,0.053620,0.060662,0.071750,0.086746,0.097857,0.103111,0.104796};
    Double_t xErrCumulant2nd0005ALICE_v3[20] = {0.};
    Double_t yErrCumulant2nd0005ALICE_v3[] = {0.000000,0.000000,0.000194,0.000192,0.000204,0.000221,0.000241,0.000265,0.000293,0.000323,
    0.000266,0.000323,0.000397,0.000486,0.000601,0.000545,0.000947,0.001541,0.002328,0.002777};
    Int_t nPointsCumulant2nd0005ALICE_v3 = sizeof(xCumulant2nd0005ALICE_v3)/sizeof(Double_t);                                      
    TGraphErrors *Cumulant2nd0005ALICE_v3 = new TGraphErrors(nPointsCumulant2nd0005ALICE_v3,xCumulant2nd0005ALICE_v3,yCumulant2nd0005ALICE_v3,
                                                          xErrCumulant2nd0005ALICE_v3,yErrCumulant2nd0005ALICE_v3);
    Cumulant2nd0005ALICE_v3->SetMarkerStyle(kFullCircle);
    Cumulant2nd0005ALICE_v3->SetMarkerColor(kBlue);

    return Cumulant2nd0005ALICE_v3;
  }
  
  if (n == 4)
  {
    // v4{2}(pt) for 0-5%:
    Double_t xCumulant2nd0005ALICE_v4[] = {0.050000,0.150000,0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,0.850000,0.950000,
    1.100000,1.300000,1.500000,1.700000,1.900000,2.250000,2.750000,3.250000,3.750000,4.500000};
    Double_t yCumulant2nd0005ALICE_v4[] = {0.000000,0.000000,0.006121,0.006137,0.005598,0.005956,0.007262,0.007991,0.009159,0.012062,
    0.015085,0.019225,0.024782,0.030092,0.035708,0.046542,0.060077,0.076088,0.082964,0.085405};
    Double_t xErrCumulant2nd0005ALICE_v4[20] = {0.};
    Double_t yErrCumulant2nd0005ALICE_v4[] = {0.000000,0.000000,0.000275,0.000278,0.000294,0.000319,0.000346,0.000380,0.000419,0.000459,
    0.000378,0.000460,0.000570,0.000700,0.000865,0.000789,0.001370,0.002227,0.003370,0.004018};
    Int_t nPointsCumulant2nd0005ALICE_v4 = sizeof(xCumulant2nd0005ALICE_v4)/sizeof(Double_t);                                      
    TGraphErrors *Cumulant2nd0005ALICE_v4 = new TGraphErrors(nPointsCumulant2nd0005ALICE_v4,xCumulant2nd0005ALICE_v4,yCumulant2nd0005ALICE_v4,
                                                          xErrCumulant2nd0005ALICE_v4,yErrCumulant2nd0005ALICE_v4);
    Cumulant2nd0005ALICE_v4->SetMarkerStyle(kFullCircle);
    Cumulant2nd0005ALICE_v4->SetMarkerColor(kBlue); 
    
    return Cumulant2nd0005ALICE_v4;
  }
  
  return 0;
}

TGraphErrors* GetFlow510(Int_t n)
{
  // private communication 02.02.11, Ante B. / Raimond

  if (n == 2)
  {
    // v2{2}(pt) for 5-10%:
    Double_t xCumulant2nd0510ALICE_v2[] = {0.050000,0.150000,0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,0.850000,0.950000,
    1.100000,1.300000,1.500000,1.700000,1.900000,2.250000,2.750000,3.250000,3.750000,4.500000};
    Double_t yCumulant2nd0510ALICE_v2[] = {0.000000,0.000000,0.019872,0.026451,0.032901,0.039085,0.044930,0.050200,0.054887,0.060253,
    0.066587,0.075080,0.083303,0.090298,0.098782,0.109632,0.124486,0.129621,0.132076,0.120697};
    Double_t xErrCumulant2nd0510ALICE_v2[20] = {0.};
    Double_t yErrCumulant2nd0510ALICE_v2[] = {0.000000,0.000000,0.000150,0.000152,0.000163,0.000178,0.000196,0.000215,0.000237,0.000261,
    0.000213,0.000256,0.000313,0.000381,0.000468,0.000423,0.000727,0.001157,0.001741,0.002064};
    Int_t nPointsCumulant2nd0510ALICE_v2 = sizeof(xCumulant2nd0510ALICE_v2)/sizeof(Double_t);                                      
    TGraphErrors *Cumulant2nd0510ALICE_v2 = new TGraphErrors(nPointsCumulant2nd0510ALICE_v2,xCumulant2nd0510ALICE_v2,yCumulant2nd0510ALICE_v2,
                                                          xErrCumulant2nd0510ALICE_v2,yErrCumulant2nd0510ALICE_v2);
    Cumulant2nd0510ALICE_v2->SetMarkerStyle(kOpenCircle);
    Cumulant2nd0510ALICE_v2->SetMarkerColor(kBlue);   
     
    return Cumulant2nd0510ALICE_v2;
  }
  
  if (n == 3)
  {
    // v3{2}(pt) for 5-10%:
    Double_t xCumulant2nd0510ALICE_v3[] = {0.050000,0.150000,0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,0.850000,0.950000,
    1.100000,1.300000,1.500000,1.700000,1.900000,2.250000,2.750000,3.250000,3.750000,4.500000};
    Double_t yCumulant2nd0510ALICE_v3[] = {0.000000,0.000000,0.008596,0.010700,0.013820,0.017524,0.021507,0.024316,0.028491,0.032880,
    0.038741,0.045830,0.052486,0.059560,0.067990,0.081006,0.097402,0.107050,0.111743,0.116434};
    Double_t xErrCumulant2nd0510ALICE_v3[20] = {0.};
    Double_t yErrCumulant2nd0510ALICE_v3[] = {0.000000,0.000000,0.000208,0.000207,0.000218,0.000235,0.000258,0.000284,0.000314,0.000347,
    0.000285,0.000345,0.000426,0.000521,0.000642,0.000586,0.001008,0.001611,0.002421,0.002853};
    Int_t nPointsCumulant2nd0510ALICE_v3 = sizeof(xCumulant2nd0510ALICE_v3)/sizeof(Double_t);                                      
    TGraphErrors *Cumulant2nd0510ALICE_v3 = new TGraphErrors(nPointsCumulant2nd0510ALICE_v3,xCumulant2nd0510ALICE_v3,yCumulant2nd0510ALICE_v3,
                                                          xErrCumulant2nd0510ALICE_v3,yErrCumulant2nd0510ALICE_v3);
    Cumulant2nd0510ALICE_v3->SetMarkerStyle(kOpenCircle);
    Cumulant2nd0510ALICE_v3->SetMarkerColor(kBlue);
    
    return Cumulant2nd0510ALICE_v3;
  }
  
  if (n == 4)
  {
    // v4{2}(pt) for 5-10%:
    Double_t xCumulant2nd0510ALICE_v4[] = {0.050000,0.150000,0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,0.850000,0.950000,
    1.100000,1.300000,1.500000,1.700000,1.900000,2.250000,2.750000,3.250000,3.750000,4.500000};
    Double_t yCumulant2nd0510ALICE_v4[] = {0.000000,0.000000,0.006466,0.006731,0.006344,0.007374,0.008775,0.010324,0.012323,0.014533,
    0.017261,0.022507,0.028776,0.035403,0.041936,0.051491,0.070340,0.080081,0.095077,0.088526};
    Double_t xErrCumulant2nd0510ALICE_v4[20] = {0.};
    Double_t yErrCumulant2nd0510ALICE_v4[] = {0.000000,0.000000,0.000292,0.000295,0.000312,0.000336,0.000366,0.000403,0.000443,0.000485,
    0.000399,0.000486,0.000603,0.000738,0.000914,0.000836,0.001443,0.002303,0.003448,0.004078};
    Int_t nPointsCumulant2nd0510ALICE_v4 = sizeof(xCumulant2nd0510ALICE_v4)/sizeof(Double_t);                                      
    TGraphErrors *Cumulant2nd0510ALICE_v4 = new TGraphErrors(nPointsCumulant2nd0510ALICE_v4,xCumulant2nd0510ALICE_v4,yCumulant2nd0510ALICE_v4,
                                                          xErrCumulant2nd0510ALICE_v4,yErrCumulant2nd0510ALICE_v4);
    Cumulant2nd0510ALICE_v4->SetMarkerStyle(kOpenCircle);
    Cumulant2nd0510ALICE_v4->SetMarkerColor(kBlue);    
    
    return Cumulant2nd0510ALICE_v4;
  }
  
  return 0;
}

TGraphErrors* GetFlow1020(Int_t n)
{
  // private communication 02.02.11, Ante B. / Raimond

  if (n == 2)
  {
    // v2{2}(pt) for 10-20%:
    Double_t xCumulant2nd1020ALICE_v2[] = {0.050000,0.150000,0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,0.850000,0.950000,
    1.100000,1.300000,1.500000,1.700000,1.900000,2.250000,2.750000,3.250000,3.750000,4.500000};
    Double_t yCumulant2nd1020ALICE_v2[] = {0.000000,0.000000,0.027683,0.037083,0.046511,0.055519,0.063979,0.071626,0.078537,0.085975,
    0.095001,0.106979,0.118456,0.129721,0.140641,0.155161,0.173402,0.179870,0.180616,0.168921};
    Double_t xErrCumulant2nd1020ALICE_v2[20] = {0.};
    Double_t yErrCumulant2nd1020ALICE_v2[] = {0.000000,0.000000,0.000121,0.000124,0.000134,0.000147,0.000163,0.000179,0.000198,0.000217,
    0.000177,0.000212,0.000257,0.000311,0.000380,0.000341,0.000569,0.000882,0.001309,0.001537};
    Int_t nPointsCumulant2nd1020ALICE_v2 = sizeof(xCumulant2nd1020ALICE_v2)/sizeof(Double_t);                                      
    TGraphErrors *Cumulant2nd1020ALICE_v2 = new TGraphErrors(nPointsCumulant2nd1020ALICE_v2,xCumulant2nd1020ALICE_v2,yCumulant2nd1020ALICE_v2,
                                                          xErrCumulant2nd1020ALICE_v2,yErrCumulant2nd1020ALICE_v2);
    Cumulant2nd1020ALICE_v2->SetMarkerStyle(kFullSquare);
    Cumulant2nd1020ALICE_v2->SetMarkerColor(kRed);     
    
    return Cumulant2nd1020ALICE_v2;
  }
  
  if (n == 3)
  {
    // v3{2}(pt) for 10-20%:
    Double_t xCumulant2nd1020ALICE_v3[] = {0.050000,0.150000,0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,0.850000,0.950000,
    1.100000,1.300000,1.500000,1.700000,1.900000,2.250000,2.750000,3.250000,3.750000,4.500000};
    Double_t yCumulant2nd1020ALICE_v3[] = {0.000000,0.000000,0.009830,0.012858,0.016111,0.020120,0.023948,0.028349,0.032741,0.037244,
    0.043385,0.051803,0.059374,0.068686,0.076763,0.090151,0.106530,0.117448,0.121383,0.118247};
    Double_t xErrCumulant2nd1020ALICE_v3[20] = {0.};
    Double_t yErrCumulant2nd1020ALICE_v3[] = {0.000000,0.000000,0.000171,0.000170,0.000180,0.000195,0.000215,0.000236,0.000261,0.000287,
    0.000236,0.000288,0.000353,0.000434,0.000536,0.000488,0.000823,0.001277,0.001892,0.002224};
    Int_t nPointsCumulant2nd1020ALICE_v3 = sizeof(xCumulant2nd1020ALICE_v3)/sizeof(Double_t);                                      
    TGraphErrors *Cumulant2nd1020ALICE_v3 = new TGraphErrors(nPointsCumulant2nd1020ALICE_v3,xCumulant2nd1020ALICE_v3,yCumulant2nd1020ALICE_v3,
                                                          xErrCumulant2nd1020ALICE_v3,yErrCumulant2nd1020ALICE_v3);
    Cumulant2nd1020ALICE_v3->SetMarkerStyle(kFullSquare);
    Cumulant2nd1020ALICE_v3->SetMarkerColor(kRed); 
    
    return Cumulant2nd1020ALICE_v3;
  }
  
  if (n == 4)
  {
    // v4{2}(pt) for 10-20%:
    Double_t xCumulant2nd1020ALICE_v4[] = {0.050000,0.150000,0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,0.850000,0.950000,
    1.100000,1.300000,1.500000,1.700000,1.900000,2.250000,2.750000,3.250000,3.750000,4.500000};
    Double_t yCumulant2nd1020ALICE_v4[] = {0.000000,0.000000,0.007423,0.007647,0.008189,0.008592,0.009489,0.011671,0.013739,0.017199,
    0.020084,0.026004,0.031843,0.038388,0.047357,0.057251,0.072433,0.086326,0.094282,0.097432};
    Double_t xErrCumulant2nd1020ALICE_v4[20] = {0.};
    Double_t yErrCumulant2nd1020ALICE_v4[] = {0.000000,0.000000,0.000243,0.000244,0.000257,0.000279,0.000306,0.000335,0.000368,0.000405,
    0.000333,0.000406,0.000502,0.000618,0.000770,0.000701,0.001185,0.001845,0.002730,0.003193};
    Int_t nPointsCumulant2nd1020ALICE_v4 = sizeof(xCumulant2nd1020ALICE_v4)/sizeof(Double_t);                                      
    TGraphErrors *Cumulant2nd1020ALICE_v4 = new TGraphErrors(nPointsCumulant2nd1020ALICE_v4,xCumulant2nd1020ALICE_v4,yCumulant2nd1020ALICE_v4,
                                                          xErrCumulant2nd1020ALICE_v4,yErrCumulant2nd1020ALICE_v4);
    Cumulant2nd1020ALICE_v4->SetMarkerStyle(kFullSquare);
    Cumulant2nd1020ALICE_v4->SetMarkerColor(kRed); 
        
    return Cumulant2nd1020ALICE_v4;
  }
  
  return 0;
}

TGraphErrors* GetFlow2030(Int_t n)
{
  // private communication 02.02.11, Ante B. / Raimond

  if (n == 2)
  {
    Double_t xCumulant2nd2030ALICE_v2[] = {0.050000,0.150000,0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,0.850000,0.950000,
    1.100000,1.300000,1.500000,1.700000,1.900000,2.250000,2.750000,3.250000,3.750000,4.500000};
    Double_t yCumulant2nd2030ALICE_v2[] = {0.000000,0.000000,0.035557,0.048064,0.060768,0.072585,0.083808,0.093772,0.103310,0.112602,
    0.124846,0.140603,0.155345,0.169450,0.183077,0.200173,0.219693,0.225741,0.223318,0.207356};
    Double_t xErrCumulant2nd2030ALICE_v2[20] = {0.};
    Double_t yErrCumulant2nd2030ALICE_v2[] = {0.000000,0.000000,0.000144,0.000147,0.000159,0.000175,0.000194,0.000214,0.000235,0.000259,
    0.000211,0.000254,0.000310,0.000377,0.000464,0.000418,0.000677,0.001027,0.001513,0.001761};
    Int_t nPointsCumulant2nd2030ALICE_v2 = sizeof(xCumulant2nd2030ALICE_v2)/sizeof(Double_t);                                      
    TGraphErrors *Cumulant2nd2030ALICE_v2 = new TGraphErrors(nPointsCumulant2nd2030ALICE_v2,xCumulant2nd2030ALICE_v2,yCumulant2nd2030ALICE_v2,
                                                          xErrCumulant2nd2030ALICE_v2,yErrCumulant2nd2030ALICE_v2);
    Cumulant2nd2030ALICE_v2->SetMarkerStyle(kOpenSquare);
    Cumulant2nd2030ALICE_v2->SetMarkerColor(kRed); 
        
    return Cumulant2nd2030ALICE_v2;
  }
  
  if (n == 3)
  {
    // v3{2}(pt) for 20-30%:
    Double_t xCumulant2nd2030ALICE_v3[] = {0.050000,0.150000,0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,0.850000,0.950000,
    1.100000,1.300000,1.500000,1.700000,1.900000,2.250000,2.750000,3.250000,3.750000,4.500000};
    Double_t yCumulant2nd2030ALICE_v3[] = {0.000000,0.000000,0.011287,0.014937,0.018479,0.022962,0.027854,0.031938,0.037038,0.041443,
    0.049872,0.058720,0.068353,0.077692,0.087859,0.099863,0.116141,0.124486,0.124786,0.119520};
    Double_t xErrCumulant2nd2030ALICE_v3[20] = {0.};
    Double_t yErrCumulant2nd2030ALICE_v3[] = {0.000000,0.000000,0.000215,0.000216,0.000228,0.000248,0.000272,0.000300,0.000332,0.000366,
    0.000301,0.000369,0.000456,0.000561,0.000700,0.000636,0.001038,0.001574,0.002307,0.002676};
    Int_t nPointsCumulant2nd2030ALICE_v3 = sizeof(xCumulant2nd2030ALICE_v3)/sizeof(Double_t);                                      
    TGraphErrors *Cumulant2nd2030ALICE_v3 = new TGraphErrors(nPointsCumulant2nd2030ALICE_v3,xCumulant2nd2030ALICE_v3,yCumulant2nd2030ALICE_v3,
                                                          xErrCumulant2nd2030ALICE_v3,yErrCumulant2nd2030ALICE_v3);
    Cumulant2nd2030ALICE_v3->SetMarkerStyle(kOpenSquare);
    Cumulant2nd2030ALICE_v3->SetMarkerColor(kRed); 
    
    return Cumulant2nd2030ALICE_v3;
  }
  
  if (n == 4)
  {
    // v4{2}(pt) for 20-30%:
    Double_t xCumulant2nd2030ALICE_v4[] = {0.050000,0.150000,0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,0.850000,0.950000,
    1.100000,1.300000,1.500000,1.700000,1.900000,2.250000,2.750000,3.250000,3.750000,4.500000};
    Double_t yCumulant2nd2030ALICE_v4[] = {0.000000,0.000000,0.008435,0.009824,0.009387,0.010408,0.011458,0.013770,0.016335,0.019165,
    0.023770,0.030071,0.036231,0.043950,0.051803,0.065222,0.082155,0.091864,0.105061,0.105167};
    Double_t xErrCumulant2nd2030ALICE_v4[20] = {0.};
    Double_t yErrCumulant2nd2030ALICE_v4[] = {0.000000,0.000000,0.000308,0.000308,0.000328,0.000355,0.000391,0.000430,0.000473,0.000524,
    0.000430,0.000524,0.000652,0.000807,0.001006,0.000919,0.001502,0.002277,0.003339,0.003871};
    Int_t nPointsCumulant2nd2030ALICE_v4 = sizeof(xCumulant2nd2030ALICE_v4)/sizeof(Double_t);                                      
    TGraphErrors *Cumulant2nd2030ALICE_v4 = new TGraphErrors(nPointsCumulant2nd2030ALICE_v4,xCumulant2nd2030ALICE_v4,yCumulant2nd2030ALICE_v4,
                                                          xErrCumulant2nd2030ALICE_v4,yErrCumulant2nd2030ALICE_v4);
    Cumulant2nd2030ALICE_v4->SetMarkerStyle(kOpenSquare);
    Cumulant2nd2030ALICE_v4->SetMarkerColor(kRed); 
        
    return Cumulant2nd2030ALICE_v4;
  }
  
  return 0;
}

TGraphErrors* GetFlow3040(Int_t n)
{
  // private communication 02.02.11, Ante B. / Raimond

  if (n == 2)
  {
    // v2{2}(pt) for 30-40%:
    Double_t xCumulant2nd3040ALICE_v2[] = {0.050000,0.150000,0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,0.850000,0.950000,
    1.100000,1.300000,1.500000,1.700000,1.900000,2.250000,2.750000,3.250000,3.750000,4.500000};
    Double_t yCumulant2nd3040ALICE_v2[] = {0.000000,0.000000,0.040797,0.055427,0.070543,0.084024,0.096562,0.109411,0.119853,0.130964,
    0.145377,0.163806,0.179825,0.196178,0.210377,0.226556,0.245686,0.247898,0.240058,0.225011};
    Double_t xErrCumulant2nd3040ALICE_v2[20] = {0.};
    Double_t yErrCumulant2nd3040ALICE_v2[] = {0.000000,0.000000,0.000177,0.000182,0.000197,0.000216,0.000239,0.000265,0.000293,0.000325,
    0.000266,0.000321,0.000395,0.000486,0.000603,0.000536,0.000840,0.001258,0.001843,0.002118};
    Int_t nPointsCumulant2nd3040ALICE_v2 = sizeof(xCumulant2nd3040ALICE_v2)/sizeof(Double_t);                                      
    TGraphErrors *Cumulant2nd3040ALICE_v2 = new TGraphErrors(nPointsCumulant2nd3040ALICE_v2,xCumulant2nd3040ALICE_v2,yCumulant2nd3040ALICE_v2,
                                                          xErrCumulant2nd3040ALICE_v2,yErrCumulant2nd3040ALICE_v2);
    Cumulant2nd3040ALICE_v2->SetMarkerStyle(kFullTriangleUp);
    Cumulant2nd3040ALICE_v2->SetMarkerSize(1.2);
    Cumulant2nd3040ALICE_v2->SetMarkerColor(kGreen+2);     
    
    return Cumulant2nd3040ALICE_v2;
  }
  
  if (n == 3)
  {
    // v3{2}(pt) for 30-40%:
    Double_t xCumulant2nd3040ALICE_v3[] = {0.050000,0.150000,0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,0.850000,0.950000,
    1.100000,1.300000,1.500000,1.700000,1.900000,2.250000,2.750000,3.250000,3.750000,4.500000};
    Double_t yCumulant2nd3040ALICE_v3[] = {0.000000,0.000000,0.012226,0.016391,0.020792,0.026208,0.030380,0.035710,0.041025,0.047062,
    0.053075,0.064201,0.074116,0.085314,0.094391,0.106819,0.124012,0.129388,0.134315,0.132330};
    Double_t xErrCumulant2nd3040ALICE_v3[20] = {0.};
    Double_t yErrCumulant2nd3040ALICE_v3[] = {0.000000,0.000000,0.000284,0.000286,0.000303,0.000329,0.000364,0.000403,0.000443,0.000492,
    0.000408,0.000500,0.000627,0.000778,0.000973,0.000874,0.001375,0.002050,0.002992,0.003438};
    Int_t nPointsCumulant2nd3040ALICE_v3 = sizeof(xCumulant2nd3040ALICE_v3)/sizeof(Double_t);                                      
    TGraphErrors *Cumulant2nd3040ALICE_v3 = new TGraphErrors(nPointsCumulant2nd3040ALICE_v3,xCumulant2nd3040ALICE_v3,yCumulant2nd3040ALICE_v3,
                                                          xErrCumulant2nd3040ALICE_v3,yErrCumulant2nd3040ALICE_v3);
    Cumulant2nd3040ALICE_v3->SetMarkerStyle(kFullTriangleUp);
    Cumulant2nd3040ALICE_v3->SetMarkerSize(1.2);
    Cumulant2nd3040ALICE_v3->SetMarkerColor(kGreen+2); 
    
    return Cumulant2nd3040ALICE_v3;
  }
  
  if (n == 4)
  {
    // v4{2}(pt) for 30-40%:
    Double_t xCumulant2nd3040ALICE_v4[] = {0.050000,0.150000,0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,0.850000,0.950000,
    1.100000,1.300000,1.500000,1.700000,1.900000,2.250000,2.750000,3.250000,3.750000,4.500000};
    Double_t yCumulant2nd3040ALICE_v4[] = {0.000000,0.000000,0.010718,0.011490,0.010994,0.012527,0.013748,0.016425,0.018857,0.021622,
    0.026853,0.034636,0.042651,0.049892,0.057795,0.070865,0.088486,0.101656,0.113886,0.118202};
    Double_t xErrCumulant2nd3040ALICE_v4[20] = {0.};
    Double_t yErrCumulant2nd3040ALICE_v4[] = {0.000000,0.000000,0.000401,0.000406,0.000433,0.000472,0.000521,0.000575,0.000634,0.000704,
    0.000580,0.000714,0.000890,0.001114,0.001398,0.001253,0.001974,0.002945,0.004290,0.004909};
    Int_t nPointsCumulant2nd3040ALICE_v4 = sizeof(xCumulant2nd3040ALICE_v4)/sizeof(Double_t);                                      
    TGraphErrors *Cumulant2nd3040ALICE_v4 = new TGraphErrors(nPointsCumulant2nd3040ALICE_v4,xCumulant2nd3040ALICE_v4,yCumulant2nd3040ALICE_v4,
                                                          xErrCumulant2nd3040ALICE_v4,yErrCumulant2nd3040ALICE_v4);
    Cumulant2nd3040ALICE_v4->SetMarkerStyle(kFullTriangleUp);
    Cumulant2nd3040ALICE_v4->SetMarkerSize(1.2);
    Cumulant2nd3040ALICE_v4->SetMarkerColor(kGreen+2); 
        
    return Cumulant2nd3040ALICE_v4;
  }
  
  return 0;
}

TGraphErrors* GetFlow4050(Int_t n)
{
  // private communication 02.02.11, Ante B. / Raimond

  if (n == 2)
  {
    // v2{2}(pt) for 40-50%:
    Double_t xCumulant2nd4050ALICE_v2[] = {0.050000,0.150000,0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,0.850000,0.950000,
    1.100000,1.300000,1.500000,1.700000,1.900000,2.250000,2.750000,3.250000,3.750000,4.500000};
    Double_t yCumulant2nd4050ALICE_v2[] = {0.000000,0.000000,0.043593,0.059522,0.075558,0.090681,0.104453,0.117913,0.129349,0.141743,
    0.157076,0.176121,0.193651,0.208844,0.222402,0.238407,0.252512,0.253592,0.243571,0.233018};
    Double_t xErrCumulant2nd4050ALICE_v2[20] = {0.};
    Double_t yErrCumulant2nd4050ALICE_v2[] = {0.000000,0.000000,0.000234,0.000241,0.000261,0.000288,0.000322,0.000357,0.000395,0.000438,
    0.000362,0.000447,0.000555,0.000685,0.000846,0.000733,0.001121,0.001665,0.002433,0.002768};
    Int_t nPointsCumulant2nd4050ALICE_v2 = sizeof(xCumulant2nd4050ALICE_v2)/sizeof(Double_t);                                      
    TGraphErrors *Cumulant2nd4050ALICE_v2 = new TGraphErrors(nPointsCumulant2nd4050ALICE_v2,xCumulant2nd4050ALICE_v2,yCumulant2nd4050ALICE_v2,
                                                          xErrCumulant2nd4050ALICE_v2,yErrCumulant2nd4050ALICE_v2);
    Cumulant2nd4050ALICE_v2->SetMarkerStyle(kOpenTriangleUp);
    Cumulant2nd4050ALICE_v2->SetMarkerSize(1.2);
    Cumulant2nd4050ALICE_v2->SetMarkerColor(kGreen+2);      
    
    return Cumulant2nd4050ALICE_v2;
  }
  
  return 0;
}
*/

// results up to high pT

TGraphErrors* GetFlow05(Int_t n)
{
  // private communication 09.03.11, Ante B. / Raimond

  if (n == 2)
  {
    // v2{2}(pt) for 0-5%:
    Double_t xCumulant2nd0005ALICE[] = {0.050000,0.150000,0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,0.850000,0.950000,
    1.050000,1.150000,1.250000,1.350000,1.450000,1.550000,1.650000,1.750000,1.850000,1.950000,
    2.100000,2.300000,2.500000,2.700000,2.900000,3.200000,3.600000,4.000000,4.400000,4.800000,
    5.500000,6.500000,7.500000,8.500000,9.500000,11.250000,13.750000,16.250000,18.750000,60.000000};
    Double_t yCumulant2nd0005ALICE[] = {0.000000,0.000000,0.018169,0.020988,0.024899,0.028309,0.031974,0.035188,0.038218,0.041202,
    0.043910,0.046560,0.048183,0.050640,0.052874,0.055334,0.056948,0.059427,0.061540,0.063218,
    0.065452,0.069222,0.072578,0.074723,0.077749,0.077178,0.080514,0.075325,0.077692,0.079710,
    0.073280,0.063849,0.068274,0.066045,0.071496,0.104352,0.091646,0.050220,0.124185,0.088535};
    Double_t xErrCumulant2nd0005ALICE[40] = {0.};
    Double_t yErrCumulant2nd0005ALICE[] = {0.000000,0.000000,0.000190,0.000191,0.000199,0.000212,0.000230,0.000248,0.000268,0.000293,
    0.000319,0.000346,0.000379,0.000418,0.000459,0.000502,0.000554,0.000614,0.000674,0.000749,
    0.000620,0.000769,0.000958,0.001182,0.001446,0.001331,0.001858,0.002552,0.003453,0.004606,
    0.004289,0.007006,0.010046,0.013853,0.017709,0.016630,0.025728,0.036763,0.045056,0.029011};
    Int_t nPointsCumulant2nd0005ALICE = sizeof(xCumulant2nd0005ALICE)/sizeof(Double_t);                                      
    TGraphErrors *Cumulant2nd0005ALICE = new TGraphErrors(nPointsCumulant2nd0005ALICE,xCumulant2nd0005ALICE,yCumulant2nd0005ALICE,
                                                          xErrCumulant2nd0005ALICE,yErrCumulant2nd0005ALICE);
    Cumulant2nd0005ALICE->SetMarkerStyle(kFullCircle);
    Cumulant2nd0005ALICE->SetMarkerColor(kBlue);
    
    return Cumulant2nd0005ALICE;
  }
  
  if (n == 3)
  {
  }
  
  if (n == 4)
  {
  }
  
  return 0;
}

TGraphErrors* GetFlow05_Rap10(Int_t n)
{
  // private communication 17.05.11, Ante B. / should correspond to machcone paper draft 

  if (n == 2)
  {
    //  v2{SP}(pt) for 00-05%:
    const Int_t nPointsSP_0005ALICE_v2_etaGap10 = 17;
    Double_t xSP_0005ALICE_v2_etaGap10[nPointsSP_0005ALICE_v2_etaGap10] = {0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,0.900000,1.100000,
        1.300000,1.500000,1.700000,1.900000,2.200000,2.600000,3.000000,3.650000,4.550000};
    Double_t ySP_0005ALICE_v2_etaGap10[nPointsSP_0005ALICE_v2_etaGap10] = {0.013672,0.015745,0.019944,0.024169,0.026921,0.030296,0.034916,0.038472,
        0.043103,0.046626,0.052386,0.057132,0.060070,0.068419,0.061459,0.056816,0.050311};
    Double_t xErrSP_0005ALICE_v2_etaGap10[nPointsSP_0005ALICE_v2_etaGap10] = {0.};
    Double_t yErrSP_0005ALICE_v2_etaGap10[nPointsSP_0005ALICE_v2_etaGap10] = {0.000813,0.000804,0.000852,0.000930,0.001029,0.001140,0.000939,0.001155,
        0.001410,0.001736,0.002131,0.002620,0.002502,0.003759,0.005615,0.006617,0.014242};
    TGraphErrors *GrSP_0005ALICE_v2_etaGap10 = new TGraphErrors(nPointsSP_0005ALICE_v2_etaGap10,xSP_0005ALICE_v2_etaGap10,ySP_0005ALICE_v2_etaGap10,
                                                              xErrSP_0005ALICE_v2_etaGap10,yErrSP_0005ALICE_v2_etaGap10);
							      
    return GrSP_0005ALICE_v2_etaGap10;
  }
  
  if (n == 3)
  {
    //  v3{SP}(pt) for 00-05%:
    const Int_t nPointsSP_0005ALICE_v3_etaGap10 = 16;
    Double_t xSP_0005ALICE_v3_etaGap10[nPointsSP_0005ALICE_v3_etaGap10] = {0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,0.900000,1.100000,
        1.300000,1.500000,1.700000,1.900000,2.300000,2.900000,3.650000,4.550000};
    Double_t ySP_0005ALICE_v3_etaGap10[nPointsSP_0005ALICE_v3_etaGap10] = {0.006303,0.009800,0.011143,0.014246,0.017628,0.019437,0.028412,0.030580,
        0.038730,0.045653,0.052469,0.062303,0.071522,0.082223,0.083373,0.076951};
    Double_t xErrSP_0005ALICE_v3_etaGap10[nPointsSP_0005ALICE_v3_etaGap10] = {0.};
    Double_t yErrSP_0005ALICE_v3_etaGap10[nPointsSP_0005ALICE_v3_etaGap10] = {0.001012,0.000996,0.001062,0.001158,0.001277,0.001415,0.001158,0.001422,
        0.001734,0.002124,0.002610,0.003206,0.002724,0.004977,0.008123,0.017180};
    TGraphErrors *GrSP_0005ALICE_v3_etaGap10 = new TGraphErrors(nPointsSP_0005ALICE_v3_etaGap10,xSP_0005ALICE_v3_etaGap10,ySP_0005ALICE_v3_etaGap10,
                                                              xErrSP_0005ALICE_v3_etaGap10,yErrSP_0005ALICE_v3_etaGap10);

							      
    return GrSP_0005ALICE_v3_etaGap10;
  }
  
  if (n == 4)
  {
    //  v4{SP}(pt) for 00-05%:
    const Int_t nPointsSP_0005ALICE_v4_etaGap10 = 11;
    Double_t xSP_0005ALICE_v4_etaGap10[nPointsSP_0005ALICE_v4_etaGap10] = {0.300000,0.500000,0.700000,0.950000,1.250000,1.550000,1.850000,2.300000,2.900000,3.650000,4.550000};
    Double_t ySP_0005ALICE_v4_etaGap10[nPointsSP_0005ALICE_v4_etaGap10] = {0.002042,0.002556,0.009693,0.013286,0.016780,0.027865,0.031797,0.051101,0.060164,
        0.095985,0.094607};
    Double_t xErrSP_0005ALICE_v4_etaGap10[nPointsSP_0005ALICE_v4_etaGap10] = {0.};
    Double_t yErrSP_0005ALICE_v4_etaGap10[nPointsSP_0005ALICE_v4_etaGap10] = {0.001460,0.001624,0.001930,0.002021,0.002737,0.003717,0.005042,0.005564,0.010160,
        0.016472,0.035083};
    TGraphErrors *GrSP_0005ALICE_v4_etaGap10 = new TGraphErrors(nPointsSP_0005ALICE_v4_etaGap10,xSP_0005ALICE_v4_etaGap10,ySP_0005ALICE_v4_etaGap10,
                                                              xErrSP_0005ALICE_v4_etaGap10,yErrSP_0005ALICE_v4_etaGap10);
   
    return GrSP_0005ALICE_v4_etaGap10;
  }
  
  if (n == 5)
  {
    //  v5{SP}(pt) for 00-05%:
    const Int_t nPointsSP_0005ALICE_v5_etaGap10 = 12;
    Double_t xSP_0005ALICE_v5_etaGap10[nPointsSP_0005ALICE_v5_etaGap10] = {0.300000,0.500000,0.700000,0.900000,1.100000,1.300000,1.600000,2.000000,2.400000,
        2.800000,3.500000,4.500000};
    Double_t ySP_0005ALICE_v5_etaGap10[nPointsSP_0005ALICE_v5_etaGap10] = {0.002016,0.003409,0.004029,0.002665,0.002765,0.003042,0.013241,0.015430,0.031845,
        0.031373,0.068504,0.017964};
    Double_t xErrSP_0005ALICE_v5_etaGap10[nPointsSP_0005ALICE_v5_etaGap10] = {0.};
    Double_t yErrSP_0005ALICE_v5_etaGap10[nPointsSP_0005ALICE_v5_etaGap10] = {0.001260,0.001386,0.001696,0.002101,0.002560,0.003119,0.002970,0.004472,0.006802,
        0.010073,0.011899,0.027756};
    TGraphErrors *GrSP_0005ALICE_v5_etaGap10 = new TGraphErrors(nPointsSP_0005ALICE_v5_etaGap10,xSP_0005ALICE_v5_etaGap10,ySP_0005ALICE_v5_etaGap10,
                                                              xErrSP_0005ALICE_v5_etaGap10,yErrSP_0005ALICE_v5_etaGap10);
    
    return GrSP_0005ALICE_v5_etaGap10;
  }
  
  return 0;
}

TGraphErrors* GetFlow510(Int_t n)
{
  // private communication 09.03.11, Ante B. / Raimond

  if (n == 2)
  {
    // v2{2}(pt) for 5-10%:
    Double_t xCumulant2nd0510ALICE[] = {0.050000,0.150000,0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,0.850000,0.950000,
    1.050000,1.150000,1.250000,1.350000,1.450000,1.550000,1.650000,1.750000,1.850000,1.950000,
    2.100000,2.300000,2.500000,2.700000,2.900000,3.200000,3.600000,4.000000,4.400000,4.800000,
    5.500000,6.500000,7.500000,8.500000,9.500000,11.250000,13.750000,16.250000,18.750000,60.000000};
    Double_t yCumulant2nd0510ALICE[] = {0.000000,0.000000,0.022354,0.028064,0.034341,0.040769,0.046265,0.051160,0.056712,0.061354,
    0.066290,0.070340,0.074453,0.078444,0.082367,0.085785,0.088086,0.092676,0.096169,0.100366,
    0.104227,0.109710,0.117032,0.121784,0.122227,0.126537,0.127966,0.117790,0.117884,0.111436,
    0.104945,0.097640,0.086481,0.091663,0.091404,0.080132,0.015864,0.101500,0.033347,0.205130};
    Double_t xErrCumulant2nd0510ALICE[40] = {0.};
    Double_t yErrCumulant2nd0510ALICE[] = {0.000000,0.000000,0.000173,0.000176,0.000186,0.000199,0.000219,0.000236,0.000260,0.000283,
    0.000312,0.000340,0.000373,0.000410,0.000454,0.000501,0.000552,0.000610,0.000675,0.000753,
    0.000620,0.000774,0.000961,0.001183,0.001431,0.001309,0.001814,0.002481,0.003342,0.004379,
    0.004122,0.006716,0.009851,0.013626,0.017517,0.015790,0.025680,0.038041,0.054144,0.038987};
    Int_t nPointsCumulant2nd0510ALICE = sizeof(xCumulant2nd0510ALICE)/sizeof(Double_t);                                      
    TGraphErrors *Cumulant2nd0510ALICE = new TGraphErrors(nPointsCumulant2nd0510ALICE,xCumulant2nd0510ALICE,yCumulant2nd0510ALICE,
                                                          xErrCumulant2nd0510ALICE,yErrCumulant2nd0510ALICE);
    Cumulant2nd0510ALICE->SetMarkerStyle(kOpenCircle);
    Cumulant2nd0510ALICE->SetMarkerColor(kBlue);
    
    return Cumulant2nd0510ALICE;
  }
  
  if (n == 3)
  {
  }
  
  if (n == 4)
  {
  }
  
  return 0;
}

TGraphErrors* GetFlow510_Rap10(Int_t n)
{
  // private communication 18.05.11, Ante B.

  if (n == 2)
  {
    //  v2{SP}(pt) for 05-10%, rapidity gap = 1.0:
    const Int_t nPointsSP_0510ALICE_v2_etaGap10 = 19;
    Double_t xSP_0510ALICE_v2_etaGap10[nPointsSP_0510ALICE_v2_etaGap10] = {0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,0.850000,0.950000,
1.100000,1.300000,1.500000,1.700000,1.900000,2.250000,2.750000,3.250000,3.750000,4.250000,
4.750000};
    Double_t ySP_0510ALICE_v2_etaGap10[nPointsSP_0510ALICE_v2_etaGap10] = {0.018634,0.025538,0.032157,0.038272,0.044020,0.049252,0.055627,0.059272,
0.066516,0.073891,0.080945,0.090386,0.094505,0.106393,0.120303,0.122586,0.121731,0.107343,
0.104059};
    Double_t xErrSP_0510ALICE_v2_etaGap10[nPointsSP_0510ALICE_v2_etaGap10] = {0.};
    Double_t yErrSP_0510ALICE_v2_etaGap10[nPointsSP_0510ALICE_v2_etaGap10] = {0.000419,0.000411,0.000435,0.000474,0.000523,0.000580,0.000643,0.000712,
0.000586,0.000717,0.000880,0.001082,0.001333,0.001189,0.001969,0.003199,0.005011,0.007602,
0.010906};
    TGraphErrors *GrSP_0510ALICE_v2_etaGap10 = new TGraphErrors(nPointsSP_0510ALICE_v2_etaGap10,xSP_0510ALICE_v2_etaGap10,
                                                              ySP_0510ALICE_v2_etaGap10,xErrSP_0510ALICE_v2_etaGap10,
                                                              yErrSP_0510ALICE_v2_etaGap10);

    return GrSP_0510ALICE_v2_etaGap10;
  }
  
  return 0;
}

TGraphErrors* GetFlow1020(Int_t n)
{
  // private communication 09.03.11, Ante B. / Raimond

  if (n == 2)
  {
    // v2{2}(pt) for 10-20%:
    Double_t xCumulant2nd1020ALICE[] = {0.050000,0.150000,0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,0.850000,0.950000,
    1.050000,1.150000,1.250000,1.350000,1.450000,1.550000,1.650000,1.750000,1.850000,1.950000,
    2.100000,2.300000,2.500000,2.700000,2.900000,3.200000,3.600000,4.000000,4.400000,4.800000,
    5.500000,6.500000,7.500000,8.500000,9.500000,11.250000,13.750000,16.250000,18.750000,60.000000};
    Double_t yCumulant2nd1020ALICE[] = {0.000000,0.000000,0.028508,0.037698,0.046914,0.055750,0.063948,0.072298,0.079396,0.086829,
    0.093382,0.099676,0.105412,0.110792,0.116715,0.122327,0.127079,0.132302,0.137183,0.142346,
    0.149059,0.156729,0.164672,0.171840,0.175346,0.176824,0.178878,0.174979,0.169229,0.153109,
    0.139676,0.125686,0.120554,0.096537,0.084736,0.118152,0.116079,0.060032,0.093764,0.206506};
    Double_t xErrCumulant2nd1020ALICE[40] = {0.};
    Double_t yErrCumulant2nd1020ALICE[] = {0.000000,0.000000,0.000134,0.000136,0.000145,0.000158,0.000173,0.000190,0.000209,0.000229,
    0.000252,0.000275,0.000301,0.000331,0.000365,0.000401,0.000442,0.000488,0.000541,0.000599,
    0.000496,0.000614,0.000757,0.000915,0.001095,0.000989,0.001360,0.001847,0.002477,0.003234,
    0.003012,0.004867,0.007123,0.009774,0.012740,0.011682,0.019629,0.028568,0.039931,0.024776};
    Int_t nPointsCumulant2nd1020ALICE = sizeof(xCumulant2nd1020ALICE)/sizeof(Double_t);                                      
    TGraphErrors *Cumulant2nd1020ALICE = new TGraphErrors(nPointsCumulant2nd1020ALICE,xCumulant2nd1020ALICE,yCumulant2nd1020ALICE,
                                                          xErrCumulant2nd1020ALICE,yErrCumulant2nd1020ALICE);
    Cumulant2nd1020ALICE->SetMarkerStyle(kFullSquare);
    Cumulant2nd1020ALICE->SetMarkerColor(kRed); 
    
    return Cumulant2nd1020ALICE;
  }
  
  if (n == 3)
  {
  }
  
  if (n == 4)
  {
  }
  
  return 0;
}

TGraphErrors* GetFlow1020_Rap10(Int_t n)
{
  // private communication 18.05.11, Ante B.

  if (n == 2)
  {
   //  v2{SP}(pt) for 10-20%, rapidity gap = 1.0:
    const Int_t nPointsSP_1020ALICE_v2_etaGap10 = 19;
    Double_t xSP_1020ALICE_v2_etaGap10[nPointsSP_1020ALICE_v2_etaGap10] = {0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,0.850000,0.950000,
1.100000,1.300000,1.500000,1.700000,1.900000,2.250000,2.750000,3.250000,3.750000,4.250000,
4.750000};
    Double_t ySP_1020ALICE_v2_etaGap10[nPointsSP_1020ALICE_v2_etaGap10] = {0.026592,0.036955,0.046103,0.055537,0.063461,0.070993,0.078751,0.085723,
0.094701,0.105631,0.117906,0.128147,0.138505,0.153494,0.166651,0.172691,0.177337,0.155068,
0.131586};
    Double_t xErrSP_1020ALICE_v2_etaGap10[nPointsSP_1020ALICE_v2_etaGap10] = {0.};
    Double_t yErrSP_1020ALICE_v2_etaGap10[nPointsSP_1020ALICE_v2_etaGap10] = {0.000302,0.000296,0.000314,0.000342,0.000377,0.000418,0.000465,0.000515,
0.000423,0.000517,0.000634,0.000779,0.000959,0.000856,0.001406,0.002266,0.003528,0.005281,
0.007561};
    TGraphErrors *GrSP_1020ALICE_v2_etaGap10 = new TGraphErrors(nPointsSP_1020ALICE_v2_etaGap10,xSP_1020ALICE_v2_etaGap10,
                                                              ySP_1020ALICE_v2_etaGap10,xErrSP_1020ALICE_v2_etaGap10,
                                                              yErrSP_1020ALICE_v2_etaGap10);

    return GrSP_1020ALICE_v2_etaGap10;
  }
  
  return 0;
}
    
TGraphErrors* GetFlow2030(Int_t n)
{
  // private communication 09.03.11, Ante B. / Raimond

  if (n == 2)
  {
    // v2{2}(pt) for 20-30%:
    Double_t xCumulant2nd2030ALICE[] = {0.050000,0.150000,0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,0.850000,0.950000,
    1.050000,1.150000,1.250000,1.350000,1.450000,1.550000,1.650000,1.750000,1.850000,1.950000,
    2.100000,2.300000,2.500000,2.700000,2.900000,3.200000,3.600000,4.000000,4.400000,4.800000,
    5.500000,6.500000,7.500000,8.500000,9.500000,11.250000,13.750000,16.250000,18.750000,60.000000};
    Double_t yCumulant2nd2030ALICE[] = {0.000000,0.000000,0.035366,0.047465,0.060083,0.072090,0.083418,0.093576,0.103842,0.113110,
    0.122193,0.130168,0.138158,0.145627,0.152851,0.159129,0.166073,0.173144,0.178698,0.186188,
    0.192045,0.202199,0.210535,0.215004,0.220457,0.223339,0.224050,0.211567,0.203955,0.189716,
    0.165994,0.147185,0.131953,0.139331,0.151293,0.127406,0.153764,0.089628,0.161247,0.511418};
    Double_t xErrCumulant2nd2030ALICE[40] = {0.};
    Double_t yErrCumulant2nd2030ALICE[] = {0.000000,0.000000,0.000155,0.000157,0.000169,0.000184,0.000202,0.000222,0.000244,0.000269,
    0.000296,0.000325,0.000357,0.000394,0.000435,0.000481,0.000532,0.000589,0.000655,0.000731,
    0.000605,0.000743,0.000904,0.001081,0.001277,0.001145,0.001568,0.002119,0.002806,0.003635,
    0.003383,0.005346,0.007935,0.010739,0.014682,0.013434,0.021531,0.032352,0.040396,0.028472};
    Int_t nPointsCumulant2nd2030ALICE = sizeof(xCumulant2nd2030ALICE)/sizeof(Double_t);                                      
    TGraphErrors *Cumulant2nd2030ALICE = new TGraphErrors(nPointsCumulant2nd2030ALICE,xCumulant2nd2030ALICE,yCumulant2nd2030ALICE,
                                                          xErrCumulant2nd2030ALICE,yErrCumulant2nd2030ALICE);
    Cumulant2nd2030ALICE->SetMarkerStyle(kOpenSquare);
    Cumulant2nd2030ALICE->SetMarkerColor(kRed); 
    
    return Cumulant2nd2030ALICE;
  }
  
  if (n == 3)
  {
  }
  
  if (n == 4)
  {
  }
  
  return 0;
}

TGraphErrors* GetFlow3040(Int_t n)
{
  // private communication 09.03.11, Ante B. / Raimond

  if (n == 2)
  {
    // v2{2}(pt) for 30-40%:
    Double_t xCumulant2nd3040ALICE[] = {0.050000,0.150000,0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,0.850000,0.950000,
    1.050000,1.150000,1.250000,1.350000,1.450000,1.550000,1.650000,1.750000,1.850000,1.950000,
    2.100000,2.300000,2.500000,2.700000,2.900000,3.200000,3.600000,4.000000,4.400000,4.800000,
    5.500000,6.500000,7.500000,8.500000,9.500000,11.250000,13.750000,16.250000,18.750000,60.000000};
    Double_t yCumulant2nd3040ALICE[] = {0.000000,0.000000,0.039584,0.054196,0.069291,0.082976,0.095622,0.108940,0.119830,0.131587,
    0.141091,0.150899,0.160573,0.168676,0.176301,0.183823,0.191672,0.199506,0.206854,0.212830,
    0.219526,0.229376,0.236513,0.240863,0.245961,0.245891,0.242608,0.234302,0.219580,0.212848,
    0.194666,0.190184,0.171036,0.159173,0.156932,0.141324,0.132809,0.182683,0.023272,0.032825};
    Double_t xErrCumulant2nd3040ALICE[40] = {0.};
    Double_t yErrCumulant2nd3040ALICE[] = {0.000000,0.000000,0.000189,0.000192,0.000205,0.000224,0.000247,0.000273,0.000302,0.000334,
    0.000369,0.000407,0.000448,0.000496,0.000552,0.000612,0.000682,0.000758,0.000844,0.000941,
    0.000774,0.000937,0.001123,0.001329,0.001565,0.001397,0.001911,0.002549,0.003378,0.004306,
    0.003987,0.006353,0.009128,0.013032,0.016891,0.015806,0.025150,0.035119,0.044487,0.050083};
    Int_t nPointsCumulant2nd3040ALICE = sizeof(xCumulant2nd3040ALICE)/sizeof(Double_t);                                      
    TGraphErrors *Cumulant2nd3040ALICE = new TGraphErrors(nPointsCumulant2nd3040ALICE,xCumulant2nd3040ALICE,yCumulant2nd3040ALICE,
                                                          xErrCumulant2nd3040ALICE,yErrCumulant2nd3040ALICE);
    Cumulant2nd3040ALICE->SetMarkerStyle(kFullTriangleUp);
    Cumulant2nd3040ALICE->SetMarkerColor(kGreen+2);
    
    return Cumulant2nd3040ALICE;
  }
  
  if (n == 3)
  {
  }
  
  if (n == 4)
  {
  }
  
  return 0;
}

TGraphErrors* GetFlow4050(Int_t n)
{
  // private communication 09.03.11, Ante B. / Raimond

  if (n == 2)
  {
    // v2{2}(pt) for 40-50%:
    Double_t xCumulant2nd4050ALICE[] = {0.050000,0.150000,0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,0.850000,0.950000,
    1.050000,1.150000,1.250000,1.350000,1.450000,1.550000,1.650000,1.750000,1.850000,1.950000,
    2.100000,2.300000,2.500000,2.700000,2.900000,3.200000,3.600000,4.000000,4.400000,4.800000,
    5.500000,6.500000,7.500000,8.500000,9.500000,11.250000,13.750000,16.250000,18.750000,60.000000};
    Double_t yCumulant2nd4050ALICE[] = {0.000000,0.000000,0.041872,0.058090,0.074444,0.089181,0.103780,0.117279,0.129769,0.142051,
    0.153185,0.163147,0.173309,0.181668,0.190998,0.197703,0.205011,0.211063,0.219587,0.223287,
    0.231163,0.239606,0.246533,0.251457,0.250034,0.252989,0.240823,0.236489,0.230268,0.204321,
    0.213476,0.200247,0.167065,0.190655,0.173573,0.166173,0.153232,0.173112,-0.123540,0.211999};
    Double_t xErrCumulant2nd4050ALICE[40] = {0.};
    Double_t yErrCumulant2nd4050ALICE[] = {0.000000,0.000000,0.000248,0.000251,0.000270,0.000296,0.000328,0.000363,0.000403,0.000447,
    0.000500,0.000556,0.000617,0.000687,0.000770,0.000857,0.000952,0.001057,0.001176,0.001297,
    0.001054,0.001265,0.001498,0.001767,0.002076,0.001849,0.002527,0.003358,0.004372,0.005573,
    0.005091,0.007981,0.011746,0.015732,0.021883,0.019866,0.032443,0.046521,0.051631,0.083717};
    Int_t nPointsCumulant2nd4050ALICE = sizeof(xCumulant2nd4050ALICE)/sizeof(Double_t);                                      
    TGraphErrors *Cumulant2nd4050ALICE = new TGraphErrors(nPointsCumulant2nd4050ALICE,xCumulant2nd4050ALICE,yCumulant2nd4050ALICE,
                                                          xErrCumulant2nd4050ALICE,yErrCumulant2nd4050ALICE);
    Cumulant2nd4050ALICE->SetMarkerStyle(kOpenTriangleUp);
    Cumulant2nd4050ALICE->SetMarkerColor(kGreen+2);
    
    return Cumulant2nd4050ALICE;
  }
  
  return 0;
}

TGraphErrors* GetFlow6070_Rap10(Int_t n)
{
  // private communication 18.05.11, Ante B. 

  if (n == 2)
  {
    //  v2{SP}(pt) for 60-70%, rapidity gap = 1.0:
    const Int_t nPointsSP_6070ALICE_v2_etaGap10 = 9;
    Double_t xSP_6070ALICE_v2_etaGap10[nPointsSP_6070ALICE_v2_etaGap10] = {0.300000,0.500000,0.700000,0.900000,1.250000,1.750000,2.500000,3.500000,4.500000};
    Double_t ySP_6070ALICE_v2_etaGap10[nPointsSP_6070ALICE_v2_etaGap10] = {0.044958,0.073313,0.105726,0.120423,0.147537,0.186749,0.205423,0.208575,0.185938};
    Double_t xErrSP_6070ALICE_v2_etaGap10[nPointsSP_6070ALICE_v2_etaGap10] = {0.};
    Double_t yErrSP_6070ALICE_v2_etaGap10[nPointsSP_6070ALICE_v2_etaGap10] = {0.001520,0.001772,0.002245,0.002842,0.002600,0.004443,0.006240,0.014665,0.028810};
    TGraphErrors *GrSP_6070ALICE_v2_etaGap10 = new TGraphErrors(nPointsSP_6070ALICE_v2_etaGap10,xSP_6070ALICE_v2_etaGap10,
                                                              ySP_6070ALICE_v2_etaGap10,xErrSP_6070ALICE_v2_etaGap10,
                                                              yErrSP_6070ALICE_v2_etaGap10);

    return GrSP_6070ALICE_v2_etaGap10;
  }
  
  return 0;
}

Float_t CalculateFlow(TH1* ptDist, Float_t ptMin, Float_t ptMax, Int_t n, Int_t centralityBegin, Int_t centralityEnd)
{
  if (centralityBegin == 0 && centralityEnd == 1)
    flow = GetFlow01_Rap10(n);
  else if (centralityBegin == 0 && centralityEnd == 2)
    flow = GetFlow02_Rap10(n);
  else if (centralityBegin == 0 && centralityEnd == 5)
    flow = GetFlow05_Rap10(n);
  else if (centralityBegin == 5 && centralityEnd == 10)
    flow = GetFlow510_Rap10(n);
  else if (centralityBegin == 20 && centralityEnd == 30)
    flow = GetFlow2030(n);
  else if (centralityBegin == 30 && centralityEnd == 40)
    flow = GetFlow3040(n);
  else if (centralityBegin == 40 && centralityEnd == 50)
    flow = GetFlow4050(n);
  else if (centralityBegin > 50)
    flow = GetFlow6070_Rap10(n);
  else if (centralityBegin == 0 && centralityEnd == 20)
  {
    flow1 = GetFlow05_Rap10(n);
    flow2 = GetFlow510_Rap10(n);
    flow3 = GetFlow1020_Rap10(n);
    
    flow = (TGraphErrors*) flow2->Clone();
    
    // centrality width * dn/deta from http://arxiv.org/PS_cache/arxiv/pdf/1012/1012.1657v2.pdf
    Float_t mult[] = { 5 * 1601, 5 * 1294, 10 * 966 };
    
    for (Int_t i=0; i<flow->GetN(); i++)
    {
      Float_t x= flow->GetX()[i];

//       Printf("%f: %f %f %f", x, flow1->Eval(x), flow2->Eval(x), flow3->Eval(x));
      flow->GetY()[i] = flow1->Eval(x) * mult[0] + flow2->Eval(x) * mult[1] + flow3->Eval(x) * mult[2];
      flow->GetY()[i] /= mult[0] + mult[1] + mult[2];
//       Printf(" --> %f", flow->GetY()[i]);
    }
  }
  else if (centralityBegin == 20 && centralityEnd == 60)
  {
    flow = GetFlow2030(n);
    flow2 = GetFlow3040(n);
    flow3 = GetFlow4050(n);
    
    // centrality width * dn/deta from http://arxiv.org/PS_cache/arxiv/pdf/1012/1012.1657v2.pdf
    Float_t mult[] = { 10 * 649, 10 * 426, 10 * (261+149) };
    
    if (flow->GetN() != flow2->GetN() || flow2->GetN() != flow3->GetN())
      AliFatal("Incompatible graphs");
    
    for (Int_t i=0; i<flow->GetN(); i++)
    {
//       Printf("%f %f %f", flow->GetY()[i], flow2->GetY()[i], flow3->GetY()[i]);
      flow->GetY()[i] = flow->GetY()[i] * mult[0] + flow2->GetY()[i] * mult[1] + flow3->GetY()[i] * mult[2];
      flow->GetY()[i] /= mult[0] + mult[1] + mult[2];
//       Printf(" --> %f", flow->GetY()[i]);
    }
  }
  else
  {
    Printf("Flow range %d %d not available", centralityBegin, centralityEnd);
    AliFatal("");
  }

  Float_t vn = 0;
  Float_t sum = 0;
  for (Int_t bin = ptDist->FindBin(ptMin + 0.01); bin <= ptDist->FindBin(ptMax - 0.01); bin++)
  {
    if (ptDist->GetBinCenter(bin) > flow->GetX()[flow->GetN()-1])
      vn += ptDist->GetBinContent(bin) * flow->GetY()[flow->GetN()-1];
    else
      vn += ptDist->GetBinContent(bin) * flow->Eval(ptDist->GetBinCenter(bin));
    sum += ptDist->GetBinContent(bin);
  }
  
  if (sum > 0)
    vn /= sum;
  
  Printf("v_{%d} = %f for %f < pT < %f", n, vn, ptMin, ptMax);
  
  return vn;
}

void CalculateFlow(const char* fileName, Int_t centralityBegin, Int_t centralityEnd)
{
  Float_t ptTrigMin = 2;
  Float_t ptTrigMax = 3;
  
  Float_t ptMin = 1;
  Float_t ptMax = 2;

  loadlibs();

  AliUEHistograms* h = (AliUEHistograms*) GetUEHistogram(fileName);
  
  cont = h->GetUEHist(2)->GetTrackHist(0);
  cont->GetGrid(6)->GetGrid()->GetAxis(3)->SetRangeUser(0.01 + centralityBegin, -0.01 + centralityEnd);
  cont->GetGrid(6)->GetGrid()->GetAxis(2)->SetRangeUser(ptTrigMin + 0.01, ptTrigMax - 0.01);
  
  ptDist = cont->ShowProjection(1, 6);
  ptDist->Scale(1.0 / ptDist->Integral());
  
  cont = h->GetUEHist(2)->GetEventHist();
  cont->GetGrid(6)->GetGrid()->GetAxis(1)->SetRangeUser(0.01 + centralityBegin, -0.01 + centralityEnd);
  
  ptDist2 = cont->ShowProjection(0, 6);
  ptDist2->Scale(1.0 / ptDist2->Integral());
  
  TString str;
  
  for (Int_t n=2; n<=4; n++)
  {
    Float_t v2A = CalculateFlow(ptDist, ptMin, ptMax, n, centralityBegin, centralityEnd);
    Float_t v2T = CalculateFlow(ptDist2, ptTrigMin, ptTrigMax, n, centralityBegin, centralityEnd);
    
    str += Form("%f * %f, ", v2A, v2T);
  }
  
  ptDist->Draw();
  ptDist2->SetLineColor(2);
  ptDist2->Draw("SAME");
  
  Printf("%s", str.Data());
}

// four dimensions for: nearside/awayside/normalization unc, trigger, centrality, case (same, same/mixed, same w/ v2 subtraction, same/mixed w/ v2 subtraction, some more cases...)
TGraphErrors***** yields = 0;
TString currentYieldFile;

void WriteYields()
{
  TFile::Open("yields.root", "RECREATE");
  for (Int_t i=0; i<3; i++)
  {
    for (Int_t j=0; j<6; j++)
    {
      for (Int_t k=0; k<4; k++)
      {
        // CINT limitation here
        TGraphErrors** tmp = yields[i][j][k];
        for (Int_t l=0; l<31; l++)
        {
          //Printf("%d %d %d %d", i, j, k, l);
          tmp[l]->Write(Form("yield_%d_%d_%d_%d", i, j, k, l));
        }
      }
    }
  }
  gFile->Close();
}

void ReadYields(const char* fileName = "yields.root")
{
  if (currentYieldFile == fileName)
    return;
    
  currentYieldFile = fileName;

  CreateYieldStructure();
  TFile::Open(fileName);
  for (Int_t i=0; i<3; i++)
  {
    for (Int_t j=0; j<6; j++)
    {
      for (Int_t k=0; k<4; k++)
      {
        // CINT limitation here
        TGraphErrors** tmp = yields[i][j][k];
        for (Int_t l=0; l<31; l++)
        {
          //Printf("%d %d %d %d", i, j, k, l);
          tmp[l] = gFile->Get(Form("yield_%d_%d_%d_%d", i, j, k, l));
        }
      }
    }
  }
}

void CreateYieldStructure()
{
  if (!yields)
  {
    yields = new TGraphErrors****[3];
    for (Int_t i=0; i<3; i++)
    {
      yields[i] = new TGraphErrors***[6];
      for (Int_t j=0; j<6; j++)
      {
        yields[i][j] = new TGraphErrors**[4];
        for (Int_t k=0; k<4; k++)
        {
          yields[i][j][k] = new TGraphErrors*[31];
          // CINT limitation here
          TGraphErrors** tmp = yields[i][j][k];
          for (Int_t l=0; l<31; l++)
          {
            //Printf("%d %d %d %d", i, j, k, l);
            TGraphErrors* graph = new TGraphErrors;
            tmp[l] = graph;
          }
        }
      }
    }
  }
}

void GraphShiftX(TGraphErrors* graph, Float_t offset)
{
  for (Int_t i=0; i<graph->GetN(); i++)
    graph->GetX()[i] += offset;
}

// Float_t kPythiaScalingFactor = 0.935;
Float_t kPythiaScalingFactor = 1;

void DrawYields(const char* fileName = "yields.root")
{
  ReadYields(fileName);
  
  c = new TCanvas("c", "c", 1800, 1200);
  c->Divide(6, 4);
  
  Int_t markers1[] = { 24, 25, 26, 30 };
  Int_t markers2[] = { 20, 21, 22, 29 };
  Int_t colors[] = { 1, 2, 3, 4 };
  //Int_t caseList[] = { 0, 10+1, 10+1, 8, 10, 11, 12, 13, 1, 9 };
  Int_t caseList[] = { 0, 18, 18, 23, 18, 19, 20, 21, 22, 9 };
  const char* caseString[] = { "", "baseline sub", "baseline vs v2 sub comp", 0, "baseline sub comp", "Same/Mixed", "Same/Mixed - v2 subtr" };
  
  Bool_t iaa = kTRUE;
  if (!iaa)
  {
    Int_t centralityList[] =  {  0,  1,  2, 0, 0, 1 };
    Int_t centralityList2[] = { -1, -1, -1, 1, 2, 2 };
    Float_t factors[]       = {  1,  1,  1, 1, 1, 1 };
    const char* centralityString[] = { "0-5%", "20-40%", "60-90%", "0-5% vs 20-40%", "0-5% vs 60-90%", "20-40% vs 60-90%" };
  }
  else
  {
    Int_t centralityList[] =  {  0,  2,  3, 0, 0, 2 };
    Int_t centralityList2[] = { -1, -1, -1, 2, 3, 3 };
    //Float_t factors[]       = {  1,  1,  1, 1, 1.0/0.9, 1.0/0.9 };
    //Float_t factors[]       = {  1,  1,  1, 1, 1.0/0.804, 1.0/0.804 };
    Float_t factors[]       = {  1,  1,  1, 1, 1.0/kPythiaScalingFactor, 1.0/kPythiaScalingFactor };
    //Float_t factors[]       = {  1,  1,  1, 1, 1.0, 1.0 };
    const char* centralityString[] = { "0-5%", "60-90%", "pp", "0-5% vs 60-90%", "0-5% vs pp", "60-90% vs pp" };
  }
  
  l = new TLegend(0.5, 0.5, 0.9, 0.9);
  l->SetFillColor(0);
  
  Float_t max[] = { 0.4, 2, 2, 2, 2, 10, 10 };
  Float_t max2[] = { 5, 2, 2, 2, 2, 2, 2, 2, 2 };
  
  dummy = new TH2F("dummy", ";p_{T,assoc};yield", 100, 0, 20, 20000, 0, 10);
  dummy->SetStats(0);
  
  TGraphErrors* prevNear[6];
  TGraphErrors* prevAway[6];
  for (Int_t caseId=0; caseId<9; caseId++) // case 0->9
  {
    for (Int_t k=0; k<6; k++) // centrality
    {
      if (caseId != 3 && caseId < 5)
      {
        if (caseId < 4)
          c->cd(1 + k + caseId*6);
        else
          c->cd(1 + k + (caseId-1)*6);
        currentDummy = dummy->DrawCopy();
        currentDummy->GetYaxis()->SetRangeUser(0, max[caseId]);
        gPad->SetGridx();
        gPad->SetGridy();
        
        latex = new TLatex(0.15, 0.95, Form("%s - %s", caseString[caseId], centralityString[k]));
        latex->SetNDC();
        latex->SetTextSize(0.05);
        latex->Draw();
      }
      else if (caseId == 3)
        c->cd(1 + k + 2*6);
      else
        c->cd(1 + k + 3*6);
      
      for (Int_t j=0; j<3; j++) // trigger pt
      {
        if (caseId >= 2 && j != 1)
          continue;
        Printf("%d %d %d %d", caseId, k, j, caseList[caseId]);
      
        if (k < 3)
        {
          // CINT limitation here
          TGraphErrors** tmp = yields[0][j][centralityList[k]];
          nearSide = tmp[caseList[caseId]];
          //new TCanvas; nearSide->Draw("AP"); return;
          
          nearSide->SetMarkerStyle(markers1[j]);
          nearSide->SetMarkerColor(colors[j]);
          nearSide->SetLineColor(colors[j]);
          
          if (caseId == 0 && k == 0)
            l->AddEntry(nearSide, Form("trig pt bin %d", j), "P");
            
          TGraphErrors** tmp = yields[1][j][centralityList[k]];
          awaySide = tmp[caseList[caseId]];
          awaySide->SetMarkerStyle(markers2[j]);
          awaySide->SetMarkerColor(colors[j]);
          awaySide->SetLineColor(colors[j]);
            
          GraphShiftX(nearSide, -0.1 + j*0.1);
          GraphShiftX(awaySide, -0.1 + j*0.1);
          
          if (caseId == 3)
          {
            nearSide->SetLineColor(caseId - 2);
            nearSide->SetMarkerColor(caseId - 2);
            awaySide->SetLineColor(caseId - 2);
            awaySide->SetMarkerColor(caseId - 2);
          }
          
          if (caseId >= 4)
          {
            nearSide->SetLineColor(caseId - 3);
            nearSide->SetMarkerColor(caseId - 3);
            awaySide->SetLineColor(caseId - 3);
            awaySide->SetMarkerColor(caseId - 3);
          }
          
          awaySide->Print();
          
          prevNear[k] = (TGraphErrors*) nearSide->DrawClone("PSAME");
          prevAway[k] = (TGraphErrors*) awaySide->DrawClone("PSAME");
          
          if (caseId == 4)
          {
            Printf("%d", k);
            for (Int_t i=0; i<nearSide->GetN(); i++)
            {
              Printf("%f: %.0f%%", nearSide->GetX()[i], 100.0 * nearSide->GetEY()[i] / nearSide->GetY()[i]);
              Printf("%f: %.0f%%", awaySide->GetX()[i], 100.0 * awaySide->GetEY()[i] / awaySide->GetY()[i]);
            }
          }
        }
        else
        {
          // CINT limitation here
          tmp = yields[0][j][centralityList[k]]; nearSideCentral = (TGraphErrors*) tmp[caseList[caseId]]->Clone();
          tmp = yields[1][j][centralityList[k]]; awaySideCentral = (TGraphErrors*) tmp[caseList[caseId]]->Clone();
          
          tmp = yields[0][j][centralityList2[k]]; 
          nearSidePeripheral = tmp[caseList[caseId]];
          if (iaa && caseId == 3 && k >= 4)
            nearSidePeripheral = tmp[caseList[caseId-1]];
          
          tmp = yields[1][j][centralityList2[k]]; 
          awaySidePeripheral = tmp[caseList[caseId]];
          if (iaa && caseId == 3 && k >= 4)
            awaySidePeripheral = tmp[caseList[caseId-1]];
            
//           nearSideCentral->Print();
//           nearSidePeripheral->Print();
          
          if (k == 3)
            currentDummy->GetYaxis()->SetTitle("I_{CP}");
          else
            currentDummy->GetYaxis()->SetTitle("I_{AA}");
          currentDummy->GetYaxis()->SetRangeUser(0, max2[caseId]);
          
          for (Int_t i=0; i<nearSideCentral->GetN(); i++)
          {
            if (i >= nearSidePeripheral->GetN())
            {
              nearSideCentral->RemovePoint(i);
              i--;
              continue;
            }
          
            //Printf("near %d %f %f", i, nearSideCentral->GetY()[i], nearSidePeripheral->GetY()[i]);
            if (nearSidePeripheral->GetY()[i] <= 1e-5 || nearSideCentral->GetY()[i] <= 1e-5)
            {
              nearSideCentral->RemovePoint(i);
              nearSidePeripheral->RemovePoint(i);
              i--;
              continue;
            }
          
            nearSideCentral->GetEY()[i] = TMath::Sqrt(
              TMath::Power(nearSideCentral->GetEY()[i] / nearSideCentral->GetY()[i], 2) + 
              TMath::Power(nearSidePeripheral->GetEY()[i] / nearSidePeripheral->GetY()[i], 2) );
              
            nearSideCentral->GetY()[i] /= nearSidePeripheral->GetY()[i];
            nearSideCentral->GetY()[i] *= factors[k];
            
            nearSideCentral->GetEY()[i] *= nearSideCentral->GetY()[i];
          }
          //Printf("done");
          
          for (Int_t i=0; i<awaySideCentral->GetN(); i++)
          {
            if (i >= awaySidePeripheral->GetN())
            {
              awaySideCentral->RemovePoint(i);
              i--;
              continue;
            }
            
            //Printf("away %d", i);
            if (awaySidePeripheral->GetY()[i] <= 1e-5 || awaySideCentral->GetY()[i] <= 1e-5)
            {
              awaySideCentral->RemovePoint(i);
              awaySidePeripheral->RemovePoint(i);
              i--;
              continue;
            }
            
            awaySideCentral->GetEY()[i] = TMath::Sqrt(
              TMath::Power(awaySideCentral->GetEY()[i] / awaySideCentral->GetY()[i], 2) + 
              TMath::Power(awaySidePeripheral->GetEY()[i] / awaySidePeripheral->GetY()[i], 2) );
            
            awaySideCentral->GetY()[i] /= awaySidePeripheral->GetY()[i];
            awaySideCentral->GetY()[i] *= factors[k];
          
            awaySideCentral->GetEY()[i] *= awaySideCentral->GetY()[i];
          }
          
          if (caseId == 3)
          {
            nearSideCentral = (TGraphErrors*) nearSideCentral->Clone();
            nearSideCentral->SetLineColor(1);
            nearSideCentral->SetMarkerColor(1);
          
            awaySideCentral = (TGraphErrors*) awaySideCentral->Clone();
            awaySideCentral->SetLineColor(1);
            awaySideCentral->SetMarkerColor(1);
          }
          
          if (caseId >= 4)
          {
            nearSideCentral = (TGraphErrors*) nearSideCentral->Clone();
            nearSideCentral->SetLineColor(caseId - 3);
            nearSideCentral->SetMarkerColor(caseId - 3);
          
            awaySideCentral = (TGraphErrors*) awaySideCentral->Clone();
            awaySideCentral->SetLineColor(caseId - 3);
            awaySideCentral->SetMarkerColor(caseId - 3);
          }
          
          //Printf("%d", caseList[caseId]);
          
//           nearSideCentral->Print();
          
          nearSideCentral->Draw("PSAME");
          awaySideCentral->Draw("PSAME");
        
          if (caseId == 3)
          {
            for (Int_t i=0; i<nearSideCentral->GetN(); i++)
              Printf("Near, bin %d pt = %f, difference %.1f%%", i, nearSideCentral->GetX()[i], 100.0 - 100.0 * nearSideCentral->GetY()[i] / prevNear[k]->Eval(nearSideCentral->GetX()[i]));
            for (Int_t i=0; i<awaySideCentral->GetN(); i++)
              Printf("Away, bin %d pt = %f, difference %.1f%%", i, awaySideCentral->GetX()[i], 100.0 - 100.0 * awaySideCentral->GetY()[i] / prevAway[k]->Eval(awaySideCentral->GetX()[i]));
          }
        
          prevNear[k] = nearSideCentral;
          prevAway[k] = awaySideCentral;
        }
      }
      if (caseId == 0 && k == 0)
        l->Draw();
    }
  }
  
  c->SaveAs("yields.png");
}

void FitGaussians(const char* fileName, Bool_t flat)
{
  CreateYieldStructure();

  aliceFile = TFile::Open(fileName);

  Int_t maxLeadingPt = 3;
  Int_t maxAssocPt = 3;
  
  TCanvas* canvas = new TCanvas("FitGaussians", "FitGaussians", 1000, 700);
  canvas->Divide(maxAssocPt, maxLeadingPt);
  
  for (Int_t i=0; i<maxLeadingPt; i++)
    for (Int_t j=0; j<maxAssocPt; j++)
    {
      TH1* first = 0;
      for (Int_t aliceCentrality=0; aliceCentrality<4; aliceCentrality++)
      {
        Printf("%d %d %d", i, j, aliceCentrality);
        
        canvas->cd(j+1 + (i) * maxAssocPt);
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.2);
        //gPad->SetTopMargin(0.01);
        gPad->SetRightMargin(0.01);
        
        hist = (TH1*) aliceFile->Get(Form("dphi_%d_%d_%d%s", i, j, aliceCentrality, (flat) ? "_fit_flat" : ""));
        if (!hist)
          continue;
	
// 	hist->Rebin(2);
        
        // two Gauss fits
        gausFit = new TF1("gausFit", "[0] + gaus(1) + gaus(4)", -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
        gausFit->SetParameters(0, 1, 0, 1, 1, TMath::Pi(), 1);
        gausFit->SetParLimits(0, -1, 10000);
        gausFit->SetParLimits(1, 0.001, 10000);
        gausFit->FixParameter(2, 0);
        gausFit->SetParLimits(3, 0.05, 10);
        gausFit->SetParLimits(4, 0.001, 10000);
        gausFit->FixParameter(5, TMath::Pi());
        gausFit->SetParLimits(6, 0.1, 10);
        
        gausFit->SetLineWidth(1);
        gausFit->SetLineColor(hist->GetLineColor());
        
        hist->Fit(gausFit, "0RIQ");
        gausFit->FixParameter(0, gausFit->GetParameter(0));
        hist->Fit(gausFit, "RI", (aliceCentrality == 0) ? "" : "SAME");
        //gausFit->DrawCopy("SAME");
        
        // get pT,assoc
        TPRegexp reg("(\\d+\\.\\d+) \\< p_\\{T,assoc\\} \\< (\\d+\\.\\d+)");
        arr = reg.MatchS(hist->GetTitle(), "", 0, 10);
        if (arr->GetEntries() != 3)
          continue;
          
        Float_t pt1 = atof(arr->At(1)->GetName());
        Float_t pt2 = atof(arr->At(2)->GetName());

        FillYield(i, aliceCentrality, (pt1 + pt2) / 2, (pt2 - pt1) / 2, 18, gausFit->GetParameter(3), gausFit->GetParError(3), gausFit->GetParameter(6), gausFit->GetParError(6));
        
        if (!first)
          first = hist;
        else
          first->GetYaxis()->SetRangeUser(TMath::Min(first->GetMinimum(), hist->GetMinimum()), TMath::Max(first->GetMaximum(), hist->GetMaximum()));
      }
      first->GetYaxis()->SetRangeUser(first->GetMinimum(), 2 * first->GetMaximum());
        
//        break;
    }
    
    
  // draw
  
  for (Int_t j=0; j<maxLeadingPt; j++)
  {
    new TCanvas;
    
    Int_t markers1[] = { 24, 25, 26, 30 };
    Int_t markers2[] = { 20, 21, 22, 29 };
    Int_t colors[] = { 1, 2, 3, 4 };
    
    dummy = new TH2F("dummy", ";p_{T,assoc};", 100, 0, 10, 1000, 0, 1);
    dummy->SetStats(0);
    dummy->Draw();
    
    for (Int_t k=0; k<4; k++) // centrality
    {
      // CINT limitation here
      TGraphErrors** tmp = yields[0][j][k];
      nearSide = tmp[18];
      nearSide->SetMarkerStyle(markers1[k]);
      nearSide->SetMarkerColor(colors[k]);
      nearSide->SetLineColor(colors[k]);  
    
      TGraphErrors** tmp = yields[1][j][k];
      awaySide = tmp[18];
      awaySide->SetMarkerStyle(markers2[k]);
      awaySide->SetMarkerColor(colors[k]);
      awaySide->SetLineColor(colors[k]);  
      
      nearSide->DrawClone("PSAME");
      awaySide->DrawClone("PSAME");
    }
  }
}

void TsallisExamples()
{
  tsallis = new TF1("tsallis", "[0] * (1-[2]*(1-[1])*x*x)**(1/(1-[1]))", -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
  
  Float_t q[] = { 1.0000001, 1.5, 2 };
  Float_t beta[] = { 0.1, 1, 10, 100 };
  
  dummy = new TH2F("dummy", "", 100, -0.5 * TMath::Pi(), 1.5 * TMath::Pi(), 100, 0, 1.2);
  dummy->SetStats(0);
  dummy->Draw();
  
  legend = new TLegend(0.65, 0.5, 1, 1);
  legend->SetFillColor(0);
  
  for (Int_t i=0; i<3; i++)
    for (Int_t j=0; j<4; j++)
    {
      tsallis->SetParameters(1, q[i], beta[j]);
      tsallis->SetLineColor(i+1);
      tsallis->SetLineStyle(j+1);
      
      legend->AddEntry(tsallis->DrawCopy("SAME"), Form("q = %.1f, #beta = %.1f", q[i], beta[j]), "L");
    }
    
  legend->Draw();
}

void DrawTsallisParams()
{
  ReadYields();
  
  c = new TCanvas("c", "c", 900, 600);
  c->Divide(3, 2);
  
  Int_t markers1[] = { 24, 25, 26, 30 };
  Int_t markers2[] = { 20, 21, 22, 29 };
  Int_t colors[] = { 1, 2, 3, 4 };
  Int_t caseList[] = { 12, 13, 8, 6, 8, 9, 10, 11, 1, 7 };
  const char* caseString[] = { "q", "beta" };
  const char* centralityString[] = { "0-5%", "20-40%", "60-90%", "0-5% vs 20-40%", "0-5% vs 60-90%", "20-40% vs 60-90%" };
  
  l = new TLegend(0.5, 0.5, 0.9, 0.9);
  l->SetFillColor(0);
  
  Float_t max[] = { 3, 1000 };
  Float_t max2[] = { 5, 2 };
  
  dummy = new TH2F("dummy", ";p_{T,assoc};", 100, 0, 20, 1000, 0, 1000);
  dummy->SetStats(0);
  
  for (Int_t caseId=0; caseId<2; caseId++) // case
  {
    for (Int_t k=0; k<3; k++) // centrality
    {
      c->cd(1 + k + caseId*3);
        
      currentDummy = dummy->DrawCopy();
      currentDummy->GetYaxis()->SetRangeUser(0, max[caseId]);
      if (caseId == 1 && k < 3)
        gPad->SetLogy();
      
      latex = new TLatex(0.15, 0.95, Form("%s - %s", caseString[caseId], centralityString[k]));
      latex->SetNDC();
      latex->SetTextSize(0.05);
      latex->Draw();
      
      for (Int_t j=0; j<2; j++) // trigger pt
      {
        if (k < 3)
        {
          // CINT limitation here
          TGraphErrors** tmp = yields[0][j][k];
          nearSide = tmp[caseList[caseId]];
          nearSide->SetMarkerStyle(markers1[j]);
          nearSide->SetMarkerColor(colors[j]);
          nearSide->SetLineColor(colors[j]);
          
          if (caseId == 0 && k == 0)
            l->AddEntry(nearSide, Form("trig pt bin %d", j), "P");
            
          TGraphErrors** tmp = yields[1][j][k];
          awaySide = tmp[caseList[caseId]];
          awaySide->SetMarkerStyle(markers2[j]);
          awaySide->SetMarkerColor(colors[j]);
          awaySide->SetLineColor(colors[j]);
            
          GraphShiftX(nearSide, -0.1 + j*0.1);
          GraphShiftX(awaySide, -0.1 + j*0.1);
          
          Printf("%d %d %d", caseId, j, k);
          nearSide->Print();
          awaySide->Print();
          
          nearSide->DrawClone("PSAME");
          awaySide->DrawClone("PSAME");
        }
        else
        {
          // CINT limitation here
          if (k == 3)
          {
            tmp = yields[0][j][0]; nearSideCentral = (TGraphErrors*) tmp[caseList[caseId]]->Clone();
            tmp = yields[1][j][0]; awaySideCentral = (TGraphErrors*) tmp[caseList[caseId]]->Clone();
            tmp = yields[0][j][1]; nearSidePeripheral = tmp[caseList[caseId]];
            tmp = yields[1][j][1]; awaySidePeripheral = tmp[caseList[caseId]];
          }
          else if (k == 4)
          {
            tmp = yields[0][j][0]; nearSideCentral = (TGraphErrors*) tmp[caseList[caseId]]->Clone();
            tmp = yields[1][j][0]; awaySideCentral = (TGraphErrors*) tmp[caseList[caseId]]->Clone();
            tmp = yields[0][j][2]; nearSidePeripheral = tmp[caseList[caseId]];
            tmp = yields[1][j][2]; awaySidePeripheral = tmp[caseList[caseId]];
          }
          else if (k == 5)
          {
            tmp = yields[0][j][1]; nearSideCentral = (TGraphErrors*) tmp[caseList[caseId]]->Clone();
            tmp = yields[1][j][1]; awaySideCentral = (TGraphErrors*) tmp[caseList[caseId]]->Clone();
            tmp = yields[0][j][2]; nearSidePeripheral = tmp[caseList[caseId]];
            tmp = yields[1][j][2]; awaySidePeripheral = tmp[caseList[caseId]];
          }
          
          currentDummy->GetYaxis()->SetTitle("I_{CP}");
          currentDummy->GetYaxis()->SetRangeUser(0, max2[caseId]);
        
          for (Int_t i=0; i<nearSideCentral->GetN(); i++)
          {
            //Printf("near %d %f %f", i, nearSideCentral->GetY()[i], nearSidePeripheral->GetY()[i]);
            if (nearSidePeripheral->GetY()[i] <= 1e-5 || nearSideCentral->GetY()[i] <= 1e-5)
            {
              nearSideCentral->RemovePoint(i);
              nearSidePeripheral->RemovePoint(i);
              i--;
              continue;
            }
          
            nearSideCentral->GetEY()[i] = TMath::Sqrt(
              TMath::Power(nearSideCentral->GetEY()[i] / nearSideCentral->GetY()[i], 2) + 
              TMath::Power(nearSidePeripheral->GetEY()[i] / nearSidePeripheral->GetY()[i], 2) );
              
            nearSideCentral->GetY()[i] /= nearSidePeripheral->GetY()[i];
            
            nearSideCentral->GetEY()[i] *= nearSideCentral->GetY()[i];
          }
          //Printf("done");
          
          for (Int_t i=0; i<awaySideCentral->GetN(); i++)
          {
            //Printf("away %d", i);
            if (awaySidePeripheral->GetY()[i] <= 1e-5 || awaySideCentral->GetY()[i] <= 1e-5)
            {
              awaySideCentral->RemovePoint(i);
              awaySidePeripheral->RemovePoint(i);
              i--;
              continue;
            }
            
            awaySideCentral->GetEY()[i] = TMath::Sqrt(
              TMath::Power(awaySideCentral->GetEY()[i] / awaySideCentral->GetY()[i], 2) + 
              TMath::Power(awaySidePeripheral->GetEY()[i] / awaySidePeripheral->GetY()[i], 2) );
            
            awaySideCentral->GetY()[i] /= awaySidePeripheral->GetY()[i];
          
            awaySideCentral->GetEY()[i] *= awaySideCentral->GetY()[i];
          }
          
          if (caseId == 3)
          {
            nearSideCentral = (TGraphErrors*) nearSideCentral->Clone();
            nearSideCentral->SetLineColor(1);
            nearSideCentral->SetMarkerColor(1);
          
            awaySideCentral = (TGraphErrors*) awaySideCentral->Clone();
            awaySideCentral->SetLineColor(1);
            awaySideCentral->SetMarkerColor(1);
          }
          
          if (caseId >= 4)
          {
            nearSideCentral = (TGraphErrors*) nearSideCentral->Clone();
            nearSideCentral->SetLineColor(caseId - 3);
            nearSideCentral->SetMarkerColor(caseId - 3);
          
            awaySideCentral = (TGraphErrors*) awaySideCentral->Clone();
            awaySideCentral->SetLineColor(caseId - 3);
            awaySideCentral->SetMarkerColor(caseId - 3);
          }
          
          //Printf("%d", caseList[caseId]);
          
          nearSideCentral->Draw("PSAME");
          awaySideCentral->Draw("PSAME");
        }
      }
      if (caseId == 0 && k == 0)
        l->Draw();
    }
  }
  
  c->SaveAs("tsallisparams.png");
}

void DrawYieldLHCRHIC_ICP(const char* fileName = "yields.root")
{
  Int_t caseId = 0;
  Int_t caseList[] = { 18 };
  Int_t j = 1;
  
  ReadYields(fileName);
  
  c = new TCanvas("lhc_rhic", "lhc_rhic", 800, 600);
  gPad->SetGridx();
  gPad->SetGridy();
  
  dummy = new TH2F("dummy", ";p_{T,assoc};I_{CP}", 100, 0, 12, 1000, 0, 10);
  dummy->SetStats(0);
  currentDummy = dummy->DrawCopy();
  currentDummy->GetYaxis()->SetRangeUser(0, 2);
        
  TGraphErrors** tmp = yields[0][j][0]; nearSideCentral = (TGraphErrors*) tmp[caseList[caseId]]->Clone();
  tmp = yields[1][j][0]; awaySideCentral = (TGraphErrors*) tmp[caseList[caseId]]->Clone();
  tmp = yields[0][j][2]; nearSidePeripheral = tmp[caseList[caseId]];
  tmp = yields[1][j][2]; awaySidePeripheral = tmp[caseList[caseId]];
        
  for (Int_t i=0; i<nearSideCentral->GetN(); i++)
  {
    //Printf("near %d %f %f", i, nearSideCentral->GetY()[i], nearSidePeripheral->GetY()[i]);
    if (nearSidePeripheral->GetY()[i] <= 1e-5 || nearSideCentral->GetY()[i] <= 1e-5)
    {
      nearSideCentral->RemovePoint(i);
      nearSidePeripheral->RemovePoint(i);
      i--;
      continue;
    }
  
    nearSideCentral->GetEY()[i] = TMath::Sqrt(
      TMath::Power(nearSideCentral->GetEY()[i] / nearSideCentral->GetY()[i], 2) + 
      TMath::Power(nearSidePeripheral->GetEY()[i] / nearSidePeripheral->GetY()[i], 2) );
      
    nearSideCentral->GetY()[i] /= nearSidePeripheral->GetY()[i];
    
    nearSideCentral->GetEY()[i] *= nearSideCentral->GetY()[i];
  }
  //Printf("done");
  
  for (Int_t i=0; i<awaySideCentral->GetN(); i++)
  {
    //Printf("away %d", i);
    if (awaySidePeripheral->GetY()[i] <= 1e-5 || awaySideCentral->GetY()[i] <= 1e-5)
    {
      awaySideCentral->RemovePoint(i);
      awaySidePeripheral->RemovePoint(i);
      i--;
      continue;
    }
    
    awaySideCentral->GetEY()[i] = TMath::Sqrt(
      TMath::Power(awaySideCentral->GetEY()[i] / awaySideCentral->GetY()[i], 2) + 
      TMath::Power(awaySidePeripheral->GetEY()[i] / awaySidePeripheral->GetY()[i], 2) );
    
    awaySideCentral->GetY()[i] /= awaySidePeripheral->GetY()[i];
  
    awaySideCentral->GetEY()[i] *= awaySideCentral->GetY()[i];
  }
          
  nearSideCentral->SetMarkerStyle(21);
  nearSideCentral->SetLineColor(1);
  nearSideCentral->SetMarkerColor(1);
          
  awaySideCentral->SetMarkerStyle(25);
  awaySideCentral->SetLineColor(1);
  awaySideCentral->SetMarkerColor(1);
  
  nearSideCentral->Draw("PSAME");
  awaySideCentral->Draw("PSAME");        

  legend = new TLegend(0.4, 0.8, 0.99, 0.99);
  legend->SetFillColor(0);
  legend->AddEntry(nearSideCentral, "ALICE I_{CP} Near 8-15 GeV 0-5%/60-90%", "P");
  legend->AddEntry(awaySideCentral, "ALICE I_{CP} Away 8-15 GeV 0-5%/60-90%", "P");
  
  if (0)
  {
    nearSideCentral->Fit("pol1", "", "SAME", 4.5, 11.5);
    awaySideCentral->Fit("pol1", "", "SAME", 4.5, 11.5);
    nearSideCentral->GetFunction("pol1")->SetLineWidth(1);
    awaySideCentral->GetFunction("pol1")->SetLineWidth(1);
  }

  if (0)
  {
    // RCP, ALICE
    dndpt_central = ReadHepdata("raa_dndpt_central.txt", kFALSE, 3);
    dndpt_peripheral = ReadHepdata("raa_dndpt_peripheral.txt", kFALSE, 3);
    
    //dndpt_central->Print();
    //dndpt_peripheral->Print();
    
    // ratio of Ncoll (from paper)
    // TODO there is also an uncertainty on that
    Double_t scale = 1690. / 15.7;
    Double_t scaleUnc = TMath::Sqrt(TMath::Power(131. / 1690., 2) + TMath::Power(0.7 / 15.7, 2));
    
    for (Int_t i=0; i<dndpt_central->GetN(); i++)
    {
      dndpt_central->GetEY()[i] = TMath::Sqrt(TMath::Power(dndpt_central->GetEY()[i] / dndpt_central->GetY()[i], 2) + 
        TMath::Power(dndpt_peripheral->GetEY()[i] / dndpt_peripheral->GetY()[i], 2));
      dndpt_central->GetY()[i] /= dndpt_peripheral->GetY()[i] * scale;
      dndpt_central->GetEY()[i] *= dndpt_central->GetY()[i];
    }
    
    //new TCanvas;
    //dndpt_central->Print();
    //dndpt_central->Draw("AP");
    
    dndpt_central->Draw("*SAME");
    
    latex = new TLatex(0.34, 0.09, "9% N_{coll} scale uncertainty on R_{CP}");
    latex->SetTextSize(0.03);
    latex->Draw();
    
    // RAA, ALICE
    raa_central = ReadHepdata("raa_alice_central.txt", kFALSE, 3);
    raa_central->SetMarkerStyle(20);
    //raa_central->Draw("PSAME");
  
    legend->AddEntry(dndpt_central, "ALICE R_{CP} 0-5%/70-80%", "P");
  }
  
  if (0)
  {
    // IAA, RHIC
    // calculate ICP from IAA(0-20) / IAA (20-60)
    // TODO only using uncertainties from IAA(0-20)
    // TODO systematic uncertainty stored in graph with _sys appended
    TFile::Open("rhic/pi0h_graphs.root");
    
    for (Int_t i=0; i<3; i+=2)
    {
      // gIAA_<ptTrigBin>_<centBin>_<angularRegionBin>
      central = (TGraphErrors*) gFile->Get(Form("gIAA_3_0_%d", i));
      peripheral = (TGraphErrors*) gFile->Get(Form("gIAA_3_1_%d", i));
      
      for (Int_t j=0; j<central->GetN(); j++)
      {
        central->GetY()[j] /= peripheral->GetY()[j];
        central->GetEY()[j] /= peripheral->GetY()[j];
      }
      
      central->SetMarkerStyle((i == 0) ? 29 : 30);
      central->SetMarkerColor(2);
      central->SetLineColor(2);
      central->Draw("PSAME");
      legend->AddEntry(central, Form("PHENIX I_{CP} %s 9-12 GeV 0-20%/20-60%", (i == 0) ? "Near" : "Away"), "P");
    }
  }
  
  legend->Draw();
  
  c->SaveAs("icp.png");
}

TGraphErrors* GetRatio(const char* fileName, Int_t centrality1, Int_t centrality2, Int_t triggerId, Int_t caseId, Int_t side, const char* fileName2 = 0)
{
  // 0 = near side; 1 = away side

  ReadYields(fileName);
  
  TGraphErrors** tmp = yields[side][triggerId][centrality1]; 
  nearSideCentral = (TGraphErrors*) tmp[caseId]->Clone();
  
  if (fileName2)
    ReadYields(fileName2);
  
  tmp = yields[side][triggerId][centrality2]; 
  nearSidePeripheral = (TGraphErrors*) tmp[caseId]->Clone();
  
//   Printf("%d %d", nearSideCentral->GetN(), nearSidePeripheral->GetN());
        
  for (Int_t i=0; i<nearSideCentral->GetN(); i++)
  {
    if (nearSidePeripheral->GetY()[i] <= 1e-5 || nearSideCentral->GetY()[i] <= 1e-5)
    {
      nearSideCentral->RemovePoint(i);
      nearSidePeripheral->RemovePoint(i);
      i--;
      continue;
    }
    
    if (nearSideCentral->GetX()[i] != nearSidePeripheral->GetX()[i])
    {
      Printf("Inconsistent x values %f %f", nearSideCentral->GetX()[i], nearSidePeripheral->GetX()[i]);
      AliFatal("");
    }
  
//     printf("%f %f %f", nearSideCentral->GetX()[i], nearSideCentral->GetY()[i], nearSidePeripheral->GetY()[i]);
    
    nearSideCentral->GetEY()[i] = TMath::Sqrt(
      TMath::Power(nearSideCentral->GetEY()[i] / nearSideCentral->GetY()[i], 2) + 
      TMath::Power(nearSidePeripheral->GetEY()[i] / nearSidePeripheral->GetY()[i], 2) );
      
    nearSideCentral->GetY()[i] /= nearSidePeripheral->GetY()[i];
    // scaling for IAA
    if (centrality2 == 3)
      nearSideCentral->GetY()[i] /= kPythiaScalingFactor;
    
//     Printf(" %f", nearSideCentral->GetY()[i]);
    
    nearSideCentral->GetEY()[i] *= nearSideCentral->GetY()[i];
  }
  
  //nearSideCentral->GetXaxis()->SetTitle("p_{T, assoc}");
  //nearSideCentral->GetYaxis()->SetTitle((centrality2 == 3) ? "I_{AA}" : "I_{CP}");

  if (0 && side == 1 && centrality2 == 3)
  {
    Printf("\n\n\nWARNING !!! Fudging different away side acceptance (only for second file)\n\n\n");
    ScaleGraph(nearSideCentral, 1.56);
  }
  
  return nearSideCentral;
}

TGraphErrors* GetIAA(const char* fileName, Int_t centrality, Int_t triggerId, Int_t caseId, Int_t side)
{
  return GetRatio(fileName, centrality, 3, triggerId, caseId, side);
}

TGraphErrors* GetICP(const char* fileName, Int_t triggerId, Int_t caseId, Int_t side)
{
  return GetRatio(fileName, 1, 2, triggerId, caseId, side);
}

void DrawYieldLHCRHIC_IAA(const char* fileName, Bool_t central, Bool_t rhic)
{
  Int_t caseId = 0;
  Int_t caseList[] = { 18 };
  Int_t j = 1;
  
  ReadYields(fileName);
  
  c = new TCanvas("lhc_rhic", "lhc_rhic", 800, 600);
  gPad->SetGridx();
  gPad->SetGridy();
  
  dummy = new TH2F("dummy", ";p_{T,assoc};I_{AA}", 100, 0, 12, 1000, 0, 10);
  dummy->SetStats(0);
  currentDummy = dummy->DrawCopy();
  currentDummy->GetYaxis()->SetRangeUser(0, 3);
        
  Int_t nominatorBin = 0;
  if (!central)
    nominatorBin = 2;
        
  nearSideCentral = GetIAA(fileName, nominatorBin, j, caseList[caseId], 0);
  awaySideCentral = GetIAA(fileName, nominatorBin, j, caseList[caseId], 1);
  
  nearSideCentral->SetMarkerStyle(21);
  nearSideCentral->SetLineColor(1);
  nearSideCentral->SetMarkerColor(1);
          
  awaySideCentral->SetMarkerStyle(25);
  awaySideCentral->SetLineColor(1);
  awaySideCentral->SetMarkerColor(1);
  
  nearSideCentral->Draw("PSAME");
  awaySideCentral->Draw("PSAME");        
  
  legend = new TLegend(0.4, 0.8, 0.99, 0.99);
  legend->SetFillColor(0);
  legend->AddEntry(nearSideCentral, Form("ALICE I_{AA} Near 8-15 GeV %s%%/Pythia6", (central) ? "0-5" : "60-90"), "P");
  legend->AddEntry(awaySideCentral, Form("ALICE I_{AA} Away 8-15 GeV %s%%/Pythia6", (central) ? "0-5" : "60-90"), "P");
  
  if (rhic)
  {
    // IAA, RHIC
    // systematic uncertainty stored in graph with _sys appended
    TFile::Open("rhic/pi0h_graphs.root");
    
    for (Int_t i=0; i<3; i+=2)
    {
      // gIAA_<ptTrigBin>_<centBin>_<angularRegionBin>
      rhic_iaa = (TGraphErrors*) gFile->Get(Form("gIAA_2_%d_%d", (central) ? 0 : 1, i));
      rhic_iaa_sys = (TGraphErrors*) gFile->Get(Form("gIAA_2_%d_%d_sys", (central) ? 0 : 1, i));
  
      rhic_iaa->SetMarkerStyle((i == 0) ? 29 : 30);
      rhic_iaa->SetMarkerColor(2);
      rhic_iaa->SetLineColor(2);
      rhic_iaa->Draw("PSAME");
      rhic_iaa_sys->SetLineColor(2);
      rhic_iaa_sys->SetMarkerColor(2);
      rhic_iaa_sys->Draw("PSAME");
      legend->AddEntry(rhic_iaa, Form("PHENIX I_{AA} %s 7-9 GeV %s%%/pp", (i == 0) ? "Near" : "Away", (central) ? "0-20" : "20-60"), "P");
    }
    for (Int_t i=0; i<3; i+=2)
    {
      // gIAA_<ptTrigBin>_<centBin>_<angularRegionBin>
      rhic_iaa = (TGraphErrors*) gFile->Get(Form("gIAA_3_%d_%d", (central) ? 0 : 1, i));
      rhic_iaa_sys = (TGraphErrors*) gFile->Get(Form("gIAA_3_%d_%d_sys", (central) ? 0 : 1, i));
  
      rhic_iaa->SetMarkerStyle((i == 0) ? 20 : 24);
      rhic_iaa->SetMarkerColor(4);
      rhic_iaa->SetLineColor(4);
      rhic_iaa->Draw("PSAME");
      rhic_iaa_sys->SetLineColor(4);
      rhic_iaa_sys->SetMarkerColor(4);
      rhic_iaa_sys->Draw("PSAME");
      legend->AddEntry(rhic_iaa, Form("PHENIX I_{AA} %s 9-12 GeV %s%%/pp", (i == 0) ? "Near" : "Away", (central) ? "0-20" : "20-60"), "P");
    }
  }
  
  legend->Draw();
  
  c->SaveAs("iaa.png");
}

void CompareIAAICP(const char* fileName1, const char* fileName2, Int_t nominatorBin, Int_t denominatorBin, Bool_t skipAway = kFALSE)
{
  Int_t j = 1;
  Int_t caseId = 18;
//   caseId = 23; Printf("WARNING: Comparing case 23");
  //Int_t caseId = 12;

  nearSide1 = GetRatio(fileName1, nominatorBin, denominatorBin, j, caseId, 0);
  if (!skipAway)
    awaySide1 = GetRatio(fileName1, nominatorBin, denominatorBin, j, caseId, 1);
  else
    awaySide1 = new TGraphErrors;
  
  kPythiaScalingFactor = 1;
  nearSide2 = GetRatio(fileName2, nominatorBin, denominatorBin, j, caseId, 0);
  if (!skipAway)
    awaySide2 = GetRatio(fileName2, nominatorBin, denominatorBin, j, caseId, 1);
  else
    awaySide2 = new TGraphErrors;
  
//   Printf("\n\n\nWARNING !!! Fudging different away side acceptance (only for second file)\n\n\n");
//   ScaleGraph(awaySide2, 1.56);
  
/*  ScaleGraph(nearSide2, 1.33);
  ScaleGraph(awaySide2, 1.33);*/
  
  c = new TCanvas;
  dummy = new TH2F("dummy", Form(";p_{T,assoc};%s", (denominatorBin == 3) ? "I_{AA}" : "I_{CP}"), 100, 2, 12, 1000, 0, 10);
  dummy->SetStats(0);
  currentDummy = dummy->DrawCopy();
  currentDummy->GetYaxis()->SetRangeUser(0, 2.5);
  
  nearSide1->SetMarkerStyle(21);
  nearSide1->SetLineColor(1);
  nearSide1->SetMarkerColor(1);
          
  nearSide2->SetMarkerStyle(21);
  nearSide2->SetLineColor(2);
  nearSide2->SetMarkerColor(2);
  
  awaySide1->SetMarkerStyle(25);
  awaySide1->SetLineColor(1);
  awaySide1->SetMarkerColor(1);
  
  awaySide2->SetMarkerStyle(25);
  awaySide2->SetLineColor(2);
  awaySide2->SetMarkerColor(2);
  
  RemovePointsBelowX(nearSide1, 3);
  RemovePointsBelowX(nearSide2, 3);
  RemovePointsBelowX(awaySide1, 3);
  RemovePointsBelowX(awaySide2, 3);
  
  gPad->SetGridx();
  gPad->SetGridy();
  
  nearSide1->DrawClone("PSAME");
  awaySide1->DrawClone("PSAME");        
  nearSide2->DrawClone("PSAME");
  awaySide2->DrawClone("PSAME");        
  
  new TCanvas;
  currentDummy = dummy->DrawCopy();
  currentDummy->GetYaxis()->SetRangeUser(0.5, 1.5);
  gPad->SetGridx();
  gPad->SetGridy();
  
  DivideGraphs(nearSide1, nearSide2);
  DivideGraphs(awaySide1, awaySide2);
  
  nearSide1->Draw("PSAME");
  awaySide1->Draw("PSAME");        
  
  nearSide1->Fit("pol0", "", "SAME", 3.01, 9.99);
  awaySide1->Fit("pol0", "", "SAME", 3.01, 9.99);
}

void ScaleGraph(TGraphErrors* graph, Float_t factor)
{
	for (Int_t i=0; i<graph->GetN(); i++)
	{
    graph->GetY()[i] *= factor;
    graph->GetEY()[i] *= factor;
  }
}

TGraphErrors* DrawBaselines(const char* fileName, Int_t numeratorBin, Int_t denominatorBin, Int_t side, Bool_t verbose = kTRUE)
{
  Int_t j = 1;
  
  TGraphErrors* first = 0;
  TGraphErrors* max = 0;
  
  if (verbose)
  {
    c = new TCanvas;
    dummy = new TH2F("dummy", Form(";p_{T,assoc};%s", (denominatorBin == 3) ? "I_{AA}" : "I_{CP}"), 100, 0, 12, 1000, 0, 10);
    dummy->SetStats(0);
    currentDummy = dummy->DrawCopy();
    currentDummy->GetYaxis()->SetRangeUser(0, 3);
    
    c2 = new TCanvas;
    currentDummy = dummy->DrawCopy();
    currentDummy->GetYaxis()->SetRangeUser(0, 2);
  }
  
//   Int_t caseList[] = { 18, 19, 21, 22 };
  //Int_t caseList[] = { 23, 24, 26, 27 };
  Int_t caseList[] = { 18, 19, 22 };
  
  for (Int_t caseId = 0; caseId < 3; caseId++)
  {
    graph = GetRatio(fileName, numeratorBin, denominatorBin, j, caseList[caseId], side);
    
    graph->SetLineColor(caseId + 1);
    graph->SetMarkerColor(caseId + 1);
    
    RemovePointsBelowX(graph, 3);
    
    if (verbose)
    {
      c->cd();
      graph->DrawClone("*SAME");
    }
    
    if (!first)
    {
      first = graph;
      max = (TGraphErrors*) graph->Clone();
      for (Int_t i=0; i<max->GetN(); i++)
      {
        max->GetY()[i] = 0;
        max->GetEY()[i] = 0;
      }
    }
    else
    {
      DivideGraphs(graph, first);
      
      if (verbose)
      {
        c2->cd();
        graph->DrawClone("*SAME");
      }
      
      for (Int_t i=0; i<graph->GetN(); i++)
      {
        if (max->GetY()[i] < TMath::Abs(graph->GetY()[i] - 1))
          max->GetY()[i] = TMath::Abs(graph->GetY()[i] - 1);
      }
    }
  }
  
  if (verbose)
  {
    new TCanvas;
    max->Draw("A*");
  }
  
  max->GetXaxis()->SetTitle("p_{T, assoc}");
  max->GetYaxis()->SetTitle((denominatorBin == 3) ? "Effect on I_{AA}" : "Effect on I_{CP}");
  max->GetYaxis()->SetRangeUser(0, 0.4);
  
  return max;
}

void DrawBaselinesAll(const char* fileName)
{
  c = new TCanvas("c", "c", 1000, 600);
  c->Divide(4, 2);
  
  for (Int_t i=0; i<2; i++)
  {
    c->cd(1+i*4); DrawBaselines(fileName, 0, 3, i, kFALSE)->Draw("A*"); gPad->SetGridy();
    c->cd(2+i*4); DrawBaselines(fileName, 1, 3, i, kFALSE)->Draw("A*"); gPad->SetGridy();
    c->cd(3+i*4); DrawBaselines(fileName, 2, 3, i, kFALSE)->Draw("A*"); gPad->SetGridy();
    c->cd(4+i*4); DrawBaselines(fileName, 0, 2, i, kFALSE)->Draw("A*"); gPad->SetGridy();
  }
}

TGraphErrors* DrawSystUncIAAPYTHIA(TGraphErrors* graph, Int_t iaa_icp, Int_t awaySide)
{
  // iaa = 0, icp = 1

  Float_t baseline = 1;
  if (awaySide == 0)
    baseline = 0.07;
  else if (awaySide == 1)
    baseline = 0.2;
    
  /*
  Float_t reference = 0;
  if (iaa_icp == 0)
    reference = 0.13;
  */
    
  Float_t efficiency = 0.08;
  
  Float_t centrality = 1;
  if (iaa_icp == 0)
    centrality = (awaySide == 0) ? 0.02 : 0.06;
  else if (iaa_icp == 1)
    centrality = (awaySide == 0) ? 0.04 : 0.08;
  
  Float_t ptResolution = 0;
  
  Float_t total = TMath::Sqrt(baseline * baseline + efficiency * efficiency + centrality * centrality + ptResolution * ptResolution);
  
  Printf("Total syst: %f", total);
  
  systUnc = (TGraphErrors*) graph->Clone();
  for (Int_t i=0; i<systUnc->GetN(); i++)
  {
    systUnc->GetEY()[i] = systUnc->GetY()[i] * total;
    systUnc->GetEX()[i] = 0;
  }
    
  //systUnc->Print();
  
  systUnc->SetLineColor(kGray + 1);
  systUnc->SetLineWidth(6);
  
  systUnc->Draw("PSAME");
    
  return systUnc;
}

TGraphErrors* DrawSystUncPreliminaries(TGraphErrors* graph, Int_t iaa_icp, Int_t awaySide, Int_t iaa, Bool_t central)
{
  // iaa = 0, icp = 1

  Float_t baseline = 1;
  if (awaySide == 0)
    baseline = 0.05;
  else if (awaySide == 1)
  {
    if (iaa == 0 && central || iaa == 1)
      baseline = 0.2;
    else if (iaa == 0 && !central)
      baseline = 0.05;
    else if (iaa == 2 && central)
      baseline = 0.1;
  }
    
  /*
  Float_t reference = 0;
  if (iaa_icp == 0)
    reference = 0.13;
  */
    
  Float_t efficiency = 0.08;
  if (iaa_icp == 1)
    efficiency = 0.05;
  
  Float_t centrality = 1;
  if (iaa_icp == 0)
    centrality = 0.04;
  else if (iaa_icp == 1)
    centrality = 0.06;
  
  Float_t ptResolution = 0.03;
  
  Float_t integrationWindow = 0;
  if (awaySide == 1)
    integrationWindow = 0.03;
  
  Float_t total = TMath::Sqrt(baseline * baseline + efficiency * efficiency + centrality * centrality + ptResolution * ptResolution + integrationWindow * integrationWindow);
  
  Printf("Total syst: %f", total);
  
  systUnc = (TGraphErrors*) graph->Clone();
  for (Int_t i=0; i<systUnc->GetN(); i++)
  {
    systUnc->GetEY()[i] = systUnc->GetY()[i] * total;
    systUnc->GetEX()[i] = 0;
  }
    
  //systUnc->Print();
  
  systUnc->SetLineColor(kGray + 1);
  systUnc->SetLineWidth(6);
  
  systUnc->Draw("PSAME");
    
  return systUnc;
}

TGraphErrors* DrawSystUnc(TGraphErrors* graph, Int_t iaa_icp, Int_t awaySide, Int_t iaa, Bool_t central, Float_t shift = 0.0)
{
  // iaa_icp: iaa = 0, icp = 1

  Float_t baseline = 1;
  if (awaySide == 0)
    baseline = 0.05;
  else if (awaySide == 1)
  {
    if (iaa == 0 && central || iaa == 1)
      baseline = 0.2;
    else if (iaa == 0 && !central)
      baseline = 0.05;
    else if (iaa == 2 && central)
      baseline = 0.1;
  }
    
  /*
  Float_t reference = 0;
  if (iaa_icp == 0)
    reference = 0.13;
  */
    
  Float_t efficiency = 0.04;
  
  Float_t centrality = 1;
  if (iaa_icp == 0)
    centrality = 0.02;
  else if (iaa_icp == 1)
    centrality = 0.03;
  
  Float_t ptResolution = 0;
  
  Float_t integrationWindow = 0;
  if (awaySide == 1)
    integrationWindow = 0.03;
  
  Float_t corrections = 1;
  if (iaa_icp == 0)
    corrections = 0.02;
  else if (iaa_icp == 1)
    corrections = 0.01;
  
  Float_t twoTrack = 0.01;
  
  Float_t total = TMath::Sqrt(baseline * baseline + efficiency * efficiency + centrality * centrality + ptResolution * ptResolution + integrationWindow * integrationWindow + corrections * corrections + twoTrack * twoTrack);
  
  Printf("%f %f %f %f %f %f %f --> Total syst: %f", baseline, efficiency, centrality, ptResolution, integrationWindow, corrections, twoTrack, total);
  
  systUnc = (TGraphErrors*) graph->Clone();
  for (Int_t i=0; i<systUnc->GetN(); i++)
  {
    systUnc->GetEY()[i] = systUnc->GetY()[i] * total;
    systUnc->GetEX()[i] = 0;
    Printf("pt = %.2f  y = %.2f +- %.2f (stat.) +- %.2f (syst.)", graph->GetX()[i] - shift, graph->GetY()[i], graph->GetEY()[i], systUnc->GetEY()[i]);
  }
    
  //systUnc->Print();
  
  systUnc->SetLineColor(kGray + 1);
  systUnc->SetLineWidth(4);
  
  systUnc->Draw("PSAME e2");
    
  return systUnc;
}

void IAA(const char* fileName, Int_t iaa, const char* fileNameEtaGap = 0)
{
  // iaa
  // 0 = IAA LHC
  // 1 = ICP
  // 2 = IAA RHIC
  // 3 = IAA LHC with theory
  
  Bool_t showTheory = 0;
  Bool_t showSTAR = 0;

  if (iaa == 3)
  {
    showTheory = kTRUE;
    iaa = 0;
  }
  
  if (iaa == 4)
  {
    showSTAR = kTRUE;
    iaa = 0;
  }
  
  if (kPythiaScalingFactor != 1)
    Printf("Using reference data scaling factor: %f", kPythiaScalingFactor);

  style();
  
  Int_t j = 1;
  Int_t caseId[] = { 18, 23 };

  c = new TCanvas((iaa != 1) ? ((iaa == 0) ? "iaa" : "iaarhic") : "icp", (iaa != 1) ? ((iaa == 0) ? "iaa" : "iaarhic") : "icp", 900, 450);
  c->Range(0, 0, 1, 1);

  TPad* pad1 = new TPad(Form("%s_2", fileName), "", 0, 0, 0.5, 1);
  pad1->Draw();

  TPad* pad2 = new TPad(Form("%s_1", fileName), "", 0.5, 0, 1, 1);
  pad2->Draw();
  
  pad1->cd();
  gPad->SetRightMargin(0);
  gPad->SetTopMargin(0.05);
  gPad->SetLeftMargin(0.13);
  //gPad->SetGridx();
//   gPad->SetGridy();
  
  pad2->cd();
  gPad->SetLeftMargin(0);
  gPad->SetTopMargin(0.05);
  gPad->SetRightMargin(0.13);
  //gPad->SetGridx();
//   gPad->SetGridy();
  
  pad1->cd();
//   dummy = new TH2F("dummy", Form(";p_{T,assoc} (GeV/c);%s", (iaa != 1) ? "I_{AA,Pythia}" : "I_{CP}"), 100, 1.5, 10.5, 1000, 0, 2.9);
  dummy = new TH2F("dummy", Form(";p_{T,assoc} (GeV/#font[12]{c});%s", (iaa != 1) ? "I_{AA}" : "I_{CP}"), 100, 1.5, 10.5, 1000, 0, 2.4);
  dummy->SetStats(0);
//   dummy->GetYaxis()->SetTitleOffset(1.3);
  dummy->GetXaxis()->SetTitleOffset(1.1);
  currentDummy = dummy->DrawCopy();
  
  gStyle->SetTextAlign(13);
  
  latex = new TLatex(0.17, 0.90, "Near-side_{ }");
  latex->SetTextSize(0.04);
  latex->SetNDC();
  latex->Draw();
  
//   if (iaa != 1)
//   {
//     box = new TBox(2, 0.87, 2.5, 1.13);
//     box->SetFillColor(kGray + 1);
//     box->SetLineColor(kGray + 1);
//     box->Draw();
//   
//     box2 = new TBox(2, 0.87 * 1.5, 2.5, 1.13 * 1.5);
//     box2->SetFillColor(kGray + 1);
//     box2->SetLineColor(kGray + 1);
//     box2->Draw();
//   }
  
  pad2->cd();
  currentDummy = dummy->DrawCopy();
  
/*  if (iaa != 1)
    box->Draw();*/
  
  latex = new TLatex(0.05, 0.90, "Away-side");
  latex->SetTextSize(0.04);
  latex->SetNDC();
  latex->Draw();
  
  legend = new TLegend(0.27, (showTheory) ? 0.60 : 0.62, (iaa != 2 && !showSTAR) ? 0.63 : 0.95, 0.77);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  legend->SetTextSize((iaa != 2) ? 0.04 : 0.035);
  
  for (Int_t i=0; i<2; i++)
  {
    nearSideCentral    = GetRatio(fileName, (iaa != 2) ? 0 : 1, (iaa != 1) ? 3 : 2, j, caseId[i], 0);
    if (iaa == 0)
      nearSidePeripheral = GetRatio(fileName, 2, 3, j, caseId[i], 0);
    awaySideCentral    = GetRatio(fileName, (iaa != 2) ? 0 : 1, (iaa != 1) ? 3 : 2, j, caseId[i], 1);
    if (iaa == 0)
      awaySidePeripheral = GetRatio(fileName, 2, 3, j, caseId[i], 1);
    
    RemovePointsBelowX(nearSideCentral, 3);
    RemovePointsBelowX(awaySideCentral, 3);
    RemovePointsAboveX(nearSideCentral, 10);
    RemovePointsAboveX(awaySideCentral, 10);
    
    RemoveXErrors(nearSideCentral); RemoveXErrors(awaySideCentral);
    
    if (iaa == 0)
    {
      RemovePointsBelowX(nearSidePeripheral, 3);
      RemovePointsBelowX(awaySidePeripheral, 3);
      RemovePointsAboveX(nearSidePeripheral, 10);
      RemovePointsAboveX(awaySidePeripheral, 10);
      ShiftPoints(nearSidePeripheral, 0.2);
      ShiftPoints(awaySidePeripheral, 0.2);
      RemoveXErrors(nearSidePeripheral); RemoveXErrors(awaySidePeripheral);
    }
    
    nearSideCentral->SetMarkerStyle(21);
    nearSideCentral->SetLineColor(1);
    nearSideCentral->SetMarkerColor(1);
            
    awaySideCentral->SetMarkerStyle(21);
    awaySideCentral->SetLineColor(1);
    awaySideCentral->SetMarkerColor(1);
    
    if (showTheory)
    {
      nearSideCentral->SetLineColor(2);
      nearSideCentral->SetMarkerColor(2);
      awaySideCentral->SetLineColor(2);
      awaySideCentral->SetMarkerColor(2);
    }
    
    if (iaa == 0)
    {
      nearSidePeripheral->SetMarkerStyle(22);
      nearSidePeripheral->SetLineColor(2);
      nearSidePeripheral->SetMarkerColor(2);
      
      awaySidePeripheral->SetMarkerStyle(22);
      awaySidePeripheral->SetLineColor(2);
      awaySidePeripheral->SetMarkerColor(2);
    }
    
    if (i == 0)
    {
//       TString denominatorStr("Pythia");
      TString denominatorStr("pp");
      if (iaa == 0)
      {
        legend->AddEntry(nearSideCentral, Form("%s0-5% / %s", (showSTAR) ? "ALICE Pb-Pb 2.76 TeV " : "", denominatorStr.Data()), "P");
	if (!showTheory)
	  legend->AddEntry(nearSidePeripheral, Form("%s60-90% / %s", (showSTAR) ? "ALICE Pb-Pb 2.76 TeV " : "", denominatorStr.Data()), "P");
      }
      else if (iaa == 2)
      {
//         legend->AddEntry(nearSideCentral, Form("0-20% / %s", denominatorStr.Data()), "P");
        legend->AddEntry(nearSideCentral, Form("ALICE 8 GeV/#font[12]{c} < p_{T,trig} < 15 GeV/#font[12]{c}", denominatorStr.Data()), "P");
      }
    }
    
    
    const char* style = "PSAME";
    if (i == 1)
      style = "LSAMEX";
    
    awaySideCentral->Print();
    
    pad1->cd();
    if (i == 0)
      DrawSystUnc(nearSideCentral, (iaa == 1), 0, iaa, kTRUE);
    nearSideCentral->DrawClone(style);
    if (iaa == 0 && !showTheory)
    {
      if (i == 0)
        DrawSystUnc(nearSidePeripheral, (iaa == 1), 0, iaa, kFALSE);
      nearSidePeripheral->DrawClone(style);        
    }
    
    pad2->cd();
    if (i == 0)
      DrawSystUnc(awaySideCentral, (iaa == 1), 1, iaa, kTRUE);
    awaySideCentral->DrawClone(style);
    if (iaa == 0 && !showTheory)
    {
      if (i == 0)
        DrawSystUnc(awaySidePeripheral, (iaa == 1), 1, iaa, kFALSE);
      awaySidePeripheral->DrawClone(style);        
    }
  }
  
  if (fileNameEtaGap)
  {
    pad1->cd();
    
    nearSideEtaGapCentral = GetRatio(fileNameEtaGap, (iaa != 2) ? 0 : 1, (iaa != 1) ? 3 : 2, j, caseId[0], 0);
    RemovePointsBelowX(nearSideEtaGapCentral, 3);
    RemovePointsAboveX(nearSideEtaGapCentral, 10);
    RemoveXErrors(nearSideEtaGapCentral); 
    nearSideEtaGapCentral->SetMarkerStyle(25);
    nearSideEtaGapCentral->SetLineColor(1);
    nearSideEtaGapCentral->SetMarkerColor(1);
    nearSideEtaGapCentral->Draw("PSAME");

    if (iaa != 1)
    {
      nearSideEtaGapPeripheral = GetRatio(fileNameEtaGap, 2, 3, j, caseId[0], 0);
      RemovePointsBelowX(nearSideEtaGapPeripheral, 3);
      RemovePointsAboveX(nearSideEtaGapPeripheral, 10);
      RemoveXErrors(nearSideEtaGapPeripheral); 
      ShiftPoints(nearSideEtaGapPeripheral, 0.2);
      nearSideEtaGapPeripheral->SetMarkerStyle(26);
      nearSideEtaGapPeripheral->SetLineColor(2);
      nearSideEtaGapPeripheral->SetMarkerColor(2);
      nearSideEtaGapPeripheral->Draw("PSAME");
    }
  }
  
  if (iaa == 2)
  {
    // IAA, RHIC, PHENIX
    // systematic uncertainty stored in graph with _sys appended
    TFile::Open("rhic/pi0h_graphs.root");
    
    legend2 = new TLegend(0.5, 0.16, 0.96, 0.27);
    legend2->SetFillColor(0);
    
    for (Int_t ptTrigBin=2; ptTrigBin<4; ptTrigBin++)
    {
      Bool_t central = kTRUE;
      for (Int_t i=0; i<3; i+=2)
      {
        // gIAA_<ptTrigBin>_<centBin>_<angularRegionBin>
        rhic_iaa = (TGraphErrors*) gFile->Get(Form("gIAA_%d_%d_%d", ptTrigBin, (central) ? 0 : 1, i));
        rhic_iaa_sys = (TGraphErrors*) gFile->Get(Form("gIAA_%d_%d_%d_sys", ptTrigBin, (central) ? 0 : 1, i));
    
        rhic_iaa->SetMarkerStyle((ptTrigBin == 2) ? 20 : 33);
        rhic_iaa->SetMarkerColor(2);
        rhic_iaa->SetLineColor(2);
        rhic_iaa_sys->SetLineColor(2);
        rhic_iaa_sys->SetMarkerColor(2);
	
	RemovePointsBelowX(rhic_iaa, 2);
	RemovePointsBelowX(rhic_iaa_sys, 2);
	
	Printf("RHIC:");
	rhic_iaa->Print();
	rhic_iaa_sys->Print();
        
        ShiftPoints(rhic_iaa, -0.05 + 0.05 * (ptTrigBin*2-4));
        ShiftPoints(rhic_iaa_sys, -0.05 + 0.05 * (ptTrigBin*2-4));
        
        if (i == 0)
          pad1->cd();
        else if (i == 2)
          pad2->cd();
        rhic_iaa->Draw("PSAME");
        rhic_iaa_sys->Draw("PSAME");
        
        if (i == 0)
          legend->AddEntry(rhic_iaa, Form("PHENIX %s GeV/c < p_{T,trig} < %s GeV/c", (ptTrigBin == 2) ? "7" : "9", (ptTrigBin == 2) ? "9" : "12", (central) ? "0-20" : "20-60"), "P");
//           legend->AddEntry(rhic_iaa, Form("PHENIX %s GeV/c < p_{T,trig} < %s GeV/c %s%% / pp", (ptTrigBin == 2) ? "7" : "9", (ptTrigBin == 2) ? "9" : "12", (central) ? "0-20" : "20-60"), "P");
      }
    }
    
//     pad1->cd();
//     legend2->Draw();
  }

  if (iaa == 0 && showSTAR)
  {
    for (Int_t i=0; i<2; i++)
    {
      const char* centralityStr = "05";
      if (i == 1)
	centralityStr = "4080";
      
      // IAA, RHIC, STAR
      // systematic uncertainty stored in graph with _sys appended
      nearSide = ReadHepdata(Form("rhic/star_iaa_%s_near.txt", centralityStr));
      awaySide = ReadHepdata(Form("rhic/star_iaa_%s_away.txt", centralityStr));
      
      nearSide->SetMarkerStyle(20 + i * 13);
      nearSide->SetMarkerColor(4);
      nearSide->SetLineColor(4);
      
      awaySide->SetMarkerStyle(20 + i * 13);
      awaySide->SetMarkerColor(4);
      awaySide->SetLineColor(4);

      ShiftPoints(nearSide, -0.1 + 0.2 * i);
      ShiftPoints(awaySide, -0.1 + 0.2 * i);
	  
      pad1->cd();
      nearSide->Draw("PSAME");
      pad2->cd();
      awaySide->Draw("PSAME");
        
      legend->AddEntry(nearSide, Form("STAR Au-Au 0.2 TeV %s / dAu", (i == 0) ? "0-5%" : "40-80%"), "P");
    }
  }

  // Theory predictions
  if (showTheory && iaa == 0)
  {
    const char* theoryList[] = { "AdS", "ASW", "YaJEM", "YaJEM-D", "XinNian" };
    Int_t nTheory = 5;
    Int_t markers[] = { 27, 28, 29, 30, 34 };
    
    for (Int_t i=0; i<nTheory; i++)
    {
      nearSide = ReadHepdata(Form("theory/IAA_near_%s.dat", theoryList[i]));
      awaySide = ReadHepdata(Form("theory/IAA_away_%s.dat", theoryList[i]));
      
      nearSide->SetMarkerStyle(markers[i]);
      awaySide->SetMarkerStyle(markers[i]);
      
      RemovePointsBelowX(nearSide, 3);
      RemovePointsBelowX(awaySide, 3);
      
      Float_t shiftBy = (i % 2 == 0) ? -0.2 : 0.2;
      ShiftPoints(nearSide, shiftBy);
      ShiftPoints(awaySide, shiftBy);
      
      pad1->cd();
      nearSide->Draw("PSAME");
      
      pad2->cd();
      awaySide->Draw("PSAME");
      
      if (i == 0)
	  theoryList[i] = "AdS/CFT";
      if (i == 4)
	  theoryList[i] = "X-N Wang";
      legend->AddEntry(nearSide, theoryList[i], "P");
    }
  }
  
  for (Int_t i=0; i<2; i++)
  {
    if (i == 0)
      pad1->cd();
    else
      pad2->cd();
      
    Float_t xC = 0.05;
    if (i == 0)
      xC = 0.17;

    if (iaa == 2)
      latex = new TLatex(xC, 0.84, "0-20% / pp");
    else
      latex = new TLatex(0.35+xC, 0.9, "8 GeV/#font[12]{c} < p_{T,trig} < 15 GeV/#font[12]{c}");
    latex->SetTextSize(0.04);
    latex->SetNDC();
    if (i == 0)
      latex->Draw();
    
    line = new TLine(1.5, 1, 10.5, 1);
    line->SetLineStyle(2);
    line->Draw();
  
    if (iaa == 1)
    {
      latex = new TLatex(xC, 0.84, "0-5% / 60-90%");
      latex->SetTextSize(0.04);
      latex->SetNDC();
      latex->Draw();
    }
  
    latex = new TLatex(0.65+xC, 0.90, "ALICE");
//     latex = new TLatex(0.5+xC, 0.90, "-- work in progress --");
    latex->SetTextSize(0.04);
    latex->SetNDC();
    if (i == 1)
      latex->Draw();
    
    latex = new TLatex(0.35+xC, 0.84, "p_{T,assoc} < p_{T,trig}");
    latex->SetTextSize(0.04);
    latex->SetNDC();
    if (i == 0)
      latex->Draw();
    
    if (iaa == 1 || iaa == 0)
    {
      if (iaa == 0 && showSTAR)
	latex = new TLatex(0.3+xC, 0.90, "|#eta| < 1.0");
      else
	latex = new TLatex(xC + 0.65, 0.20, "|#eta| < 1.0");
      latex->SetTextSize(0.04);
      latex->SetNDC();
      if (i == 0)
        latex->Draw();
    }
    
    if (iaa == 1 || (iaa == 0 && !showSTAR))
    {
      latex = new TLatex(0.5+xC, 0.72, "Pb-Pb 2.76 TeV");
      latex->SetTextSize(0.04);
      latex->SetNDC();
//       latex->Draw();
    }
    
    if (iaa == 2)
    {
      if (i == 0)
	latex = new TLatex(0.5, 0.24, "ALICE: Pb-Pb 2.76 TeV |#eta| < 1.0");
      else
	latex = new TLatex(xC, 0.6, "ALICE: Pb-Pb 2.76 TeV |#eta| < 1.0");
      latex->SetTextSize(0.035);
      latex->SetNDC();
      latex->Draw();
      
      if (i == 0)
	latex = new TLatex(0.5, 0.18, "PHENIX: Au-Au 0.2 TeV |#eta| < 0.35");
      else
	latex = new TLatex(xC, 0.54, "PHENIX: Au-Au 0.2 TeV |#eta| < 0.35");
      latex->SetTextSize(0.035);
      latex->SetNDC();
      latex->Draw();
    }
    
    if (showTheory && iaa == 0 && i == 1)
      latex = new TLatex(xC, 0.18, "Points: flat pedestal");
    else
      latex = new TLatex(xC + 0.5 * i, 0.26, "Points: flat pedestal");
    latex->SetTextSize(0.04);
    latex->SetNDC();
    if (i == 0)
      latex->Draw();
  
    latex = new TLatex(xC + 0.5 * i, 0.20, "Line: v_{2} subtracted");
    latex->SetTextSize(0.04);
    latex->SetNDC();
    if (i == 0)
      latex->Draw();
    
    if (iaa != 1 && i == 1)
    {
      legendClone = (TLegend*) legend->DrawClone();
      if (showTheory && iaa == 0 && i == 0)
	legendClone->GetListOfPrimitives()->RemoveLast();
//       else
// 	legendClone->SetY1(legendClone->GetY1()-0.05);
	
      legend->SetX1(legend->GetX1()-0.12);
      legend->SetX2(legend->GetX2()-0.12);
    }

/*    if (i == 0 && iaa != 2)
      DrawALICELogo(0.83,0.15,0.98,0.30);
    else
      DrawALICELogo(0.62 + xC,0.48,0.77 + xC,0.63);*/
  }
  
  c->SaveAs(Form("%s.eps", c->GetTitle()));
  c->SaveAs(Form("%s.png", c->GetTitle()));
}

void IAAPaper(const char* fileName, Int_t iaa, const char* fileNameEtaGap = 0)
{
  // iaa
  // 0 = IAA LHC
  // 1 = ICP
  
  Bool_t showTheory = 0;
  Bool_t showSTAR = 0;

  if (iaa == 3)
  {
    showTheory = kTRUE;
    iaa = 0;
  }
  
  if (iaa == 4)
  {
    showSTAR = kTRUE;
    iaa = 0;
  }
  
  if (kPythiaScalingFactor != 1)
    Printf("Using reference data scaling factor: %f", kPythiaScalingFactor);

  style();
  
  Int_t j = 1;
  Int_t caseId[] = { 18, 23 };

  c = new TCanvas((iaa != 1) ? ((iaa == 0) ? "iaa" : "iaarhic") : "icp", (iaa != 1) ? ((iaa == 0) ? "iaa" : "iaarhic") : "icp", 900, 450);
  c->Range(0, 0, 1, 1);

  TPad* pad1 = new TPad(Form("%s_2", fileName), "", 0, 0, 0.5, 1);
  pad1->Draw();

  TPad* pad2 = new TPad(Form("%s_1", fileName), "", 0.5, 0, 1, 1);
  pad2->Draw();
  
  pad1->cd();
  gPad->SetRightMargin(0);
  gPad->SetTopMargin(0.05);
  gPad->SetLeftMargin(0.13);
  //gPad->SetGridx();
//   gPad->SetGridy();
  
  pad2->cd();
  gPad->SetLeftMargin(0);
  gPad->SetTopMargin(0.05);
  gPad->SetRightMargin(0.13);
  //gPad->SetGridx();
//   gPad->SetGridy();
  
  pad1->cd();
//   dummy = new TH2F("dummy", Form(";p_{T,assoc} (GeV/c);%s", (iaa != 1) ? "I_{AA,Pythia}" : "I_{CP}"), 100, 1.5, 10.5, 1000, 0, 2.9);
  dummy = new TH2F("dummy", Form(";#font[12]{p}_{t,assoc} (GeV/#font[12]{c});%s", (iaa != 1) ? "#font[12]{I}_{AA}" : "#font[12]{I}_{CP} (0-5% / 60-90%)"), 100, 1.5, 10.5, 1000, 0, 2.4);
  dummy->SetStats(0);
//   dummy->GetYaxis()->SetTitleOffset(1.3);
  dummy->GetXaxis()->SetTitleOffset(1.1);
  currentDummy = dummy->DrawCopy();
  
  gStyle->SetTextAlign(13);
  
  latex = new TLatex(0.17, 0.90, "Near-side_{ }");
  latex->SetTextSize(0.04);
  latex->SetNDC();
  latex->Draw();
  
  latex = new TLatex(0.17, /*(iaa == 1) ? 0.78 :*/ 0.84, "#sqrt{#font[12]{s}_{NN}} = 2.76 TeV");
  latex->SetTextSize(0.04);
  latex->SetNDC();
  latex->Draw();

  latex = new TLatex(0.17, 0.2, (iaa == 1) ? "b)" : "a)");
  latex->SetTextSize(0.04);
  latex->SetNDC();
  latex->Draw();
  
//   if (iaa != 1)
//   {
//     box = new TBox(2, 0.87, 2.5, 1.13);
//     box->SetFillColor(kGray + 1);
//     box->SetLineColor(kGray + 1);
//     box->Draw();
//   
//     box2 = new TBox(2, 0.87 * 1.5, 2.5, 1.13 * 1.5);
//     box2->SetFillColor(kGray + 1);
//     box2->SetLineColor(kGray + 1);
//     box2->Draw();
//   }
  
  pad2->cd();
  currentDummy = dummy->DrawCopy();
  
/*  if (iaa != 1)
    box->Draw();*/
  
  latex = new TLatex(0.05, 0.90, "Away-side");
  latex->SetTextSize(0.04);
  latex->SetNDC();
  latex->Draw();

  if (iaa == 1)
    legend = new TLegend(0.22, 0.60, 0.67, 0.73);
  else
    legend = new TLegend(0.12, 0.60, 0.76, 0.79);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  legend->SetTextSize((iaa != 2) ? 0.04 : 0.035);
  
  for (Int_t i=0; i<2; i++)
  {
    nearSideCentral    = GetRatio(fileName, (iaa != 2) ? 0 : 1, (iaa != 1) ? 3 : 2, j, caseId[i], 0);
    if (iaa == 0)
      nearSidePeripheral = GetRatio(fileName, 2, 3, j, caseId[i], 0);
    awaySideCentral    = GetRatio(fileName, (iaa != 2) ? 0 : 1, (iaa != 1) ? 3 : 2, j, caseId[i], 1);
    if (iaa == 0)
      awaySidePeripheral = GetRatio(fileName, 2, 3, j, caseId[i], 1);
    
    RemovePointsBelowX(nearSideCentral, 3);
    RemovePointsBelowX(awaySideCentral, 3);
    RemovePointsAboveX(nearSideCentral, 10);
    RemovePointsAboveX(awaySideCentral, 10);

    ShiftPoints(nearSideCentral, -0.05);
    ShiftPoints(awaySideCentral, -0.05);

    if (i == 1)
    {
      ShiftPoints(nearSideCentral, -0.2);
      ShiftPoints(awaySideCentral, -0.2);
      
      if (iaa == 0)
      {
	ShiftPoints(nearSidePeripheral, 0.2);
	ShiftPoints(awaySidePeripheral, 0.2);
      }
    }
    
    RemoveXErrors(nearSideCentral); RemoveXErrors(awaySideCentral);
    
    if (iaa == 0)
    {
      RemovePointsBelowX(nearSidePeripheral, 3);
      RemovePointsBelowX(awaySidePeripheral, 3);
      RemovePointsAboveX(nearSidePeripheral, 10);
      RemovePointsAboveX(awaySidePeripheral, 10);
      ShiftPoints(nearSidePeripheral, 0.05);
      ShiftPoints(awaySidePeripheral, 0.05);
      RemoveXErrors(nearSidePeripheral); RemoveXErrors(awaySidePeripheral);
    }
    
    nearSideCentral->SetMarkerStyle((i == 0) ? 25 : 27);
    nearSideCentral->SetLineColor(1);
    nearSideCentral->SetMarkerColor(1);
            
    awaySideCentral->SetMarkerStyle((i == 0) ? 25 : 27);
    awaySideCentral->SetLineColor(1);
    awaySideCentral->SetMarkerColor(1);
    
    if (showTheory)
    {
      nearSideCentral->SetLineColor(2);
      nearSideCentral->SetMarkerColor(2);
      awaySideCentral->SetLineColor(2);
      awaySideCentral->SetMarkerColor(2);
    }
    
    if (iaa == 0)
    {
      nearSidePeripheral->SetMarkerStyle((i == 0) ? 21 : 33);
      nearSidePeripheral->SetLineColor(2);
      nearSidePeripheral->SetMarkerColor(2);
      
      awaySidePeripheral->SetMarkerStyle((i == 0) ? 21 : 33);
      awaySidePeripheral->SetLineColor(2);
      awaySidePeripheral->SetMarkerColor(2);
    }
    
//     if (i == 0)
    {
//       TString denominatorStr("Pythia");
      TString denominatorStr("pp");
      if (iaa == 0)
      {
	legend->SetHeader("0-5% Pb-Pb/pp     60-90% Pb-Pb/pp");
	legend->SetNColumns(2);
      }
      legend->AddEntry(nearSideCentral, (i == 0) ? "Flat bkg" : "#font[12]{v}_{2} bkg", "P");
      if (iaa == 0)
	legend->AddEntry(nearSidePeripheral, (i == 0) ? "Flat bkg" : "#font[12]{v}_{2} bkg", "P");
	
//         legend->AddEntry(nearSideCentral, Form("%s0-5% / %s", (showSTAR) ? "ALICE Pb-Pb 2.76 TeV " : "", denominatorStr.Data()), "P");
// 	if (!showTheory)
// 	  legend->AddEntry(nearSidePeripheral, Form("%s60-90% / %s", (showSTAR) ? "ALICE Pb-Pb 2.76 TeV " : "", denominatorStr.Data()), "P");
      if (iaa == 2)
      {
//         legend->AddEntry(nearSideCentral, Form("0-20% / %s", denominatorStr.Data()), "P");
        legend->AddEntry(nearSideCentral, Form("ALICE 8 GeV/#font[12]{c} < p_{T,trig} < 15 GeV/#font[12]{c}", denominatorStr.Data()), "P");
      }
    }
    
    
    const char* style = "PSAME";
//     if (i == 1)
//       style = "LSAMEX";
    
    pad1->cd();
//     if (i == 0)
      DrawSystUnc(nearSideCentral, (iaa == 1), 0, iaa, kTRUE, -0.05 - 0.2 * i);
    nearSideCentral->DrawClone(style);
    if (iaa == 0 && !showTheory)
    {
//       if (i == 0)
        DrawSystUnc(nearSidePeripheral, (iaa == 1), 0, iaa, kFALSE, 0.05 + ((iaa == 0) ? i * 0.2 : 0));
      nearSidePeripheral->DrawClone(style);        
    }
    
    pad2->cd();
//     if (i == 0)
      DrawSystUnc(awaySideCentral, (iaa == 1), 1, iaa, kTRUE, -0.05 - 0.2 * i);
    awaySideCentral->DrawClone(style);
    if (iaa == 0 && !showTheory)
    {
//       if (i == 0)
        DrawSystUnc(awaySidePeripheral, (iaa == 1), 1, iaa, kFALSE, 0.05 + ((iaa == 0) ? i * 0.2 : 0));
      awaySidePeripheral->DrawClone(style);        
    }
  }
  
  if (fileNameEtaGap)
  {
    pad1->cd();
    
    nearSideEtaGapCentral = GetRatio(fileNameEtaGap, (iaa != 2) ? 0 : 1, (iaa != 1) ? 3 : 2, j, caseId[0], 0);
    RemovePointsBelowX(nearSideEtaGapCentral, 3);
    RemovePointsAboveX(nearSideEtaGapCentral, 10);
    RemoveXErrors(nearSideEtaGapCentral); 
    nearSideEtaGapCentral->SetMarkerStyle(24);
    nearSideEtaGapCentral->SetLineColor(1);
    nearSideEtaGapCentral->SetMarkerColor(1);
    ShiftPoints(nearSideEtaGapCentral, -0.45);
    DrawSystUnc(nearSideEtaGapCentral, (iaa == 1), 0, iaa, kTRUE, -0.45);
    nearSideEtaGapCentral->Draw("PSAME");

    legend->AddEntry(nearSideEtaGapCentral, "#font[12]{#eta}-gap", "P");


    if (iaa != 1)
    {
      nearSideEtaGapPeripheral = GetRatio(fileNameEtaGap, 2, 3, j, caseId[0], 0);
      RemovePointsBelowX(nearSideEtaGapPeripheral, 3);
      RemovePointsAboveX(nearSideEtaGapPeripheral, 10);
      RemoveXErrors(nearSideEtaGapPeripheral); 
      nearSideEtaGapPeripheral->SetMarkerStyle(20);
      nearSideEtaGapPeripheral->SetLineColor(2);
      nearSideEtaGapPeripheral->SetMarkerColor(2);
      ShiftPoints(nearSideEtaGapPeripheral, 0.45);
      DrawSystUnc(nearSideEtaGapPeripheral, (iaa == 1), 0, iaa, kFALSE, 0.45);
      nearSideEtaGapPeripheral->Draw("PSAME");

      legend->AddEntry(nearSideEtaGapPeripheral, "#font[12]{#eta}-gap", "P");
    }
  }
  
  if (iaa == 0)
  {
    // add pythia line showing gluon filtering
    // from Jan Rak, by mail, 22.09.11
    // pp CF calculated without any adjustment
    // AA CF also no adjustment. Particles (h+-) in acceptance are identified with their parent string parton (status 71)
    // and the flavor of parent hard outgoing parton (status 23) which is  the parent of string parton 71 is checked.
    // If it is a gluon the  particle is ignored.
    
    Float_t xVal[] = { 3.5, 4.5, 5.5, 7, 9};
    Float_t yVal[] = { 0.779418, 0.794234, 0.82759, 0.856457, 0.880526 };
    gluonfiltering = new TGraph(5, xVal, yVal);
    pad1->cd();
    gluonfiltering->Draw("L");
    
    legend3 = new TLegend(0.51, 0.15, 0.95, 0.25);
    legend3->SetFillColor(0);
    legend3->SetBorderSize(0);
    legend3->SetTextSize(0.04);
    legend3->AddEntry(gluonfiltering, "Pythia gluon filtering", "L");
    legend3->Draw();
  }
  
  if (iaa == 2)
  {
    // IAA, RHIC, PHENIX
    // systematic uncertainty stored in graph with _sys appended
    TFile::Open("rhic/pi0h_graphs.root");
    
    legend2 = new TLegend(0.5, 0.16, 0.96, 0.27);
    legend2->SetFillColor(0);
    
    for (Int_t ptTrigBin=2; ptTrigBin<4; ptTrigBin++)
    {
      Bool_t central = kTRUE;
      for (Int_t i=0; i<3; i+=2)
      {
        // gIAA_<ptTrigBin>_<centBin>_<angularRegionBin>
        rhic_iaa = (TGraphErrors*) gFile->Get(Form("gIAA_%d_%d_%d", ptTrigBin, (central) ? 0 : 1, i));
        rhic_iaa_sys = (TGraphErrors*) gFile->Get(Form("gIAA_%d_%d_%d_sys", ptTrigBin, (central) ? 0 : 1, i));
    
        rhic_iaa->SetMarkerStyle((ptTrigBin == 2) ? 20 : 33);
        rhic_iaa->SetMarkerColor(2);
        rhic_iaa->SetLineColor(2);
        rhic_iaa_sys->SetLineColor(2);
        rhic_iaa_sys->SetMarkerColor(2);
        
        ShiftPoints(rhic_iaa, -0.05 + 0.05 * (ptTrigBin*2-4));
        ShiftPoints(rhic_iaa_sys, -0.05 + 0.05 * (ptTrigBin*2-4));
        
        if (i == 0)
          pad1->cd();
        else if (i == 2)
          pad2->cd();
        rhic_iaa->Draw("PSAME");
        rhic_iaa_sys->Draw("PSAME");
        
        if (i == 0)
          legend->AddEntry(rhic_iaa, Form("PHENIX %s GeV/c < p_{T,trig} < %s GeV/c", (ptTrigBin == 2) ? "7" : "9", (ptTrigBin == 2) ? "9" : "12", (central) ? "0-20" : "20-60"), "P");
//           legend->AddEntry(rhic_iaa, Form("PHENIX %s GeV/c < p_{T,trig} < %s GeV/c %s%% / pp", (ptTrigBin == 2) ? "7" : "9", (ptTrigBin == 2) ? "9" : "12", (central) ? "0-20" : "20-60"), "P");
      }
    }
    
//     pad1->cd();
//     legend2->Draw();
  }

  if (iaa == 0 && showSTAR)
  {
    for (Int_t i=0; i<2; i++)
    {
      const char* centralityStr = "05";
      if (i == 1)
	centralityStr = "4080";
      
      // IAA, RHIC, STAR
      // systematic uncertainty stored in graph with _sys appended
      nearSide = ReadHepdata(Form("rhic/star_iaa_%s_near.txt", centralityStr));
      awaySide = ReadHepdata(Form("rhic/star_iaa_%s_away.txt", centralityStr));
      
      nearSide->SetMarkerStyle(20 + i * 13);
      nearSide->SetMarkerColor(4);
      nearSide->SetLineColor(4);
      
      awaySide->SetMarkerStyle(20 + i * 13);
      awaySide->SetMarkerColor(4);
      awaySide->SetLineColor(4);

      ShiftPoints(nearSide, -0.1 + 0.2 * i);
      ShiftPoints(awaySide, -0.1 + 0.2 * i);
	  
      pad1->cd();
      nearSide->Draw("PSAME");
      pad2->cd();
      awaySide->Draw("PSAME");
        
      legend->AddEntry(nearSide, Form("STAR Au-Au 0.2 TeV %s / dAu", (i == 0) ? "0-5%" : "40-80%"), "P");
    }
  }

  // Theory predictions
  if (showTheory && iaa == 0)
  {
    const char* theoryList[] = { "AdS", "ASW", "YaJEM", "YaJEM-D", "XinNian" };
    Int_t nTheory = 5;
    Int_t markers[] = { 27, 28, 29, 30, 34 };
    
    for (Int_t i=0; i<nTheory; i++)
    {
      nearSide = ReadHepdata(Form("theory/IAA_near_%s.dat", theoryList[i]));
      awaySide = ReadHepdata(Form("theory/IAA_away_%s.dat", theoryList[i]));
      
      nearSide->SetMarkerStyle(markers[i]);
      awaySide->SetMarkerStyle(markers[i]);
      
      RemovePointsBelowX(nearSide, 3);
      RemovePointsBelowX(awaySide, 3);
      
      Float_t shiftBy = (i % 2 == 0) ? -0.2 : 0.2;
      ShiftPoints(nearSide, shiftBy);
      ShiftPoints(awaySide, shiftBy);
      
      pad1->cd();
      nearSide->Draw("PSAME");
      
      pad2->cd();
      awaySide->Draw("PSAME");
      
      if (i == 0)
	  theoryList[i] = "AdS/CFT";
      if (i == 4)
	  theoryList[i] = "X-N Wang";
      legend->AddEntry(nearSide, theoryList[i], "P");
    }
  }
  
  for (Int_t i=0; i<2; i++)
  {
    if (i == 0)
      pad1->cd();
    else
      pad2->cd();
      
    Float_t xC = 0.05;
    if (i == 0)
      xC = 0.17;

    if (iaa == 2)
      latex = new TLatex(xC, 0.84, "0-20% / pp");
    else
      latex = new TLatex(0.35+xC, 0.9, "8 GeV/#font[12]{c} < #font[12]{p}_{t,trig} < 15 GeV/#font[12]{c}");
    latex->SetTextSize(0.04);
    latex->SetNDC();
    if (i == 0)
      latex->Draw();
    
    line = new TLine(1.5, 1, 10.5, 1);
    line->SetLineStyle(2);
    line->Draw();
  
//     if (iaa == 1 && i == 0)
//     {
//       latex = new TLatex(xC, 0.84, "0-5% / 60-90%");
//       latex->SetTextSize(0.04);
//       latex->SetNDC();
//       latex->Draw();
//     }
  
    latex = new TLatex(0.65+xC, 0.90, "ALICE");
//     latex = new TLatex(0.5+xC, 0.90, "-- work in progress --");
    latex->SetTextSize(0.04);
    latex->SetNDC();
    if (i == 1)
      latex->Draw();
    
    latex = new TLatex(0.35+xC, 0.84, "#font[12]{p}_{t,assoc} < #font[12]{p}_{t,trig}");
    latex->SetTextSize(0.04);
    latex->SetNDC();
    if (i == 0)
      latex->Draw();
    
    if (iaa == 1 || iaa == 0)
    {
      if (iaa == 0 && showSTAR)
	latex = new TLatex(0.3+xC, 0.90, "|#font[12]{#eta}| < 1.0");
      else
	latex = new TLatex(0.82, 0.84, "|#font[12]{#eta}| < 1.0");
      latex->SetTextSize(0.04);
      latex->SetNDC();
      if (i == 0)
        latex->Draw();
    }
    
    if (iaa == 1 || (iaa == 0 && !showSTAR))
    {
      latex = new TLatex(0.5+xC, 0.72, "Pb-Pb 2.76 TeV");
      latex->SetTextSize(0.04);
      latex->SetNDC();
//       latex->Draw();
    }
    
    if (iaa == 2)
    {
      if (i == 0)
	latex = new TLatex(0.5, 0.24, "ALICE: Pb-Pb 2.76 TeV |#eta| < 1.0");
      else
	latex = new TLatex(xC, 0.6, "ALICE: Pb-Pb 2.76 TeV |#eta| < 1.0");
      latex->SetTextSize(0.035);
      latex->SetNDC();
      latex->Draw();
      
      if (i == 0)
	latex = new TLatex(0.5, 0.18, "PHENIX: Au-Au 0.2 TeV |#eta| < 0.35");
      else
	latex = new TLatex(xC, 0.54, "PHENIX: Au-Au 0.2 TeV |#eta| < 0.35");
      latex->SetTextSize(0.035);
      latex->SetNDC();
      latex->Draw();
    }
    
    if (i == 1)
    {
      legendClone = (TLegend*) legend->DrawClone();
      if (showTheory && iaa == 0 && i == 0)
	legendClone->GetListOfPrimitives()->RemoveLast();
//       else
// 	legendClone->SetY1(legendClone->GetY1()-0.05);
	
      legend->SetX1(legend->GetX1()-0.12);
      legend->SetX2(legend->GetX2()-0.12);
    }
    
    

/*    if (i == 0 && iaa != 2)
      DrawALICELogo(0.83,0.15,0.98,0.30);
    else
      DrawALICELogo(0.62 + xC,0.48,0.77 + xC,0.63);*/
  }
  
  c->SaveAs(Form("%s.eps", c->GetTitle()));
  c->SaveAs(Form("%s.png", c->GetTitle()));
}


void DrawALICELogo(Float_t x1, Float_t y1, Float_t x2, Float_t y2, Bool_t debug = kFALSE)
{
  TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo", x1, y1, x2, y2);
  if (debug)
    myPadLogo->SetFillColor(2); // color to first figure out where is the pad then comment !
  myPadLogo->SetLeftMargin(0);
  myPadLogo->SetTopMargin(0);
  myPadLogo->SetRightMargin(0);
  myPadLogo->SetBottomMargin(0);
  myPadLogo->Draw();
  myPadLogo->cd();
  TASImage *myAliceLogo = new TASImage("~/alice_logo_transparent.png");
  myAliceLogo->Draw();
}

void RemoveXErrors(TGraphErrors* graph)
{
	for (Int_t i=0; i<graph->GetN(); i++)
		graph->GetEX()[i] = 0;
}

void RemovePointsBelowX(TGraphErrors* graph, Float_t minX)
{
  Int_t i=0;
  while (i<graph->GetN())
  {
    if (graph->GetX()[i] < minX)
      graph->RemovePoint(i);
    else
      i++;
  }
}

void RemovePointsAboveX(TGraphErrors* graph, Float_t maxX)
{
  Int_t i=0;
  while (i<graph->GetN())
  {
    if (graph->GetX()[i] > maxX)
      graph->RemovePoint(i);
    else
      i++;
  }
}

void ShiftPoints(TGraphErrors* graph, Float_t dx)
{
  Int_t i=0;
  while (i<graph->GetN())
  {
    graph->GetX()[i] += dx;
    i++;
  }
}

void DivideGraphs(TGraphErrors* graph1, TGraphErrors* graph2)
{
  graph1->Print();
  graph2->Print();

  for (Int_t bin1 = 0; bin1 < graph1->GetN(); bin1++)
  {
    Float_t x = graph1->GetX()[bin1];

    Int_t bin2 = 0;
    for (bin2 = 0; bin2<graph2->GetN(); bin2++)
      if (graph2->GetX()[bin2] >= x)
        break;

    if (bin2 == graph2->GetN())
            bin2--;

    if (bin2 > 0)
      if (TMath::Abs(graph2->GetX()[bin2-1] - x) < TMath::Abs(graph2->GetX()[bin2] - x))
        bin2--;

    if (graph1->GetY()[bin1] == 0 || graph2->GetY()[bin2] == 0 || bin2 == graph2->GetN())
    {
      Printf("%d %d removed", bin1, bin2);
      graph1->RemovePoint(bin1--);
      continue;
    }

    Float_t graph2Extrapolated = graph2->GetY()[bin2];
    if (TMath::Abs(x - graph2->GetX()[bin2]) > 0.0001)
    {
      Printf("%f %f %d %d not found", x, graph2->GetX()[bin2], bin1, bin2);
      graph1->RemovePoint(bin1--);
      continue;
    }

    Float_t value = graph1->GetY()[bin1] / graph2Extrapolated;
    Float_t error = value * TMath::Sqrt(TMath::Power(graph1->GetEY()[bin1] / graph1->GetY()[bin1], 2) + TMath::Power(graph2->GetEY()[bin2] / graph2->GetY()[bin2], 2));

    graph1->GetY()[bin1] = value;
    graph1->GetEY()[bin1] = error;

    Printf("%d %d %f %f %f %f", bin1, bin2, x, graph2Extrapolated, value, error);
  }
}

Float_t ExtractYields(TH1* hist, Int_t trigId, Int_t centralityId, Float_t ptA, Float_t ptW, Int_t caseId, Double_t normUnc)
{
  CreateYieldStructure();

  Float_t width = 0.69;
//   Float_t width = 0.5; Printf("WARNING. Using integration width 0.5");
//   Float_t width = 0.85; Printf("WARNING. Using integration width 0.85");
  
  //new TCanvas; hist->DrawCopy();
  
  Int_t binBegin = hist->FindBin(-width);
  Int_t binEnd = hist->FindBin(width);
  Double_t nearSideError = 0;
  Float_t nearSideYield = hist->IntegralAndError(binBegin, binEnd, nearSideError);
  nearSideYield *= hist->GetXaxis()->GetBinWidth(1);
  nearSideError *= hist->GetXaxis()->GetBinWidth(1);

  Int_t binBegin = hist->FindBin(TMath::Pi() - width);
  Int_t binEnd = hist->FindBin(TMath::Pi() + width);
  Double_t awaySideError = 0;
  Float_t awaySideYield = hist->IntegralAndError(binBegin, binEnd, awaySideError);
  awaySideYield *= hist->GetXaxis()->GetBinWidth(1);
  awaySideError *= hist->GetXaxis()->GetBinWidth(1);
  
  Printf("%d %d %f %d: %f+-%f %f+-%f (%f)", trigId, centralityId, ptA, caseId, nearSideYield, nearSideError, awaySideYield, awaySideError, normUnc);
  
  FillYield(trigId, centralityId, ptA, ptW, caseId, nearSideYield, TMath::Sqrt(nearSideError * nearSideError + normUnc * normUnc), awaySideYield, TMath::Sqrt(awaySideError * awaySideError + normUnc * normUnc));
  
  // store normUnc for error propagation later
  TGraphErrors** tmp = yields[2][trigId][centralityId];
  tmp[caseId]->SetPoint(tmp[caseId]->GetN(), ptA, normUnc);
  
  return awaySideYield;
}

void FillYield(Int_t trigId, Int_t centralityId, Float_t ptA, Float_t ptW, Int_t caseId, Double_t nearSideYield, Double_t nearSideError, Double_t awaySideYield, Double_t awaySideError)
{
  // CINT limitation here
  TGraphErrors** tmp = yields[0][trigId][centralityId];
  tmp[caseId]->SetPoint(tmp[caseId]->GetN(), ptA, nearSideYield);
  tmp[caseId]->SetPointError(tmp[caseId]->GetN()-1, ptW, nearSideError);
  //tmp[caseId]->Print();
  
  tmp = yields[1][trigId][centralityId];
  tmp[caseId]->SetPoint(tmp[caseId]->GetN(), ptA, awaySideYield);
  tmp[caseId]->SetPointError(tmp[caseId]->GetN()-1, ptW, awaySideError);
}

Double_t hypgeo(Double_t a, Double_t b, Double_t c, Double_t z)
{
  // ROOT::Math::hyperg is only defined for |z| < 1
  if (z > -1)
    return ROOT::Math::hyperg(a, b, c, z);
    
  // for z < -1 we can use an identity from http://en.wikipedia.org/wiki/Hypergeometric_function
  return TMath::Power(1-z, -b) * ROOT::Math::hyperg(c-a, b, c, z / (z-1));
}

Double_t hypgeoder2(Double_t a, Double_t b, Double_t c, Double_t z)
{
  // derivate of hypgeo with parameter b: d/db hypgeo
  
  static TF1* func = 0;
  if (!func)
    func = new TF1("hypgeoder2_func", "hypgeo([0], x, [1], [2])");
  
  func->SetParameters(a, c, z);
  return func->Derivative(b);
}

void TsallisYieldUncertainty(Double_t n, Double_t q, Double_t beta, TMatrixDSym covMatrix, Int_t offset, Double_t& yield, Double_t& unc)
{
  gSystem->Load("libMathMore");

  TF1 tsallis("tsallis", "[0] * (1-[2]*(1-[1])*x*x)**(1/(1-[1]))", -0.5 * TMath::Pi(), 1.5 * TMath::Pi());

  // Integrate[N (1 - b (1 - q) x^2)^(1/(1 - q)), x] --> N x Hypergeometric2F1[1/2, 1/(-1 + q), 3/2, -b (-1 + q) x^2]
  // [0] = n, [1] = q, [2] = beta
  TF1 tsallisIntegral("tsallisIntegral",               "[0]*x*hypgeo(0.5, 1/([1]-1), 1.5, -[2]*([1]-1)*x*x)");
  TF1 tsallisIntegralDerN("tsallisIntegralDerN",           "x*hypgeo(0.5, 1/([1]-1), 1.5, -[2]*([1]-1)*x*x)");
  
  TF1 tsallisIntegralDerq("tsallisIntegralDerq",       "[0]*x*(((1+[2]*([1]-1)*x*x)**(-1/([1]-1))-hypgeo(0.5, 1/([1]-1), 1.5, -[2]*([1]-1)*x*x))/2/([1]-1) - hypgeoder2(0.5, 1/([1]-1), 1.5, -[2]*([1]-1)*x*x)/([1]-1)**2)");
  
  TF1 tsallisIntegralDerbeta("tsallisIntegralDerbeta", "[0]*x*((1+[2]*([1]-1)*x*x)**(-1/([1]-1))-hypgeo(0.5, 1/([1]-1), 1.5, -[2]*([1]-1)*x*x))/2/[2]");
  
  const Float_t width = 0.7;
  
  TF1* derivates[4];
  derivates[0] = new TF1("const", "1");
  derivates[1] = &tsallisIntegralDerN;
  derivates[2] = &tsallisIntegralDerq;
  derivates[3] = &tsallisIntegralDerbeta;
  
  Printf("%f %f %f", n, q, beta);
  
  tsallis.SetParameters(n, q, beta);
  tsallisIntegral.SetParameters(n, q, beta);
  for (Int_t i=1; i<4; i++)
    derivates[i]->SetParameters(n, q, beta);
  
  yield = tsallisIntegral.Eval(width) - tsallisIntegral.Eval(-width);
  Double_t yield2 = tsallis.Integral(-width, width);
  
  unc = 0;
  for (Int_t i=0; i<4; i++)
    for (Int_t j=0; j<4; j++)
    {
      Int_t covI = i;
      Int_t covJ = j;
      if (i > 0)
        covI += offset;
      if (j > 0)
        covJ += offset;
        
      Double_t der1 = (derivates[i]->Eval(width) - derivates[i]->Eval(-width));
      Double_t der2 = (derivates[j]->Eval(width) - derivates[j]->Eval(-width));
        
      unc += der1 * der2 * covMatrix(covI, covJ);
      
      Printf("%d %d % .3e % .3e % .3e % .3e", i, j, der1, der2, covMatrix(covI, covJ), unc);
    }
    
  unc = TMath::Sqrt(unc);
  
  yield /= width * 2;
  yield2 /= width * 2;
  unc /= width * 2;
  
  Printf("%f (%f) +- %f", yield, yield2, unc);
}

double NLowestBinAverage(TH1 &h, int N)
{
  // Returns mean content of N lowest bins in h.

  double avg = 0;
  int counter = 0;
  std::list<double> content;
  std::list<double>::iterator i;

  for (int jbin=1; jbin<=h.GetNbinsX(); jbin++) {
    content.push_back(h.GetBinContent(jbin));
  }
  content.sort();

  for (i = content.begin(); i != content.end(); ++i) {

    if (counter == N)
      break;

    avg += *i/N;
    counter++;
  }

  if (counter < N)
    Warning("NLowestBinAverage()",
        "Avg. of %d bins requested, %d used", N, counter);

  return avg;
} 

void DrawFlow(Float_t v2, TH1* hist, Float_t ptT, Float_t ptA, TH1* histMixed, Int_t trigId, Int_t centralityId, Int_t caseOffset, Float_t ptACenter, Float_t ptAWidth, Bool_t flatBaseLine = kFALSE, Int_t baseLineDetermination = 0, Float_t* vn = 0)
{
  // caseOffsets
  // 0 same with v2 (14 tsallis for baseline [28/29 tsallis params]; 15 yield from tsallis function)
  // [disabled] 1 same/mixed with v2 (15 tsallis for baseline)
  // 2 same no v2 (16 tsallis for baseline; 17 yield from tsallis function)
  // 4 same no v2 (18 flat fit 1)
  // 5 same no v2 (19 flat fit 2)
  // 6 same no v2 (20 flat avg over 4)
  // 7 same no v2 (21 flat avg over 8)
  // 8 same no v2 (22 flat avg over 16)
  // 9 same with v2 (23 flat fit 1)
  // 10 same with v2 (24 flat fit 2)
  // 11 same with v2 (25 flat avg over 4)
  // 12 same with v2 (26 flat avg over 8)
  // 13 same with v2 (27 flat avg over 16)
  
  // 30 gaussian fit width

  // same stuff for same/mixed
  if (0 && caseOffset == 0 && histMixed)
  {
    TH1* clone = (TH1*) hist->Clone();
    clone->Divide(histMixed);
    DrawFlow(v2, clone, ptT, ptA, histMixed, trigId, centralityId, 1, ptACenter, ptAWidth, kFALSE, 0, vn);
  }
  
  // same for other baseline subtractions
  if (caseOffset == 0)
  {
    file = TFile::Open("dphi_corr.root", "UPDATE");
    hist->Write();
    file->Close();
    
    TH1* clone = (TH1*) hist->Clone();
    DrawFlow(v2, clone, ptT, ptA, histMixed, trigId, centralityId, 2, ptACenter, ptAWidth, kTRUE, 0, vn);
  
    clone = (TH1*) hist->Clone();
    DrawFlow(v2, clone, ptT, ptA, histMixed, trigId, centralityId, 4, ptACenter, ptAWidth, kTRUE, 1, vn);
  
    clone = (TH1*) hist->Clone();
    DrawFlow(v2, clone, ptT, ptA, histMixed, trigId, centralityId, 5, ptACenter, ptAWidth, kTRUE, 2, vn);
    
    clone = (TH1*) hist->Clone();
    DrawFlow(v2, clone, ptT, ptA, histMixed, trigId, centralityId, 6, ptACenter, ptAWidth, kTRUE, 3, vn);
    
    clone = (TH1*) hist->Clone();
    DrawFlow(v2, clone, ptT, ptA, histMixed, trigId, centralityId, 7, ptACenter, ptAWidth, kTRUE, 4, vn);
    
    clone = (TH1*) hist->Clone();
    DrawFlow(v2, clone, ptT, ptA, histMixed, trigId, centralityId, 8, ptACenter, ptAWidth, kTRUE, 5, vn);
    
    clone = (TH1*) hist->Clone();
    DrawFlow(v2, clone, ptT, ptA, histMixed, trigId, centralityId, 9, ptACenter, ptAWidth, kFALSE, 1, vn);
  
    clone = (TH1*) hist->Clone();
    DrawFlow(v2, clone, ptT, ptA, histMixed, trigId, centralityId, 10, ptACenter, ptAWidth, kFALSE, 2, vn);
  
    clone = (TH1*) hist->Clone();
    DrawFlow(v2, clone, ptT, ptA, histMixed, trigId, centralityId, 11, ptACenter, ptAWidth, kFALSE, 3, vn);
  
    clone = (TH1*) hist->Clone();
    DrawFlow(v2, clone, ptT, ptA, histMixed, trigId, centralityId, 12, ptACenter, ptAWidth, kFALSE, 4, vn);
  
    clone = (TH1*) hist->Clone();
    DrawFlow(v2, clone, ptT, ptA, histMixed, trigId, centralityId, 13, ptACenter, ptAWidth, kFALSE, 5, vn);
  }

  Float_t awaySideYield = ExtractYields(hist, trigId, centralityId, ptACenter, ptAWidth, 0 + caseOffset, 0);
  
  if (flatBaseLine)
  {
    v2 = 0;
    vn = 0;
  }
    
/*    Float_t v2Trig  = flowGraph->Eval(ptT);
    Float_t v2Assoc = flowGraph->Eval(ptA);
    Float_t v2 = v2Trig * v2Assoc;
    Printf("%f %f: %.2f %.2f --> %.4f", ptT, ptA, v2Trig, v2Assoc, v2);*/
  
  TF1* flowFunc = new TF1("flowFunc", "[0] * (1+2*[1]*cos(2*x)+2*[2]*cos(3*x)+2*[3]*cos(4*x)+2*[4]*cos(5*x))", -0.5 * TMath::Pi(), 1.5 * TMath::Pi());

  // find minimum
  if (0)
  {
    Float_t minBin = hist->FindBin(0.6);
    for (Int_t bin=hist->FindBin(0.6); bin<=hist->FindBin(1.8); bin++)
      if (hist->GetBinContent(bin) < hist->GetBinContent(minBin))
        minBin = bin;
            
    Float_t norm = hist->GetBinContent(minBin) / (1 + 2 * v2 * cos(2 * hist->GetBinCenter(minBin)));
    Double_t normUnc = 0;
  }
  else
  {
    if (baseLineDetermination == 0)
    {
      return;
      /*
      if (centralityId == 0)
        hist->Fit("pol0", "0", "", 0.8., 1.2);
      else
        hist->Fit("pol0", "0", "", 1.2, 1.6);
      */
      //tsallis = new TF1("tsallis", "(1-[1]*(1-[0])*x*x)**(1/(1-[0]))", 0, 10)
      // (1-b(1-q)x*x)**(1/(1-q))
      
      // combination of two tsallis functions
      total = new TF1("total", "[0] + [1] * (1-[3]*(1-[2])*x*x)**(1/(1-[2])) + [4] * (1-[6]*(1-[5])*(x-TMath::Pi())*(x-TMath::Pi()))**(1/(1-[5]))", -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
      //total = new TF1("total", "[0] + [1] * (1-[3]*(1-[2])*x*x)**(1/(1-[2])) + [4] * (1-[6]*(1-[5])*(x-TMath::Pi())*(x-TMath::Pi()))**(1/(1-[5])) + [1] * (1-[3]*(1-[2])*(x-TMath::TwoPi())*(x-TMath::TwoPi()))**(1/(1-[2])) + [4] * (1-[6]*(1-[5])*(x+TMath::Pi())*(x+TMath::Pi()))**(1/(1-[5]))", -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
      // the latter is the correct function, but a very wide tsallis on the away side can fake a baseline, and then the baseline is significantly off
      total->SetParLimits(0, 0, 10000);
      total->SetParLimits(1, 0, 10000);
      total->SetParLimits(4, 0, 10000);
      
      total->SetParLimits(2, 1.0000001, 5);
      total->SetParLimits(5, 1.0000001, 5);
      
      total->SetParLimits(3, 0.1, 1000);
      total->SetParLimits(6, 0.1, 1000);
      
      total->SetParameters(1, 1, 1.5, 3, 1, 1.5, 1.1);
      
      // limit q>1
      hist->Fit(total, "0RI");
      fitResult = hist->Fit(total, "0RSI");
      
      if (0 && fitResult->CovMatrixStatus() > 0)
      {
        TMatrixDSym cov(fitResult->GetCovarianceMatrix());
        cov.Print();
        
        Double_t nearSideYield, nearSideError;
        TsallisYieldUncertainty(total->GetParameter(1), total->GetParameter(2), total->GetParameter(3), cov, 0, nearSideYield, nearSideError);
        
        Double_t awaySideYield2, awaySideError;
        TsallisYieldUncertainty(total->GetParameter(4), total->GetParameter(5), total->GetParameter(6), cov, 3, awaySideYield2, awaySideError);
        
        FillYield(trigId, centralityId, ptACenter, ptAWidth, 14 + caseOffset + 1, nearSideYield, nearSideError, awaySideYield2, awaySideError);
        
        Printf("Chi2 %f / ndf %d = %f", total->GetChisquare(), total->GetNDF(), (total->GetNDF() > 0) ? total->GetChisquare() / total->GetNDF() : -1);
        
        total->SetLineColor(hist->GetLineColor());
        total->SetLineWidth(1);
        if (caseOffset == 2)
        {
          if (0)
            total->DrawCopy("SAME");
          //total->SetParameter(5, total->GetParameter(5) - TMath::Sqrt(cov[5][5]));
          //total->DrawCopy("SAME");
          //total->SetParameter(5, total->GetParameter(5) + 2 * TMath::Sqrt(cov[5][5]));
          //Printf("%f", total->Integral(-0.7,0.7));
          //total->DrawCopy("SAME");
        }
        
        if (caseOffset == 0)
        {
          // store fit parameters
          
          // q
          FillYield(trigId, centralityId, ptACenter, ptAWidth, 28, total->GetParameter(2), total->GetParError(2), total->GetParameter(5), total->GetParError(5));
          // beta
          FillYield(trigId, centralityId, ptACenter, ptAWidth, 29, total->GetParameter(3), total->GetParError(3), total->GetParameter(6), total->GetParError(6));
          
          if (0)
          {
            // two Gauss fits
            gausFit = new TF1("gausFit", "[0] + gaus(1) + gaus(4)", -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
            gausFit->SetParameters(1, 1, 0, 1, 1, TMath::Pi(), 0);
            gausFit->SetParLimits(0, 0, 10000);
            gausFit->SetParLimits(1, 0, 10000);
            gausFit->FixParameter(2, 0);
            gausFit->SetParLimits(3, 0.01, 10);
            gausFit->SetParLimits(4, 0, 10000);
            gausFit->FixParameter(5, TMath::Pi());
            gausFit->SetParLimits(6, 0.01, 10);
            
            hist->Fit(gausFit, "0RI");
            gausFit->SetLineWidth(1);
            gausFit->DrawCopy("SAME");
          
            FillYield(trigId, centralityId, ptACenter, ptAWidth, 30, gausFit->GetParameter(3), total->GetParError(3), total->GetParameter(6), total->GetParError(6));
          }
        }
        
        if (0 && cov[0][0] > 0)
        {
          for (Int_t i=0; i<cov.GetNcols(); i++)
            for (Int_t j=0; j<cov.GetNrows(); j++)
              cov[j][i] /= TMath::Sqrt(cov[i][i]) * TMath::Sqrt(cov[j][j]);
          cov.Print();
        }
      }
      
      func = total;
      Float_t norm = func->GetParameter(0) / (1 - 2 * v2);
      Double_t normUnc = func->GetParError(0) / (1 - 2 * v2);
      //return;
    }
    else if (baseLineDetermination <= 2)
    {
      Float_t values[3];
      Float_t errors[3];
      
      if (baseLineDetermination == 1)
      {
        Float_t regionBegin[3] = { -TMath::Pi() / 2,        TMath::Pi() / 2 - 0.4, 1.5 * TMath::Pi() - 0.4 };
        Float_t regionEnd[3] =   { -TMath::Pi() / 2 + 0.4,  TMath::Pi() / 2 + 0.4, 1.5 * TMath::Pi() };
      }
      else if (baseLineDetermination == 2)
      {
        Float_t regionBegin[3] = { -TMath::Pi() / 2,              TMath::Pi() / 2 - 0.4 - 0.2, 1.5 * TMath::Pi() - 0.4 - 0.2};
        Float_t regionEnd[3] =   { -TMath::Pi() / 2 + 0.4 - 0.2,  TMath::Pi() / 2 + 0.4 - 0.2, 1.5 * TMath::Pi() };
      }
        
      // weighted mean
      Float_t sum = 0;
      Float_t weight = 0;
      for (Int_t i=0; i<3; i++)
      {
        hist->Fit("pol0", "0Q", "", regionBegin[i], regionEnd[i]);
        func = hist->GetFunction("pol0");
        if (!func)
          continue;
        sum += func->GetParameter(0) / func->GetParError(0) / func->GetParError(0);
        weight += 1. / func->GetParError(0) / func->GetParError(0);
      }
      
      if (weight == 0)
        return;
      
      sum /= weight;
      weight = TMath::Sqrt(1. / weight);
      
      Float_t norm = sum / (1 - 2 * v2);
      Double_t normUnc = weight / (1 - 2 * v2);
    }
    else
    {
      Int_t bins = 2;
      if (baseLineDetermination == 4)
        bins = 4;
      if (baseLineDetermination == 5)
        bins = 8;
      Float_t norm = NLowestBinAverage(*hist, bins);
      Double_t normUnc = 0;
      
      Printf("NLowestBinAverage %d --> %f", bins, norm);
      
      norm /= (1 - 2 * v2);
    }
    
    Printf("Baseline: %f +- %f", norm, normUnc);
  }
  
  if (caseOffset == 0)
  {
    Float_t awaySideYieldNoBaseline = awaySideYield - norm * (1 - 2 * v2);
    Float_t v2YieldNoBaseline = norm * 2 * v2 * (1 + 1.4) / 1.4;
    Printf("Relative v2 contribution at centrality = %d, pT,t = %.1f, pT,a = %.1f: %.1f%% (%f %f)", centralityId, ptT, ptA, 100.0 * v2YieldNoBaseline / awaySideYieldNoBaseline, awaySideYieldNoBaseline, v2YieldNoBaseline);
  }

  //flowFunc->SetParameters(hist->GetBinContent(hist->FindBin(1.4)) / (1.0 - 2.0 * v2), v2);
  flowFunc->SetParameters(norm, v2, 0, 0, 0);
  if (vn)
    flowFunc->SetParameters(norm, vn[1], vn[2], vn[3], vn[4]);
  flowFunc->SetLineWidth(2);
  if (caseOffset >= 4 && caseOffset <= 9 && caseOffset != 6 && caseOffset != 7)
  {
    if (caseOffset != 4 && caseOffset != 9)
      flowFunc->SetLineStyle(2);
    flowFunc->DrawCopy("SAME"); //->SetLineColor(caseOffset - 3);
    if (vn)
    {
      flowFuncTmp = (TF1*) flowFunc->Clone("flowFuncTmp");
      for (Int_t i=1; i<=4; i++)
      {
	flowFuncTmp->SetParameters(flowFuncTmp->GetParameter(0), 0, 0, 0, 0);
	flowFuncTmp->SetParameter(i, vn[i]);
	flowFuncTmp->SetLineStyle(2);
	flowFuncTmp->Print();
	flowFuncTmp->DrawCopy("SAME");
      }
    }
  }
  
  hist->Add(flowFunc, -1);
  if (caseOffset == 0)
  {
    file = TFile::Open("dphi_corr.root", "UPDATE");
    hist->SetName(TString(hist->GetName()) + "_tsallis_v2");
    hist->Write();
    file->Close();
    //  hist->DrawCopy("SAME");
  }
  if (caseOffset == 2 && flatBaseLine)
  {
    file = TFile::Open("dphi_corr.root", "UPDATE");
    hist->SetName(TString(hist->GetName()) + "_tsallis_flat");
    hist->Write();
    file->Close();
    //  hist->DrawCopy("SAME");
  }
  if (caseOffset == 4 && flatBaseLine)
  {
    file = TFile::Open("dphi_corr.root", "UPDATE");
    hist->SetName(TString(hist->GetName()) + "_fit_flat");
    hist->Write();
    file->Close();
    //  hist->DrawCopy("SAME");
  }
  if (caseOffset == 9)
  {
    file = TFile::Open("dphi_corr.root", "UPDATE");
    hist->SetName(TString(hist->GetName()) + "_fit_v2");
    hist->Write();
    file->Close();
    //  hist->DrawCopy("SAME");
  }
  
  ExtractYields(hist, trigId, centralityId, ptACenter, ptAWidth, 14 + caseOffset, normUnc);
  
  if (0 && histMixed)
  {
    Printf("%f", histMixed->Integral() / histMixed->GetNbinsX());
    flowFunc->SetParameters(histMixed->Integral() / histMixed->GetNbinsX(), v2);
    flowFunc->SetLineColor(2);
    //if (caseOffset == 0)
    //  flowFunc->DrawCopy("SAME");
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
    
  return;
  
  for (Int_t i=1; i<=hist->GetNbinsX(); i++)
    hist->SetBinContent(i, hist->GetBinContent(i) - zyam);
}
 
void* cacheIds[10];
TH2* cacheMixed[10];

Int_t gHistCount = 0;
void GetDistAndFlow(void* hVoid, void* hMixedVoid, TH1** hist, Float_t* v2, Int_t step, Int_t centralityBegin, Int_t centralityEnd, Float_t ptBegin, Float_t ptEnd, Int_t twoD = 0, Bool_t equivMixedBin = kFALSE, Float_t* vn = 0, Bool_t scaleToPairs = kTRUE)
{
  h = (AliUEHistograms*) hVoid;
  hMixed = (AliUEHistograms*) hMixedVoid;

  Int_t centralityBeginBin = 0;
  Int_t centralityEndBin = -1;
  
  if (centralityEnd >= centralityBegin)
  {
    centralityBeginBin = h->GetUEHist(2)->GetEventHist()->GetGrid(step)->GetGrid()->GetAxis(1)->FindBin(0.01 + centralityBegin);
    centralityEndBin = h->GetUEHist(2)->GetEventHist()->GetGrid(step)->GetGrid()->GetAxis(1)->FindBin(-0.01 + centralityEnd);
  }
  
  // 2d same and mixed event
  TH2* sameTwoD  = h->GetUEHist(2)->GetUEHist(step, 0, ptBegin, ptEnd, centralityBeginBin, centralityEndBin, 1, kFALSE);
  
  if (hMixed)
  {
    if (!equivMixedBin)
    {
      // No centrality, nor pT associated dep of the mixed event observed. Use a larger range to get more statistics
      
      // TODO HACK. mixed event is not propagated to step0
      Int_t stepMixed = 6;
      
      Int_t cacheId = -1;
      
      for (Int_t i=0; i<10; i++)
	if (cacheIds[i] == hMixed)
	{
	  cacheId = i;
	  break;
	}
	
      // not found
      if (cacheId == -1)
	for (Int_t i=0; i<10; i++)
	  if (cacheIds[i] == 0)
	  {
	    cacheId = i;
	    break;
	  }

      if (!cacheIds[cacheId])
      {
	hMixed->SetPtRange(1.0, 10);
	cacheMixed[cacheId] = (TH2*) hMixed->GetUEHist(2)->GetUEHist(stepMixed, 0, 1.0, 10.0, 1, 15, 1, kFALSE);
	cacheIds[cacheId] = hMixed;
	Printf("GetDistAndFlow: Cached for %p on slot %d", hMixed, cacheId);
      }
      
      TH2* mixedTwoD = cacheMixed[cacheId];
    }
    else
    {
      // use same bin for mixing
      
      // TODO HACK. mixed event is not propagated to step0
      Int_t stepMixed = 6;
      TH2* mixedTwoD = (TH2*) hMixed->GetUEHist(2)->GetUEHist(stepMixed, 0, ptBegin, ptEnd, centralityBeginBin, centralityEndBin, 1, kFALSE);
    }
    
    if (0)
    {
      // asssume flat in dphi, gain in statistics
      Printf("NOTE: Assuming flat acceptance in phi!");
      
      TH1* histMixedproj = mixedTwoD->ProjectionY();
      histMixedproj->Scale(1.0 / mixedTwoD->GetNbinsX());
      
      for (Int_t x=1; x<=mixedTwoD->GetNbinsX(); x++)
	for (Int_t y=1; y<=mixedTwoD->GetNbinsY(); y++)
	  mixedTwoD->SetBinContent(x, y, histMixedproj->GetBinContent(y));
    }
    
    // get mixed event normalization by assuming full acceptance at deta of 0 (only works for flat dphi)
    if (scaleToPairs)
    {
      Double_t mixedNorm = mixedTwoD->Integral(1, mixedTwoD->GetNbinsX(), mixedTwoD->GetYaxis()->FindBin(-0.01), mixedTwoD->GetYaxis()->FindBin(0.01));
      mixedNorm /= mixedTwoD->GetNbinsX() * (mixedTwoD->GetYaxis()->FindBin(0.01) - mixedTwoD->GetYaxis()->FindBin(-0.01) + 1);
    }
    else
      Double_t mixedNorm = mixedTwoD->Integral() / sameTwoD->Integral();
    
    // divide and scale
    sameTwoD->Divide(mixedTwoD);
    sameTwoD->Scale(mixedNorm);
    
/*    new TCanvas;
    sameTwoD->Draw("SURF1");
    dfdsafd;*/
  }
  
  TString histName;
  histName.Form("GetDistAndFlow_%d", gHistCount++);
  
  // extract dphi distribution if requested
  if (twoD == 1)
  {
    *hist = sameTwoD;
  }

  //  Float_t etaLimit = 0.8;
  Float_t etaLimit = 1.0;

  // 20: return corr in |delta eta| < 1 from which 1 < |delta eta| < 2 is subtracted
  if (twoD == 0 || twoD == 10 || twoD == 20)
  {
    Int_t etaBegin = 1;
    Int_t etaEnd = sameTwoD->GetNbinsY();
    
    if (twoD == 10 || twoD == 20)
    {
      etaBegin = sameTwoD->GetYaxis()->FindBin(-etaLimit + 0.01);
      etaEnd   = sameTwoD->GetYaxis()->FindBin(etaLimit - 0.01);
    }

    *hist = sameTwoD->ProjectionX(histName, etaBegin, etaEnd);
    
    if (!scaleToPairs)
      (*hist)->Scale(1.0 / (etaEnd - etaBegin + 1));
  }
  
  if (twoD == 11 || twoD == 20)
  {
    // errors --> are ok
    
//       Float_t outerLimit = 2.0;
    Float_t outerLimit = 1.8;
    Printf("Phi dist: Using outer limit %.2f", outerLimit);
//       Float_t outerLimit = etaLimit * 2;
    
    histTmp = sameTwoD->ProjectionX(histName + "1", TMath::Max(1, sameTwoD->GetYaxis()->FindBin(-outerLimit + 0.01)), sameTwoD->GetYaxis()->FindBin(-etaLimit - 0.01));
    Int_t etaBins = sameTwoD->GetYaxis()->FindBin(-etaLimit - 0.01) - TMath::Max(1, sameTwoD->GetYaxis()->FindBin(-outerLimit + 0.01)) + 1;

    TH1D* tracksTmp = sameTwoD->ProjectionX(histName + "2", sameTwoD->GetYaxis()->FindBin(etaLimit + 0.01), TMath::Min(sameTwoD->GetYaxis()->GetNbins(), sameTwoD->GetYaxis()->FindBin(outerLimit - 0.01)));
    etaBins += TMath::Min(sameTwoD->GetYaxis()->GetNbins(), sameTwoD->GetYaxis()->FindBin(outerLimit - 0.01)) - sameTwoD->GetYaxis()->FindBin(etaLimit + 0.01) + 1;

//       printf("%f +- %f  %f +- %f ", (*hist)->GetBinContent(1), (*hist)->GetBinError(1), tracksTmp->GetBinContent(1), tracksTmp->GetBinError(1));
    histTmp->Add(tracksTmp);
//       Printf(" --> %f +- %f", (*hist)->GetBinContent(1), (*hist)->GetBinError(1));
    
    if (!scaleToPairs)
      histTmp->Scale(1.0 / etaBins);
    
    if (twoD == 11)
      *hist = histTmp;
    else if (twoD == 20)
    {
      // calculate acc with 2 * (deta - 0.5 * deta*deta / 1.6)
      if (!hMixedVoid)
	histTmp->Scale(0.75 / 0.25);
      
      histTmp->Scale(1.0 / 0.8);
      
      (*hist)->Add(histTmp, -1);
    }
  }
  
//   (*hist)->Rebin(2); (*hist)->Scale(0.5);
  
  //*hist = h->GetUEHist(2)->GetUEHist(step, 0, ptBegin, ptEnd, h->GetUEHist(2)->GetEventHist()->GetGrid(step)->GetGrid()->GetAxis(1)->FindBin(0.01 + centralityBegin), h->GetUEHist(2)->GetEventHist()->GetGrid(step)->GetGrid()->GetAxis(1)->FindBin(-0.01 + centralityEnd), twoD);
  
  TString str;
  str.Form("%.1f < p_{T,trig} < %.1f", ptBegin - 0.01, ptEnd + 0.01);
  
  TString str2;
  str2.Form("%.2f < p_{T,assoc} < %.2f", gpTMin - 0.01, gpTMax + 0.01);
    
  TString newTitle;
  newTitle.Form("%s - %s - %d-%d%%", str.Data(), str2.Data(), centralityBegin, centralityEnd);
  (*hist)->SetTitle(newTitle);
  
  if (0 && hMixed)
  {
    histMixed = hMixed->GetUEHist(2)->GetUEHist(step, 0, ptBegin, ptEnd, hMixed->GetUEHist(2)->GetEventHist()->GetGrid(step)->GetGrid()->GetAxis(1)->FindBin(0.01 + centralityBegin), hMixed->GetUEHist(2)->GetEventHist()->GetGrid(step)->GetGrid()->GetAxis(1)->FindBin(-0.01 + centralityEnd));
    
    //new TCanvas; (*hist)->DrawCopy(); histMixed->DrawCopy("SAME")->SetLineColor(2);
    
    Float_t totalPairs = (*hist)->Integral();
    
    (*hist)->Divide(histMixed);
    (*hist)->Scale(totalPairs / (*hist)->Integral());
    
    //(*hist)->DrawCopy("SAME")->SetLineColor(4);
  }
  
  if (v2 || vn)
  {
    // calculate v2trigger
    h->GetUEHist(2)->GetEventHist()->GetGrid(step)->GetGrid()->GetAxis(1)->SetRangeUser(0.01 + centralityBegin, -0.01 + centralityEnd);
    ptDist = h->GetUEHist(2)->GetEventHist()->Project(step, 0);
    Float_t vTrig[5];
    for (Int_t i=2; i<=((vn) ? 5 : 2); i++)
      vTrig[i-1] = CalculateFlow(ptDist, ptBegin, ptEnd, i, centralityBegin, centralityEnd);
    delete ptDist;
    
    // calculate v2 assoc
    cont = h->GetUEHist(2)->GetTrackHist(0);
    h->GetUEHist(2)->GetTrackHist(0)->GetGrid(step)->GetGrid()->GetAxis(3)->SetRangeUser(0.01 + centralityBegin, -0.01 + centralityEnd);
    h->GetUEHist(2)->GetTrackHist(0)->GetGrid(step)->GetGrid()->GetAxis(2)->SetRangeUser(ptBegin, ptEnd);
    ptDist = h->GetUEHist(2)->GetTrackHist(0)->Project(step, 1);
    Float_t vAssoc[5];
    for (Int_t i=2; i<=((vn) ? 5 : 2); i++)
      vAssoc[i-1] = CalculateFlow(ptDist, gpTMin, gpTMax, i, centralityBegin, centralityEnd);
    delete ptDist;
  
    if (v2)
      *v2 = vTrig[2-1] * vAssoc[2-1];
    if (vn)
      for (Int_t i=2; i<=5; i++)
	vn[i-1] = vTrig[i-1] * vAssoc[i-1];
  }
}
 
void GetSumOfRatios(void* hVoid, void* hMixedVoid, TH1** hist, Int_t step, Int_t centralityBegin, Int_t centralityEnd, Float_t ptBegin, Float_t ptEnd, Bool_t useVertexBins)
{
  h = (AliUEHistograms*) hVoid;
  hMixed = (AliUEHistograms*) hMixedVoid;

  Int_t centralityBeginBin = 0;
  Int_t centralityEndBin = -1;
  
  if (centralityEnd >= centralityBegin)
  {
    centralityBeginBin = h->GetUEHist(2)->GetEventHist()->GetGrid(step)->GetGrid()->GetAxis(1)->FindBin(0.01 + centralityBegin);
    centralityEndBin = h->GetUEHist(2)->GetEventHist()->GetGrid(step)->GetGrid()->GetAxis(1)->FindBin(-0.01 + centralityEnd);
  }
  
  // 2d same and mixed event
  *hist  = h->GetUEHist(2)->GetSumOfRatios(hMixed->GetUEHist(2), step, 0, ptBegin, ptEnd, centralityBeginBin, centralityEndBin, kFALSE, useVertexBins);
  
  TString str;
  str.Form("%.1f < p_{T,trig} < %.1f", ptBegin - 0.01, ptEnd + 0.01);
  
  TString str2;
  str2.Form("%.2f < p_{T,assoc} < %.2f", gpTMin - 0.01, gpTMax + 0.01);
    
  TString newTitle;
  newTitle.Form("%s - %s - %d-%d%%", str.Data(), str2.Data(), centralityBegin, centralityEnd);
  (*hist)->SetTitle(newTitle);
}
 
void PlotDeltaPhiDistributions(const char* fileName1, const char* fileName2, Float_t yMax = 0.1, Int_t twoD = 0, Int_t centrBegin = 1, Int_t centrEnd = 1)
{
  loadlibs();
  
  Bool_t veryCentral = 0;
  Bool_t flowComparison = 0;
  Bool_t rhicOverlay = 0;
  Bool_t highStatBinning = 0;

  file = TFile::Open("dphi_corr.root", "RECREATE");
  file->Close();
  
   Int_t leadingPtOffset = 1;
    
  if (veryCentral || flowComparison)
  {
    Int_t maxLeadingPt = 2;
    Int_t maxAssocPt = 2;
    Float_t leadingPtArr[] = { 2.0, 3.0, 4.0, 10.0, 20.0, 40.0 };
    //Float_t assocPtArr[] =   { 0.15, 0.5, 1.0, 2.0, 4.0, 6.0, 10.0, 20.0, 40.0 };
    Float_t assocPtArr[] =   { 1.0, 2.0, 3.0, 6.0, 10.0, 20.0, 40.0 };
  }
  else if (rhicOverlay) // RHIC binning
  {
    Int_t maxLeadingPt = 4;
    Int_t maxAssocPt = 5;
    Float_t leadingPtArr[] =   { 4.0, 5.0, 7.0, 9.0, 12.0 };
    Float_t assocPtArr[] =     { 0.5, 1.0, 2.0, 3.0, 5.0, 7.0 };
  }
  else if (highStatBinning) 
  {
    Int_t maxLeadingPt = 3;
    Int_t maxAssocPt = 2;
    Float_t leadingPtArr[] = { 4.0, 6.0, 8.0, 15.0 };
    Float_t assocPtArr[] =   { 1.0, 4.0, 10.0 };
  }
  else // ALICE binning
  {
    if (1) // binning from preliminaries
    {
      Int_t maxLeadingPt = 2;
      Int_t maxAssocPt = 7;
//       Float_t leadingPtArr[] = { 6.0, 8.0, 10.0, 10.0, 15.0 };
      Float_t leadingPtArr[] = { 6.0, 8.0, 10.0, 15.0, 15.0 };
      Float_t assocPtArr[] =     { 0.5, 1.5, 3.0, 4.0, 6.0, 8.0, 10.0, 12.0 };
      leadingPtOffset = 2;
    }
    else if (0)
    {
      Int_t maxLeadingPt = 1;
      Int_t maxAssocPt = 4;
      Float_t leadingPtArr[] = { 8.0, 10.0, 15.0, 15.0 };
      Float_t assocPtArr[] =     { 3.0, 4.0, 6.0, 8.0, 10.0, 12.0 };
      leadingPtOffset = 2;
    }
    else     
    {
      Int_t maxLeadingPt = 3;
      Int_t maxAssocPt = 3;
      Float_t leadingPtArr[] = { 6.0, 8.0, 10.0, 15.0, 20.0 };
      Float_t assocPtArr[] =     { 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 12.0 };
      leadingPtOffset = 2;
    }
  }
  
  Int_t nCentralityBins = 5;
  Int_t centralityBins[] = { 1, 7, 9, 11, 13, 16 };
  //Int_t centralityBins[] = { 1, 3, 5, 7, 9, 13 };
  
  AliUEHistograms* h = (AliUEHistograms*) GetUEHistogram(fileName1);
  
//   h->SetZVtxRange(-0.5, 0.5);
//   h->SetZVtxRange(1.5, 2.5);
  
  AliUEHistograms* hMixed = 0;
  AliUEHistograms* hMixed2 = 0; // GetUEHistogram(fileName2, 0, kTRUE);

  if (twoD)
  {
    hMixed = (AliUEHistograms*) GetUEHistogram(fileName1, 0, kTRUE);
    hMixed2 = (AliUEHistograms*) GetUEHistogram(fileName2, 0, kTRUE);
  }


  if (veryCentral)
  {
    Printf("WARNING: Reading mixed event from preliminaries/corrected_110317.root");
    hMixed = (AliUEHistograms*) GetUEHistogram("preliminaries/corrected_110317.root", 0, kTRUE);
  }
  
  AliUEHistograms* h2 = 0;
  if (!twoD)
    h2 = (AliUEHistograms*) GetUEHistogram(fileName2);

  
  TCanvas* canvas = new TCanvas("DeltaPhi", "DeltaPhi", 1000, 700);
  canvas->Divide(maxAssocPt, maxLeadingPt);
  
  TCanvas* canvas2 = new TCanvas("Centrality", "Centrality", 800, 600);
  centralityHist = (TH1*) h->GetCentralityDistribution();
  NormalizeToBinWidth(centralityHist);
  centralityHist->Draw();
  gPad->SetLogy();
  
  TLegend* legend = new TLegend(0.2, 0.5, 0.95, 0.90);
  TLegend* legend2 = new TLegend(0.5, 0.63, 0.95, 0.90);
  legend2->SetFillColor(0);
  legend2->SetTextSize(0.04);
  
  for (Int_t i=0; i<maxLeadingPt; i++)
    for (Int_t j=0; j<maxAssocPt; j++)
    {
      TString str;
      str.Form("%.1f < p_{T,trig} < %.1f", leadingPtArr[i], leadingPtArr[i+leadingPtOffset]);
      
      if (j == 0)
      {
        canvas2->cd();
        h->GetUEHist(2)->GetEventHist()->GetGrid(6)->SetRangeUser(0, leadingPtArr[i] + 0.01, leadingPtArr[i+leadingPtOffset] - 0.01);
        centralityHist = h->GetUEHist(2)->GetEventHist()->ShowProjection(1, 6);
        centralityHist->SetLineColor(i+2);
        NormalizeToBinWidth(centralityHist);
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
      
      gpTMin = assocPtArr[j] + 0.01;
      gpTMax = assocPtArr[j+1] - 0.01;
      
      TString str2;
      str2.Form("%.1f < p_{T,assoc} < %.1f", gpTMin - 0.01, gpTMax + 0.01);
      
      SetupRanges(h);

      if (h2)
      {
        SetupRanges(h2); // SetEtaRange(0, 0) does not need to be called for the leading track result
      }
      
      if (hMixed)
      {
        SetupRanges(hMixed);
      }
      
      // delta phi
      if (!twoD)
      {
        if (assocPtArr[j] >= leadingPtArr[i+leadingPtOffset])
          continue;
    
        // 0-5% --> 1, 5
        // 0-10% --> 1, 6
        // 0-20% --> 1, 8
        // 20-40% --> 9, 10
        // 40-80% --> 11, 14
        // > 40% --> 11, 16
        
        TString hist1Str, hist2Str, hist2bStr;
        
        Float_t v2[3];
        for (Int_t k=0; k<3; k++)
          v2[k] = 0;
	Float_t vn[3][3];
        
        if (veryCentral)
        {
          Int_t step = 0;
          TH1* hist1 = 0;
          TH1* hist2 = 0;
          TH1* hist2b = 0;
	  TH1* hist3 = 0;
	  
	  GetDistAndFlow(h, hMixed, &hist1,  v2, step, 0,  2,  leadingPtArr[i] + 0.01, leadingPtArr[i+leadingPtOffset] - 0.01); hist1Str = "0-2%";
// 	  GetDistAndFlow(h, hMixed, &hist2,  v2, step, 1,  3,  leadingPtArr[i] + 0.01, leadingPtArr[i+leadingPtOffset] - 0.01); hist2Str = "1-3%";
	  GetDistAndFlow(h, hMixed, &hist2b,  v2+2, step, 30,  40,  leadingPtArr[i] + 0.01, leadingPtArr[i+leadingPtOffset] - 0.01); hist2bStr = "30-40%";
          
	  //TH1* hist1 = h->GetUEHist(2)->GetUEHist(step, 0, leadingPtArr[i] + 0.01, leadingPtArr[i+leadingPtOffset] - 0.01, 1, 2);  hist1Str = "0-2%";
          //TH1* hist2 = h->GetUEHist(2)->GetUEHist(step, 0, leadingPtArr[i] + 0.01, leadingPtArr[i+leadingPtOffset] - 0.01, 2, 3);  hist2Str = "2-3%";
          //TH1* hist2b = h->GetUEHist(2)->GetUEHist(step, 0, leadingPtArr[i] + 0.01, leadingPtArr[i+leadingPtOffset] - 0.01, 10, 10); hist2bStr = "30-40%";
        }
        else if (flowComparison)
        {
          TH1* hist1 = h->GetUEHist(2)->GetUEHist(6, 0, leadingPtArr[i] + 0.01, leadingPtArr[i+leadingPtOffset] - 0.01, 1, 5);  hist1Str = "0-5%";
          TH1* hist2 = h->GetUEHist(2)->GetUEHist(6, 0, leadingPtArr[i] + 0.01, leadingPtArr[i+leadingPtOffset] - 0.01, 9, 10);  hist2Str = "20-40%";
          TH1* hist2b = h->GetUEHist(2)->GetUEHist(6, 0, leadingPtArr[i] + 0.01, leadingPtArr[i+leadingPtOffset] - 0.01, 10, 10); hist2bStr = "30-40%";
          TH1* hist3 = 0; // h2->GetUEHist(2)->GetUEHist(0, 0, leadingPtArr[i] + 0.01, leadingPtArr[i+leadingPtOffset] - 0.01);
        }
        else if (rhicOverlay)
        {
          TH1* hist1 = h->GetUEHist(2)->GetUEHist(6, 0, leadingPtArr[i] + 0.01, leadingPtArr[i+leadingPtOffset] - 0.01, 1, 8);   hist1Str = "0-20%";
          TH1* hist2 = h->GetUEHist(2)->GetUEHist(6, 0, leadingPtArr[i] + 0.01, leadingPtArr[i+leadingPtOffset] - 0.01, 9, 12);  hist2Str = "20-60%";
          TH1* hist2b = h->GetUEHist(2)->GetUEHist(6, 0, leadingPtArr[i] + 0.01, leadingPtArr[i+leadingPtOffset] - 0.01, 13, 15);hist2bStr = "60-90%";
          TH1* hist3 = h2->GetUEHist(2)->GetUEHist(0, 0, leadingPtArr[i] + 0.01, leadingPtArr[i+leadingPtOffset] - 0.01);
        }
        else
        {
          Int_t step = 6;
          TH1* hist1 = 0;
          TH1* hist2 = 0;
          TH1* hist2b = 0;
	  
	  Bool_t equivMixedBin = 1;
	  Int_t histType = 0;
// 	  histType = 20; Printf("WARNING: Using histogram type 20");
	  
//           GetDistAndFlow(h, hMixed, &hist1,  v2, step, 0,  2,  leadingPtArr[i] + 0.01, leadingPtArr[i+leadingPtOffset] - 0.01, histType, equivMixedBin, vn[0]); hist1Str = "0-2%";
//           GetDistAndFlow(h, hMixed, &hist1,  v2, step, 0,  5,  leadingPtArr[i] + 0.01, leadingPtArr[i+leadingPtOffset] - 0.01, histType, equivMixedBin, vn[0]); hist1Str = "0-5%";
          GetDistAndFlow(h, hMixed, &hist1,  v2, step, 0,  5,  leadingPtArr[i] + 0.01, leadingPtArr[i+leadingPtOffset] - 0.01, histType, equivMixedBin); hist1Str = "0-5%";
	  
	  /*
	  new TCanvas;
	  hist1->Draw();

	  TH1* histTmp1 = 0;
          GetDistAndFlow(h, hMixed, &histTmp1,  v2, step, 0,  5,  leadingPtArr[i] + 0.01, leadingPtArr[i+leadingPtOffset] - 0.01, 10, equivMixedBin); hist1Str = "0-5%";
	  histTmp1->SetLineColor(2);
	  histTmp1->Draw("SAME");
	  histTmp1->Scale(1 / 0.75);
	  
	  histTmp1 = 0;
          GetDistAndFlow(h, hMixed, &histTmp1,  v2, step, 0,  5,  leadingPtArr[i] + 0.01, leadingPtArr[i+leadingPtOffset] - 0.01, 20, equivMixedBin); hist1Str = "0-5%";
	  histTmp1->SetLineColor(4);
	  histTmp1->Draw("SAME");
	  histTmp1->Scale(1 / 0.75);

	  histTmp1 = 0;
          GetDistAndFlow(h, hMixed, &histTmp1,  v2, step, 0,  5,  leadingPtArr[i] + 0.01, leadingPtArr[i+leadingPtOffset] - 0.01, 11, equivMixedBin); hist1Str = "0-5%";
	  histTmp1->SetLineColor(3);
	  histTmp1->Draw("SAME");
	  histTmp1->Scale(1 / 0.25);

	  return;
	  */

	  GetDistAndFlow(h, hMixed, &hist2,  v2+1, step, 0, 20, leadingPtArr[i] + 0.01, leadingPtArr[i+leadingPtOffset] - 0.01, histType, equivMixedBin); hist2Str = "0-20%";
//           GetDistAndFlow(h, hMixed, &hist2b, v2[2], step, 60, 80, leadingPtArr[i] + 0.01, leadingPtArr[i+leadingPtOffset] - 0.01, 0, equivMixedBin); hist2bStr = "60-80%";
          GetDistAndFlow(h, hMixed, &hist2b, v2+2, step, 60, 90, leadingPtArr[i] + 0.01, leadingPtArr[i+leadingPtOffset] - 0.01, histType, equivMixedBin); hist2bStr = "60-90%";
          
          Printf("%f %f %f", v2[0], v2[1], v2[2]);
//           Printf("%f %f %f", vn[0][1], vn[0][2], vn[0][3]);
          
//           TH1* hist1 = h->GetUEHist(2)->GetUEHist(step, 0, leadingPtArr[i] + 0.01, leadingPtArr[i+leadingPtOffset] - 0.01, 1, 8);    hist1Str = "0-20%";
//           TH1* hist2 = h->GetUEHist(2)->GetUEHist(step, 0, leadingPtArr[i] + 0.01, leadingPtArr[i+leadingPtOffset] - 0.01, 9, 12);   hist2Str = "20-60%";
//           TH1* hist2b = h->GetUEHist(2)->GetUEHist(step, 0, leadingPtArr[i] + 0.01, leadingPtArr[i+leadingPtOffset] - 0.01, 13, 15); hist2bStr = "60-90%";
          
          step = 6;
//           TH1* hist3Old = h2->GetUEHist(2)->GetUEHist(step, 0, leadingPtArr[i] + 0.01, leadingPtArr[i+leadingPtOffset] - 0.01);
	  TH1* hist3 = 0;
          GetDistAndFlow(h2, hMixed2, &hist3,  0, step, 0,  -1,  leadingPtArr[i] + 0.01, leadingPtArr[i+leadingPtOffset] - 0.01, histType, equivMixedBin);
// 	  hist3->Rebin(2); hist3->Scale(0.5);
	  
/*	  new TCanvas;
	  hist3->Draw();
	  hist3Old->DrawCopy("SAME")->SetLineColor(2); */
	}
        
        /*
        RemoveBaseLine(hist1);
        RemoveBaseLine(hist2);
        RemoveBaseLine(hist2b);
        RemoveBaseLine(hist3);
        */
        
        TString newTitle;
        newTitle.Form("%s - %s", str.Data(), str2.Data());
        if (hist1)
        {
          hist1->SetName(Form("dphi_%d_%d_%d", i, j, 0));
          hist1->SetTitle(newTitle + " - " + hist1Str);
        }
        if (hist2)
        {
          hist2->SetName(Form("dphi_%d_%d_%d", i, j, 1));
          hist2->SetTitle(newTitle + " - " + hist2Str);
        }
        if (hist2b)
        {
          hist2b->SetName(Form("dphi_%d_%d_%d", i, j, 2));
          hist2b->SetTitle(newTitle + " - " + hist2bStr);
        }
        if (hist3)
        {
          hist3->SetName(Form("dphi_%d_%d_%d", i, j, 3));
          hist3->SetTitle(newTitle + " - pp");
        }
        
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
            legend->AddEntry(hist2b, "Pb+Pb 60-90%");
          if (hist3)
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
          yMin = -0.01; //TMath::Min(hist1->GetMinimum(), hist2->GetMinimum()) * 0.97;
          yMax2 = TMath::Max((hist3) ? hist3->GetMaximum() : 0.0, TMath::Max(hist1->GetMaximum(), (hist2) ? hist2->GetMaximum() : 0.0)) * 1.03;
        }
        
        yMax2 *= 1.4;
      
        TH2F* dummy = new TH2F("dummy", "", 100, hist1->GetXaxis()->GetBinLowEdge(1), hist1->GetXaxis()->GetBinUpEdge(hist1->GetNbinsX()), 1000, yMin, yMax2); //TMath::Max(hist1->GetMaximum(), hist2->GetMaximum()) * 1.1);
        dummy->SetStats(kFALSE);
        dummy->SetXTitle(hist1->GetXaxis()->GetTitle());
        dummy->SetYTitle(hist1->GetYaxis()->GetTitle());
        dummy->SetYTitle("1/N_{trig} dN/d#Delta#phi"); 
        Prepare1DPlot(dummy);
        
        dummy->GetYaxis()->SetTitleOffset(0.8);
      
        dummy->GetXaxis()->SetLabelSize(0.08);
        dummy->GetYaxis()->SetLabelSize(0.08);
        dummy->GetXaxis()->SetTitleSize(0.08);
        dummy->GetYaxis()->SetTitleSize(0.08);
        /*
        dummy->GetYaxis()->SetTitleOffset(0.8);
        */
        
        dummyTmp = dummy->DrawCopy();
        
        hist1->DrawCopy("SAME");
        
        if (hMixed)
        {
          SetupRanges(hMixed);
          // for HI file do not set range in eta anymore after it was changed to delta eta axis
          hMixed->SetEtaRange(0, 0);
        }
        TH1* hist1Mixed = 0; //hMixed->GetUEHist(2)->GetUEHist(6, 0, leadingPtArr[i] + 0.01, leadingPtArr[i+leadingPtOffset] - 0.01, 1, 5);
        //DrawFlow(v2[0], hist1, leadingPtArr[i], assocPtArr[j], hist1Mixed, i, 0, 0, (assocPtArr[j] + assocPtArr[j+1]) / 2, (assocPtArr[j+1] - assocPtArr[j]) / 2, kFALSE, 0, vn[0]);
        DrawFlow(v2[0], hist1, leadingPtArr[i], assocPtArr[j], hist1Mixed, i, 0, 0, (assocPtArr[j] + assocPtArr[j+1]) / 2, (assocPtArr[j+1] - assocPtArr[j]) / 2, kFALSE, 0);
        
        //hist1Mixed->Draw("SAME");
        if (hist2)
        {
          hist2->SetLineColor(2);
          hist2->DrawCopy("SAME");
          
          TH1* hist2Mixed = 0; //hMixed->GetUEHist(2)->GetUEHist(6, 0, leadingPtArr[i] + 0.01, leadingPtArr[i+leadingPtOffset] - 0.01, 9, 10);
          //hist2Mixed->SetLineColor(2);
          //hist2Mixed->Draw("SAME");

          DrawFlow(v2[1], hist2, leadingPtArr[i], assocPtArr[j], hist2Mixed, i, 1, 0, (assocPtArr[j] + assocPtArr[j+1]) / 2, (assocPtArr[j+1] - assocPtArr[j]) / 2);
          //DrawFlow(GetFlow05(), hist2, leadingPtArr[i], assocPtArr[j], hist2Mixed);
        }
        if (hist2b)
        {
          hist2b->SetLineColor(3);
          hist2b->DrawCopy("SAME");
        
          TH1* hist2bMixed = 0; //hMixed->GetUEHist(2)->GetUEHist(6, 0, leadingPtArr[i] + 0.01, leadingPtArr[i+leadingPtOffset] - 0.01, 13, 15);
          DrawFlow(v2[2], hist2b, leadingPtArr[i], assocPtArr[j], hist2bMixed, i, 2, 0, (assocPtArr[j] + assocPtArr[j+1]) / 2, (assocPtArr[j+1] - assocPtArr[j]) / 2);
          //DrawFlow(GetFlow05(), hist2b, leadingPtArr[i], assocPtArr[j], hist2bMixed);
        }
        if (hist3)
        {
          hist3->SetLineColor(4);
          hist3->DrawCopy("SAME");
          
          DrawFlow(0, hist3, leadingPtArr[i], assocPtArr[j], 0, i, 3, 0, (assocPtArr[j] + assocPtArr[j+1]) / 2, (assocPtArr[j+1] - assocPtArr[j]) / 2);
        }
        //dummyTmp->GetYaxis()->SetRangeUser(0, 1.1 * TMath::Max(TMath::Max(hist1->GetMaximum(), hist2->GetMaximum()), hist2b->GetMaximum()));
      }
      else // delta eta delta phi
      {
        if (twoD == 1)
        {
          if (assocPtArr[j] > leadingPtArr[i])
            continue;
        }
        else
        {
          Int_t jRef = 1;
        
          if (assocPtArr[jRef] > leadingPtArr[i])
            continue;
          
          // fix pt assoc
          gpTMin = assocPtArr[jRef] + 0.01;
          gpTMax = assocPtArr[jRef+1] - 0.01;
          
          str2.Form("%.1f < p_{T,assoc} < %.1f", gpTMin - 0.01, gpTMax + 0.01);
      
          // use j for centrality
          if (j >= nCentralityBins)
            continue;
            
          centrBegin = centralityBins[j];
          centrEnd = centralityBins[j+1] - 1;
        }
        
        SetupRanges(h);
        // for HI file do not set range in eta anymore after it was changed to delta eta axis
        h->SetEtaRange(0, 0);
        
        SetupRanges(hMixed);
        // for HI file do not set range in eta anymore after it was changed to delta eta axis
        hMixed->SetEtaRange(0, 0);
          
        TH2* histSame = (TH2*) h->GetUEHist(2)->GetUEHist(6, 0, leadingPtArr[i] + 0.01, leadingPtArr[i+leadingPtOffset] - 0.01, centrBegin, centrEnd, kTRUE);
        TH2* histMixed = (TH2*) hMixed->GetUEHist(2)->GetUEHist(6, 0, leadingPtArr[i] + 0.01, leadingPtArr[i+leadingPtOffset] - 0.01, centrBegin, centrEnd, kTRUE);
        
        // rebin 
        histSame->Rebin2D(2, 2);
        histMixed->Rebin2D(2, 2);
        
        if (1)
        {
          // fit delta eta, assuming no dependence on dphi
          
          //new TCanvas; histMixed->DrawCopy("SURF1");
          
          histMixedproj = histMixed->ProjectionY();
          histMixedproj->Scale(1.0 / histMixed->GetNbinsX());
          
          for (Int_t x=1; x<=histMixed->GetNbinsX(); x++)
            for (Int_t y=1; y<=histMixed->GetNbinsY(); y++)
              histMixed->SetBinContent(x, y, histMixedproj->GetBinContent(y));
          
          //new TCanvas; histMixed->DrawCopy("SURF1");
        }
        
        histSame->SetStats(0);
        histSame->GetYaxis()->SetRangeUser(-1.5, 1.5);
        histSame->SetTitle("");
        histSame->Divide(histMixed);
        histSame->DrawCopy("SURF1");
        
        TString str3;
        str3.Form("%d-%d%%", (Int_t) h->GetCentralityDistribution()->GetXaxis()->GetBinLowEdge(centrBegin), (Int_t) h->GetCentralityDistribution()->GetXaxis()->GetBinUpEdge(centrEnd));
        latex = new TLatex(0.15, 0.95, str3);
        latex->SetNDC();
        latex->SetTextSize(0.08);
        latex->Draw();
      }
      
      latex = new TLatex(0.55, 0.8, str);
      latex->SetNDC();
      latex->SetTextSize(0.06);
      latex->Draw();
      
      latex = new TLatex(0.55, 0.88, str2);
      latex->SetNDC();
      latex->SetTextSize(0.06);
      latex->Draw();
      
//            if (i == 0)        return;
    }

  canvas->SaveAs(Form("DeltaPhi_%.2f.png", yMax));
  
  canvas2->cd();
  legend2->Draw();

  //TString name;
  //name.Form("%s_%.2f_%.2f_%.2f_%.2f.png", TString(gSystem->BaseName(fileName1)).Tokenize(".")->First()->GetName(), leadingPtArr[i], leadingPtArr[i+1], assocPtArr[j], assocPtArr[j+1]);
}
   
void ExamplePhiEtaGap(const char* fileNamePbPb, const char* fileNamePbPbMix)
{
  loadlibs();
  
  if (!fileNamePbPbMix)
    fileNamePbPbMix = fileNamePbPb;
  
  Int_t leadingPtOffset = 1;
    
  Int_t maxLeadingPt = 4;
  Int_t maxAssocPt = 5;
  Float_t leadingPtArr[] = { 2.0, 3.0, 6.0, 6.0, 8.0, 10.0, 15.0, 20.0 };
  Float_t assocPtArr[] =     { 0.15, 0.5, 1.0, 4.0, 3.0, 4.0, 6.0, 8.0, 10.0, 12.0 };
  leadingPtOffset = 1;
  
  AliUEHistograms* h = (AliUEHistograms*) GetUEHistogram(fileNamePbPb);
  hMixed = (AliUEHistograms*) GetUEHistogram(fileNamePbPbMix, 0, kTRUE);
  
  Int_t i=1;
  Int_t j=2;
  
  gpTMin = assocPtArr[j] + 0.01;
  gpTMax = assocPtArr[j+1] - 0.01;
  
  SetupRanges(h);
  SetupRanges(hMixed);

  if (assocPtArr[j] >= leadingPtArr[i+leadingPtOffset])
    continue;

  TString hist1Str, hist2Str, hist2bStr;
  
  Int_t step = 6;
  TH1* hist1 = 0;
  TH1* hist2 = 0;
  TH1* hist3 = 0;
  
  Float_t v2[3];
  
  Bool_t equivMixedBin = kTRUE;
  
  GetDistAndFlow(h, hMixed, &hist1,  0, step, 0,  5,  leadingPtArr[i] + 0.01, leadingPtArr[i+leadingPtOffset] - 0.01, 20, equivMixedBin); 

  GetDistAndFlow(h, hMixed, &hist2,  0, step, 0,  5,  leadingPtArr[i] + 0.01, leadingPtArr[i+leadingPtOffset] - 0.01, 10, equivMixedBin); 
//   hist2->Scale(1.0 / 0.8);
  
  GetDistAndFlow(h, hMixed, &hist3,  0, step, 0,  5,  leadingPtArr[i] + 0.01, leadingPtArr[i+leadingPtOffset] - 0.01, 11, equivMixedBin); 
  hist3->Scale(1.0 / 0.8);
  
  hist1->Draw();
  hist2->SetLineColor(2);
  hist2->Draw("SAME");
  hist3->SetLineColor(4);
  hist3->Draw("SAME");
}  
   
void PlotDeltaPhiEtaGap(const char* fileNamePbPb, const char* fileNamePbPbMix, const char* fileNamepp)
{
  loadlibs();
  
  if (!fileNamePbPbMix)
    fileNamePbPbMix = fileNamePbPb;
  
  file = TFile::Open("dphi_corr.root", "RECREATE");
  file->Close();
  
  Int_t leadingPtOffset = 1;
    
  Int_t maxLeadingPt = 5;
  Int_t maxAssocPt = 6;
  if (1)
  {
    Float_t leadingPtArr[] = { 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 15.0, 20.0 };
    Float_t assocPtArr[] =     { 0.15, 0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 12.0 };
  }
  else
  {
    Float_t leadingPtArr[] = { 0.15, 10.0 };
    Float_t assocPtArr[] =     { 0.15, 10.0 };
  }
  leadingPtOffset = 1;
  
  AliUEHistograms* h = (AliUEHistograms*) GetUEHistogram(fileNamePbPb);
  hMixed = (AliUEHistograms*) GetUEHistogram(fileNamePbPbMix, 0, kTRUE);
  
  AliUEHistograms* h2 = (AliUEHistograms*) GetUEHistogram(fileNamepp);
  hMixed2 = (AliUEHistograms*) GetUEHistogram(fileNamepp, 0, kTRUE);

  if (0)
  {
    h->SetZVtxRange(-0.99, 0.99);
    hMixed->SetZVtxRange(-0.99, 0.99);
    h2->SetZVtxRange(-0.99, 0.99);
    hMixed2->SetZVtxRange(-0.99, 0.99);
  }
  
  for (Int_t i=0; i<maxLeadingPt; i++)
    for (Int_t j=0; j<maxAssocPt; j++)
    {
      gpTMin = assocPtArr[j] + 0.01;
      gpTMax = assocPtArr[j+1] - 0.01;
      
      SetupRanges(h);
      SetupRanges(hMixed);
      SetupRanges(h2);
      SetupRanges(hMixed2);

      if (assocPtArr[j] >= leadingPtArr[i+leadingPtOffset])
	continue;
  
      Int_t step = 6;
      TH1* hist1 = 0;
      TH1* hist2 = 0;
      TH1* hist3 = 0;
      TH1* hist4 = 0;
      
      Bool_t equivMixedBin = 1; //kFALSE; // TODO ?
      Bool_t scaleToPairs = kTRUE;
      
      Int_t histType = 1;
      
//       GetSumOfRatios(h, hMixed, &hist1,  step, 0,  5, leadingPtArr[i] + 0.01, leadingPtArr[i+leadingPtOffset] - 0.01, kTRUE); 

      GetDistAndFlow(h, hMixed, &hist1,  0, step, 0,  5,  leadingPtArr[i] + 0.01, leadingPtArr[i+leadingPtOffset] - 0.01, histType, equivMixedBin, 0, scaleToPairs); 
      
      GetDistAndFlow(h, hMixed, &hist4,  0, step, 20,  30,  leadingPtArr[i] + 0.01, leadingPtArr[i+leadingPtOffset] - 0.01, histType, equivMixedBin, 0, scaleToPairs); 

      GetDistAndFlow(h, hMixed, &hist2,  0, step, 60,  90,  leadingPtArr[i] + 0.01, leadingPtArr[i+leadingPtOffset] - 0.01, histType, equivMixedBin, 0, scaleToPairs);

      GetDistAndFlow(h2, hMixed2, &hist3,  0, step, 0,  -1,  leadingPtArr[i] + 0.01, leadingPtArr[i+leadingPtOffset] - 0.01, histType, equivMixedBin, 0, scaleToPairs);

      file = TFile::Open("dphi_corr.root", "UPDATE");
      
      if (hist1)
      {
	hist1->SetName(Form("dphi_%d_%d_%d", i, j, 0));
	hist1->Write();
      }
      
      if (hist2)
      {
	hist2->SetName(Form("dphi_%d_%d_%d", i, j, 1));
	hist2->Write();
      }
      
      if (hist4)
      {
	hist4->SetName(Form("dphi_%d_%d_%d", i, j, 3));
	hist4->Write();
      }

      if (hist3)
      {
	hist3->SetName(Form("dphi_%d_%d_%d", i, j, 2));
	TString title(hist3->GetTitle());
	title.ReplaceAll("0--1%", "pp");
	hist3->SetTitle(title);
	hist3->Write();
      }
      
      file->Close();
      
//       return;
    }
}

void DrawLatex(Float_t x, Float_t y, Int_t color, const char* text)
{
  latex = new TLatex(x, y, text);
  latex->SetNDC();
  latex->SetTextSize(0.06);
  latex->SetTextColor(color);
  latex->Draw();
}

void DrawChi2NDF(TF1* func, TH1* hist, Float_t x, Float_t y, Int_t color = 1)
{
  Float_t chi2 = 0;
  Int_t ndf = 0;
  for (Int_t i=hist->GetXaxis()->FindBin(-0.8); i<=hist->GetXaxis()->FindBin(0.8); i++)
  {
    chi2 += TMath::Power((hist->GetBinContent(i) - func->Integral(hist->GetXaxis()->GetBinLowEdge(i), hist->GetXaxis()->GetBinUpEdge(i)) / (hist->GetXaxis()->GetBinUpEdge(i) - hist->GetXaxis()->GetBinLowEdge(i))) / hist->GetBinError(i), 2);
    ndf++;
  }
  ndf -= func->GetNumberFreeParameters();
  
  printf("#chi^{2}/ndf = %.1f/%d = %.1f  ", func->GetChisquare(), func->GetNDF(), func->GetChisquare() / func->GetNDF());
  Printf("#chi^{2}/ndf = %.1f/%d = %.1f", chi2, ndf, chi2 / ndf);

  DrawLatex(x, y, color, Form("#chi^{2}/ndf = %.1f/%d = %.1f", func->GetChisquare(), func->GetNDF(), func->GetChisquare() / func->GetNDF()));
  DrawLatex(x, y - 0.05, color, Form("#chi^{2}/ndf = %.1f/%d = %.1f", chi2, ndf, chi2 / ndf));
}   
   
void FitDeltaPhiEtaGap(TH1* hist, Int_t color, TGraphErrors* graph, Float_t x, Float_t yPosChi2)
{
  hist->SetLineColor(color);
  hist->DrawCopy("SAME");

  Bool_t twoGauss = kFALSE;
  
  if (!twoGauss)
    func = new TF1("func", "[0]+gaus(1)");
  else
  {
    func = new TF1("func", "[0]+gaus(1)+gaus(4)");
    func->FixParameter(5, 0);
    func->SetParLimits(3, 0.1, 10);
    func->SetParLimits(6, 0.1, 10);
    func->SetParLimits(1, 0, 10);
    func->SetParLimits(4, 0, 10);
  }
  
  func->SetParameters(0, 1, 0, 0.3, 1, 0, 1);
  func->FixParameter(2, 0);
  func->SetLineColor(color);

  hist->Fit(func, "", "SAME");
//   hist->Fit(func, "IM", "SAME");
  
  if (twoGauss)
  {
    func2 = new TF1("func2", "[0]+gaus(1)", -1.5, 4.5);
    func2->SetParameters(func->GetParameter(0), func->GetParameter(1), func->GetParameter(2), func->GetParameter(3));
    func2->SetLineColor(color);
    func2->SetLineWidth(1);
    func2->SetLineStyle(2);
    func2->Draw("SAME");
    
    func2 = new TF1("func2", "[0]+gaus(1)", -1.5, 4.5);
    func2->SetParameters(func->GetParameter(0), func->GetParameter(4), func->GetParameter(5), func->GetParameter(6));
    func2->SetLineColor(color);
    func2->SetLineWidth(1);
    func2->SetLineStyle(2);
    func2->Draw("SAME");
  }
      
  if (twoGauss)
  {
    Bool_t firstIsMin = func->GetParameter(3) < func->GetParameter(6);
    
    Bool_t onlyOne = kFALSE;
    if (func->GetParameter(1) / func->GetParameter(4) < 0.1)
    {
      firstIsMin = kFALSE;
      onlyOne = kTRUE;
    }
    if (func->GetParameter(1) / func->GetParameter(4) > 10)
    {
      firstIsMin = kTRUE;
      onlyOne = kTRUE;
    }
    
    graph->SetPoint(graph->GetN(), x - 0.1, func->GetParameter((firstIsMin) ? 3 : 6));
    graph->SetPointError(graph->GetN()-1, 0, func->GetParError((firstIsMin) ? 3 : 6));

    if (!onlyOne)
    {
      graph->SetPoint(graph->GetN(), x + 0.1, TMath::Abs(func->GetParameter((!firstIsMin) ? 3 : 6)));
      graph->SetPointError(graph->GetN()-1, 0, func->GetParError((!firstIsMin) ? 3 : 6));
    }
  }
  else
  {
    graph->SetPoint(graph->GetN(), x, TMath::Abs(func->GetParameter(3)));
    graph->SetPointError(graph->GetN()-1, 0, func->GetParError(3));
  }  
    
  DrawChi2NDF(func, hist, 0.5, yPosChi2, color);
}

void AnalyzeDeltaPhiEtaGap(const char* fileName)
{
  TFile::Open(fileName);
  
  Int_t maxLeadingPt = 4;
  Int_t maxAssocPt = 5;

  TCanvas* canvas = new TCanvas("DeltaPhi", "DeltaPhi", 1000, 700);
  canvas->Divide(maxAssocPt, maxLeadingPt);
      
  TGraphErrors* width1 = new TGraphErrors;
  TGraphErrors* width2 = new TGraphErrors;
  TGraphErrors* width3 = new TGraphErrors;
  
  for (Int_t i=0; i<maxLeadingPt; i++)
    for (Int_t j=0; j<maxAssocPt; j++)
    {
      canvas->cd(j+1 + i * maxAssocPt);
      gPad->SetLeftMargin(0.15);
      gPad->SetBottomMargin(0.2);
      gPad->SetTopMargin(0.01);
      gPad->SetRightMargin(0.01);
      
      hist1 = (TH1*) gFile->Get(Form("dphi_%d_%d_%d", i, j, 0));
      hist2 = (TH1*) gFile->Get(Form("dphi_%d_%d_%d", i, j, 1));
      hist3 = (TH1*) gFile->Get(Form("dphi_%d_%d_%d", i, j, 2));
      
      if (!hist1)
	continue;

      TString tmpStr(hist1->GetTitle());
      tmpStr.ReplaceAll(" - ", "#");
      tokens = tmpStr.Tokenize("#");
      TString str(tokens->At(0)->GetName());
      TString str2(tokens->At(1)->GetName());
      
      Prepare1DPlot(hist1);
      Prepare1DPlot(hist2);
      Prepare1DPlot(hist3);

      // remove baseline
      hist1->Fit("pol0", "0", "", 1, 4);
      if (!hist1->GetFunction("pol0"))
	continue;
      hist1->GetFunction("pol0")->SetRange(-10, 10);
      hist1->Add(hist1->GetFunction("pol0"), -1);
      
      hist3->Fit("pol0", "0", "", 1, 4);
      hist3->GetFunction("pol0")->SetRange(-10, 10);
      hist3->Add(hist3->GetFunction("pol0"), -1);
      
      Double_t yMin = 0.01;
      Double_t yMax = -1;
      
      if (yMax < 0)
      {
	yMin = TMath::Min(hist1->GetMinimum(), hist2->GetMinimum()) * 0.97;
	yMax = TMath::Max(hist1->GetMaximum(), (hist2) ? hist2->GetMaximum() : 0.0) * 1.1;
      }
      
//       yMin = hist1->GetMinimum() * 0.9;
      yMax *= 1.5;
    
      TH2F* dummy = new TH2F("dummy", "", 100, hist1->GetXaxis()->GetBinLowEdge(1), hist1->GetXaxis()->GetBinUpEdge(hist1->GetNbinsX()), 1000, yMin, yMax);
      dummy->SetStats(kFALSE);
      dummy->SetXTitle(hist1->GetXaxis()->GetTitle());
      dummy->SetYTitle(hist1->GetYaxis()->GetTitle());
      dummy->SetYTitle("1/N_{trig} dN/d#Delta#phi"); 
      Prepare1DPlot(dummy);
      
      dummy->GetYaxis()->SetTitleOffset(0.8);
    
      dummy->GetXaxis()->SetLabelSize(0.08);
      dummy->GetYaxis()->SetLabelSize(0.08);
      dummy->GetXaxis()->SetTitleSize(0.08);
      dummy->GetYaxis()->SetTitleSize(0.08);
      /*
      dummy->GetYaxis()->SetTitleOffset(0.8);
      */
      
      dummyTmp = dummy->DrawCopy();
      
      // TODO plot yield? baseline problematic?

      Float_t xPos = width1->GetN();
      xPos = j*7+i;

      FitDeltaPhiEtaGap(hist1, 1, width1, xPos, 0.7);
      FitDeltaPhiEtaGap(hist2, 4, width2, xPos, 0.5);
      FitDeltaPhiEtaGap(hist3, 2, width3, xPos, 0.6);
      
      latex = new TLatex(0.3, 0.8, str);
      latex->SetNDC();
      latex->SetTextSize(0.06);
      latex->Draw();
      
      latex = new TLatex(0.3, 0.88, str2);
      latex->SetNDC();
      latex->SetTextSize(0.06);
      latex->Draw();
      
      DrawLatex(0.8, 0.9,  1, "0-5%");
      DrawLatex(0.8, 0.85, 2, "60-90%");
      DrawLatex(0.8, 0.8,  4, "pp");

//       return;
//       i = 10; j = 10;
    }
    
  new TCanvas;
  width1->SetMarkerStyle(20);
  width1->Draw("AP");

  width2->SetMarkerStyle(24);
  width2->SetMarkerColor(4);
  width2->Draw("P SAME");

  width3->SetMarkerStyle(25);
  width3->SetMarkerColor(2);
  width3->Draw("P SAME");
}

void CheckWing(const char* fileName)
{
  TFile::Open(fileName);
  
  Int_t maxLeadingPt = 4;
  Int_t maxAssocPt = 5;

  TCanvas* canvas = new TCanvas("DeltaPhi", "DeltaPhi", 1000, 700);
  canvas->Divide(maxAssocPt, maxLeadingPt);
      
  for (Int_t i=0; i<maxLeadingPt; i++)
    for (Int_t j=0; j<maxAssocPt; j++)
    {
      canvas->cd(j+1 + i * maxAssocPt);
      gPad->SetLeftMargin(0.15);
      gPad->SetBottomMargin(0.2);
//       gPad->SetTopMargin(0.01);
      gPad->SetRightMargin(0.01);
      
      hist1 = (TH1*) gFile->Get(Form("dphi_%d_%d_%d", i, j, 0));
      
      if (!hist1)
	continue;
      
//       hist1->Draw("COLZ");

      Float_t width = 0.5;

      proj = ((TH2*) hist1)->ProjectionY(Form("%s_projx", hist1->GetName()), hist1->GetXaxis()->FindBin(TMath::Pi() - width),hist1->GetXaxis()->FindBin(TMath::Pi() + width));
      
//       proj->GetXaxis()->SetRangeUser(-1.59, 1.59);
      proj->SetStats(kFALSE);
      proj->Draw();

      proj2 = ((TH2*) hist1)->ProjectionY(Form("%s_proj2x", hist1->GetName()), hist1->GetXaxis()->FindBin(TMath::Pi() / 2 - width),hist1->GetXaxis()->FindBin(TMath::Pi() / 2 + width));

//       proj2->GetXaxis()->SetRangeUser(-1.59, 1.59);
      proj2->SetLineColor(2);
      proj2->Draw("SAME");
      
      proj->SetMinimum(proj2->GetMinimum());
    }
}

void CheckWing()
{
  Int_t maxLeadingPt = 4;
  Int_t maxAssocPt = 5;

  TCanvas* canvas = new TCanvas("DeltaPhi", "DeltaPhi", 1000, 700);
  canvas->Divide(maxAssocPt, maxLeadingPt);

//   const char* fileNames[] = { "dphi_corr_allpt_zcentral_01.root", "dphi_corr_allpt_01.root" };
//   const char* fileNames[] = { "dphi_corr_allpt_zcentral.root", "dphi_corr_allpt.root" };
//   const char* fileNames[] = { "dphi_corr_allpt_zcentral.root", "dphi_corr.root" };
//   const char* fileNames[] = { "dphi_corr_allpt_01_zcentral.root", "dphi_corr_allpt_01_zsumofratios.root" };
  const char* fileNames[] = { "dphi_corr_allpt_cfct_01_zcentral.root", "dphi_corr_allpt_cfct_01_zsumofratios.root" };
//   const char* fileNames[] = { "dphi_corr_2d.root", "dphi_corr.root" };
//   const char* fileNames[] = { "dphi_corr_2d_01.root", "dphi_corr.root" };
//   const char* fileNames[] = { "dphi_corr_2d.root", "dphi_corr_2d_vtxzcentral.root" };
//   const char* fileNames[] = { "dphi_corr_2d_01.root", , "dphi_corr_2d_01centr_zvtxcentral.root" }; 
  
  for (Int_t i=0; i<maxLeadingPt; i++)
    for (Int_t j=0; j<maxAssocPt; j++)
    {
      canvas->cd(j+1 + i * maxAssocPt);
      gPad->SetLeftMargin(0.15);
      gPad->SetBottomMargin(0.2);
//       gPad->SetTopMargin(0.01);
      gPad->SetRightMargin(0.01);

      TH1* first = 0;
      for (Int_t fileId = 0; fileId < 2; fileId++)
      {
	TFile::Open(fileNames[fileId]);
  
	hist1 = (TH1*) gFile->Get(Form("dphi_%d_%d_%d", i, j, 0));
	
	if (!hist1)
	  continue;
	
  //       hist1->Draw("COLZ");

	Float_t width = 0.5;

	proj = ((TH2*) hist1)->ProjectionY(Form("%s_%d_projx", hist1->GetName(), fileId), hist1->GetXaxis()->FindBin(TMath::Pi() - width),hist1->GetXaxis()->FindBin(TMath::Pi() + width));
// 	proj->Rebin(2); proj->Scale(0.5);
	
	proj->GetXaxis()->SetRangeUser(-1.79, 1.79);
	proj->SetStats(kFALSE);
	proj->SetLineColor(fileId + 1);
	proj->Draw((fileId == 0) ? "" : "SAME");
	if (!first)
	  first = proj;
	
	proj->Scale(1. / 12);
	
	first->SetMinimum(TMath::Min(first->GetMinimum(), proj->GetMinimum()));
	first->SetMaximum(TMath::Max(first->GetMaximum(), proj->GetMaximum()));
      }
  }
}

void FitDeltaPhiEtaGap2D(TH2* hist, Bool_t scale, TVirtualPad* pad1, TVirtualPad* pad2, TVirtualPad* pad3, TGraphErrors* width1, TGraphErrors* width2, Float_t x, Float_t yPosChi2)
{
  Float_t etaLimit = 1.0;
  Float_t outerLimit = 1.8;
  
  TString histName(hist->GetName());

  TH1D* etaGap = hist->ProjectionX(histName + "_1", TMath::Max(1, hist->GetYaxis()->FindBin(-outerLimit + 0.01)), hist->GetYaxis()->FindBin(-etaLimit - 0.01));
  Int_t etaBins = hist->GetYaxis()->FindBin(-etaLimit - 0.01) - TMath::Max(1, hist->GetYaxis()->FindBin(-outerLimit + 0.01)) + 1;

  TH1D* tracksTmp = hist->ProjectionX(histName + "_2", hist->GetYaxis()->FindBin(etaLimit + 0.01), TMath::Min(hist->GetYaxis()->GetNbins(), hist->GetYaxis()->FindBin(outerLimit - 0.01)));
  etaBins += TMath::Min(hist->GetYaxis()->GetNbins(), hist->GetYaxis()->FindBin(outerLimit - 0.01)) - hist->GetYaxis()->FindBin(etaLimit + 0.01) + 1;
  
  etaGap->Add(tracksTmp);

  // get per bin result
  etaGap->Scale(1.0 / etaBins);
  
//   new TCanvas; etaGap->DrawCopy();
  
  histTmp2D = (TH2*) hist->Clone("histTmp2D");
  histTmp2D->Reset();
  
  for (Int_t xbin=1; xbin<=histTmp2D->GetNbinsX(); xbin++)
    for (Int_t y=1; y<=histTmp2D->GetNbinsY(); y++)
      histTmp2D->SetBinContent(xbin, y, etaGap->GetBinContent(xbin));
    
  if (scale)
  {
    // mixed event does not reproduce away-side perfectly
    // --> extract scaling factor on the away-side from ratios of eta gap and central region
    TH1D* centralRegion = hist->ProjectionX(histName + "_3", hist->GetYaxis()->FindBin(-etaLimit + 0.01), hist->GetYaxis()->FindBin(etaLimit - 0.01));
    etaBins = hist->GetYaxis()->FindBin(etaLimit - 0.01) - hist->GetYaxis()->FindBin(-etaLimit + 0.01) + 1;
    centralRegion->Scale(1.0 / etaBins);
    
//     new TCanvas; centralRegion->DrawCopy(); etaGap->SetLineColor(2); etaGap->DrawCopy("SAME");
    centralRegion->Divide(etaGap);
//     new TCanvas; centralRegion->Draw();
    centralRegion->Fit("pol0", "0", "", TMath::Pi() - 1, TMath::Pi() + 1);
    Float_t scalingFactor = centralRegion->GetFunction("pol0")->GetParameter(0);
    histTmp2D->Scale(scalingFactor);
  }
    
//   new TCanvas; hist->DrawCopy("SURF1");

  hist->Add(histTmp2D, -1);

//   new TCanvas; hist->DrawCopy("SURF1");

  hist->GetYaxis()->SetRangeUser(-1.59, 1.59);
  
  pad1->cd();
  hist->SetStats(0);
  hist->DrawCopy("SURF1");
  
  Float_t min = hist->GetMinimum();
  Float_t max = hist->GetMaximum();
  
  // ranges are to exclude eta gap region from fit
  func = new TF2("func", "[0]+[1]*exp(-0.5*((x/[2])**2+(y/[3])**2))", -5, 5, -1, 1);
  func->SetParameters(0, 1, 0.3, 0.3);
  func->SetParLimits(1, 0, 10);
  func->SetParLimits(2, 0.1, 10);
  func->SetParLimits(3, 0.1, 10);
  
  hist->Fit(func, "0R", "");
//   hist->Fit(func, "IM", "SAME");

  pad2->cd();
  funcHist = (TH2*) hist->Clone("funcHist");
  funcHist->Reset();
  funcHist->Add(func);
  funcHist->SetMinimum(min);
  funcHist->SetMaximum(max);
  funcHist->Draw("SURF1");
  
  pad3->cd();
  hist->Add(func, -1);
  hist->SetMinimum(min);
  hist->SetMaximum(max);
  hist->DrawCopy("SURF1");
  
  width1->SetPoint(width1->GetN(), x, TMath::Abs(func->GetParameter(2)));
  width1->SetPointError(width1->GetN()-1, 0, func->GetParError(2));
    
  width2->SetPoint(width2->GetN(), x, TMath::Abs(func->GetParameter(3)));
  width2->SetPointError(width2->GetN()-1, 0, func->GetParError(3));

  Float_t chi2 = 0;
  Int_t ndf = 0;
  for (Int_t i=hist->GetXaxis()->FindBin(-0.8); i<=hist->GetXaxis()->FindBin(0.8); i++)
    for (Int_t j=hist->GetYaxis()->FindBin(-0.8); j<=hist->GetYaxis()->FindBin(0.8); j++)
    {
      if (hist->GetBinError(i, j) > 0)
      {
	chi2 += TMath::Power(hist->GetBinContent(i, j) / hist->GetBinError(i, j), 2);
	ndf++;
      }
    }
  ndf -= func->GetNumberFreeParameters();
  
  printf("#chi^{2}/ndf = %.1f/%d = %.1f  ", func->GetChisquare(), func->GetNDF(), func->GetChisquare() / func->GetNDF());
  Printf("#chi^{2}/ndf = %.1f/%d = %.1f", chi2, ndf, chi2 / ndf);

  DrawLatex(0.5, yPosChi2, 1, Form("#chi^{2}/ndf = %.1f/%d = %.1f", func->GetChisquare(), func->GetNDF(), func->GetChisquare() / func->GetNDF()));
  DrawLatex(0.5, yPosChi2 - 0.05, 1, Form("#chi^{2}/ndf = %.1f/%d = %.1f", chi2, ndf, chi2 / ndf));
}

void AnalyzeDeltaPhiEtaGap2D(const char* fileName)
{
  TFile::Open(fileName);
  
  Int_t maxLeadingPt = 5;
  Int_t maxAssocPt = 6;

  TGraphErrors* width1[4];
  TGraphErrors* width2[4];
  
  Int_t nHists = 3;
  for (Int_t histId = 0; histId < nHists; histId++)
  {
    width1[histId] = new TGraphErrors;
    width2[histId] = new TGraphErrors;
    for (Int_t i=0; i<maxLeadingPt; i++)
    {
      TCanvas* canvas = new TCanvas(Form("DeltaPhi_%d_%d", histId, i), Form("DeltaPhi_%d_%d", histId, i), 1000, 1000);
      canvas->Divide(3, maxAssocPt);
      
      for (Int_t j=0; j<maxAssocPt; j++)
      {
	for (Int_t k=1; k<=3; k++)
	{
	  canvas->cd(3 * j + k);
	  gPad->SetLeftMargin(0.15);
	  gPad->SetBottomMargin(0.2);
	  gPad->SetTopMargin(0.01);
	  gPad->SetRightMargin(0.01);
	}
	
// 	if (i != 1 || j != 2)
// 	  continue;
    
	hist1 = (TH1*) gFile->Get(Form("dphi_%d_%d_%d", i, j, histId));
	if (!hist1)
	  continue;
	
	Float_t xPos = j*8+i;

	FitDeltaPhiEtaGap2D((TH2*) hist1, kTRUE, canvas->cd(3 * j + 1), canvas->cd(3 * j + 2), canvas->cd(3 * j + 3), width1[histId], width2[histId], xPos, 0.9);
      }
    }
    
//     break;
  }
  
  Int_t marker[] = { 20, 24, 25, 26 };
  Int_t colors[] = { 1, 2, 4, 3 };
  const char* labels[] = { "0-5%", "60-90%", "pp", "20-30%" };
    
  new TCanvas;
  for (Int_t histId = 0; histId < nHists; histId++)
  {
    width1[histId]->SetMarkerStyle(marker[histId]);
    width1[histId]->SetMarkerColor(colors[histId]);
    width1[histId]->Draw((histId == 0) ? "AP" : "PSAME");
    DrawLatex(0.7, 0.8 - 0.05 * histId, colors[histId], labels[histId]);
  }

  new TCanvas;
  for (Int_t histId = 0; histId < nHists; histId++)
  {
    width2[histId]->SetMarkerStyle(marker[histId]);
    width2[histId]->SetMarkerColor(colors[histId]);
    width2[histId]->Draw((histId == 0) ? "AP" : "PSAME");
    DrawLatex(0.7, 0.8 - 0.05 * histId, colors[histId], labels[histId]);
  }
}

  
void PlotPtDistributions(const char* fileName1, Int_t centrBegin = 1, Int_t centrEnd = 2)
{
  loadlibs();

  Int_t maxLeadingPt = 3;
  Float_t leadingPtArr[] = { 1.0, 10.0, 20.0, 40.0 };
  
  Int_t nCentralityBins = 5;
  Int_t centralityBins[] = { 1, 7, 9, 11, 13, 16 };
  
  AliUEHistograms* h = (AliUEHistograms*) GetUEHistogram(fileName1);
    
  TCanvas* canvas = new TCanvas("Pt", "Pt", 1000, 1000);
  canvas->Divide(2, maxLeadingPt+2);
  
  TLegend* legend = new TLegend(0.2, 0.2, 0.95, 0.90);
  legend->SetFillColor(0);
  legend->SetTextSize(0.08);
  
  TLegend* legendB = new TLegend(0.2, 0.2, 0.95, 0.90);
  legendB->SetFillColor(0);
  legendB->SetTextSize(0.08);
  
  TLegend* legend2 = new TLegend(0.2, 0.2, 0.95, 0.90);
  legend2->SetFillColor(0);
  legend2->SetTextSize(0.08);
  
  Int_t colors[] = { 1, 2, 4, 6 };
  Int_t markers[] = { 20, 21, 22, 23, 24, 25 };
  
  for (Int_t i=0; i<maxLeadingPt; i++)
  {
    Double_t ptMin = leadingPtArr[i] + 0.01;
    //Double_t ptMax = leadingPtArr[i+1] - 0.01;
    Double_t ptMax = 39.99;
      
    TString str;
    str.Form("%.1f < p_{T,trig} < %.1f", ptMin - 0.01, ptMax + 0.01);
    
    canvas->cd(2*i+1+2);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.2);
    gPad->SetTopMargin(0.01);
    gPad->SetRightMargin(0.01);
    gPad->SetLogy();
    gPad->SetLogx();
    gPad->SetGridx();
    gPad->SetGridy();
    
    TH2F* dummy = new TH2F("dummy", "", 100, 1, 40, 100, 1e-5, 1e3);
    dummy->SetStats(kFALSE);
    dummy->SetXTitle("p_{T,assoc}");
    dummy->SetYTitle("");
    dummy->GetYaxis()->SetTitleOffset(1);
    Prepare1DPlot(dummy);
  
    dummy->GetXaxis()->SetLabelSize(0.06);
    dummy->GetYaxis()->SetLabelSize(0.06);
    dummy->GetXaxis()->SetTitleSize(0.06);
    dummy->GetYaxis()->SetTitleSize(0.06);
    dummy->DrawCopy();
    
    Float_t phiRange[] = { 0, TMath::Pi() / 2, TMath::Pi() };
    const char* phiLabels[] = { "Towards", "Transverse", "Away" };
    
    for (Int_t j=0; j<3; j++)
    {
      TH1* centralEta = 0;
      TH1* sideEta1 = 0;
      
      if (j == 0)
      {
        centralEta = h->GetUEHist(2)->GetPtHist(6, 0, ptMin, ptMax, centrBegin, centrEnd, phiRange[j] - TMath::Pi() / 4, phiRange[j] + TMath::Pi() / 4, -0.69, 0.69);
        sideEta1 = h->GetUEHist(2)->GetPtHist(6, 0, ptMin, ptMax, centrBegin, centrEnd, phiRange[j] - TMath::Pi() / 4, phiRange[j] + TMath::Pi() / 4, -1.39, -0.71);
        TH1* sideEta2 = h->GetUEHist(2)->GetPtHist(6, 0, ptMin, ptMax, centrBegin, centrEnd, phiRange[j] - TMath::Pi() / 4, phiRange[j] + TMath::Pi() / 4, 0.71, 1.39);
        sideEta1->Add(sideEta2); // TODO can be done smarter? what about the errors?
      
        Prepare1DPlot(sideEta1);
        Prepare1DPlot(centralEta);
        
        centralEta->SetLineColor(colors[j]);
        sideEta1->SetLineColor(colors[j+1]);
      
        if (i == 0)
        {
          legend->AddEntry(centralEta->Clone(), Form("Jet, %s: |#eta| < 0.7, #phi ~ %.1f", phiLabels[j], phiRange[j]));
          legend->AddEntry(sideEta1->Clone(), Form("Ridge, %s: 0.7 < |#eta| < 1.4, #phi ~ %.1f", phiLabels[j], phiRange[j]));
        }
      }
      else
      {
        centralEta = h->GetUEHist(2)->GetPtHist(6, 0, ptMin, ptMax, centrBegin, centrEnd, phiRange[j] - TMath::Pi() / 4, phiRange[j] + TMath::Pi() / 4, -1.39, 1.39);
        centralEta->Scale(0.5);
        centralEta->SetLineColor(colors[j+1]);
      
        Prepare1DPlot(centralEta);
        
        if (i == 0)
          legend->AddEntry(centralEta->Clone(), Form("%s: |#eta| < 1.4, #phi ~ %.1f", phiLabels[j], phiRange[j]));
      }
      
      canvas->cd(2*i+1+2);
      centralEta->DrawCopy("SAME");
      if (sideEta1)
        sideEta1->DrawCopy("SAME");
      
      centralEta->SetLineColor(colors[i]);
      centralEta->Scale(100.0 / centralEta->Integral());
      if (sideEta1)
      {
        sideEta1->SetLineColor(colors[i]);
        sideEta1->Scale(100.0 / sideEta1->Integral());
      }
    
      if (j < 2)
        canvas->cd(j*4+4);
      else
        canvas->cd(j*4+2);
      if (i == 0)
      {
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.2);
        gPad->SetTopMargin(0.01);
        gPad->SetRightMargin(0.01);
        
        dummy->DrawCopy();
        
        gPad->SetLogy();
        gPad->SetLogx();
        gPad->SetGridx();
        gPad->SetGridy();
        
        TString str3;
        if (j == 0)
          str3.Form("Jet, %s: |#eta| < 0.7, #phi ~ %.1f", phiLabels[j], phiRange[j]);
        else
          str3.Form("%s: |#eta| < 1.4, #phi ~ %.1f", phiLabels[j], phiRange[j]);
        latex = new TLatex(0.2, 0.3, str3);
        latex->SetNDC();
        latex->SetTextSize(0.08);
        latex->Draw();
      }
      
      if (j == 0)
        legend2->AddEntry(centralEta->Clone(), str);
        
      centralEta->DrawCopy("SAME");
      
      if (sideEta1)
      {
        canvas->cd(j*4+4+2);
        if (i == 0)
        {
          gPad->SetLeftMargin(0.15);
          gPad->SetBottomMargin(0.2);
          gPad->SetTopMargin(0.01);
          gPad->SetRightMargin(0.01);
          
          dummy->DrawCopy();
          
          gPad->SetLogy();
          gPad->SetLogx();
          gPad->SetGridx();
          gPad->SetGridy();
          
          TString str3;
          if (j == 0)
            str3.Form("Ridge, %s: 0.7 < |#eta| < 1.4, #phi ~ %.1f", phiLabels[j], phiRange[j]);
          else
            str3.Form("%s: |#eta| < 1.4, #phi ~ %.1f", phiLabels[j], phiRange[j]);
          latex = new TLatex(0.2, 0.3, str3);
          latex->SetNDC();
          latex->SetTextSize(0.08);
          latex->Draw();
        }
        
        sideEta1->DrawCopy("SAME");
      }
    }
    
    TString str3;
    str3.Form("%d-%d%%", (Int_t) h->GetCentralityDistribution()->GetXaxis()->GetBinLowEdge(centrBegin), (Int_t) h->GetCentralityDistribution()->GetXaxis()->GetBinUpEdge(centrEnd));
    latex = new TLatex(0.2, 0.3, str3);
    latex->SetNDC();
    latex->SetTextSize(0.08);
    
    latex2 = new TLatex(0.55, 0.8, str);
    latex2->SetNDC();
    latex2->SetTextSize(0.06);
    
    canvas->cd(2*i+1+2);
    latex->Draw();
    latex2->Draw();
    
    //break;
  }
  
  canvas->cd(1);
  legend->Draw();

  canvas->cd(2);
  legend2->Draw();
}

TH1* GetNonSubtractedRHICYield(Int_t ptT, Int_t ptA, Int_t centrality)
{
  rhicFile = TFile::Open("rhic/postcorr_AuAu_iter6_ppg106final.root");
  
  corrFunc = (TH1*) rhicFile->Get(Form("pi0hdphi_%d_%d_%d_2", ptT, ptA, centrality));
  //corrFunc->Draw();
  
  // integral
  corrFunc->Scale(1.0 / corrFunc->Integral());
  // bin width
  corrFunc->Scale(TMath::TwoPi() / corrFunc->GetBinWidth(1));
  
  same     = (TH1*) rhicFile->Get(Form("pi0hdphi_%d_%d_%d_0", ptT, ptA, centrality));
  Float_t nPairs = same->Integral();
  
  triggers = (TH1*) rhicFile->Get(Form("pi0pt_%d_%d", ptT, centrality));
  Float_t nTrigs = triggers->Integral();
  
  graph = (TGraph*) rhicFile->Get(Form("gS_%d_%d", ptT, centrality));
  
  Printf("%f %f %f %f", nPairs, nTrigs, nPairs / nTrigs, graph->GetY()[ptA]);
  
  corrFunc->Scale(nPairs / nTrigs);
 
  //corrFunc->Add(new TF1("func", "1", -5, 5), -1. / 3 * (corrFunc->GetBinContent(corrFunc->FindBin(1.5)) + corrFunc->GetBinContent(corrFunc->FindBin(1.2)) + corrFunc->GetBinContent(corrFunc->FindBin(1.0))));
  
  Float_t values[3];
  Float_t errors[3];
  
  Float_t regionBegin[3] = { -TMath::Pi() / 2,        TMath::Pi() / 2 - 0.4, 1.5 * TMath::Pi() - 0.4 };
  Float_t regionEnd[3] =   { -TMath::Pi() / 2 + 0.4,  TMath::Pi() / 2 + 0.4, 1.5 * TMath::Pi() };
    
  // weighted mean
  Float_t sum = 0;
  Float_t weight = 0;
  for (Int_t i=0; i<3; i++)
  {
    corrFunc->Fit("pol0", "0Q", "", regionBegin[i], regionEnd[i]);
    func = corrFunc->GetFunction("pol0");
    if (!func)
      continue;
    sum += func->GetParameter(0) / func->GetParError(0) / func->GetParError(0);
    weight += 1. / func->GetParError(0) / func->GetParError(0);
  }
  
  if (weight == 0)
    return 0;
  
  sum /= weight;
  weight = TMath::Sqrt(1. / weight);
  
  corrFunc->Add(new TF1("func", "1", -5, 5), -sum);
  
  return corrFunc;
  
  new TCanvas;
  corrFunc->Draw();
  compHist = (TH1*) rhicFile->Get(Form("ptyMSMP_0_%d%d%d", ptT, ptA, centrality));
  compHist->DrawCopy("SAME");
  
  //compHist->Add((TF1*) compHist->GetListOfFunctions()->First());
  //compHist->DrawCopy("SAME");
}
  
void DeltaPhiVsRHIC(Int_t rhicCentrality = 0, Int_t aliceCentrality = 0, Bool_t reduced = kFALSE)
{
  aliceFile = TFile::Open("alice_dphi_corr_rhicbinning.root");
  rhicFile = TFile::Open("rhic/postcorr_AuAu_iter6_ppg106final.root");

  Int_t maxLeadingPt = 4;
  Int_t maxAssocPt = 5;
  
  if (reduced)
  {
    Int_t maxSelected = 4;
    Int_t selectedLead[] = { 0, 0, 2, 2 };
    Int_t selectedAssoc[] = { 1, 2, 3, 4 };
   
    maxLeadingPt = TMath::Sqrt(maxSelected);
    maxAssocPt = TMath::Sqrt(maxSelected);
  }
    
  TCanvas* canvas = new TCanvas("DeltaPhi", "DeltaPhi", 1000, 700);
  canvas->Divide(maxAssocPt, maxLeadingPt);
  
  for (Int_t i=0; i<maxLeadingPt; i++)
    for (Int_t j=0; j<maxAssocPt; j++)
    {
      Printf("%d %d", i, j);
    
      canvas->cd(j+1 + i * maxAssocPt);
      gPad->SetLeftMargin(0.15);
      gPad->SetBottomMargin(0.2);
      //gPad->SetTopMargin(0.01);
      gPad->SetRightMargin(0.01);
      
      Int_t iSel = i;
      Int_t jSel = j;
      
      if (reduced)
      {
        iSel = selectedLead[j + i * maxAssocPt];
        jSel = selectedAssoc[j + i * maxAssocPt];
      }
      
      //rhic = (TH1*) rhicFile->Get(Form("VptyMSMP_0_%d%d%d", iSel, jSel, rhicCentrality));
      rhic = GetNonSubtractedRHICYield(iSel, jSel, rhicCentrality);
      if (!rhic)
        continue;
      rhic->SetLineColor(2);
      rhic->SetMarkerStyle(1);
      
      alice = (TH1*) aliceFile->Get(Form("dphi_%d_%d_%d_fit_flat", iSel, jSel, aliceCentrality));
      if (!alice)
        continue;
        
      // match near side yield
      if (1)
      {
        Float_t factor = 0.5 * alice->Integral(alice->FindBin(-0.1), alice->FindBin(0.1)) / rhic->Integral(rhic->FindBin(-0.1), rhic->FindBin(0.1));
        
        Printf("%f", factor);
        rhic->Scale(factor);
      }
        
      alice->SetLineColor(1);
        
      //alice->Rebin(36); rhic->Rebin(30);
      
      clone = alice->DrawCopy("");
      rhic->DrawCopy("SAME");
      
      //Printf("chi2 test: chi2/ndf = %f", alice->Chi2Test(rhic, "WW CHI2/NDF"));
      Float_t chi2 = 0;
      Float_t n = 0;
      for (Int_t k=1; k<=rhic->GetNbinsX(); k++)
      {
        chi2 += TMath::Power(rhic->GetBinContent(k) - alice->Interpolate(rhic->GetBinCenter(k)), 2) / (TMath::Power(rhic->GetBinError(k), 2) + TMath::Power(alice->GetBinError(alice->FindBin(rhic->GetBinCenter(k))), 2));
        n++;
      }
      
      chi2 /= n;
      Printf("chi2 test: chi2/ndf = %f", chi2);      
      
      clone->GetYaxis()->SetRangeUser(TMath::Min(alice->GetMinimum(), rhic->GetMinimum()) * 1.1, TMath::Max(alice->GetMaximum(), rhic->GetMaximum()) * 1.1);
      
      for (Int_t bin=1; bin<=rhic->GetNbinsX(); bin++)
      {
        Double_t aliceValue = alice->GetBinContent(alice->FindBin(rhic->GetXaxis()->GetBinCenter(bin)));
        if (aliceValue == 0)
          continue;
        rhic->SetBinContent(bin, rhic->GetBinContent(bin) / aliceValue);
        rhic->SetBinError(bin, 0);
      }
      
      rhic->Rebin(2);
      rhic->Scale(0.5);
      
      rhic->SetLineColor(3);
      //rhic->DrawCopy("SAME");
      
      //return;
    }
    
    canvas->SaveAs("yield_comparison.png");
}
  
void DeltaPhi(const char* fileName, const char* fileName2 = 0, Bool_t reduced = kFALSE, Bool_t ppComparison = kFALSE, Int_t mode = 0, const char* dataTag = "", const char* histName = "_fit_flat")
{
  // DeltaPhi("high_stat_binning_pp7_pt8.root", "high_stat_binning_pp7_pythia_pt8.root", 0, 1, 0, "pp 7 TeV uncorrected")
  // DeltaPhi("high_stat_binning_pp900_pt8.root", "high_stat_binning_pp900_pythia_pt8.root", 0, 1, 0, "pp 0.9 TeV uncorrected")

  style(2);

  aliceFile = TFile::Open(fileName);
  if (fileName2)
    secondFile = TFile::Open(fileName2);

  if (1)
  {
    Int_t maxLeadingPt = 3;
    Int_t maxAssocPt = 3;
  }
  else
  {
    Int_t maxLeadingPt = 3;
    Int_t maxAssocPt = 7;
  }
  
  if (reduced)
  {
    Int_t maxSelected = 4;
    Int_t selectedLead[] = { 1, 1, 1, 1 };
    Int_t selectedAssoc[] = { 2, 3, 4, 5 };
   
    maxLeadingPt = TMath::Sqrt(maxSelected);
    maxAssocPt = TMath::Sqrt(maxSelected);
  }
  
  factorGraph = new TGraphErrors;
    
//   TCanvas* canvas = new TCanvas(Form("%s_%s", fileName, fileName2 ? fileName2 : ""), Form("%s_%s", fileName, fileName2 ? fileName2 : ""), 600, 900);
  TCanvas* canvas = new TCanvas(Form("%s_%s", fileName, fileName2 ? fileName2 : ""), Form("%s_%s", fileName, fileName2 ? fileName2 : ""), 1000, 1000);
  canvas->Divide(maxAssocPt, maxLeadingPt);
  
  for (Int_t i=0; i<maxLeadingPt; i++)
    for (Int_t j=0; j<maxAssocPt; j++)
    {
      TH1* first = 0;
      TH1* peripheral = 0;
      for (Int_t aliceCentrality=0; aliceCentrality<4; aliceCentrality++)
      {
        Printf("%d %d %d", i, j, aliceCentrality);
	
	if (aliceCentrality == 1)
	  continue;
      
        canvas->cd(j+1 + i * maxAssocPt);
        
        Int_t iSel = i;
        Int_t jSel = j;
        Int_t centralitySel = aliceCentrality;
        currentFile = aliceFile;
        
        if (reduced)
        {
          iSel = selectedLead[j + i * maxAssocPt];
          jSel = selectedAssoc[j + i * maxAssocPt];
        }
        
        if (fileName2)
        {
          if (ppComparison)
          {
            centralitySel = 3;
            if (aliceCentrality == 1)
              currentFile = secondFile;
            else if (aliceCentrality == 2)
              break;
          }
          else
            if (aliceCentrality == 3)
              currentFile = secondFile;
        }
        
        //alice = (TH1*) currentFile->Get(Form("dphi_%d_%d_%d%s", iSel, jSel, centralitySel, (centralitySel < 3) ? "_tsallis_flat" : ((flatOrTsallis) ? "" : "_fit_flat")));
        alice = (TH1*) currentFile->Get(Form("dphi_%d_%d_%d%s", iSel, jSel, centralitySel, histName));
        if (!alice)
          continue;

	if (0)
	{
	  Printf("WARNING: Applying some scaling and rebinning! Only for 2.76 data-MC comparison!");
	  if (aliceCentrality == 0)
	  {
	    alice->Rebin(2); alice->Scale(0.5); 
  	  alice->Scale(1.0 / 1.6);
	  }
// 	  else
// 	    alice->Scale(1.6);
	}
	
	    alice->Rebin(2); alice->Scale(0.5); 
// 	    alice->Scale(1.0 / 1.6);

	// match near side yield to peripheral
        if (1 && centralitySel == 3 && peripheral)
        {
         if (mode == 0 || mode == 1)
          {
            Double_t width = 0.5;
            Double_t error1, error2;
            Double_t integral1 = peripheral->IntegralAndError(peripheral->FindBin(-width), peripheral->FindBin(width), error1);
            Double_t integral2 = alice->IntegralAndError(alice->FindBin(-width), alice->FindBin(width), error2);
            if (mode == 1)
            {
              Double_t tmpErr = 0;
              integral1 += peripheral->IntegralAndError(peripheral->FindBin(TMath::Pi() - width), peripheral->FindBin(TMath::Pi() + width), tmpErr);
              error1 = TMath::Sqrt(error1 * error1 + tmpErr * tmpErr);
              
              integral2 += alice->IntegralAndError(alice->FindBin(TMath::Pi() - width), alice->FindBin(TMath::Pi() + width), tmpErr);
              error2 = TMath::Sqrt(error2 * error2 + tmpErr * tmpErr);
            }
          }
          else if (mode == 2)
          {
            Double_t error1, error2;
            Double_t integral1 = peripheral->IntegralAndError(1, peripheral->GetNbinsX(), error1);
            Double_t integral2 = alice->IntegralAndError(1, alice->GetNbinsX(), error2);
          }
          else if (mode == 3)
          {
            Double_t width = 1.0;
            Double_t error1, error2;
            Double_t integral1 = peripheral->IntegralAndError(peripheral->FindBin(TMath::Pi() - width), peripheral->FindBin(TMath::Pi() + width), error1);
            Double_t integral2 = alice->IntegralAndError(alice->FindBin(TMath::Pi() - width), alice->FindBin(TMath::Pi() + width), error2);
          }
            
          Double_t factor = integral1 / integral2;
          
          //factor = 0.804 * 0.9;
          Printf("%f", factor);
//           alice->Scale(factor);
          factorGraph->SetPoint(factorGraph->GetN(), factorGraph->GetN(), factor);
          factorGraph->SetPointError(factorGraph->GetN() - 1, 0, factor * TMath::Sqrt(TMath::Power(error1 / integral1, 2) + TMath::Power(error2 / integral2, 2)));
        }
        
        if (ppComparison && aliceCentrality == 0)
          peripheral = alice;
        else if (aliceCentrality == 2)
          peripheral = alice;          
          
//         alice->SetYTitle("1/(N_{trig} #Delta#eta) dN_{assoc}/d#Delta#phi (1/rad)");
        alice->SetYTitle("1/N_{trig} dN_{assoc}/d#Delta#phi (1/rad)");
        alice->SetXTitle("#Delta#phi (rad)");
        alice->SetLineColor(aliceCentrality+1);
        alice->SetLineWidth(2);
        alice->SetMarkerColor(aliceCentrality+1);
        alice->GetYaxis()->SetTitleOffset(1.7);
        clone = alice->DrawCopy((aliceCentrality > 0) ? "SAME" : "");
        clone->SetTitle("");
        
        TString str(alice->GetTitle());
        str.ReplaceAll(" - ", "#");
        tokens = str.Tokenize("#");
          
        if (aliceCentrality == 0)
        {
          for (Int_t k=0; k<2; k++)
          {
            TString str(tokens->At(k)->GetName());
            str.ReplaceAll(".0", "");
            str.ReplaceAll("< p", "GeV/c < p");
            str += " GeV/c";
            latex = new TLatex(0.48, 0.92-k*0.06, str);
            latex->SetNDC();
            latex->SetTextSize(0.04);
            latex->Draw();
          }
          
        }
        
        if (!first)
          first = clone;
        else
          first->GetYaxis()->SetRangeUser(TMath::Min(first->GetMinimum(), clone->GetMinimum()), TMath::Max(first->GetMaximum(), clone->GetMaximum()));
          
        
          
        //return;
      }
      if (first)
        //first->GetYaxis()->SetRangeUser(first->GetMinimum(), first->GetMaximum() * 1.2);
        first->GetYaxis()->SetRangeUser(first->GetMinimum(), first->GetMaximum() * 1.5);
    }
    
  if (0)
  {
    for (Int_t i=1; i<=6; i++)
    {
      canvas->cd(i);
      
      latex = new TLatex(0.58, 0.8, "ALICE preliminary");
      latex->SetTextSize(0.04);
      latex->SetNDC();
      latex->Draw();
    
      latex = new TLatex(0.58, 0.74, Form("%s", dataTag));
      latex->SetTextSize(0.04);
      latex->SetNDC();
      latex->Draw();
      
      latex = new TLatex(0.58, 0.68, "Stat. uncertainties only");
      latex->SetTextSize(0.04);
      latex->SetNDC();
      latex->Draw();
    
      latex = new TLatex(0.58, 0.62, "|#eta| < 0.8");
      latex->SetTextSize(0.04);
      latex->SetNDC();
      latex->Draw();
    
      if (ppComparison)
      {
	legend = new TLegend(0.3, 0.8, 0.47, 0.95);
	legend->SetFillColor(0);
	legend->SetTextSize(0.04);
	legend->AddEntry(peripheral, "Data");
	legend->AddEntry(alice, "Pythia");
	legend->Draw();
      }

      DrawALICELogo(0.75, 0.47, 0.95, 0.6);
    }
  }
  
  canvas->SaveAs(Form("%s.eps", canvas->GetName()));
  canvas->SaveAs(Form("%s.png", canvas->GetName()));
    
  new TCanvas;
  peripheral->Divide(alice);
  peripheral->DrawCopy();
  //alice->DrawCopy("SAME");
  
    
  graphCanvas = (TCanvas*) gROOT->GetListOfCanvases()->FindObject("graphCanvas");
  factorGraph->SetMarkerStyle(20);
  if (!graphCanvas)
  {
    graphCanvas = new TCanvas("graphCanvas", "graphCanvas", 800, 600);
    factorGraph->Draw("AP");
  }
  else
  {
    graphCanvas->cd();
    factorGraph->SetLineColor(2);
    factorGraph->SetMarkerColor(2);
    factorGraph->Draw("SAMEP");
  }
  factorGraph->GetYaxis()->SetRangeUser(0.6, 1.4);
  
  factorGraph->Print();
}

void ComparePPYields(const char* histName = "_fit_flat")
{
  const char* files[] = { "high_stat_binning_pp900_pt8.root", "high_stat_binning_pp900_pythia_pt8.root", "high_stat_binning_pp276_pythia_pt8.root", "high_stat_binning_pp7_pythia_pt8.root", "high_stat_binning_pp7_pt8.root" };
  const char* titles[] = { "0.9 data", "0.9 Pythia", "2.76 Pythia", "7 Pythia", "7 data" };
  
  Int_t colors[] = { 1, 2, 3, 4, 6 };
  
  Int_t maxLeadingPt = 3;
  Int_t maxAssocPt = 2;
  
  TCanvas* canvas = new TCanvas("ComparePPYields", "ComparePPYields", 600, 900);
  canvas->Divide(maxAssocPt, maxLeadingPt);
  
  legend = new TLegend(0.5, 0.5, 0.8, 0.8);
  legend->SetFillColor(0);
  
  TGraphErrors** graphs = new TGraphErrors*[5];
  TGraphErrors** graphs2 = new TGraphErrors*[5];
      
  for (Int_t nFiles = 0; nFiles < 5; nFiles++)
  {
    graphs[nFiles] = new TGraphErrors;
    graphs2[nFiles] = new TGraphErrors;
  }
  
  for (Int_t i=0; i<maxLeadingPt; i++)
    for (Int_t j=0; j<maxAssocPt; j++)
    {
      canvas->cd(j+1 + i * maxAssocPt);
      gPad->SetLeftMargin(0.18);
      gPad->SetBottomMargin(0.1);
      gPad->SetTopMargin(0.01);
      gPad->SetRightMargin(0.01);
        
      TH1* first = 0;
        
      for (Int_t nFiles = 0; nFiles < 5; nFiles++)
      {
        currentFile = TFile::Open(files[nFiles]);

        Int_t aliceCentrality = 3;
          
        Printf("%d %d %d", i, j, nFiles);
      
        alice = (TH1*) currentFile->Get(Form("dphi_%d_%d_%d%s", i, j, aliceCentrality, histName));
        if (!alice)
          continue;
          
        alice->SetYTitle("1/N_{trig} dN/dp_{T,assoc}");
        alice->SetLineColor(colors[nFiles]);
        alice->SetMarkerColor(colors[nFiles]);
        alice->GetYaxis()->SetTitleOffset(1.7);
        clone = alice->DrawCopy((nFiles > 0) ? "SAME" : "");
        clone->SetTitle("");
        
        Double_t width = 0.7;
        Double_t error, error2;
        Double_t integral = alice->IntegralAndError(alice->FindBin(-width), alice->FindBin(width), error);
        Double_t integral2 = alice->IntegralAndError(alice->FindBin(TMath::Pi() - width), alice->FindBin(TMath::Pi() + width), error2);
        
        graphs[nFiles]->SetPoint(graphs[nFiles]->GetN(), j + i * maxAssocPt - 0.2, integral);
        graphs[nFiles]->SetPointError(graphs[nFiles]->GetN()-1, 0, error);
        
        graphs2[nFiles]->SetPoint(graphs2[nFiles]->GetN(), j + i * maxAssocPt + 0.2, integral2);
        graphs2[nFiles]->SetPointError(graphs2[nFiles]->GetN()-1, 0, error2);
        
        if (!first)
          first = clone;
        else
          first->GetYaxis()->SetRangeUser(TMath::Min(first->GetMinimum(), clone->GetMinimum()), TMath::Max(first->GetMaximum(), clone->GetMaximum()));
          
        if (nFiles == 0)
        {
          TString str(alice->GetTitle());
          str.ReplaceAll(" - ", "#");
          tokens = str.Tokenize("#");
            
          for (Int_t k=0; k<2; k++)
          {
            latex = new TLatex(0.6, 0.88-k*0.07, tokens->At(k)->GetName());
            latex->SetNDC();
            latex->SetTextSize(0.04);
            latex->Draw();
          }
        }
          
        if (i == 0 && j == 0)
          legend->AddEntry(clone, titles[nFiles]);
      }
      if (first)
        first->GetYaxis()->SetRangeUser(first->GetMinimum(), first->GetMaximum() * 1.1);
    }
    
  canvas->cd(1);
  
  legend->Draw();
  
  new TCanvas;
  for (Int_t nFiles = 0; nFiles < 5; nFiles++)
  {
    graphs[nFiles]->SetMarkerStyle(25);
    graphs[nFiles]->SetLineColor(colors[nFiles]);
    graphs[nFiles]->SetMarkerColor(colors[nFiles]);
    clone2 = (TGraphErrors*) graphs[nFiles]->DrawClone((nFiles == 0) ? "AP" : "PSAME");
    clone2->GetYaxis()->SetRangeUser(0, 5);
  
    graphs2[nFiles]->SetMarkerStyle(26);
    graphs2[nFiles]->SetLineColor(colors[nFiles]);
    graphs2[nFiles]->SetMarkerColor(colors[nFiles]);
    graphs2[nFiles]->DrawClone("PSAME");
  }  

  new TCanvas;
  for (Int_t nFiles = 0; nFiles < 5; nFiles++)
  {
    DivideGraphs(graphs[nFiles], graphs[4]);
    DivideGraphs(graphs2[nFiles], graphs2[4]);
    
    graphs[nFiles]->SetMarkerStyle(25);
    graphs[nFiles]->SetLineColor(colors[nFiles]);
    graphs[nFiles]->SetMarkerColor(colors[nFiles]);
    clone2 = (TGraphErrors*) graphs[nFiles]->DrawClone((nFiles == 0) ? "AP" : "PSAME");
    clone2->GetYaxis()->SetRangeUser(0, 1.5);
  
    graphs2[nFiles]->SetMarkerStyle(26);
    graphs2[nFiles]->SetLineColor(colors[nFiles]);
    graphs2[nFiles]->SetMarkerColor(colors[nFiles]);
    graphs2[nFiles]->DrawClone("PSAME");
  }  
}

void DeltaPhiPreliminary(const char* fileName, const char* histName)
{
  style(2);

  currentFile = TFile::Open(fileName);

  TCanvas* canvas = new TCanvas(Form("dphi%s", histName), "dphi", 800, 800);
  canvas->Divide(2, 2);
  
  TLegend* legend = new TLegend(0.55, 0.47, 0.85, 0.65);
  legend->SetFillColor(0);
  legend->SetTextSize(0.04);
  
  Int_t colors[] = { 1, 2, 4, 6 };
  
  Int_t pad = 1;
  for (Int_t i=1; i<2; i++)
    for (Int_t j=2; j<6; j++)
    {
      canvas->cd(pad++);
  
      for (Int_t aliceCentrality=0; aliceCentrality<3; aliceCentrality+=2)
      {
        Printf("%d %d %d", i, j, aliceCentrality);
      
        hist = (TH1*) currentFile->Get(Form("dphi_%d_%d_%d%s", i, j, aliceCentrality, histName));
        if (!hist)
          continue;
          
        if (aliceCentrality == 3)
          hist->Scale(kPythiaScalingFactor);
          
        hist->GetYaxis()->SetTitleOffset(1.9);
        hist->SetLineColor(colors[aliceCentrality]);
        hist->SetMarkerColor(colors[aliceCentrality]);
        hist->SetLineWidth(2);
        hist->SetYTitle("1/N_{trig} dN_{assoc}/d#Delta#phi");
        hist->SetXTitle("#Delta#phi (rad)");
        
        clone = hist->DrawCopy((aliceCentrality > 0) ? "SAME" : "");
        clone->SetTitle("");
        
        TString str(hist->GetTitle());
        str.ReplaceAll(" - ", "#");
        tokens = str.Tokenize("#");
        hist->SetTitle("");
          
        if (aliceCentrality == 0)
        {
          for (Int_t k=0; k<2; k++)
          {
            TString str(tokens->At(k)->GetName());
            str.ReplaceAll(".0", "");
            str.ReplaceAll("< p", "GeV/c < p");
            str += " GeV/c";
          
            latex = new TLatex(0.5, 0.88-k*0.06, str);
            latex->SetNDC();
            latex->SetTextSize(0.04);
            latex->Draw();
          }
          
          latex = new TLatex(0.5, 0.76, "ALICE preliminary");
          latex->SetTextSize(0.04);
          latex->SetNDC();
          latex->Draw();
        
          latex = new TLatex(0.5, 0.70, "Stat. uncertainties only");
          latex->SetTextSize(0.04);
          latex->SetNDC();
          latex->Draw();

          clone->GetYaxis()->SetRangeUser(0, clone->GetMaximum() * 1.2);
        }
        
        if (pad == 2)
        {
          if (aliceCentrality == 3)
            legend->AddEntry(hist, "Pythia");
          else
            legend->AddEntry(hist, tokens->At(2)->GetName());
        }
        
        legend->Draw();
        
      }
    }
    
  canvas->cd(3);
  
  canvas->SaveAs(Form("%s.eps", canvas->GetName()));
  canvas->SaveAs(Form("%s.png", canvas->GetName()));
}

void DeltaPhiBaseLinePreliminary(const char* fileName)
{
  style(2);
  
  loadlibs();

  Float_t leadingPtArr[] = { 6.0, 8.0, 10.0, 15.0, 15.0 };
  Float_t assocPtArr[] =     { 0.5, 1.5, 3.0, 4.0, 6.0, 8.0, 10.0, 12.0 };
  AliUEHistograms* h = (AliUEHistograms*) GetUEHistogram(fileName);

  Int_t i = 1;
  Int_t j = 3;
  Int_t step = 0;

  gpTMin = assocPtArr[j] + 0.01;
  gpTMax = assocPtArr[j+1] - 0.01;
  SetupRanges(h);
      
  TString str;
  str.Form("%.0f GeV/c < p_{T,trig} < %.0f GeV/c", leadingPtArr[i], leadingPtArr[i+2]);
  
  TString str2;
  str2.Form("%.0f GeV/c < p_{T,assoc} < %.0f GeV/c", gpTMin - 0.01, gpTMax + 0.01);
          
  Float_t v2[3];
  TH1* hist1 = 0;
  TH1* hist2b = 0;
  
  GetDistAndFlow(h, 0, &hist1,  v2, step, 0,  5,  leadingPtArr[i] + 0.01, leadingPtArr[i+2] - 0.01); 
  GetDistAndFlow(h, 0, &hist2b, v2+2, step, 60, 90, leadingPtArr[i] + 0.01, leadingPtArr[i+2] - 0.01); 
  
  hist1->SetLineWidth(2);
  hist2b->SetLineWidth(2);

  canvas = new TCanvas("dphi_baseline", "dphi_baseline", 1200, 1200);
  canvas->Divide(2, 2);
  
  canvas->cd(1);
//   gPad->SetLeftMargin(0.20);
  
  hist1->GetYaxis()->SetTitleOffset(1.9);
  hist1->SetYTitle("1/N_{trig} dN_{assoc}/d#Delta#phi");
  hist1->SetXTitle("#Delta#phi (rad)");
  hist1->SetTitle("");
  hist1->DrawCopy()->GetYaxis()->SetRangeUser(0, 0.59);
  
  DrawFlow(v2[0], (TH1*) hist1->Clone(), leadingPtArr[i], assocPtArr[j], 0, i, 0, 0, (assocPtArr[j] + assocPtArr[j+1]) / 2, (assocPtArr[j+1] - assocPtArr[j]) / 2);
          
  latex = new TLatex(0.5, 0.78, "0-5%");
  latex->SetNDC();
  latex->SetTextSize(0.04);
  latex->Draw();
  
  if (1)
  {
    canvas->cd(2);
  //   gPad->SetLeftMargin(0.20);
    
    hist1->DrawCopy()->GetYaxis()->SetRangeUser(0.18, 0.295);
    DrawFlow(v2[0], hist1, leadingPtArr[i], assocPtArr[j], 0, i, 0, 0, (assocPtArr[j] + assocPtArr[j+1]) / 2, (assocPtArr[j+1] - assocPtArr[j]) / 2);
    
    latex->Draw();
    
    canvas->cd(3);
  //   gPad->SetLeftMargin(0.20);
    
    hist2b->GetYaxis()->SetTitleOffset(1.9);
    hist2b->SetYTitle("1/N_{trig} dN_{assoc}/d#Delta#phi");
    hist2b->SetXTitle("#Delta#phi (rad)");
    hist2b->SetTitle("");
    hist2b->DrawCopy();
            
    DrawFlow(v2[2], (TH1*) hist2b->Clone(), leadingPtArr[i], assocPtArr[j], 0, i, 2, 0, (assocPtArr[j] + assocPtArr[j+1]) / 2, (assocPtArr[j+1] - assocPtArr[j]) / 2);
    
    latex = new TLatex(0.5, 0.78, "60-90%");
    latex->SetNDC();
    latex->SetTextSize(0.04);
    latex->Draw();
    
    canvas->cd(4);
  //   gPad->SetLeftMargin(0.20);
        
    hist2b->DrawCopy()->GetYaxis()->SetRangeUser(0.01, 0.13);
    DrawFlow(v2[2], hist2b, leadingPtArr[i], assocPtArr[j], 0, i, 2, 0, (assocPtArr[j] + assocPtArr[j+1]) / 2, (assocPtArr[j+1] - assocPtArr[j]) / 2);
        
    latex->Draw();
  }
  
  for (Int_t i=1; i<=4; i++)
  {
    canvas->cd(i);
    
    latex = new TLatex(0.5, 0.84, str);
    latex->SetNDC();
    latex->SetTextSize(0.04);
    latex->Draw();
    
    latex = new TLatex(0.5, 0.90, str2);
    latex->SetNDC();
    latex->SetTextSize(0.04);
    latex->Draw();
  
    latex = new TLatex(0.5, 0.72, "ALICE preliminary");
    latex->SetTextSize(0.04);
    latex->SetNDC();
    latex->Draw();
  
    latex = new TLatex(0.5, 0.66, "Stat. uncertainties only");
    latex->SetTextSize(0.04);
    latex->SetNDC();
    latex->Draw();
  }
    
  canvas->SaveAs(Form("%s.eps", canvas->GetName()));
  canvas->SaveAs(Form("%s.png", canvas->GetName()));
}

/*
void DeltaPhiBaseLinePaperPlot(const char* fileName, const char* yieldFile)
{
  style(2);
  
  Float_t fontSize = 0.08;
  Float_t titleOffset = 0.85;
  
  gStyle->SetTextSize(fontSize);
  gStyle->SetLabelSize(fontSize, "xy");
  gStyle->SetTitleSize(fontSize, "xy");
  gStyle->SetTitleOffset(titleOffset, "y");
  gStyle->SetHistLineWidth(2);
  gROOT->ForceStyle();
  
  loadlibs();

  Float_t leadingPtArr[] = { 6.0, 8.0, 10.0, 15.0, 15.0 };
  Float_t assocPtArr[] =     { 0.5, 1.5, 3.0, 4.0, 6.0, 8.0, 10.0, 12.0 };
  AliUEHistograms* h = (AliUEHistograms*) GetUEHistogram(fileName);

  Int_t i = 1;
  Int_t j = 3;
  Int_t step = 0;

  gpTMin = assocPtArr[j] + 0.01;
  gpTMax = assocPtArr[j+1] - 0.01;
  SetupRanges(h);
      
  TString str;
  str.Form("%.0f GeV/#font[12]{c} < p_{t,trig} < %.0f GeV/#font[12]{c}", leadingPtArr[i], leadingPtArr[i+2]);
  
  TString str2;
  str2.Form("%.0f GeV/#font[12]{c} < p_{t,assoc} < %.0f GeV/#font[12]{c}", gpTMin - 0.01, gpTMax + 0.01);
          
  Float_t v2[3];
  TH1* hist1 = 0;
  TH1* hist2b = 0;
  
  GetDistAndFlow(h, 0, &hist1,  v2, step, 0,  5,  leadingPtArr[i] + 0.01, leadingPtArr[i+2] - 0.01); 
  
  hist1->SetLineWidth(2);

  canvas = new TCanvas("dphi_baseline", "dphi_baseline", 600, 840);
  canvas->Range(0, 0, 1, 1);
  
  pad1 = new TPad("pad1", "pad1", 0, 0.69, 1, 1);
  pad1->Draw();
  
  pad2 = new TPad("pad2", "pad2", 0, 0.38, 1, 0.69);
  pad2->Draw();

  pad3 = new TPad("pad3", "pad3", 0, 0   , 1, 0.38);
  pad3->Draw();
  
  pad1->cd();
  gPad->SetMargin(0.15, 0.02, 0, 0.01);

//   hist1->GetYaxis()->SetTitleOffset(0.85);
  hist1->SetYTitle("1/N_{trig} dN_{assoc}/d#Delta#varphi");
  hist1->SetXTitle("#Delta#varphi (rad)");
  hist1->SetTitle("");
  
  hist1->SetLabelSize(fontSize * 0.38 / 0.31, "xy");
  hist1->SetTitleSize(fontSize * 0.38 / 0.31, "xy");
  hist1->SetTitleOffset(titleOffset / 0.38 * 0.31, "y");
  
  hist1->DrawCopy(); //->GetYaxis()->SetRangeUser(0, 0.59);
  
  //DrawFlow(v2[0], (TH1*) hist1->Clone(), leadingPtArr[i], assocPtArr[j], 0, i, 0, 0, (assocPtArr[j] + assocPtArr[j+1]) / 2, (assocPtArr[j+1] - assocPtArr[j]) / 2);
          
  pad2->cd();
  gPad->SetMargin(0.15, 0.02, 0, 0);

  hist1->DrawCopy()->GetYaxis()->SetRangeUser(0.381, 0.457);
  DrawFlow(v2[0], hist1, leadingPtArr[i], assocPtArr[j], 0, i, 0, 0, (assocPtArr[j] + assocPtArr[j+1]) / 2, (assocPtArr[j+1] - assocPtArr[j]) / 2);
  
  pad3->cd();
  gPad->SetMargin(0.15, 0.02, 0.2, 0);
  
  // get zero subtracted plots
  TFile::Open(yieldFile);
  hist1sub = (TH1*) gFile->Get(Form("dphi_%d_%d_%d_fit_flat", i, j, 0));
  hist1sub->SetYTitle("1/N_{trig} dN_{assoc}/d#Delta#varphi");
  hist1sub->SetXTitle("#Delta#varphi (rad)");
  hist1sub->SetTitle("");
  hist1sub->DrawCopy();
  
  hist2sub = (TH1*) gFile->Get(Form("dphi_%d_%d_%d_fit_flat", i, j, 2));
  hist2sub->SetMarkerStyle(24);
  hist2sub->SetLineColor(2);
  hist2sub->SetMarkerColor(2);
//   hist2sub->SetMarkerSize(0.8);
  hist2sub->DrawCopy("SAMEP E X0");
  
  hist3sub = (TH1*) gFile->Get(Form("dphi_%d_%d_%d_fit_flat", i, j, 3));
  hist3sub->SetMarkerStyle(25);
  hist3sub->SetLineColor(4);
  hist3sub->SetMarkerColor(4);
//   hist3sub->SetMarkerSize(0.8);
  hist3sub->DrawCopy("SAMEP E X0");

  legend = new TLegend(0.4, 0.55, 0.9, 0.9);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  legend->SetTextSize(0.08);
  legend->AddEntry(hist1sub, "Pb-Pb 0-5% centrality", "L");
  legend->AddEntry(hist2sub, "Pb-Pb 60-90% centrality", "P");
  legend->AddEntry(hist3sub, "pp", "P");
  legend->Draw();

  pad1->cd();
  
  latex = new TLatex(0.47, 0.85, str);
  latex->SetTextSize(fontSize * 0.38 / 0.31);
  latex->SetNDC();
  latex->Draw();
    
  latex = new TLatex(0.47, 0.73, str2);
  latex->SetTextSize(fontSize * 0.38 / 0.31);
  latex->SetNDC();
  latex->Draw();
  
  latex = new TLatex(0.5, 0.42, "0-5% centrality");
  latex->SetNDC();
//   latex->Draw();
  
  latex = new TLatex(0.45, 0.62, "Stat. uncertainties only");
  latex->SetNDC();
//   latex->Draw();
    
  canvas->SaveAs(Form("%s.eps", canvas->GetName()));
  canvas->SaveAs(Form("%s.png", canvas->GetName()));
}
*/

void DeltaPhiBaseLinePaperPlot(const char* fileName, const char* yieldFile)
{
  style(2);
  
  Float_t fontSize = 0.08;
  Float_t titleOffset = 0.85;
  
  gStyle->SetTextSize(fontSize);
  gStyle->SetLabelSize(fontSize, "xy");
  gStyle->SetTitleSize(fontSize, "xy");
  gStyle->SetTitleOffset(titleOffset, "y");
  gStyle->SetHistLineWidth(2);
  gROOT->ForceStyle();
  
  loadlibs();

  Float_t leadingPtArr[] = { 6.0, 8.0, 10.0, 15.0, 15.0 };
  Float_t assocPtArr[] =     { 0.5, 1.5, 3.0, 4.0, 6.0, 8.0, 10.0, 12.0 };
  AliUEHistograms* h = (AliUEHistograms*) GetUEHistogram(fileName);

  Int_t i = 1;
  Int_t j = 3;
  Int_t step = 0;

  gpTMin = assocPtArr[j] + 0.01;
  gpTMax = assocPtArr[j+1] - 0.01;
  SetupRanges(h);
      
  TString str;
  str.Form("%.0f < #font[12]{p}_{T,trig} < %.0f GeV/#font[12]{c}", leadingPtArr[i], leadingPtArr[i+2]);
  
  TString str2;
  str2.Form("%.0f < #font[12]{p}_{T,assoc} < %.0f GeV/#font[12]{c}", gpTMin - 0.01, gpTMax + 0.01);
          
  Float_t v2[3];
  TH1* hist1 = 0;
  TH1* hist2b = 0;
  
  GetDistAndFlow(h, 0, &hist1,  v2, step, 0,  5,  leadingPtArr[i] + 0.01, leadingPtArr[i+2] - 0.01); 
  
  hist1->SetLineWidth(2);

  canvas = new TCanvas("dphi_baseline", "dphi_baseline", 600, 840);
  canvas->Range(0, 0, 1, 1);
  
  pad1 = new TPad("pad1", "pad1", 0, 0.69, 1, 1);
  pad1->Draw();
  
  pad2 = new TPad("pad2", "pad2", 0, 0.38, 1, 0.69);
  pad2->Draw();

  pad3 = new TPad("pad3", "pad3", 0, 0   , 1, 0.38);
  pad3->Draw();
  
  pad1->cd();
  gPad->SetMargin(0.15, 0.02, 0, 0.01);

//   hist1->GetYaxis()->SetTitleOffset(0.85);
  hist1->SetYTitle("");
  hist1->SetXTitle("#Delta#font[12]{#varphi} (rad)");
  hist1->SetTitle("");
  
  hist1->SetLabelSize(fontSize * 0.38 / 0.31, "xy");
  hist1->SetTitleSize(fontSize * 0.38 / 0.31, "xy");
  hist1->SetTitleOffset(titleOffset / 0.38 * 0.31, "y");
  
  hist1->DrawCopy("HISTE"); //->GetYaxis()->SetRangeUser(0, 0.59);
  
  //DrawFlow(v2[0], (TH1*) hist1->Clone(), leadingPtArr[i], assocPtArr[j], 0, i, 0, 0, (assocPtArr[j] + assocPtArr[j+1]) / 2, (assocPtArr[j+1] - assocPtArr[j]) / 2);
          
  pad2->cd();
  gPad->SetMargin(0.15, 0.02, 0, 0);

  hist1->SetYTitle("1/#font[12]{N}_{trig} d#font[12]{N}_{assoc}/d#Delta#font[12]{#varphi} (rad^{-1})");
  hist1->DrawCopy("HISTE")->GetYaxis()->SetRangeUser(0.381, 0.457);
  DrawFlow(v2[0], hist1, leadingPtArr[i], assocPtArr[j], 0, i, 0, 0, (assocPtArr[j] + assocPtArr[j+1]) / 2, (assocPtArr[j+1] - assocPtArr[j]) / 2);
  
  pad3->cd();
  gPad->SetMargin(0.15, 0.02, 0.2, 0);
  
  // get zero subtracted plots
  TFile::Open(yieldFile);

  hist1sub = (TH1*) gFile->Get(Form("dphi_%d_%d_%d_fit_flat", i, j, 0));
  hist1sub->SetYTitle("");
  hist1sub->SetXTitle("#Delta#varphi (rad)");
  hist1sub->SetTitle("");
//   hist1sub->SetMarkerStyle(25);
//   hist1sub->SetLineColor(4);
//   hist1sub->SetMarkerColor(4);
  hist1sub->SetMaximum(0.79);
  hist1sub->DrawCopy("HIST E");
  
  hist3sub = (TH1*) gFile->Get(Form("dphi_%d_%d_%d_fit_flat", i, j, 3));
//   hist3sub->SetFillColor(4); 
//   hist3sub->SetFillStyle(3354);
//   hist3sub->SetLineWidth(0);
//   hist3sub->SetLineColor(0);
  hist3sub->SetMarkerStyle(25);
  hist3sub->SetLineColor(4);
  hist3sub->SetMarkerColor(4);
  hist3sub->DrawCopy("SAME E X0");

  hist1sub->DrawCopy("SAME HIST E");
  
  hist2sub = (TH1*) gFile->Get(Form("dphi_%d_%d_%d_fit_flat", i, j, 2));
  hist2sub->SetMarkerStyle(24);
  hist2sub->SetLineColor(2);
  hist2sub->SetMarkerColor(2);
  hist2sub->DrawCopy("SAMEP E X0");
  
  legend = new TLegend(0.4, 0.45, 0.9, 0.8);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  legend->SetTextSize(fontSize);
  legend->AddEntry(hist1sub, "Pb-Pb 0-5% centrality", "L");
  legend->AddEntry(hist2sub, "Pb-Pb 60-90% centrality", "P");
  legend->AddEntry(hist3sub, "pp", "P");
  legend->Draw();

  pad1->cd();
  
  latex = new TLatex(0.59, 0.70, str);
  latex->SetTextSize(fontSize * 0.38 / 0.31);
  latex->SetNDC();
  latex->Draw();
    
  latex = new TLatex(0.59, 0.58, str2);
  latex->SetTextSize(fontSize * 0.38 / 0.31);
  latex->SetNDC();
  latex->Draw();
  
  latex = new TLatex(0.59, 0.46, "#sqrt{#font[12]{s}_{NN}} = 2.76 TeV");
  latex->SetTextSize(fontSize * 0.38 / 0.31);
  latex->SetNDC();
  latex->Draw();

  latex = new TLatex(0.43, 0.90, "a) not background subtracted");
  latex->SetTextSize(fontSize * 0.38 / 0.31);
  latex->SetNDC();
  latex->Draw();

  latex = new TLatex(0.5, 0.42, "0-5% centrality");
  latex->SetNDC();
//   latex->Draw();
  
  latex = new TLatex(0.45, 0.62, "Stat. uncertainties only");
  latex->SetNDC();
//   latex->Draw();

  pad2->cd();

  latex = new TLatex(0.43, 0.90, "b) zoomed");
  latex->SetTextSize(fontSize * 0.38 / 0.31);
  latex->SetNDC();
  latex->Draw();
  
  pad3->cd();

  latex = new TLatex(0.43, 0.90, "c) background subtracted");
  latex->SetNDC();
  latex->Draw();
  
    
  canvas->SaveAs(Form("%s.eps", canvas->GetName()));
  canvas->SaveAs(Form("%s.png", canvas->GetName()));
}


void DeltaPhiRidgePreliminary(const char* fileName)
{
  style();
  
  loadlibs();

  Float_t leadingPtArr[] = { 2.0, 2.0, 3.0, 4.0, 10.0, 20.0, 40.0 };
  Float_t assocPtArr[] =   { 1.0, 2.0, 3.0, 6.0, 10.0, 20.0, 40.0 };
//   Float_t leadingPtArr[] = { 6.0, 8.0, 10.0, 15.0, 15.0 };
//   Float_t assocPtArr[] =     { 0.5, 1.5, 3.0, 4.0, 6.0, 8.0, 10.0, 12.0 };
  AliUEHistograms* h = (AliUEHistograms*) GetUEHistogram(fileName, 0, kFALSE);
  AliUEHistograms* hMixed = (AliUEHistograms*) GetUEHistogram(fileName, 0, kTRUE);

//   Int_t i = 1;
//   Int_t j = 3;
  Int_t i = 0;
  Int_t j = 0;
  Int_t step = 6;

  gpTMin = assocPtArr[j] + 0.01;
  gpTMax = assocPtArr[j+1] - 0.01;
  SetupRanges(h);
  SetupRanges(hMixed);
      
  TString str;
  str.Form("%.1f < p_{T,trig} < %.1f", leadingPtArr[i], leadingPtArr[i+2]);
  
  TString str2;
  str2.Form("%.1f < p_{T,assoc} < %.1f", gpTMin - 0.01, gpTMax + 0.01);
          
  Float_t v2[3];
  TH1* hist1 = 0;
  TH1* hist2b = 0;
  TH1* hist1Peak = 0;
  TH1* hist2bPeak = 0;
  TH1* hist1Ridge = 0;
  TH1* hist2bRidge = 0;
  
  TH1* hist2d = 0;
  TH1* hist2dMixed = 0;
  
//   GetDistAndFlow(h, 0, &hist2d,  v2, step, 0,  5,  leadingPtArr[i] + 0.01, leadingPtArr[i+2] - 0.01, 1); 
//   GetDistAndFlow(hMixed, 0, &hist2dMixed,  v2, step, 0,  5,  leadingPtArr[i] + 0.01, leadingPtArr[i+2] - 0.01, 1); 
//   new TCanvas;
//   //hist2d->Divide(hist2dMixed);
//   hist2dMixed->Draw("SURF1");
  
  GetDistAndFlow(h, hMixed, &hist1,  0, step, 0,  5,  leadingPtArr[i] + 0.01, leadingPtArr[i+2] - 0.01); 
//   GetDistAndFlow(h, 0, &hist2b, v2[2], step, 60, 90, leadingPtArr[i] + 0.01, leadingPtArr[i+2] - 0.01); 
  GetDistAndFlow(h, hMixed, &hist1Peak,  0, step, 0,  5,  leadingPtArr[i] + 0.01, leadingPtArr[i+2] - 0.01, 10); 

  GetDistAndFlow(h, hMixed, &hist1Ridge,  0, step, 0,  5,  leadingPtArr[i] + 0.01, leadingPtArr[i+2] - 0.01, 11); 

  // TODO normalize

  canvas = new TCanvas("dphi_baseline", "dphi_baseline", 800, 800);
  canvas->Divide(2, 2);
  
  canvas->cd(1);
  gPad->SetLeftMargin(0.20);
  
  hist1->GetYaxis()->SetTitleOffset(1.9);
  hist1->SetYTitle("1/N_{trig} dN/dp_{T,assoc}");
  hist1->SetTitle("");
  hist1->DrawCopy();
  
  hist1Peak->SetLineColor(2);
  hist1Peak->DrawCopy("SAME");
  
  hist1Ridge->SetLineColor(4);
  hist1Ridge->DrawCopy("SAME");
  
//   hist1Peak->Add(hist1Ridge);
//   hist1Peak->DrawCopy("SAME")->SetLineColor(3);
  
  return;
  
  DrawFlow(v2[0], (TH1*) hist1->Clone(), leadingPtArr[i], assocPtArr[j], 0, i, 0, 0, (assocPtArr[j] + assocPtArr[j+1]) / 2, (assocPtArr[j+1] - assocPtArr[j]) / 2);
          
  latex = new TLatex(0.65, 0.84, str);
  latex->SetNDC();
  latex->SetTextSize(0.04);
  latex->Draw();
  
  latex = new TLatex(0.65, 0.90, str2);
  latex->SetNDC();
  latex->SetTextSize(0.04);
  latex->Draw();

  latex = new TLatex(0.84, 0.78, "0-5%");
  latex->SetNDC();
  latex->SetTextSize(0.04);
  latex->Draw();
  
  canvas->cd(2);
  gPad->SetLeftMargin(0.20);
  
  hist1->DrawCopy()->GetYaxis()->SetRangeUser(0.18, 0.28);
  DrawFlow(v2[0], hist1, leadingPtArr[i], assocPtArr[j], 0, i, 0, 0, (assocPtArr[j] + assocPtArr[j+1]) / 2, (assocPtArr[j+1] - assocPtArr[j]) / 2);
  
  latex->Draw();
  
  latex = new TLatex(0.5, 0.90, "ALICE preliminary");
  latex->SetTextSize(0.04);
  latex->SetNDC();
  latex->Draw();

  latex = new TLatex(0.5, 0.84, "Statistical uncertainties only");
  latex->SetTextSize(0.04);
  latex->SetNDC();
  latex->Draw();
  
  canvas->cd(3);
  gPad->SetLeftMargin(0.20);
  
  hist2b->GetYaxis()->SetTitleOffset(1.9);
  hist2b->SetYTitle("1/N_{trig} dN/dp_{T,assoc}");
  hist2b->SetTitle("");
  hist2b->DrawCopy();
           
  DrawFlow(v2[2], (TH1*) hist2b->Clone(), leadingPtArr[i], assocPtArr[j], 0, i, 2, 0, (assocPtArr[j] + assocPtArr[j+1]) / 2, (assocPtArr[j+1] - assocPtArr[j]) / 2);
  
  latex = new TLatex(0.84, 0.78, "60-90%");
  latex->SetNDC();
  latex->SetTextSize(0.04);
  latex->Draw();
  
  canvas->cd(4);
  gPad->SetLeftMargin(0.20);
      
  hist2b->DrawCopy()->GetYaxis()->SetRangeUser(0.01, 0.1);
  DrawFlow(v2[2], hist2b, leadingPtArr[i], assocPtArr[j], 0, i, 2, 0, (assocPtArr[j] + assocPtArr[j+1]) / 2, (assocPtArr[j+1] - assocPtArr[j]) / 2);
      
  latex->Draw();
  
  canvas->SaveAs(Form("%s.eps", canvas->GetName()));
  canvas->SaveAs(Form("%s.png", canvas->GetName()));
}

void MCScalingFactor(Int_t caseId)
{
  // 900 GeV
  // x[5]=5, y[5]=0.813743, ex[5]=0, ey[5]=0.158446
  // x[5]=5, y[5]=1.09804, ex[5]=0, ey[5]=0.265434
  // 7 TeV
  // x[5]=5, y[5]=1.02276, ex[5]=0, ey[5]=0.0199991
  // x[5]=5, y[5]=0.870247, ex[5]=0, ey[5]=0.0294444
  
  Float_t factors[] = { 0.814, 1.098, 1.022, 0.870 };
  Float_t errors[]  = { 0.158, 0.265, 0.020, 0.029 };
  
  Float_t avgs[2];
  
  for (Int_t i=0; i<2; i++)
  {
    switch (caseId)
    {
      case 0:
        avgs[i] = factors[2*i] / errors[2*i] / errors[2*i] + factors[2*i+1] / errors[2*i+1] / errors[2*i+1];
        avgs[i] /= 1. / errors[2*i] / errors[2*i] + 1. / errors[2*i+1] / errors[2*i+1];
        break;
        
      case 1:
        avgs[i] = factors[2*i];
        break;
    
      case 2:
        avgs[i] = factors[2*i+1];
        break;
    
      case 3:
        avgs[0] = factors[0];
        avgs[1] = factors[3];
        break;
    
      case 4:
        avgs[0] = factors[1];
        avgs[1] = factors[2];
        break;
    }
  }
  
  Float_t a = avgs[0];
  Float_t b = avgs[1];
  
  Printf("%f %f", a, b);
  
  Float_t linEx = a + (b-a)/(7.-0.9) * (2.76-0.9);
  Float_t logEx = a + (b-a)/(log(7.)-log(0.9)) * (log(2.76)-log(0.9));
  
  Printf("%f %f", linEx, logEx);
}
  
void CompareDeltaPhi(const char* fileName1, const char* fileName2, Int_t centralityID, const char* dist = "_fit_flat", Bool_t reduced = kFALSE)
{
  firstFile = TFile::Open(fileName1);
  secondFile = TFile::Open(fileName2);

  Int_t maxLeadingPt = 2;
  Int_t maxAssocPt = 7;
  Int_t minLeadingPt = 0;
  Int_t minAssocPt = 0;
  
  if (reduced)
  {
    maxLeadingPt = 2;
    maxAssocPt = 5;
    minLeadingPt = 1;
    minAssocPt = 2;
  }
  
  TCanvas* canvas = new TCanvas(Form("%s %s", fileName1, fileName2 ? fileName2 : ""), Form("%s %s", fileName1, fileName2 ? fileName2 : ""), 1000, 700);
  canvas->Divide(maxAssocPt-minAssocPt, maxLeadingPt-minLeadingPt);
  
  for (Int_t i=minLeadingPt; i<maxLeadingPt; i++)
    for (Int_t j=minAssocPt; j<maxAssocPt; j++)
    {
      Printf("%d %d", i, j);
    
      canvas->cd(j+1-minAssocPt + (i-minLeadingPt) * maxAssocPt);
      gPad->SetLeftMargin(0.15);
      gPad->SetBottomMargin(0.2);
      //gPad->SetTopMargin(0.01);
      gPad->SetRightMargin(0.01);
      
      hist1 = (TH1*) firstFile->Get(Form("dphi_%d_%d_%d%s", i, j, centralityID, dist));
      hist2 = (TH1*) secondFile->Get(Form("dphi_%d_%d_%d%s", i, j, centralityID, dist));
      
      if (!hist1)
        continue;
      
/*      hist1->Rebin(2); hist1->Scale(0.5); 
      hist2->Rebin(2); hist2->Scale(0.5); */
//       hist1->Scale(1.6);
        
      hist1->SetLineColor(1);
      hist2->SetLineColor(2);
      
      hist1->DrawCopy()->GetYaxis()->SetRangeUser(hist1->GetMinimum(), hist1->GetMaximum() * 1.2);
      hist2->DrawCopy("SAME");
      
/*      if (strlen(dist) == 0)
        Printf("chi2 test: %f", hist1->Chi2Test(hist2, "UU NORM CHI2/NDF"));*/
      
      //DrawRatio(hist1, hist2, hist1->GetTitle());
      
      //return;
    }
}
  
TH1* MACHCone(const char* fileName, Int_t centrality, const char* drawingOption, Bool_t plotContributions, Int_t twoD, TF1** flowFuncP = 0)
{
  if (0)
  {
    TFile::Open("machcone_input.root");
    
    switch (centrality)
    {
      case 0:
        hist = (TH1*) gFile->Get("dphi_0_0_0");
        break;
      
      case 1:
        hist = (TH1*) gFile->Get("dphi_1_1_0");
        break;
        
      case 2:
        hist = (TH1*) gFile->Get("dphi_0_0_2");
        break;
    
      case 3:
        hist = (TH1*) gFile->Get("dphi_1_1_2");
        break;
    }
  }
  else
  {
    loadlibs();
    
    AliUEHistograms* h = (AliUEHistograms*) GetUEHistogram(fileName);
    AliUEHistograms* hMixed = (AliUEHistograms*) GetUEHistogram(fileName, 0, kTRUE);
    
    Float_t leadingPtArr[] = { 2.0, 3.0, 4.0, 10.0, 20.0, 40.0 };
    Float_t assocPtArr[] =   { 1.0, 2.0, 3.0, 6.0, 10.0, 20.0, 40.0 };
//     Float_t assocPtArr[] =   { 2.0, 3.0, 6.0, 10.0, 20.0, 40.0 };
  
    Int_t i = centrality % 2;
    Int_t step = 6;
    
    Int_t j = 0;
    
    gpTMin = assocPtArr[i] + 0.01;
    gpTMax = assocPtArr[i+1] - 0.01;
    
    Int_t centralityBegin = 0;
    Int_t centralityEnd = 1;
    
    if (centrality >= 2)
    {
      centralityBegin = 30;
      centralityEnd = 40;
    }
      
    SetupRanges(h);
    SetupRanges(hMixed);
    
    TH1* hist = 0;
    Float_t vn[5];
    vn[0] = 0;
    vn[1] = 0;
    vn[2] = 0;
    vn[3] = 0;
    vn[4] = 0;
    
    Bool_t scaleToPairs = 0;
    
    GetDistAndFlow(h, hMixed, &hist, 0, step, centralityBegin, centralityEnd, leadingPtArr[j] + 0.01, leadingPtArr[j+1] - 0.01, twoD, kTRUE, vn, scaleToPairs); 
    if (scaleToPairs && twoD == 11)
      hist->Scale(2);  
    
    Printf("%f %f %f", vn[2], vn[3], vn[4]);
    
    // mirror delta phi to improve stats
    if (1)
    {
      for (Int_t bin=1; bin<=hist->GetNbinsX(); bin++)
      {
	if (hist->GetBinCenter(bin) < 0 || hist->GetBinCenter(bin) > TMath::Pi())
	{
	  if (hist->GetBinCenter(bin) < 0)
	    Int_t bin2 = hist->FindBin(-1.0 * hist->GetBinCenter(bin));
	  else
	    Int_t bin2 = hist->FindBin(TMath::TwoPi() - hist->GetBinCenter(bin));
	  
	  Float_t combValue = hist->GetBinContent(bin) / hist->GetBinError(bin) / hist->GetBinError(bin) + hist->GetBinContent(bin2) / hist->GetBinError(bin2) / hist->GetBinError(bin2);
	  
	  Float_t weight = 1. / hist->GetBinError(bin) / hist->GetBinError(bin) + 1. / hist->GetBinError(bin2) / hist->GetBinError(bin2);
	  combValue /= weight;
	  
	  Float_t combError = TMath::Sqrt(1.0 / weight);
	  
	  hist->SetBinContent(bin2, combValue);
	  hist->SetBinError(bin2, combError);

	  hist->SetBinContent(bin, combValue);
	  hist->SetBinError(bin, combError);
	  
	  Printf("%d %d %f %f", bin, bin2, combValue, combError);
	}
      }
    }
  }
  
  TString str(hist->GetTitle());
  str.ReplaceAll(" - ", "#");
  tokens = str.Tokenize("#");
  hist->SetTitle("");
  
  hist->SetLineColor(1);
  hist->SetYTitle("1/N_{trig} dN_{assoc}/d#Delta#phi");
  hist->SetXTitle("#Delta#phi (rad)");
  hist = (TH1*) hist->Clone();
  if (!drawingOption)
    return hist;
  
  hist = hist->DrawCopy(drawingOption);
  
  for (Int_t i=0; i<tokens->GetEntries()-1; i++)
  {
    if (centrality == 1)
      latex = new TLatex(0.6, 0.9 - i*0.05 - gPad->GetTopMargin(), tokens->At(i)->GetName());
    else
      latex = new TLatex(0.6, 0.12 - i*0.05, tokens->At(i)->GetName());
    
    latex->SetNDC();
    latex->SetTextSize(0.04);
    latex->Draw();
  }
  
  TF1* flowFunc = new TF1("flowFunc", "[0] * (1+2*[1]*cos(2*x)+2*[2]*cos(3*x)+2*[3]*cos(4*x)+2*[4]*cos(5*x))", -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
  //flowFunc->SetLineWidth(1);
  flowFunc->SetLineColor(2);
  flowFunc->SetLineStyle(2);
  for (Int_t i=0; i<4; i++)
    flowFunc->FixParameter(i+1, vn[i+1]);
  //flowFunc->SetParameter(0, 84.95); flowFunc->DrawClone("SAME");
  
  //flowFunc->SetParameter(0, 85.4); flowFunc->DrawClone("SAME");
  //flowFunc->SetParameter(0, 85.6); flowFunc->DrawClone("SAME");
  
//   hist->Fit(flowFunc, "", "SAME", 1.8, 4.5);
  hist->Fit(flowFunc, "N", "", -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
  hist->GetYaxis()->SetRangeUser(hist->GetMinimum() * 0.95, hist->GetMaximum() * 1.05);
  flowFunc->SetRange(-0.5 * TMath::Pi(), 1.5 * TMath::Pi());
  flowFunc->Draw("SAME");
  
  if (flowFuncP)
    *flowFuncP = flowFunc;
  
  //hist->GetYaxis()->SetRangeUser(84, 89);
  //flowFunc->SetParameter(4, 1);
  
  Float_t chi2 = 0;
  Int_t n=0;
  for (Int_t bin=hist->FindBin(1.8); bin<=hist->FindBin(4.5); bin++)
  {
    chi2 += TMath::Power((hist->GetBinContent(bin) - flowFunc->Eval(hist->GetBinCenter(bin))) / hist->GetBinError(bin), 2);
    n++;
    //Printf("%f", TMath::Power((hist->GetBinContent(bin) - flowFunc->Eval(hist->GetBinCenter(bin))) / hist->GetBinError(bin), 2));
  }
    
  chi2 /= n;
  
  Printf("chi2/ndf = %f", chi2);
  
  if (!plotContributions)
    return hist;
  
  for (Int_t n=1; n<=5; n++)
  {
    flowFuncPart = new TF1("flowFuncPart", "[0] * (1+2*[1]*cos([2]*x))", -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
    flowFuncPart->SetParameters(flowFunc->GetParameter(0), vn[n-1], n);
    flowFuncPart->SetLineWidth(1);
    flowFuncPart->SetLineStyle(n);
    flowFuncPart->Draw("SAME");
  }
  
  return hist;
}

void MACHConeAll(const char* fileName)
{
  loadlibs();

  SetFlowStyle();

  c = new TCanvas("c", "c", 300, 800);
  c->SetTopMargin(0);
  c->SetLeftMargin(0);
  c->SetRightMargin(0);
  c->SetBottomMargin(0);
  c->Divide(1, 2, 0, 0);
  
  Float_t rangesMin[] = { 136.1, 17.45, 18, 2.4 };
  Float_t rangesMax[] = { 139.4, 18.95, 23.5, 3.6 };
  
  for (Int_t i=0; i<2; i++)
  {
    c->cd(i+1);

    gPad->SetRightMargin(0.05);
    gPad->SetLeftMargin(0.20);
    gPad->SetTopMargin(0.15);
    gPad->SetBottomMargin(0.15);

    if (i == 0)
    {
      gPad->SetBottomMargin(0.01);
    }
    else if (i == 1)
    {
      gPad->SetTopMargin(0.01);
      gPad->SetBottomMargin(0.7);
    }
  }

  Int_t i = 0;

  TF1* flowFunc = 0;
  
  c->cd(i+1);
//   hist = MACHCone(fileName, i, "", kTRUE, 0, &flowFunc);
  hist = (TH1*) MACHCone(fileName, i, 0, kFALSE, 0)->Clone("hist");
  hist->GetYaxis()->SetRangeUser(rangesMin[i], rangesMax[i]);
  hist->GetYaxis()->SetTitleOffset(1.7);
  hist->Draw();
//     continue;
  
  c->cd(i+1);
  //c->cd(i+2);
//   hist2 = MACHCone(fileName, i, 0, kFALSE, 11);
  hist2 = MACHCone(fileName, i, "SAME", kTRUE, 11, &flowFunc);
  //hist2->GetYaxis()->SetRangeUser(rangesMin[i], rangesMax[i]);
//   hist2->Scale(2);
  hist2->SetLineColor(4);
  //hist2->Draw("SAME");
//    return;

  // sample flowFunc with same binning
  flowFuncHist = (TH1*) hist->Clone("flowFuncHist");
  flowFuncHist->Reset();
  
  for (Int_t i=1; i<=flowFuncHist->GetNbinsX(); i++)
  {
    flowFuncHist->SetBinContent(i, flowFunc->Integral(flowFuncHist->GetXaxis()->GetBinLowEdge(i), flowFuncHist->GetXaxis()->GetBinUpEdge(i)));
    flowFuncHist->SetBinError(i, 0);
  }
  
  i = 0;
  
  flowFuncHist->Scale(1.0 / flowFuncHist->GetBinWidth(1));
  
//   new TCanvas;
//   flowFuncHist->Draw();
//   flowFunc->Draw("SAME");
//   return;

  // residuals
  residuals = (TH1*) hist->Clone("res1");
  residuals->SetYTitle("Residuals");
  residuals->Divide(flowFuncHist);
  c->cd(i+2);
  residuals->Draw();
  residuals->GetYaxis()->SetNdivisions(505);
  residuals->GetYaxis()->SetRangeUser(0.995, 1.005);

  residuals = (TH1*) hist2->Clone("res2");
  residuals->Divide(flowFuncHist);
  c->cd(i+2);
  residuals->SetLineColor(4);
  residuals->Draw("SAME");
  
  gPad->SetGridy();

  if (i == 0)
  {
    legend = new TLegend(0.43, 0.66, 0.91, 0.82);
    legend->SetFillColor(0);
    legend->SetTextSize(0.034);
    legend->SetTextAlign(22);
    
    legend->SetHeader("Centrality 0-1%");
    legend->AddEntry(hist, "|#eta| < 0.8 & All #Delta#eta", "L");
    legend->AddEntry(hist2, "|#eta| < 0.8 & |#Delta#eta| > 0.8 (2x)", "L");
    
    c->cd(i+1);
    legend->Draw();
  }
  
  c->SaveAs(Form("machcone_%d.png", 0));
  c->SaveAs(Form("machcone_%d.eps", 0));
  c->SaveAs(Form("machcone_%d.C", 0));
  
  return;
  
  c->cd(1);
  MACHCone(3, "", kTRUE)->SetTitle("");
  
  c->cd(2);
  hist1 = MACHCone(0, "", kFALSE);
  hist2 = MACHCone(1, "SAME", kFALSE);
  hist3 = MACHCone(2, "SAME", kFALSE);
  
  hist1->SetTitle("");
  hist1->GetYaxis()->SetRangeUser(1, 7);
}
  
void SetFlowStyle()
{
 // Set style which will affect all plots.
 
 gStyle->Reset();
 // gStyle->SetOptitle(0);
 // gStyle->SetOptStat(0);
 //gStyle->SetOptDate(1);
 // gStyle->SetPalette(8,0);  // (1,0)
 gStyle->SetPalette(1);  // (1,0)
 gStyle->SetDrawBorder(0);
 // gStyle->SetFillColor(0);  // kills palete ???
 gStyle->SetCanvasColor(0);
 gStyle->SetPadColor(0);
 // gStyle->SetFillColor(0); // otherwize it affects Fill colors later
 gStyle->SetFrameFillColor(0);
 gStyle->SetCanvasBorderMode(0);
 gStyle->SetFrameLineWidth(2);
 // gStyle->SetFrameFillStyle(4000);
 gStyle->SetPadBorderMode(0);
 gStyle->SetPadTickX(1);
 gStyle->SetPadTickY(1);
 gStyle->SetPadBottomMargin(0.15);
 gStyle->SetPadLeftMargin(0.15);
 gStyle->SetHistLineWidth(2);
 gStyle->SetFuncWidth(2);
 gStyle->SetLineWidth(2);
 gStyle->SetLabelSize(0.05,"xyz");
 gStyle->SetLabelOffset(0.01,"y");
 gStyle->SetLabelColor(kBlack,"xyz");
 gStyle->SetTitleSize(0.06,"xyz");
 gStyle->SetTitleOffset(1.3,"y");
 gStyle->SetTitleFillColor(0);
 gStyle->SetLineWidth(2);  
 gStyle->SetHistLineColor(1);
 gStyle->SetTextColor(1);
 gStyle->SetTitleTextColor(1);
 TGaxis::SetMaxDigits(4);
 gStyle->SetOptStat(0); // removes stat. box from all histos
 gROOT->ForceStyle();

} // end of void SetFlowStyle()

void PlotPtDistributionsSubtracted(const char* fileName1, Int_t centrBegin = 1, Int_t centrEnd = 6, Int_t fit = 0)
{
  loadlibs();

  Int_t nCentralityBins = 5;
  Int_t centralityBins[] = { 1, 7, 9, 11, 13, 16 };
  
  AliUEHistograms* h = (AliUEHistograms*) GetUEHistogram(fileName1);
    
  TCanvas* canvas = new TCanvas("Pt", "Pt", 1000, 1000);
  canvas->Divide(2, 2);
  
  TLegend* legend = new TLegend(0.3, 0.8, 0.99, 0.99);
  legend->SetFillColor(0);
  
  Int_t colors[] = { 1, 2, 4, 6, 3, 5 };
  Int_t markers[] = { 20, 21, 22, 23, 24, 25 };
  
  Double_t ptMin = 10.01;
  Double_t ptMax = 39.99;
    
  TString str;
  str.Form("%.1f < p_{T,trig} < %.1f", ptMin - 0.01, ptMax + 0.01);
  
  canvas->cd(1);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.2);
  gPad->SetTopMargin(0.01);
  gPad->SetRightMargin(0.01);
  gPad->SetLogy();
  gPad->SetLogx();
  gPad->SetGridx();
  gPad->SetGridy();
  
  TH2F* dummy = new TH2F("dummy", "", 100, 1, 40, 100, 1e-5, 1e3);
  dummy->SetStats(kFALSE);
  dummy->SetXTitle("p_{T,assoc}");
  dummy->SetYTitle("1/N_{trig} dN/dp_{T,assoc}");
  dummy->GetYaxis()->SetTitleOffset(1);
  Prepare1DPlot(dummy);

  dummy->GetXaxis()->SetLabelSize(0.05);
  dummy->GetYaxis()->SetLabelSize(0.05);
  dummy->GetXaxis()->SetTitleSize(0.05);
  dummy->GetYaxis()->SetTitleSize(0.05);
  dummy->DrawCopy();
  
  canvas->cd(2);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.2);
  gPad->SetTopMargin(0.01);
  gPad->SetRightMargin(0.01);
  gPad->SetLogy();
  gPad->SetLogx();
  gPad->SetGridx();
  gPad->SetGridy();
  dummy->DrawCopy();
  
  canvas->cd(3);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.2);
  gPad->SetTopMargin(0.01);
  gPad->SetRightMargin(0.01);
  gPad->SetLogx();
  gPad->SetGridx();
  gPad->SetGridy();
  
  dummy = new TH2F("dummy2", "", 100, 1, 40, 100, 0, 100);
  dummy->SetStats(kFALSE);
  dummy->SetXTitle("p_{T,assoc}");
  dummy->SetYTitle("Ratio: Distribution / Inclusive");
  dummy->GetYaxis()->SetTitleOffset(1);
  Prepare1DPlot(dummy);
  dummy->GetXaxis()->SetLabelSize(0.05);
  dummy->GetYaxis()->SetLabelSize(0.05);
  dummy->GetXaxis()->SetTitleSize(0.05);
  dummy->GetYaxis()->SetTitleSize(0.05);
  dummy->DrawCopy();
  
  canvas->cd(4);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.2);
  gPad->SetTopMargin(0.01);
  gPad->SetRightMargin(0.01);
  gPad->SetLogx();
  gPad->SetGridx();
  gPad->SetGridy();
  dummy->DrawCopy();
  
  
  // see GetAcceptanceScalingFactor
  Float_t scalingFactor = 0.328396; // 0.8
  //Float_t scalingFactor = 0.433268; // 0.7
  //Float_t scalingFactor = 0.895593; // 0.5 --> 1.6
      
  Float_t phiRange[] = { 0, TMath::Pi() / 2, TMath::Pi() };
  const char* phiLabels[] = { "Towards", "Transverse", "Away" };
  
  //Float_t phiSizeTowards = TMath::Pi() / 3;
  Float_t phiSizeTowards = 0.75;
  //Float_t etaLimit = 0.5;
  Float_t etaLimit = 0.8;
  
  towardsCentralEta = h->GetUEHist(2)->GetPtHist(6, 0, ptMin, ptMax, centrBegin, centrEnd, phiRange[0] - phiSizeTowards, phiRange[0] + phiSizeTowards, -etaLimit + 0.01, etaLimit - 0.01);
  towardsSideEta1 = h->GetUEHist(2)->GetPtHist(6, 0, ptMin, ptMax, centrBegin, centrEnd, phiRange[0] - phiSizeTowards, phiRange[0] + phiSizeTowards, -1.59, -etaLimit - 0.01);
  TH1* towardsSideEta2 = h->GetUEHist(2)->GetPtHist(6, 0, ptMin, ptMax, centrBegin, centrEnd, phiRange[0] - phiSizeTowards, phiRange[0] + phiSizeTowards, etaLimit + 0.01, 1.59);
  towardsSideEta1->Add(towardsSideEta2); // TODO can be done smarter? what about the errors?

  Prepare1DPlot(towardsSideEta1);
  Prepare1DPlot(towardsCentralEta);
  
  towardsCentralEta->SetLineColor(colors[0]);
  towardsSideEta1->SetLineColor(colors[3]);
  //towardsSideEta1->SetLineStyle(2);

  legend->AddEntry(towardsCentralEta->Clone(), Form("Jet, %s: |#Delta#eta| < 0.8, #phi ~ %.1f", phiLabels[0], phiRange[0]));
  legend->AddEntry(towardsSideEta1->Clone(), Form("Ridge, %s: 0.8 < |#Delta#eta| < 1.6, #phi ~ %.1f", phiLabels[0], phiRange[0]));
  
  // TODO update when data with new phi binning is available
  transverse = h->GetUEHist(2)->GetPtHist(6, 0, ptMin, ptMax, centrBegin, centrEnd, phiRange[1] - TMath::Pi() / 3, phiRange[1] + TMath::Pi() / 3, -1.59, 1.59);
  /*
  transverse2 = h->GetUEHist(2)->GetPtHist(6, 0, ptMin, ptMax, centrBegin, centrEnd, -0.5 * TMath::Pi(), -TMath::Pi() / 3, -1.59, 1.59, kTRUE);
  transverse3 = h->GetUEHist(2)->GetPtHist(6, 0, ptMin, ptMax, centrBegin, centrEnd, TMath::Pi() + TMath::Pi() / 3, 1.5 * TMath::Pi(), -1.59, 1.59, kTRUE);
  transverse->Add(transverse2);
  transverse->Add(transverse3);
  transverse->Scale(1.0 / (2.0 / 3 * TMath::Pi()));
  */
    
  //transverseSideEta1 = h->GetUEHist(2)->GetPtHist(6, 0, ptMin, ptMax, centrBegin, centrEnd, phiRange[1] - TMath::Pi() / 4, phiRange[1] + TMath::Pi() / 4, -1.59, -0.81);
  //TH1* transverseSideEta2 = h->GetUEHist(2)->GetPtHist(6, 0, ptMin, ptMax, centrBegin, centrEnd, phiRange[1] - TMath::Pi() / 4, phiRange[1] + TMath::Pi() / 4, 0.81, 1.59);
  //transverseSideEta1->Add(transverseSideEta2); // TODO can be done smarter? what about the errors?
  
  transverse->SetLineColor(colors[2]);
  //transverseSideEta1->SetLineColor(colors[2]);
  //transverseSideEta1->SetLineStyle(2);

  Prepare1DPlot(transverse);
  //Prepare1DPlot(transverseSideEta1);
  legend->AddEntry(transverse->Clone(), Form("%s: |#Delta#eta| < 1.6, #phi ~ %.1f", phiLabels[1], phiRange[1]));
  //legend->AddEntry(transverseSideEta1->Clone(), Form("%s: 0.8 < |#Delta #eta| < 1.6, #phi ~ %.1f", phiLabels[1], phiRange[1]));
  
  away = h->GetUEHist(2)->GetPtHist(6, 0, ptMin, ptMax, centrBegin, centrEnd, phiRange[2] - TMath::Pi() / 3, phiRange[2] + TMath::Pi() / 3, -1.59, 1.59);

  /*
  away2 = h->GetUEHist(2)->GetPtHist(6, 0, ptMin, ptMax, centrBegin, centrEnd, phiRange[2] - TMath::Pi() / 4, phiRange[2] + TMath::Pi() / 4, -0.79, 0.79);
  away3 = h->GetUEHist(2)->GetPtHist(6, 0, ptMin, ptMax, centrBegin, centrEnd, phiRange[2] - TMath::Pi() / 4, phiRange[2] + TMath::Pi() / 4, -1.59, -0.81);
  away4 = h->GetUEHist(2)->GetPtHist(6, 0, ptMin, ptMax, centrBegin, centrEnd, phiRange[2] - TMath::Pi() / 4, phiRange[2] + TMath::Pi() / 4,  0.81, 1.59);
  away3->Add(away4);
  */
  
  away->SetLineColor(colors[1]);
  //away3->SetLineColor(colors[1]);
  //away3->SetLineStyle(2);
  
  Prepare1DPlot(away);
  legend->AddEntry(away->Clone(), Form("%s: |#Delta#eta| < 1.6, #phi ~ %.1f", phiLabels[2], phiRange[2]));
  
  //inclusive = h->GetUEHist(2)->GetPtHist(6, 0, ptMin, ptMax, centrBegin, centrEnd, -0.5 * TMath::Pi(), 1.5 * TMath::Pi(), -1.59, 1.59);
  
  /*
  h->GetUEHist(2)->GetEventHist()->GetGrid(6)->GetGrid()->GetAxis(1)->SetRange(centrBegin, centrEnd);
  inclusive = h->GetUEHist(2)->GetEventHist()->ShowProjection(0, 6);
  h->GetUEHist(2)->GetEventHist()->GetGrid(6)->GetGrid()->GetAxis(1)->SetRange(0, -1);
  */
  
  axis = h->GetUEHist(2)->GetEventHist()->GetGrid(6)->GetGrid()->GetAxis(1);
  
  inclusiveTmp = h->GetCorrelationpT()->ProjectionY("inclusiveTmp", h->GetCorrelationpT()->GetXaxis()->FindBin(axis->GetBinLowEdge(centrBegin) + 0.01), h->GetCorrelationpT()->GetXaxis()->FindBin(axis->GetBinUpEdge(centrEnd) - 0.01));
  
  //rebin to match other histograms
  //inclusive2 = inclusive->Rebin(away->GetNbinsX(), "inclusive2", away->GetXaxis()->GetXbins()->GetArray());
  
  // manually, as usual the ROOT function makes a SEGV
  inclusive = (TH1*) away->Clone("inclusive");
  inclusive->Reset();
  inclusive->Sumw2();
  for (Int_t bin=1; bin<=inclusiveTmp->GetNbinsX(); bin++)
    inclusive->Fill(inclusiveTmp->GetBinCenter(bin), inclusiveTmp->GetBinContent(bin));
  
  for (Int_t bin=1; bin<=inclusive->GetNbinsX(); bin++)
    inclusive->SetBinError(bin, TMath::Sqrt(inclusive->GetBinContent(bin)));
    
  // normalization: events, phase space
  inclusive->Scale(1.0 / h->GetCentralityDistribution()->Integral(centrBegin, centrEnd));
  inclusive->Scale(1.0 / (TMath::TwoPi() * 0.8 * 2));
  
  // bin width
  for (Int_t i=1; i<=inclusive->GetNbinsX(); i++)
  {
    inclusive->SetBinContent(i, inclusive->GetBinContent(i) / inclusive->GetXaxis()->GetBinWidth(i));
    inclusive->SetBinError  (i, inclusive->GetBinError(i) / inclusive->GetXaxis()->GetBinWidth(i));
  }
  
  Prepare1DPlot(inclusive);
  inclusive->SetLineColor(colors[4]);
  legend->AddEntry(inclusive->Clone(), "Inclusive");
  
  //Prepare1DPlot(away3);
  //legend->AddEntry(away3->Clone(), Form("%s: 0.8 < |#Delta #eta| < 1.6, #phi ~ %.1f", phiLabels[2], phiRange[2]));
  
  // scale for acceptance to match acceptance of towardsCentralEta
  towardsSideEta1->Scale(1.0 / scalingFactor); // contains only "side-eta"
  transverse->Scale(1.0 / (1.0 + scalingFactor)); // full eta
  away->Scale(1.0 / (1.0 + scalingFactor)); // full eta
  inclusive->Scale(1.0 / (1.0 + scalingFactor)); // full eta (phi scaling implicit)
  
  canvas->cd(1);
  towardsCentralEtaClone = towardsCentralEta->DrawCopy("SAME");
  towardsSideEta1Clone = towardsSideEta1->DrawCopy("SAME");
  transverseClone = transverse->DrawCopy("SAME");
  //transverseSideEta1->DrawCopy("SAME");
  awayClone = away->DrawCopy("SAME");
  inclusiveClone = inclusive->DrawCopy("SAME");
  //inclusive2->DrawCopy("SAME");
  //away3->DrawCopy("SAME");
  
  canvas->cd(3);
  
  towardsCentralEtaR = (TH1*) towardsCentralEta->Clone("towardsCentralEtaR");
  towardsSideEta1R = (TH1*) towardsSideEta1->Clone("towardsSideEta1R");
  transverseR = (TH1*) transverse->Clone("transverseR");
  awayR = (TH1*) away->Clone("awayR");
  
  towardsCentralEtaR->Divide(inclusive);
  towardsSideEta1R->Divide(inclusive);
  transverseR->Divide(inclusive);
  awayR->Divide(inclusive);
  
  towardsCentralEtaR->DrawCopy("SAME");
  towardsSideEta1R->DrawCopy("SAME");
  transverseR->DrawCopy("SAME");
  awayR->DrawCopy("SAME");
  
  if (fit == 1)
  {
    canvas->cd(1);
    
    awayClone->Fit("expo", "", "SAME", 1.1, 2.9);
    awayClone->GetFunction("expo")->SetRange(1, 10);
    gPad->SetLogx(0);
  }
  else if (fit == 2)
  {
    canvas->cd(1);
    
    func = new TF1("func", "[0] * x**[1]");
    func->SetParameters(1, 1);
    
    Float_t limitLow = 2.1;
    Float_t limitHigh = 4.9;
    
    func2 = (TF1*) func->Clone();
    towardsCentralEtaClone->Fit(func2, "", "SAME", limitLow, limitHigh);
    towardsCentralEtaClone->GetFunction("func")->SetRange(1, 10);
    towardsCentralEtaClone->GetFunction("func")->SetLineColor(towardsCentralEtaClone->GetLineColor());
    
    func2 = (TF1*) func->Clone();
    towardsSideEta1Clone->Fit(func2, "", "SAME", limitLow, limitHigh);
    towardsSideEta1Clone->GetFunction("func")->SetRange(1, 10);
    towardsSideEta1Clone->GetFunction("func")->SetLineColor(towardsSideEta1Clone->GetLineColor());
  
    func2 = (TF1*) func->Clone();
    transverseClone->Fit(func2, "", "SAME", limitLow, limitHigh);
    transverseClone->GetFunction("func")->SetRange(1, 10);
    transverseClone->GetFunction("func")->SetLineColor(transverseClone->GetLineColor());
  
    func2 = (TF1*) func->Clone();
    awayClone->Fit(func2, "", "SAME", limitLow, limitHigh);
    awayClone->GetFunction("func")->SetRange(1, 10);
    awayClone->GetFunction("func")->SetLineColor(awayClone->GetLineColor());
  
    func2 = (TF1*) func->Clone();
    inclusiveClone->Fit(func2, "", "SAME", limitLow, limitHigh);
    inclusiveClone->GetFunction("func")->SetRange(1, 10);
    inclusiveClone->GetFunction("func")->SetLineColor(inclusiveClone->GetLineColor());
  }
  
  // subtract transverse part
  away->Add(transverse, -1);
  //away3->Add(transverseSideEta1, -1);
  towardsCentralEta->Add(transverse, -1);
  towardsSideEta1->Add(transverse, -1);
  //inclusive->Add(transverse, -1);
  
  canvas->cd(2);
  towardsCentralEta->DrawCopy("SAME");
  towardsSideEta1->DrawCopy("SAME");
  away->DrawCopy("SAME");
  inclusive->DrawCopy("SAME");
  //away3->DrawCopy("SAME");
  
  canvas->cd(4);
  
  towardsCentralEtaR = (TH1*) towardsCentralEta->Clone("towardsCentralEtaR2");
  towardsSideEta1R = (TH1*) towardsSideEta1->Clone("towardsSideEta1R2");
  awayR = (TH1*) away->Clone("awayR2");
  
  towardsCentralEtaR->Divide(inclusive);
  towardsSideEta1R->Divide(inclusive);
  awayR->Divide(inclusive);
  
  towardsCentralEtaR->DrawCopy("SAME");
  towardsSideEta1R->DrawCopy("SAME");
  awayR->DrawCopy("SAME");
  
  if (fit == 1)
  {
    canvas->cd(2);
    
    away->Fit("expo", "", "SAME", 1.1, 2.9);
    away->GetFunction("expo")->SetRange(1, 10);
    gPad->SetLogx(0);
  }
  else if (fit == 2)
  {
    canvas->cd(2);
    
    func = new TF1("func", "[0] * x**[1]");
    func->SetParameters(1, 1);
    
    func2 = (TF1*) func->Clone();
    towardsCentralEta->Fit(func2, "", "SAME", limitLow, limitHigh);
    func2->SetRange(1, 10);
    
    func2 = (TF1*) func->Clone();
    towardsSideEta1->Fit(func2, "", "SAME", limitLow, limitHigh);
    func2->SetRange(1, 10);
  
    func2 = (TF1*) func->Clone();
    away->Fit(func2, "", "SAME", limitLow, limitHigh);
    func2->SetRange(1, 10);
  
    func2 = (TF1*) func->Clone();
    inclusive->Fit(func2, "", "SAME", limitLow, limitHigh);
    func2->SetRange(1, 10);
  }
  
  for (Int_t i=1; i<=2; i++)
  {
    canvas->cd(i);
    
    TString str3;
    str3.Form("%d-%d%%", (Int_t) h->GetCentralityDistribution()->GetXaxis()->GetBinLowEdge(centrBegin), (Int_t) h->GetCentralityDistribution()->GetXaxis()->GetBinUpEdge(centrEnd));
    latex = new TLatex(0.2, 0.3, str3);
    latex->SetNDC();
    latex->SetTextSize(0.06);
    
    latex2 = new TLatex(0.55, 0.6, str);
    latex2->SetNDC();
    latex2->SetTextSize(0.06);
    
    latex->Draw();
    latex2->Draw();
    
    legend->Draw();
  }
}

void CombineDeltaPhiWithWeighting()
{
  // From Hermes:
  // For LHC11a2a:
  // xsection: 11.879829,  ntrials: 8132994.000000
  // 
  // For LHC11a2b:
  // xsection: 0.623421,  ntrials: 2293420.000000
  // 
  // For LHC11a2c:
  // xsection: 0.043815,  ntrials: 1314525.375000
  
  // TODO is merging only same event distribution at present
  
  loadlibs();
  
  Int_t nInput = 3;
  const char* inputList[] = { "mergejob/LHC11a2a_110124.root", "mergejob/LHC11a2b_110131.root", "mergejob/LHC11a2c_110131.root" };
  
  Float_t xSections[] = { 11.879829 * 1e6 / 8132994, 0.623421 * 1e6 / 2293420, 0.043815 * 1e6 / 1314525 };
  xSections[2] /= xSections[0];
  xSections[1] /= xSections[0];
  xSections[0] /= xSections[0];
  
  AliUEHistograms* files[3];
  TList* finalList = 0;
  for (Int_t i=0; i<nInput; i++)
  {
    files[i] = (AliUEHistograms*) GetUEHistogram(inputList[i], (i == 0) ? &finalList : 0);
    if (i > 0) 
    {
      files[i]->Scale(xSections[i]);
      
      TList* list2 = new TList;
      list2->Add(files[i]);
      files[0]->Merge(list2);
    }
  }
  
  TFile* file3 = TFile::Open("out.root", "RECREATE");
  file3->mkdir("PWG4_PhiCorrelations");
  file3->cd("PWG4_PhiCorrelations");
  finalList->Write(0, TObject::kSingleKey);
  file3->Close();       
}

void NormalizeToBinWidth(TH1* hist)
{
  //
  // normalizes a 1-d histogram to its bin width
  //

  for (Int_t i=1; i<=hist->GetNbinsX(); ++i)
  {
    hist->SetBinContent(i, hist->GetBinContent(i) / hist->GetBinWidth(i));
    hist->SetBinError(i, hist->GetBinError(i) / hist->GetBinWidth(i));
  }
}

void GetAcceptanceScalingFactor(const char* fileName1, Float_t eta1 = 0.8, Float_t eta2 = 1.6)
{
  loadlibs();

  AliUEHistograms* h = (AliUEHistograms*) GetUEHistogram(fileName1, 0, kTRUE);
  h->SetPtRange(1.01, 39.99);
 
  new TCanvas;
  mixed = ((TH2*)h->GetUEHist(2)->GetUEHist(6, 0, 1.01, 39.99, 0, 16, kTRUE))->ProjectionY();
  mixed->DrawCopy();
  
  Float_t left = mixed->Integral(mixed->FindBin(-eta2 + 0.01), mixed->FindBin(-eta1 - 0.01));
  Float_t center = mixed->Integral(mixed->FindBin(-eta1 + 0.01), mixed->FindBin(eta1 - 0.01));
  Float_t right = mixed->Integral(mixed->FindBin(eta1 + 0.01), mixed->FindBin(eta2 - 0.01));
  
  Printf("%f %f %f", left, center, right);
  Printf("%f", (left + right) / center);
}

void TrackingEfficiencyCentralityDependence(const char* fileName, Int_t step1 = 2, Int_t step2 = 4)
{
  Int_t nCentralityBins = 3;
  Float_t centralityBins[] = { 0, 20, 40, 90 };
  
  Int_t colors[] = { 1, 2, 4, 6 };

  loadlibs();

  AliUEHistograms* h = (AliUEHistograms*) GetUEHistogram(fileName);
  
  new TCanvas;
  dummy = new TH2F("dummy", ";p_{T};tracking efficiency", 100, 0, 20, 100, 0.7, 1.3);
  dummy->SetStats(0);
  dummy->Draw();
  
  legend = new TLegend(0.13, 0.67, 0.31, 0.87);
  legend->SetTextSize(0.04);
  legend->SetFillColor(0);
  
  for (Int_t i=0; i<nCentralityBins; i++)
  {
    eventHist = h->GetUEHist(2)->GetEventHist();
    
    eventHist->GetGrid(step1)->SetRangeUser(1, centralityBins[i] + 0.1, centralityBins[i+1] - 0.1);
    eventHist->GetGrid(step2)->SetRangeUser(1, centralityBins[i] + 0.1, centralityBins[i+1] - 0.1);
  
    TH1* hist1 = eventHist->ShowProjection(0, step1);
    TH1* hist2 = eventHist->ShowProjection(0, step2);
    
    hist2->Divide(hist1);
    
    hist2->SetLineColor(colors[i]);
    hist2->GetYaxis()->SetTitle("tracking efficiency");
    hist2->Draw("SAME");
    
    legend->AddEntry(hist2, Form("%.0f-%.0f%%", centralityBins[i], centralityBins[i+1]));
  }
  
  legend->Draw();
}


TGraphErrors* ReadHepdata(const char* fileName, Bool_t errorsAreAdded = kFALSE, Int_t skipYErrors = 0, Int_t skipXerrors = 1)
{
  // expected format: x [x2] y [ye] [ye2] [xe]
  //
  // skipYErrors:   0 --> ye present
  //                1 --> no errors ye
  //                2 --> y and ye are lower and upper error, i.e. y' = (y + ye) / 2 and ye = (ye - y) / 2
  //                3 --> ye and ye2 are stat and syst error, will be added in quadrature
  // 
  // skipXerrors:   0 --> xe present
  //                1 --> no errors xe
  //                2 --> x2 present, xe not present and is calculated from x2 - x
  
  ifstream fin(fileName);

  graph = new TGraphErrors(0);

  Double_t sum = 0;

  while (fin.good())
  {
    char buffer[2000];
    if (fin.peek() == '#')
    {
      fin.getline(buffer, 2000);
      continue;
    }
  
    Double_t x = -1;
    Double_t x2 = -1;
    Double_t y = -1;
    Double_t ye = 0;
    Double_t xe = 0;

    fin >> x;
    
    if (skipXerrors == 2)
    {
      fin >> x2;
      xe = (x2 - x + 1) / 2;
      x = x + (x2 - x) / 2;
    }
    
    fin >> y;

    if (y == -1)
      continue;

    if (skipYErrors == 0)
    {
      ye = -1;
      fin >> ye;
      if (ye == -1)
        continue;
    }
    else if (skipYErrors == 2)
    {
      ye = -1;
      fin >> ye;
      if (ye == -1)
        continue;
      
      Double_t newy = (y + ye) / 2;
      ye = (ye - y) / 2;
      y = newy;
    }
    else if (skipYErrors == 3)
    {
      ye = -1;
      fin >> ye;
      if (ye == -1)
        continue;
      
      Double_t ye2 = -1;
      fin >> ye2;
      if (ye2 == -1)
        continue;

      ye = TMath::Sqrt(ye*ye + ye2*ye2);
    }

    if (skipXerrors == 0)
    {
      xe = -1;
      fin >> xe;
      if (xe == -1)
        continue;
    }

    //Printf("%f %f %f %f", x, y, xe, ye);

    if (errorsAreAdded)
      ye -= y;

    graph->SetPoint(graph->GetN(), x, y);
    graph->SetPointError(graph->GetN()-1, xe, ye);

    sum += y;
    
    // read rest until end of line...
    fin.getline(buffer, 2000);
  }
  fin.close();

  Printf("%s: %f", fileName, sum);

  return graph;
}

void EvaluateParticleEfficiency(const char* fileName)
{
  Int_t centralityBegin = 1;
  Int_t centralityEnd = 15;
  
  if (1)
  {
    Int_t step1 = 4;
    Int_t step2 = 5;
  }
  else
  {
    Int_t step1 = 2;
    Int_t step2 = 4;
  }
  
  Float_t ptTriggerBegin = 4.01;
  Float_t ptTriggerEnd = 19.99;
  
  Bool_t verbose = 1;
  
  loadlibs();
  
  AliUEHistograms* h = (AliUEHistograms*) GetUEHistogram(fileName);
  AliUEHist* cont = h->GetUEHist(2);
  
  cont->GetEventHist()->GetGrid(step1)->GetGrid()->GetAxis(1)->SetRange(centralityBegin, centralityEnd);
  cont->GetEventHist()->GetGrid(step2)->GetGrid()->GetAxis(1)->SetRange(centralityBegin, centralityEnd);
  
  Int_t stepEff1 = step1 - (Int_t) AliUEHist::kCFStepAnaTopology;
  if (stepEff1 == -1)
    stepEff1 = 0;
  Int_t stepEff2 = step2 - (Int_t) AliUEHist::kCFStepAnaTopology;
  if (stepEff2 == -1)
    stepEff2 = 0;
  
  cont->GetTrackHistEfficiency()->GetGrid(stepEff1)->GetGrid()->GetAxis(3)->SetRange(centralityBegin, centralityEnd);
  cont->GetTrackHistEfficiency()->GetGrid(stepEff2)->GetGrid()->GetAxis(3)->SetRange(centralityBegin, centralityEnd);
  
  TH1* triggerParticle = cont->GetEventHist()->Project(step2, 0);
  triggerParticle->Divide(cont->GetEventHist()->Project(step1, 0));
  
  TH1* singleParticle = cont->GetTrackHistEfficiency()->Project(stepEff2, 1);
  singleParticle->Divide(cont->GetTrackHistEfficiency()->Project(stepEff1, 1));
  
  //singleParticle = singleParticle->Rebin(triggerParticle->GetNbinsX(), "singleRebin", triggerParticle->GetXaxis()->GetXbins()->GetArray());
  
  TGraphErrors* effectiveEffect = new TGraphErrors;
  TGraphErrors* triggerCorrelatedEffect = new TGraphErrors;
  TGraphErrors* afterBaseLine = new TGraphErrors;
  TGraphErrors* triggerCorrelatedEffect2D = new TGraphErrors;
  
  //TVirtualFitter::SetDefaultFitter("Linear");
  
  for (Int_t bin=3; bin<=triggerParticle->GetNbinsX(); bin++)
  {
    if (triggerParticle->GetBinCenter(bin) > 10)
      continue;
    
    cont->SetPtRange(triggerParticle->GetXaxis()->GetBinLowEdge(bin) + 0.01, triggerParticle->GetXaxis()->GetBinUpEdge(bin) - 0.01);
      
    TH1* hist1 = cont->GetUEHist(step1, 0, ptTriggerBegin, ptTriggerEnd, centralityBegin, centralityEnd);
    TH1* hist2 = cont->GetUEHist(step2, 0, ptTriggerBegin, ptTriggerEnd, centralityBegin, centralityEnd);
    
    // TODO the uncertainties on the ratios should be properly calculated! (if that helps ;))
    clone = (TH1*) hist2->Clone("clone");
    clone->Divide(hist1);
    
    //func = new TF1("func", "[0]", -10, 10);
    
    if ((Int_t) clone->Fit("pol0", "0", "", -10, 10) == 0)
    {
      func = clone->GetFunction("pol0");
      effectiveEffect->SetPoint(effectiveEffect->GetN(), triggerParticle->GetBinCenter(bin) - 0.2, func->GetParameter(0));
      effectiveEffect->SetPointError(effectiveEffect->GetN()-1, 0, func->GetParError(0));
      
      if (verbose)
      {
        new TCanvas;
        clone->GetYaxis()->SetRangeUser(0.9, 1.3);
        clone->Draw();
        func->DrawClone("SAME");
      }
    }
    
    if ((Int_t) clone->Fit("pol0", "0", "", -0.3, 0.3) == 0)
    {
      func = clone->GetFunction("pol0");
      triggerCorrelatedEffect->SetPoint(triggerCorrelatedEffect->GetN(), triggerParticle->GetBinCenter(bin) - 0.1, func->GetParameter(0));
      triggerCorrelatedEffect->SetPointError(triggerCorrelatedEffect->GetN()-1, 0, func->GetParError(0));
      
      if (verbose)
      {
        func->SetLineColor(2);
        func->DrawClone("SAME");
      }
    }
    
    if ((Int_t) hist2->Fit("pol0", "0", "", 1, 2) == 0)
    {
      func = hist2->GetFunction("pol0");
      func->SetRange(-10, 10);
      hist2->Add(func, -1);
    
      if ((Int_t) hist1->Fit("pol0", "0", "", 1, 2) == 0)
      {
        func = hist1->GetFunction("pol0");
        func->SetRange(-10, 10);
        hist1->Add(func, -1);
        
        //new TCanvas; hist1->DrawCopy();
        
        hist2->Divide(hist1);
        
        if ((Int_t) hist2->Fit("pol0", "0", "", -0.3, 0.3) == 0)
        {
          func = hist2->GetFunction("pol0");
          afterBaseLine->SetPoint(afterBaseLine->GetN(), triggerParticle->GetBinCenter(bin) + 0.2, func->GetParameter(0));
          afterBaseLine->SetPointError(afterBaseLine->GetN()-1, 0, func->GetParError(0));
        
          if (verbose)
          {
            new TCanvas;
            hist2->GetYaxis()->SetRangeUser(0.9, 1.3);
            hist2->DrawCopy();
            func->SetLineColor(4);
            func->DrawClone("SAME");
          }
        }
      }
    }
      
    // 2d
    TH2* hist1_2D = cont->GetUEHist(step1, 0, ptTriggerBegin, ptTriggerEnd, centralityBegin, centralityEnd, 1);
    TH2* hist2_2D = cont->GetUEHist(step2, 0, ptTriggerBegin, ptTriggerEnd, centralityBegin, centralityEnd, 1);
    
    //((TH2*)hist1)->Rebin2D(2, 2); ((TH2*)hist2)->Rebin2D(2, 2);
    
    hist2_2D->Divide(hist1_2D);
    
    if (verbose)
    {
      new TCanvas;
      hist2_2D->Draw("COLZ");
    }
    
    //Printf("%d %d %d %d", hist2_2D->GetXaxis()->FindBin(-0.01), hist2_2D->GetXaxis()->FindBin(0.01), hist2_2D->GetYaxis()->FindBin(-0.01), hist2_2D->GetYaxis()->FindBin(0.01));
    
    Double_t error;
    Float_t integral = hist2_2D->IntegralAndError(hist2_2D->GetXaxis()->FindBin(-0.01), hist2_2D->GetXaxis()->FindBin(0.01), hist2_2D->GetYaxis()->FindBin(-0.01), hist2_2D->GetYaxis()->FindBin(0.01), error);
    
    integral /= hist2_2D->GetXaxis()->FindBin(0.01) - hist2_2D->GetXaxis()->FindBin(-0.01) + 1;
    integral /= hist2_2D->GetYaxis()->FindBin(0.01) - hist2_2D->GetYaxis()->FindBin(-0.01) + 1;
    
    error /= hist2_2D->GetXaxis()->FindBin(0.01) - hist2_2D->GetXaxis()->FindBin(-0.01) + 1;
    error /= hist2_2D->GetYaxis()->FindBin(0.01) - hist2_2D->GetYaxis()->FindBin(-0.01) + 1;
    
    triggerCorrelatedEffect2D->SetPoint(triggerCorrelatedEffect2D->GetN(), triggerParticle->GetBinCenter(bin) + 0.1, integral);
    triggerCorrelatedEffect2D->SetPointError(triggerCorrelatedEffect2D->GetN()-1, 0, error);
    
    if (verbose)
      break;
  }
  
  new TCanvas;
  triggerParticle->Draw();
  triggerParticle->GetXaxis()->SetRangeUser(0, 9.9);
  
  singleParticle->SetLineColor(2);
  //singleParticle->Draw("SAME");
  
  effectiveEffect->SetMarkerStyle(24);
  effectiveEffect->Draw("PSAME");
  
  triggerCorrelatedEffect->SetMarkerStyle(25);
  triggerCorrelatedEffect->SetMarkerColor(2);
  triggerCorrelatedEffect->SetLineColor(2);
  triggerCorrelatedEffect->Draw("PSAME");
  
  triggerCorrelatedEffect2D->SetMarkerStyle(27);
  triggerCorrelatedEffect2D->SetMarkerColor(3);
  triggerCorrelatedEffect2D->SetLineColor(3);
  triggerCorrelatedEffect2D->Draw("PSAME");
  
  afterBaseLine->SetMarkerStyle(26);
  afterBaseLine->SetMarkerColor(4);
  afterBaseLine->SetLineColor(4);
  afterBaseLine->Draw("PSAME");

  legend = new TLegend(0.66, 0.15, 0.88, 0.38);
  legend->SetFillColor(0);
  legend->SetTextSize(0.04);
  legend->AddEntry(triggerParticle, "trigger", "L");
  //legend->AddEntry(singleParticle, "single", "L");
  legend->AddEntry(effectiveEffect, "effective", "P");
  legend->AddEntry(triggerCorrelatedEffect, "at 0", "P");
  legend->AddEntry(triggerCorrelatedEffect2D, "at 0 (2d)", "P");
  legend->AddEntry(afterBaseLine, "baseline", "P");
  legend->Draw();
}

void EvaluateParticleEfficiency2D(const char* fileName)
{
  Int_t centralityBegin = 1;
  Int_t centralityEnd = 15;
  
  if (1)
  {
    Int_t step1 = 4;
    Int_t step2 = 5;
  }
  else
  {
    Int_t step1 = 2;
    Int_t step2 = 4;
  }
  
  Bool_t verbose = 0;
  
  loadlibs();
  
  AliUEHistograms* h = (AliUEHistograms*) GetUEHistogram(fileName);
  AliUEHist* cont = h->GetUEHist(2);
  
  cont->GetEventHist()->GetGrid(step1)->GetGrid()->GetAxis(1)->SetRange(centralityBegin, centralityEnd);
  cont->GetEventHist()->GetGrid(step2)->GetGrid()->GetAxis(1)->SetRange(centralityBegin, centralityEnd);
  
  Int_t stepEff1 = step1 - (Int_t) AliUEHist::kCFStepAnaTopology;
  if (stepEff1 == -1)
    stepEff1 = 0;
  Int_t stepEff2 = step2 - (Int_t) AliUEHist::kCFStepAnaTopology;
  if (stepEff2 == -1)
    stepEff2 = 0;
  
  cont->GetTrackHistEfficiency()->GetGrid(stepEff1)->GetGrid()->GetAxis(3)->SetRange(centralityBegin, centralityEnd);
  cont->GetTrackHistEfficiency()->GetGrid(stepEff2)->GetGrid()->GetAxis(3)->SetRange(centralityBegin, centralityEnd);
  
  TH1* triggerParticle = cont->GetEventHist()->Project(step2, 0);
  triggerParticle->Divide(cont->GetEventHist()->Project(step1, 0));
  
  TH1* singleParticle = cont->GetTrackHistEfficiency()->Project(stepEff2, 1);
  singleParticle->Divide(cont->GetTrackHistEfficiency()->Project(stepEff1, 1));
  
  cont->GetTrackHist(0)->GetGrid(step1)->SetRangeUser(0, -0.01, 0.01); // delta eta
  cont->GetTrackHist(0)->GetGrid(step1)->SetRangeUser(4, -0.01, 0.01); // delta phi
  TH2* tracksStep1 = (TH2*) cont->GetTrackHist(0)->Project(step1, 1, 2);
  
  cont->GetTrackHist(0)->GetGrid(step2)->SetRangeUser(0, -0.01, 0.01); // delta eta
  cont->GetTrackHist(0)->GetGrid(step2)->SetRangeUser(4, -0.01, 0.01); // delta phi
  TH2* tracksStep2 = (TH2*) cont->GetTrackHist(0)->Project(step2, 1, 2);
  
  tracksStep2->Divide(tracksStep1);
  
  TH1* triggersStep1 = cont->GetEventHist()->Project(step1, 0);
  TH1* triggersStep2 = cont->GetEventHist()->Project(step2, 0);
  
  for (Int_t x=1; x<=tracksStep2->GetNbinsX(); x++)
    for (Int_t y=1; y<=tracksStep2->GetNbinsY(); y++)
      if (singleParticle->GetBinContent(x) > 0)
        tracksStep2->SetBinContent(x, y, tracksStep2->GetBinContent(x, y) / triggersStep2->GetBinContent(y) * triggersStep1->GetBinContent(y) / singleParticle->GetBinContent(x));
      else
        tracksStep2->SetBinContent(x, y, 0);
  
  new TCanvas;
  tracksStep2->Draw("COLZ");
  
  TGraphErrors** triggerCorrelatedEffect2D = new TGraphErrors*[triggerParticle->GetNbinsX()+1];
  TGraphErrors** triggerCorrelatedEffect2DSub = new TGraphErrors*[triggerParticle->GetNbinsX()+1];
  
  for (Int_t bin1=0; bin1<=triggerParticle->GetNbinsX(); bin1++)
  {
    triggerCorrelatedEffect2D[bin1] = 0;
    triggerCorrelatedEffect2DSub[bin1] = 0;
  }
  
  for (Int_t bin1=5; bin1<=triggerParticle->GetNbinsX(); bin1++)
  {
    if (triggerParticle->GetBinCenter(bin1) > 10)
      continue;
    
    Float_t ptTriggerBegin = triggerParticle->GetXaxis()->GetBinLowEdge(bin1) + 0.01;
    Float_t ptTriggerEnd = triggerParticle->GetXaxis()->GetBinUpEdge(bin1) - 0.01;
    
    for (Int_t bin2=3; bin2<=triggerParticle->GetNbinsX(); bin2++)
    {
      if (triggerParticle->GetBinCenter(bin2) > 8)
        continue;
        
      if (bin2 > bin1)
        continue;
    
      cont->SetPtRange(triggerParticle->GetXaxis()->GetBinLowEdge(bin2) + 0.01, triggerParticle->GetXaxis()->GetBinUpEdge(bin2) - 0.01);
  
      // 2d
      TH2* hist1_2D = cont->GetUEHist(step1, 0, ptTriggerBegin, ptTriggerEnd, centralityBegin, centralityEnd, 1);
      TH2* hist2_2D = cont->GetUEHist(step2, 0, ptTriggerBegin, ptTriggerEnd, centralityBegin, centralityEnd, 1);
      
      //((TH2*)hist1)->Rebin2D(2, 2); ((TH2*)hist2)->Rebin2D(2, 2);
      
      hist2_2D->Divide(hist1_2D);
      
      if (verbose)
      {
        new TCanvas;
        hist2_2D->Draw("COLZ");
      }
      
      //Printf("%d %d %d %d", hist2_2D->GetXaxis()->FindBin(-0.01), hist2_2D->GetXaxis()->FindBin(0.01), hist2_2D->GetYaxis()->FindBin(-0.01), hist2_2D->GetYaxis()->FindBin(0.01));
      
      Double_t error;
      Float_t integral = hist2_2D->IntegralAndError(hist2_2D->GetXaxis()->FindBin(-0.01), hist2_2D->GetXaxis()->FindBin(0.01), hist2_2D->GetYaxis()->FindBin(-0.01), hist2_2D->GetYaxis()->FindBin(0.01), error);
      
      integral /= hist2_2D->GetXaxis()->FindBin(0.01) - hist2_2D->GetXaxis()->FindBin(-0.01) + 1;
      integral /= hist2_2D->GetYaxis()->FindBin(0.01) - hist2_2D->GetYaxis()->FindBin(-0.01) + 1;
      
      error /= hist2_2D->GetXaxis()->FindBin(0.01) - hist2_2D->GetXaxis()->FindBin(-0.01) + 1;
      error /= hist2_2D->GetYaxis()->FindBin(0.01) - hist2_2D->GetYaxis()->FindBin(-0.01) + 1;
      
      if (integral <= 0)
        continue;
      
      if (!triggerCorrelatedEffect2D[bin1])
        triggerCorrelatedEffect2D[bin1] = new TGraphErrors;
        
      triggerCorrelatedEffect2D[bin1]->SetPoint(triggerCorrelatedEffect2D[bin1]->GetN(), triggerParticle->GetBinCenter(bin2) - 0.5 + 0.1 * bin1, integral);
      triggerCorrelatedEffect2D[bin1]->SetPointError(triggerCorrelatedEffect2D[bin1]->GetN()-1, 0, error);
    
      if (!triggerCorrelatedEffect2DSub[bin1])
        triggerCorrelatedEffect2DSub[bin1] = new TGraphErrors;
        
      triggerCorrelatedEffect2DSub[bin1]->SetPoint(triggerCorrelatedEffect2DSub[bin1]->GetN(), triggerParticle->GetBinCenter(bin2) - 0.5 + 0.1 * bin1, integral / triggerParticle->GetBinContent(bin2));
      triggerCorrelatedEffect2DSub[bin1]->SetPointError(triggerCorrelatedEffect2DSub[bin1]->GetN()-1, 0, error / triggerParticle->GetBinContent(bin2));
      
      if (verbose)
        break;
    }
      
    if (verbose)
      break;
  }
  
  new TCanvas;
  triggerParticle->Draw();
  triggerParticle->GetXaxis()->SetRangeUser(0, 9.9);
  
  legend = new TLegend(0.66, 0.15, 0.88, 0.38);
  legend->SetFillColor(0);
  legend->SetTextSize(0.04);
  legend->AddEntry(triggerParticle, "trigger", "L");
  
  Int_t marker = 20;
  Int_t color = 1;
  
  for (Int_t bin1=3; bin1<=triggerParticle->GetNbinsX(); bin1++)
    if (triggerCorrelatedEffect2D[bin1])
    {
      triggerCorrelatedEffect2D[bin1]->SetMarkerStyle(marker++);
      triggerCorrelatedEffect2D[bin1]->SetMarkerColor(color);
      triggerCorrelatedEffect2D[bin1]->SetLineColor(color++);
      triggerCorrelatedEffect2D[bin1]->Draw("PSAME");
      legend->AddEntry(triggerCorrelatedEffect2D[bin1], Form("pt %.2f", triggerParticle->GetBinCenter(bin1)), "P");
    }
  
  legend->Draw();

  new TCanvas;
  triggerParticle->Draw();
  triggerParticle->GetXaxis()->SetRangeUser(0, 9.9);
  
  marker = 20;
  color = 1;
  
  for (Int_t bin1=3; bin1<=triggerParticle->GetNbinsX(); bin1++)
    if (triggerCorrelatedEffect2DSub[bin1])
    {
      triggerCorrelatedEffect2DSub[bin1]->SetMarkerStyle(marker++);
      triggerCorrelatedEffect2DSub[bin1]->SetMarkerColor(color);
      triggerCorrelatedEffect2DSub[bin1]->SetLineColor(color++);
      triggerCorrelatedEffect2DSub[bin1]->Draw("PSAME");
    }
  
  legend->Draw();
}

void DrawEventCount(const char* fileName)
{
  loadlibs();
  
  AliUEHistograms* h = (AliUEHistograms*) GetUEHistogram(fileName);
  h->GetEventCount()->Draw("TEXT");
}

void FitDCADistributions(const char* fileName1)
{
  loadlibs();
  
  TFile::Open(fileName1);
  list = (TList*) gFile->Get("histosPhiCorrelationsQA");
  prim = (TH2*) list->FindObject("fDCAPrimaries");
  sec  = (TH2*) list->FindObject("fDCASecondaries");
  
  Float_t zCut = 0.5;

  TH2* histList[] = { prim, sec };

  TH1* primProj = 0;
  TH1* secProj = 0;
  
  for (Int_t i=0; i<2; i++)
  {
    new TCanvas;
    gPad->SetGridx();
    gPad->SetGridy();
    gPad->SetLogy();
  
//     hist = histList[i]->ProjectionX("proj", prim->GetYaxis()->FindBin(-zCut), prim->GetYaxis()->FindBin(zCut))->DrawCopy("");
    hist = histList[i]->ProjectionY("proj", prim->GetXaxis()->FindBin(-zCut), prim->GetXaxis()->FindBin(zCut))->DrawCopy("");
    hist->SetStats(0);
    
    func = new TF1("func", "gaus(0)+gaus(3)");
    func->SetParameters(1e6, 0, 0.02, 1e3, 0, 2);
    func->FixParameter(1, 0);
    func->FixParameter(4, 0);
    func->SetParLimits(2, 0, 0.1);
    func->SetParLimits(5, 1, 3);
    
    hist->Fit(func, "", "", -2, 2);
    
    func2 = new TF1("func2", "gaus(0)");
    func2->SetParameters(1e3, 0, 2);
    func2->FixParameter(1, 0);
    func2->SetParLimits(2, 0.5, 5);
    
    hist->Fit(func2, "+", "", -2, -1);
    func2->SetRange(-3, 3);
    func2->Draw("SAME");

//     break;
  }
}

void CompareDCADistributions(const char* fileName1, const char* fileName2)
{
  TFile::Open(fileName1);
  list = (TList*) gFile->Get("histosPhiCorrelationsQA");
  prim = (TH2*) list->FindObject("fDCAPrimaries");
  sec  = (TH2*) list->FindObject("fDCASecondaries");
  
  TFile::Open(fileName2);
  list = (TList*) gFile->Get("histosPhiCorrelationsQA");
  all = (TH2*) list->FindObject("fDCAPrimaries");
  
  Float_t zCut = 3.2;
  
  TH2* histList[] = { prim, sec, all };
  
  Float_t count = 0;
  
  TH1* primProj = 0;
  TH1* secProj = 0;
  
  new TCanvas;
  gPad->SetGridx();
  gPad->SetGridy();
  
  for (Int_t i=0; i<3; i++)
  {
    hist = histList[i]->ProjectionX("proj", prim->GetYaxis()->FindBin(-zCut), prim->GetYaxis()->FindBin(zCut))->DrawCopy(i>0?"SAME":"");
    hist->SetStats(0);
    //hist = histList[i]->ProjectionY("proj", prim->GetXaxis()->FindBin(-zCut), prim->GetXaxis()->FindBin(zCut))->DrawCopy(i>0?"SAME":"");
    hist->SetLineColor(i+1);
    if (i == 0)
      primProj = hist;
    if (i == 1)
      secProj = hist;
    if (i == 2)
      hist->Scale(1.0 / hist->Integral(hist->GetXaxis()->FindBin(-0.5), hist->GetXaxis()->FindBin(0.5)) * count);
    else
      count += hist->Integral(hist->GetXaxis()->FindBin(-0.5), hist->GetXaxis()->FindBin(0.5));
  }
  
  gPad->SetLogy();
  
  ratio = (TH1*) hist->Clone("ratio");
  ratio->Add(primProj, -1);
  
/*  for (Int_t i=hist->GetXaxis()->FindBin(-0.1); i<=hist->GetXaxis()->FindBin(0.1); i++)
    ratio->SetBinContent(i, 0);*/
  
  for (Int_t i=1; i<=hist->GetNbinsX(); i++)
    ratio->SetBinError(i, 0);
  
  ratio->SetLineColor(4);
  ratio->DrawCopy("SAME E");
  
  new TCanvas;
  gPad->SetGridx();
  gPad->SetGridy();
  
  ratio->Rebin(4);
  secProj = (TH1*) secProj->Clone();
  secProj->Rebin(4);
  
  ratio->Divide(secProj);
  
/*  for (Int_t i=1; i<=hist->GetNbinsX(); i++)
    ratio->SetBinError(i, 0);*/
  
  ratio->Draw("HIST");
  ratio->GetYaxis()->SetRangeUser(0.7, 1.3);
  
  ratio->Fit("pol0", "0+", "", -5, -1);
  ratio->Fit("pol0", "0+", "", 1, 5);
}
 
void TrackCuts_CompareParameters(const char* fileName1, const char* fileName2, const char* histName, const char* cutFolder = "cuts_quality_only")
{
  // plotWhich: 0 = only before
  //            1 = both
  //            2 = only after
  //
  // mirror: kTRUE --> project negative values on the positive side
  
  Int_t plotWhich = 0;
  Bool_t mirror = kFALSE;
  
  TFile* files[2];
  files[0] = TFile::Open(fileName1);
  files[1] = TFile::Open(fileName2);

  Int_t count = 0;
  Int_t colors[] = { 1, 2, 3, 4, 5, 6 };

  TLegend* legend = new TLegend(0.7, 0.85, 0.93, 0.98);
  legend->SetTextSize(0.04);
  legend->SetFillColor(0);

  TCanvas* c1 = new TCanvas(histName, histName, 1200, 600);
  c1->Divide(2, 1);
  //TCanvas* c2 = new TCanvas("c2", "c2", 800, 600);
  //TCanvas* c3 = new TCanvas("c3", "c3", 800, 600);

  const char* folders2[] = { "before_cuts", "after_cuts" };
  const char* names[] = { "MC", "Data" };
  Bool_t first = kTRUE;
  for (Int_t j = ((plotWhich < 2) ? 0 : 1); j < ((plotWhich > 0) ? 2 : 1); j++)
  {
    TH1* base = 0;
    TH1* base2 = 0;
    for (Int_t i = 0; i < 2; i++)
    {
      Printf("%d %d", j, i);
      TString folder;
      folder.Form("%s/%s/%s", cutFolder, folders2[j], histName);
      TH1* hist = (TH1*) files[i]->Get(folder);
      
      if (mirror)
      {
        for (Int_t bin=1; bin<=hist->GetXaxis()->FindBin(-0.0001); bin++)
        {
          Int_t newBin = hist->GetXaxis()->FindBin(-hist->GetXaxis()->GetBinCenter(bin));
          if (bin != newBin)
          {
            hist->Fill(-hist->GetXaxis()->GetBinCenter(bin), hist->GetBinContent(bin));
            hist->SetBinContent(bin, 0);
          }
        }
      }
      
      legend->AddEntry(hist, Form("%s %s", names[i], (plotWhich == 1) ? folders2[j] : " "));

      c1->cd(1);
      hist->SetLineColor(colors[count]);
      hist->Scale(1.0 / hist->Integral());
      hist->SetStats(0);
      hist->DrawCopy((count == 0) ? "" : "SAME");

      switch (i)
      {
        case 0: base = hist; break;
        case 1: base2 = hist; break;
      }

      count++;
    }
    
    TH1* ratio = base;
    ratio->Divide(base2);

    ratio->GetYaxis()->SetRangeUser(0, 2);
    
    c1->cd(2);
    ratio->DrawCopy((first) ? "" : "SAME");
    first = kFALSE;
  }

  c1->cd(1)->SetLogy();
  c1->cd(1)->SetGridx();
  c1->cd(1)->SetGridy();
  legend->Draw();
  
  c1->cd(2)->SetGridx();
  c1->cd(2)->SetGridy();
  
  c1->SaveAs(Form("%s.png", histName));
}

void PlotQA(const char* fileName)
{
  loadlibs();
  
  if (!gGrid && TString(fileName).BeginsWith("alien://"))
    TGrid::Connect("alien://");
  
  TFile::Open(fileName);
  
  // phys sel
  Int_t runNumber = 0;
  list = (TList*) gFile->Get("PhysSel");
  if (list)
  {
    physSel = (AliPhysicsSelection*) list->FindObject("AliPhysicsSelection");
//     runNumber = physSel->GetCurrentRun();
  }
  
  TString tmp(fileName);
  tmp.ReplaceAll("alien:///alice/cern.ch/user/j/jgrosseo/gridjob/dir_", "");
  tmp.ReplaceAll(".root", ".png");
  tmp.ReplaceAll("/", "-");
  TString title;
  title.Form("QA_%d_%s", runNumber, tmp.Data());
  c = new TCanvas(title, title, 1200, 800);
  c->Divide(3, 3);

  // QA task
  list = (TList*) gFile->Get("histosPhiCorrelationsQA");
  if (list)
  {
    prim = (TH2*) list->FindObject("fDCAPrimaries");
    dcaxy = prim->ProjectionX("dcaxy", prim->GetYaxis()->FindBin(-3.2), prim->GetYaxis()->FindBin(3.2));
    dcaz = prim->ProjectionY("dcaz", prim->GetXaxis()->FindBin(-2.4), prim->GetXaxis()->FindBin(2.4));
    centrCorr = (TH2*) list->FindObject("fCentralityCorrelation");
    
    c->cd(1); dcaxy->Draw(); dcaz->SetLineColor(2); dcaz->Draw("SAME");  gPad->SetLogy(); 
    c->cd(6); centrCorr->Draw("COLZ"); gPad->SetLogz();
  
    cuts = (AliESDtrackCuts*) list->FindObject("cuts_quality_dca");
    if (cuts)
    {
      cluster = cuts->GetNClustersTPC(1);
      c->cd(3); cluster->Draw();
    
      ptall = (TH1F*) cuts->GetPtHist(1)->Clone("ptall");
      if (ptall->Integral() > 0)
	ptall->Scale(1.0 / ptall->Integral());
    
      c->cd(7); 
      ptall->Draw(); 

      TH1* ptIts = 0;
      check_its = (AliESDtrackCuts*) list->FindObject("check_its");
      if (check_its)
      {
	ptIts = (TH1F*) check_its->GetPtHist(1)->Clone("ptIts");
	if (ptIts->Integral() > 0)
	  ptIts->Scale(1.0 / ptIts->Integral());
      }
      
      TH1* ptItsAcc = 0;
      global_cuts = (AliESDtrackCuts*) list->FindObject("global_cuts");
      if (global_cuts)
      {
	ptItsAcc = (TH1*) global_cuts->GetPtHist(1)->Clone("ptItsAcc");
	if (ptItsAcc->Integral() > 0)
	  ptItsAcc->Scale(1.0 / ptItsAcc->Integral());
      }
    
      if (ptIts)
      {
	ptIts->SetLineColor(2); 
	ptIts->Draw("SAME");
      }
      
      if (ptItsAcc)
      {
	ptItsAcc->SetLineColor(4); 
	ptItsAcc->Draw("SAME"); 
      }
      
      gPad->SetLogy();
    }
  }
  
  // centrality task
  list = (TList*) gFile->Get("CentralityStat");
  TH1* centrQuality = 0;
  if (list)
  {
    centrQuality = (TH1*) list->FindObject("fHOutQuality");
    
    c->cd(4); if (centrQuality) centrQuality->Draw();
  }
  
  // phi corr task
  AliUEHistograms* h = (AliUEHistograms*) GetUEHistogram(fileName);
  
  if (h->GetUEHist(2)->GetTrackHist(0)->GetGrid(6)->GetGrid()->GetNbins() == 0)
  {
    Printf("We have %d axes", ((AliTHn*) h->GetUEHist(2)->GetTrackHist(0)->GetNVar()));
    
//     ((AliTHn*) h->GetUEHist(2)->GetTrackHist(0))->ReduceAxis();
    ((AliTHn*) h->GetUEHist(2)->GetTrackHist(0))->FillParent();
    ((AliTHn*) h->GetUEHist(2)->GetTrackHist(0))->DeleteContainers();
  }
  centr = h->GetCentralityDistribution();
  NormalizeToBinWidth(centr);
  Int_t events = (Int_t) h->GetEventCount()->ProjectionX()->GetBinContent(3);
  
  h->SetPtRange(1.01, 3.99);
  dphi_corr = h->GetUEHist(2)->GetUEHist(AliUEHist::kCFStepReconstructed, 0, 4.01, 14.99, 1, 8);
  
  AliUEHistograms* hMixed = (AliUEHistograms*) GetUEHistogram(fileName, 0, kTRUE);
  if (hMixed->GetUEHist(2)->GetTrackHist(0)->GetGrid(6)->GetGrid()->GetNbins() == 0)
  {
//     ((AliTHn*) hMixed->GetUEHist(2)->GetTrackHist(0))->ReduceAxis();
    ((AliTHn*) hMixed->GetUEHist(2)->GetTrackHist(0))->FillParent();
  }

  hMixed->SetPtRange(1.01, 3.99);
  dphi_corr_mixed = hMixed->GetUEHist(2)->GetUEHist(AliUEHist::kCFStepReconstructed, 0, 4.01, 14.99, 1, 8);
  
  if (runNumber != 0 && runNumber != h->GetRunNumber())
    AliFatal("Run numbers inconsistent");

  Printf("%d", h->GetRunNumber());
  if (runNumber == 0)
    runNumber = h->GetRunNumber();
 
  c->cd(2); dphi_corr->Draw(); dphi_corr->GetYaxis()->SetRangeUser(dphi_corr->GetMinimum() * 0.9, dphi_corr->GetMaximum() * 1.1); dphi_corr_mixed->SetLineColor(2); dphi_corr_mixed->Draw("SAME");
  c->cd(5); centr->Draw("HIST");
  
  c->cd(1);
  latex = new TLatex(0.15, 0.8, Form("%d events", events));
  latex->SetTextSize(0.075);
  latex->SetNDC();
  latex->Draw();
  
  c->cd(8);
  h->GetEventCount()->Draw("TEXT");
  
  title.Form("QA_%d_%s", runNumber, tmp.Data());
  c->SetTitle(title);

  c->SaveAs(Form("qa/%s", c->GetTitle()));
}

void CompareStepsOnePlot(const char* fileName, Int_t caseId)
{
  loadlibs();

  AliUEHistograms* h = (AliUEHistograms*) GetUEHistogram(fileName);
  
  Int_t stepList[] = { 0, 1, 2, 4, 5, 6 };
  Int_t stepListNames[] = { 0, 1, 2, 3, 4, 5 };
  Int_t nSteps = 6;
  //const char* names[] = { "All", "PhysSel+Vertex", "Reco Primaries", "Reco "
//   //Int_t stepList[] = { 2, 4, 5, 6 };
//   
//   
//   TH1* hist1 = h->GetUEHist(2)->GetUEHist(0, 0, 6.01, 9.99, 6, 14);
//   TH1* hist2 = h->GetUEHist(2)->GetUEHist(0, 0, 6.01, 7.99, 6, 14);
//   TH1* hist3 = h->GetUEHist(2)->GetUEHist(0, 0, 8.01, 9.99, 6, 14);
//   
//   hist1->Draw();
//   hist2->SetLineColor(2);
//   hist2->Draw("SAME");
//   hist3->SetLineColor(4);
//   hist3->Draw("SAME");
//   return;
//   
  
  legend = new TLegend(0.65, 0.45, 0.85, 0.85);
  legend->SetTextSize(0.04);
  legend->SetFillColor(0);
  
  TH1* histList[10];
  
  for (Int_t i=0; i<nSteps; i++)
  {
    //TH1* hist1 = h->GetUEHist(2)->GetEventHist()->Project(stepList[i], 0);
    
/*    h->GetUEHist(2)->GetTrackHist(0)->GetGrid(stepList[i])->GetGrid()->GetAxis(1)->SetRangeUser(1.01, 3.99);
    h->GetUEHist(2)->GetTrackHist(0)->GetGrid(stepList[i])->GetGrid()->GetAxis(2)->SetRangeUser(4.01, 19.99);
    TH1* hist1 = h->GetUEHist(2)->GetTrackHist(0)->Project(stepList[i], 4);*/
    
    if (caseId == 0)
      TH1* hist1 = h->GetUEHist(2)->GetUEHist(stepList[i], 0, 1.01, 19.99);
    else if (caseId == 1)
      TH1* hist1 = h->GetUEHist(2)->GetEventHist()->Project(stepList[i], 0);
    
    //hist1->Rebin(2);
    //hist1->Rebin(4);
    hist1->SetMarkerStyle(24+i);
    hist1->SetTitle("");
    hist1->SetStats(0);
    hist1->DrawCopy((i==0)?"":"SAME");
    hist1->Fit("pol0", "0");
    //hist1->Fit("pol0");
    histList[i] = hist1;
    legend->AddEntry(hist1, Form("Step %d", stepListNames[i]), "P");
  }
  
  legend->Draw();
  
  new TCanvas;
  legend = new TLegend(0.65, 0.45, 0.85, 0.85);
  legend->SetTextSize(0.04);
  legend->SetFillColor(0);
  
  for (Int_t i=1; i<nSteps; i++)
  {
    histList[i]->DrawCopy((i==1)?"":"SAME")->Divide(histList[i-1]);
    legend->AddEntry(histList[i], Form("Step %d/%d", stepListNames[i], stepListNames[i-1]), "P");
  }
  
  legend->Draw();
  
  new TCanvas;
  for (Int_t i=0; i<nSteps; i++)
  {
    hist1 = histList[i];
    func = new TF1("func", "[0]", -10, 10);
    hist1->Fit(func, "0", "", 1, 2);
    hist1->Add(func, -1);
    hist1->DrawCopy((i==0)?"":"SAME");
  }
  
  legend->Draw();
  
  new TCanvas;
  for (Int_t i=1; i<nSteps; i++)
  {
    histList[i]->DrawCopy((i==1)?"":"SAME")->Divide(histList[i-1]);
  }
  
  legend->Draw();
}

void PtShift(const char* fileName)
{
  loadlibs();

  AliUEHistograms* h = (AliUEHistograms*) GetUEHistogram(fileName);
  
  new TCanvas;
  gPad->SetLogy();
  
  ptHist3 = h->GetCorrelationpT()->ProjectionY();
  ptHist3->GetXaxis()->SetRangeUser(0, 20);
  //ptHist3->Sumw2();
  
  /*
  func = new TF1("func", "[0]+[1]*x**[2]", 4, 20);
  func->SetParLimits(2, -10, -1);
  func->SetParameters(0, 1, -4);
  ptHist3->Fit(func, "", "", 4, 20);
  
  return;
  */
  
  ptHist3Shifted = (TH1*) ptHist3->Clone("ptHist3Shifted");
  for (Int_t i=1; i<ptHist3Shifted->GetNbinsX(); i++)
    //ptHist3Shifted->SetBinContent(i, ptHist3Shifted->GetBinContent(i+3));
    ptHist3Shifted->SetBinContent(i, ptHist3Shifted->GetBinContent(1.1*i));
  
  // 10% or 750 MeV
    
  ptHist3Shifted->SetLineColor(2);
  
  ptHist3->Rebin(2);
  ptHist3Shifted->Rebin(2);
  
  ptHist3->DrawCopy();
  ptHist3Shifted->DrawCopy("SAME");
  
  new TCanvas;
  ptHist3->Divide(ptHist3Shifted);
  ptHist3->Draw();
  
  //ptHist3->Rebin(2); ptHist3->Scale(1.0 / 2);
  
  gPad->SetGridx();
  gPad->SetGridy();
}

void RAA(const char* fileName, const char* fileName2)
{
  loadlibs();
  
  TH1* hist[2];
  
  for (Int_t j=0; j<2; j++)
  {
    if (j == 0)
      AliUEHistograms* h = (AliUEHistograms*) GetUEHistogram(fileName);
    else
    {
      if (!fileName2)
        break;
        
      AliUEHistograms* h = (AliUEHistograms*) GetUEHistogram(fileName2);
    }
  
    ptHist3 = h->GetCorrelationpT()->ProjectionY(Form("proj1_%d", j), 1, 5);
    ptHist3->GetXaxis()->SetRangeUser(0, 20);
    
    //h->GetCentralityDistribution()->Draw(); new TCanvas;
    
    Int_t nEvents = h->GetCentralityDistribution()->Integral(h->GetCentralityDistribution()->FindBin(0.01), h->GetCentralityDistribution()->FindBin(4.99));
    Printf("%d", nEvents);
    ptHist3->Scale(1.0 / 1.6 / TMath::TwoPi() / nEvents / ptHist3->GetBinWidth(1));
    
    hist[j] = ptHist3;
  }
  
  new TCanvas;
  hist[0]->DrawCopy("");
  hist[1]->DrawCopy("SAME")->SetLineColor(2);
  gPad->SetLogy();
  
  new TCanvas;
  hist[1]->Divide(hist[0]);
  hist[1]->DrawCopy();
  
  raa_central = ReadHepdata("raa_alice_central.txt", kFALSE, 3);
  raa_central->SetMarkerStyle(20);
  raa_central->Draw("PSAME");
}

void PtComparison(const char* fileName, const char* fileName2 = 0)
{
  loadlibs();

  c = new TCanvas;
  c2 = new TCanvas;
    
  for (Int_t j=0; j<2; j++)
  {
    if (j == 0)
      AliUEHistograms* h = (AliUEHistograms*) GetUEHistogram(fileName);
    else
    {
      if (!fileName2)
        break;
        
      AliUEHistograms* h = (AliUEHistograms*) GetUEHistogram(fileName2);
    }
    
    //new TCanvas; h->GetCorrelationpT()->Draw("COLZ");
    
    ptHist3 = h->GetCorrelationpT()->ProjectionY(Form("proj1_%d", j), 1, 5);
    ptHist3->GetXaxis()->SetRangeUser(0, 20);
    
    //h->GetCentralityDistribution()->Draw(); new TCanvas;
    
    Int_t nEvents = h->GetCentralityDistribution()->Integral(h->GetCentralityDistribution()->FindBin(0.01), h->GetCentralityDistribution()->FindBin(4.99));
    Printf("%d", nEvents);
    ptHist3->Scale(1.0 / 1.6 / TMath::TwoPi() / nEvents / ptHist3->GetBinWidth(1));
    
    //ptHist3->Scale(1690);
    //ptHist3->Scale(1.0 / ptHist3->Integral(ptHist3->GetXaxis()->FindBin(1.01), ptHist3->GetNbinsX()) / 3 / 0.6);
    
    for (Int_t i=2; i<ptHist3->GetNbinsX(); i++)
      ptHist3->SetBinContent(i, ptHist3->GetBinContent(i) / ptHist3->GetBinCenter(i));
      //ptHist3->SetBinContent(i, ptHist3->GetBinContent(i) / ptHist3->GetBinLowEdge(i));
    
    //AliPWG0Helper::NormalizeToBinWidth(ptHist3);
    
    c->cd();
    ptHist3->SetLineStyle(j+1);
    ptHist3->DrawCopy((j == 0) ? "HIST" : "HISTSAME")->SetLineColor(1);

    centralPt = (TH1*) ptHist3->Clone();;
  
    if (1)
    {
      ptHist3 = h->GetCorrelationpT()->ProjectionY(Form("proj2_%d", j), 71, 80);
      ptHist3->GetXaxis()->SetRangeUser(0, 20);
    
      Int_t nEvents = h->GetCentralityDistribution()->Integral(h->GetCentralityDistribution()->FindBin(70.01), h->GetCentralityDistribution()->FindBin(79.99));
      ptHist3->Scale(1.0 / 1.6 / TMath::TwoPi() / nEvents / ptHist3->GetBinWidth(1));
      //ptHist3->Scale(10.0 / ptHist3->Integral(ptHist3->GetXaxis()->FindBin(1.01), ptHist3->GetNbinsX()) / 3 / 0.6);
      Printf("%d", nEvents);
    
      for (Int_t i=2; i<ptHist3->GetNbinsX(); i++)
        ptHist3->SetBinContent(i, ptHist3->GetBinContent(i) / ptHist3->GetBinCenter(i));
        //ptHist3->SetBinContent(i, ptHist3->GetBinContent(i) / ptHist3->GetBinLowEdge(i));
    }
    
    ptHist3->SetLineStyle(j+1);
    ptHist3->DrawCopy("HISTSAME")->SetLineColor(2);
    
    dndpt_central = ReadHepdata("raa_dndpt_central.txt", kFALSE, 3);
    dndpt_peripheral = ReadHepdata("raa_dndpt_peripheral.txt", kFALSE, 3);
  
//   RemovePointsBelowX(dndpt_central, 1);
//   RemovePointsBelowX(dndpt_peripheral, 1);
//   
//   NormalizeTo(dndpt_central, 1);
//   NormalizeTo(dndpt_peripheral, 10);

    dndpt_central->Draw("*SAME");
    dndpt_peripheral->SetLineColor(2);
    dndpt_peripheral->SetMarkerColor(2);
    dndpt_peripheral->Draw("*SAME");
    
    gPad->SetLogy();
    
    c2->cd();
    
    for (Int_t i=1; i<ptHist3->GetNbinsX(); i++)
      ptHist3->SetBinContent(i, ptHist3->GetBinContent(i) / dndpt_peripheral->Eval(ptHist3->GetBinCenter(i)));
  
    for (Int_t i=1; i<centralPt->GetNbinsX(); i++)
      centralPt->SetBinContent(i, centralPt->GetBinContent(i) / dndpt_central->Eval(centralPt->GetBinCenter(i)));
  
    ptHist3->Rebin(2); ptHist3->Scale(0.5);
    centralPt->Rebin(2); centralPt->Scale(0.5);
  
    ptHist3->DrawCopy((j == 0) ? "" : "SAME")->SetLineColor(2);
    centralPt->DrawCopy("SAME")->SetLineColor(1);
  }

  return;
  
  ReadYields("preliminaries/yields_110303.root");
  
  ptHist3->Scale(100);
  
  TGraphErrors** tmp = yields[0][1][0];
  nearSide = tmp[18];
  nearSide->Draw("* SAME");  
 
  TGraphErrors** tmp = yields[1][1][0];
  awaySide = tmp[18];
  awaySide->SetLineColor(2);
  awaySide->SetMarkerColor(2);
  awaySide->Draw("* SAME");  
  
}

void style(Int_t styleId = 1)
{
  gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetHistFillColor(10);
  gStyle->SetHistFillStyle(0);
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetAxisColor(1, "X");
  gStyle->SetAxisColor(1, "Y");
  gStyle->SetAxisColor(1, "Z");
  gStyle->SetLabelColor(1, "X");
  gStyle->SetLabelColor(1, "Y");
  gStyle->SetLabelColor(1, "Z");
  gStyle->SetTickLength(0.03, "X");
  gStyle->SetTickLength(0.03, "Y");
  gStyle->SetTickLength(0.03, "Z");
  gStyle->SetTitleXSize(0.05);
  gStyle->SetTitleYSize(0.05);
  gStyle->SetNdivisions(506, "X");
  gStyle->SetNdivisions(506, "Y");
  gStyle->SetNdivisions(506, "Z");
  
  //gStyle->SetPadGridX(1);
  //gStyle->SetPadGridY(1);

  //gStyle->SetLabelOffset(0.02, "X");
  //gStyle->SetLabelOffset(0.02, "Y");
  //gStyle->SetLabelOffset(0.02, "Z");
  gStyle->SetLabelSize(0.05, "X");
  gStyle->SetLabelSize(0.05, "Y");
  gStyle->SetLabelSize(0.05, "Z");

  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadRightMargin(0.02);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadTopMargin(0.02);

  gStyle->SetMarkerSize(1.4); // CKB

  const int iFont = 42; // type * 10 + prec  6: helvetica bold 13: times normal 2: times bold

/*                                                               italic     weigth
*-*        1 : times-medium-i-normal      "Times New Roman"      1           4
*-*        2 : times-bold-r-normal        "Times New Roman"      0           7
*-*        3 : times-bold-i-normal        "Times New Roman"      1           7
*-*        4 : helvetica-medium-r-normal  "Arial"                0           4
*-*        5 : helvetica-medium-o-normal  "Arial"                1           4
*-*        6 : helvetica-bold-r-normal    "Arial"                0           7
*-*        7 : helvetica-bold-o-normal    "Arial"                1           7
*-*        8 : courier-medium-r-normal    "Courier New"          0           4
*-*        9 : courier-medium-o-normal    "Courier New"          1           4
*-*       10 : courier-bold-r-normal      "Courier New"          0           7
*-*       11 : courier-bold-o-normal      "Courier New"          1           7
*-*       12 : symbol-medium-r-normal     "Symbol"               0           6
*-*       13 : times-medium-r-normal      "Times New Roman"      0           4
*-*       14 :                            "Wingdings"            0           4
*/

  //gStyle->SetTitleXOffset(1); // 1.1
  //gStyle->SetTitleYOffset(1); // 1-4

  gStyle->SetLabelFont(iFont, "xyz");
  gStyle->SetStatFont(iFont);
  gStyle->SetTitleFont(iFont, "xyz");
  gStyle->SetTextFont(iFont);
  
  if (styleId == 2)
  {
    gStyle->SetLabelSize(0.07, "X");
    gStyle->SetLabelSize(0.07, "Y");
    gStyle->SetLabelSize(0.07, "Z");
    gStyle->SetTitleXSize(0.07);
    gStyle->SetTitleYSize(0.07);
  
    gStyle->SetPadLeftMargin(0.26);
    gStyle->SetPadRightMargin(0.01);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadTopMargin(0.01);
  }

  //gStyle->SetEndErrorSize(0.0);

  gROOT->ForceStyle();

}

void NormalizeTo(TGraphErrors* graph, Float_t normalizeTo)
{
	Float_t sum = 0;
	for (Int_t i=0; i<graph->GetN(); i++)
		sum += graph->GetY()[i];
	
	if (normalizeTo > 0 && sum > 0)
	{
		sum /= normalizeTo;
		for (Int_t i=0; i<graph->GetN(); i++)
		{
			graph->SetPoint(i, graph->GetX()[i],  graph->GetY()[i] / sum);
			graph->SetPointError(i, graph->GetEX()[i],  graph->GetEY()[i] / sum);
		}
	}	
}

void CompareMixedEvent(const char* fileName)
{
  loadlibs();
  
  Float_t leadingPtArr[] = { 6.0, 8.0, 10.0, 15.0, 15.0 };
  Float_t assocPtArr[] =     { 0.5, 1.5, 3.0, 4.0, 6.0, 8.0, 10.0, 12.0 };
  Int_t leadingPtOffset = 2;
  Int_t centralityBins[] = { 0, 0, 1, 6, 9, 16 };

  AliUEHistograms* hMixed = (AliUEHistograms*) GetUEHistogram(fileName, 0, kTRUE);  
  
  Int_t i = 1;
  for (Int_t j=2; j<5; j++)
  {
    gpTMin = assocPtArr[j] + 0.01;
    gpTMax = assocPtArr[j+1] - 0.01;
    
    gpTMin = 3.0;
    gpTMax = 6.0;

    SetupRanges(hMixed);
    
    TH2* mixed = hMixed->GetUEHist(2)->GetUEHist(6, 0, leadingPtArr[i] + 0.01, leadingPtArr[i+leadingPtOffset] - 0.01, centralityBins[j], centralityBins[j+1]-1, 1);
  
    // compare deta
    
    TH1* histMixedproj = mixed->ProjectionY();
    histMixedproj->Scale(1.0 / mixed->GetNbinsX());
    
    for (Int_t x=1; x<=mixed->GetNbinsX(); x++)
      for (Int_t y=1; y<=mixed->GetNbinsY(); y++)
	mixed->SetBinContent(x, y, histMixedproj->GetBinContent(y));

    histMixedproj->Scale(1.0 / (0.5 * (histMixedproj->GetBinContent(histMixedproj->GetXaxis()->FindBin(-0.01)) + histMixedproj->GetBinContent(histMixedproj->GetXaxis()->FindBin(0.01)))));
      
    histMixedproj->DrawCopy((j == 2) ? "" : "SAME")->SetLineColor(j-1);
  }
}

void FillParentTHnSparse(const char* fileName, Bool_t reduce = kTRUE)
{
  loadlibs();

  TList* list = 0;
  
  AliUEHistograms* h = (AliUEHistograms*) GetUEHistogram(fileName, &list);
  Printf("We have %d axes", ((AliTHn*) h->GetUEHist(2)->GetTrackHist(0)->GetNVar()));
  
  if (reduce)
    ((AliTHn*) h->GetUEHist(2)->GetTrackHist(0))->ReduceAxis();
  ((AliTHn*) h->GetUEHist(2)->GetTrackHist(0))->FillParent();
  ((AliTHn*) h->GetUEHist(2)->GetTrackHist(0))->DeleteContainers();
  
  AliUEHistograms* hMixed = (AliUEHistograms*) GetUEHistogram(fileName, 0, kTRUE);  
  if (reduce)
    ((AliTHn*) hMixed->GetUEHist(2)->GetTrackHist(0))->ReduceAxis();
  ((AliTHn*) hMixed->GetUEHist(2)->GetTrackHist(0))->FillParent();
  ((AliTHn*) hMixed->GetUEHist(2)->GetTrackHist(0))->DeleteContainers();
  
  TString newFileName(fileName);
  newFileName.ReplaceAll(".root", "");
  if (reduce)
    newFileName += "_.root";
  else
    newFileName += "_zvtx.root";

  file3 = TFile::Open(newFileName, "RECREATE");
  file3->mkdir("PWG4_PhiCorrelations");
  file3->cd("PWG4_PhiCorrelations");
  list->Write("histosPhiCorrelations", TObject::kSingleKey);
  file3->Close();
}

void CompareZVertex(const char* fileName)
{
  loadlibs();
  
  AliUEHistograms* h = (AliUEHistograms*) GetUEHistogram(fileName);
  AliUEHistograms* hMixed = (AliUEHistograms*) GetUEHistogram(fileName, 0, kTRUE);  
  
  axis = h->GetUEHist(2)->GetEventHist()->GetAxis(2, 6);
  
  gpTMin = 2.01;
  gpTMax = 7.99;
  
  SetupRanges(h);
  SetupRanges(hMixed);
  
  TFile::Open("test.root", "RECREATE");
	
  for (Int_t i=0; i<=axis->GetNbins(); i++)
  {
    TH1* hist = 0;
    if (i > 0)
    {
      Printf("%d %f %f", i, axis->GetBinLowEdge(i) + 0.01, axis->GetBinUpEdge(i) - 0.01);
      h->SetZVtxRange(axis->GetBinLowEdge(i) + 0.01, axis->GetBinUpEdge(i) - 0.01);
      hMixed->SetZVtxRange(axis->GetBinLowEdge(i) + 0.01, axis->GetBinUpEdge(i) - 0.01);
    }
      
    GetDistAndFlow(h, hMixed, &hist, 0, 6, 0, 10, 2.01, 14.99, 1, kTRUE, 0, kFALSE);
    
    new TCanvas;
    hist->DrawCopy("SURF1");
    
    hist->Write(Form("detadphi_%d", i));
    
//     if (i == 0)   break;

    continue;
    
    hist->SetLineColor(i+1);
    hist->Scale(1.0 / hist->Integral() / hist->GetBinWidth(1));
    hist->Draw((i == 0) ? "" : "SAME");
  }
  
  gFile->Close();
}

void DrawZRanges(Float_t min, Float_t max)
{
  legend = new TLegend(0.5, 0.5, 0.8, 0.8);
  
  TFile::Open("test.root");
  
  for (Int_t i=0; i<8; i++)
  {
    if (i == 1 || i == 7)
      continue;
    
    hist = (TH2*) gFile->Get(Form("detadphi_%d", i));
    hist->Rebin2D(2, 2);
    hist->Scale(0.25);
    
    proj = hist->ProjectionY("proj", hist->GetXaxis()->FindBin(min), hist->GetXaxis()->FindBin(max));
    proj->Scale(1.0 / (hist->GetXaxis()->FindBin(max) - hist->GetXaxis()->FindBin(min) + 1));
    
    proj->SetLineColor(i+1);
    proj->DrawCopy((i == 0) ? "" : "SAME HIST");
    
    legend->AddEntry(proj->Clone(), Form("%d", i));
  }
  
  legend->Draw();
}

void PlotCorrections(const char* fileName)
{
  loadlibs();
  
  AliUEHistograms* h = (AliUEHistograms*) GetUEHistogram(fileName);
  
  c = new TCanvas("c", "c", 1200, 800);
  c->Divide(3, 3);

  h->SetEtaRange(-0.89, 0.89);
  
  for (Int_t i=0; i<5; i++)
  {
    h->GetUEHist(2)->SetCentralityRange(100.0/5*i + 0.1, 100.0/5*(i+1) - 0.1);
    c->cd(i+1);
    h->GetUEHist(2)->GetTrackingEfficiency()->DrawClone("COLZ");
    
    c->cd(6);
    proj = h->GetUEHist(2)->GetTrackingEfficiency(1);
    proj->SetLineColor(i+1);
    proj->DrawClone((i == 0) ? "" : "SAME");
    
//     return;
  }

  h->GetUEHist(2)->SetCentralityRange(0, 100);

  c->cd(7);
  h->GetUEHist(2)->GetTrackingContamination()->Draw("COLZ");
  
  c->cd(8);
  hist = h->GetUEHist(2)->GetCorrelatedContamination();
  if (hist->GetEntries() > 0)
    hist->Draw("COLZ");
}
 
void ComparePPHIMixedEvent(const char* ppFile, const char* pbpbFile)
{
  loadlibs();
  
  AliUEHistograms* hpp = (AliUEHistograms*) GetUEHistogram(ppFile);
  AliUEHistograms* hpbpb = (AliUEHistograms*) GetUEHistogram(pbpbFile);
  
  new TCanvas;
  hpp->SetPtRange(2, 10);
  ppEff = hpp->GetUEHist(2)->GetTrackingEfficiency(0);
  ppEff->Draw();

  hpbpb->SetPtRange(2, 10);
  pbpbEff = hpbpb->GetUEHist(2)->GetTrackingEfficiency(0);
  pbpbEff->DrawCopy("SAME")->SetLineColor(2);
  
  new TCanvas;
  ppEff2 = hpp->GetUEHist(2)->GetTrackingEfficiency(1);
  ppEff2->Draw();

  pbpbEff2 = hpbpb->GetUEHist(2)->GetTrackingEfficiency(1);
  pbpbEff2->DrawCopy("SAME")->SetLineColor(2);

  mixed = ComparePPHIMixedEventGetMixed(ppEff);
  mixed2 = ComparePPHIMixedEventGetMixed(pbpbEff);
  
  new TCanvas;
  mixed->DrawCopy();
  mixed2->DrawCopy("SAME")->SetLineColor(2);
  
  new TCanvas;
  mixed->Divide(mixed2);
  mixed->Draw();
}

TH1* ComparePPHIMixedEventGetMixed(TH1* eff)
{
  eff->Fit("pol0", "0W");
  Float_t avgEff = eff->GetFunction("pol0")->GetParameter(0);

  eff->Fit("pol0", "0W", "", -0.89, 0.89);
  Float_t avgEffCenter = eff->GetFunction("pol0")->GetParameter(0);
  Printf("Avg is %f and avg in center is %f", avgEff, avgEffCenter);
	 
  TH1* mixed = new TH1F("mixed", "", 100, -2, 2);
  
  Float_t etaLimit = 1.0;
  
  Int_t n = 10000;
  Int_t n2 = 100;
  for (Int_t i=0; i<n; i++)
  {
    Float_t etaTrig = gRandom->Uniform(-etaLimit, etaLimit);
    
    for (Int_t j=0; j<n2; j++)
    {
      Float_t etaAssoc = gRandom->Uniform(-etaLimit, etaLimit);
      
      if (gRandom->Uniform(0, 1) > eff->GetBinContent(eff->FindBin(etaAssoc)))
	continue;
      
      mixed->Fill(etaTrig - etaAssoc);
    }
  }
  
//   mixed->Scale(1.0 / avgEffCenter);
  mixed->Scale(1.0 / avgEff);

  Printf("We have %f pairs and put in %d", mixed->Integral(), n*n2);    
  
  return mixed;
}

void MACHConeEvolution(const char* fileName, const char* fileNameMixed = 0)
{
  loadlibs();
  
  if (!fileNameMixed)
    fileNameMixed = fileName;
  
  AliUEHistograms* h = (AliUEHistograms*) GetUEHistogram(fileName);
  AliUEHistograms* hMixed = (AliUEHistograms*) GetUEHistogram(fileNameMixed, 0, kTRUE);
    
  Float_t leadingPtArr[] = { 2.0, 3.0, 4.0, 10.0, 20.0, 40.0 };
  Float_t assocPtArr[] =   { 1.0, 2.0, 3.0, 6.0, 10.0, 20.0, 40.0 };
  
  Int_t i = 1;
  Int_t step = 6;
  Int_t j = 0;
  
  gpTMin = assocPtArr[i] + 0.01;
  gpTMax = assocPtArr[i+1] - 0.01;

  for (Int_t centrBin = 0; centrBin < 5; centrBin++)
  {
    Int_t centralityBegin = centrBin;
    Int_t centralityEnd = centrBin+1;
    
    SetupRanges(h);
    SetupRanges(hMixed);
    
    TH1* hist = 0;

    Bool_t scaleToPairs = 0;
    
    GetDistAndFlow(h, hMixed, &hist, 0, step, centralityBegin, centralityEnd, leadingPtArr[j] + 0.01, leadingPtArr[j+1] - 0.01, 11, kTRUE, 0, scaleToPairs); 
    hist->Rebin(2); hist->Scale(0.5);

    copy = hist->DrawCopy((centrBin == 0) ? "" : "SAME");
    copy->SetLineColor(centrBin+1);
  }
}
    
void PlotTwoTrackEfficiencyControlPlots(const char* fileName)
{
  loadlibs();
  AliUEHistograms* h = (AliUEHistograms*) GetUEHistogram(fileName);
  
  Float_t ptRange[] = { 0.0, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0 };
  
  graph = new TGraphErrors;
  
  for (Int_t ptID = 0; ptID < 7; ptID++)
  {
    c = new TCanvas(Form("c%d", ptID), Form("%.0f < pT < %.0f", ptRange[ptID], ptRange[ptID+1]), 800, 400);
    c->Divide(2, 1);
  
    for (Int_t i=0; i<1; i++)
    {
      c->cd(i+1);
      h->GetTwoTrackDistance(i)->GetZaxis()->SetRangeUser(ptRange[ptID] + 0.01, ptRange[ptID+1] - 0.01);
      p = h->GetTwoTrackDistance(i)->Project3D(Form("yx_%d_%d", i, ptID));
//       if (ptID == 0)
	p->SetStats(0);
      
      pc = (TH2*) p->Clone("pc");
      pc->Reset();
      
      // reduce to one quadrant
      for (Int_t x=1; x<=p->GetNbinsX(); x++)
	for (Int_t y=1; y<=p->GetNbinsY(); y++)
	  pc->Fill(TMath::Abs(p->GetXaxis()->GetBinCenter(x)), TMath::Abs(p->GetYaxis()->GetBinCenter(y)), p->GetBinContent(x, y));
      
      p->DrawCopy("COLZ");
      
      c->cd(i+2);
      copy = pc->DrawCopy("COLZ");
      copy->GetXaxis()->SetRangeUser(0, 10);
      copy->GetYaxis()->SetRangeUser(0.00101, 10);
      
      // extract excess
      if (1)
      {
	Float_t center = pc->Integral(pc->GetXaxis()->FindBin(0), pc->GetXaxis()->FindBin(0.009999), 1, pc->GetNbinsY());
	Float_t outside = pc->Integral(pc->GetXaxis()->FindBin(0.01001), pc->GetXaxis()->FindBin(0.04999), 1, pc->GetNbinsY());
	
	Float_t excess1 = center - outside / 4;

	Float_t center = pc->Integral(pc->GetXaxis()->FindBin(0), pc->GetXaxis()->FindBin(0.001999), 1, pc->GetNbinsY());
	Float_t outside = pc->Integral(pc->GetXaxis()->FindBin(0.002001), pc->GetXaxis()->FindBin(0.00999), 1, pc->GetNbinsY());

	Float_t excess2 = center - outside / 4;

	Printf("%d %f %f", ptID, excess1, excess2);
      }

      // fit
      if (0 && i == 0)
      {
	p2 = ((TH2*)pc)->ProjectionX("p2", 52, 52+4);
// 	p3 = ((TH2*)p)->ProjectionX("p3", 49-4, 49);
	//p2->Add(p3);
	//new TCanvas; p2->Draw();
	//return;
	p2->Fit("pol0", "0");
	Float_t avg = p2->GetFunction("pol0")->GetParameter(0);
	p2->Fit("pol0", "0", "", -0.002, 0.002); 
	Float_t min = p2->GetFunction("pol0")->GetParameter(0);
	Float_t mine = p2->GetFunction("pol0")->GetParError(0);
	
	if (avg > 0)
	{
	  graph->SetPoint(graph->GetN(), ptRange[ptID], min / avg);
	  graph->SetPointError(graph->GetN()-1, 0, mine / avg);
	}
      }
      
    }
    
    c->SaveAs(Form("twotrack_pt_%d_%d.png", (Int_t) ptRange[ptID], (Int_t) ptRange[ptID+1]));
    c->SaveAs(Form("twotrack_pt_%d_%d.eps", (Int_t) ptRange[ptID], (Int_t) ptRange[ptID+1]));
  }
  
  new TCanvas;
  graph->Print();
  graph->Draw("A*");
}

void PlotTwoTrackEfficiencyControlPlots2(const char* fileName)
{
  loadlibs();
  AliUEHistograms* h = (AliUEHistograms*) GetUEHistogram(fileName);
  
  Float_t ptRange[] = { 0.0, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0 };
  
  graph = new TGraphErrors;
  
  for (Int_t ptID = 0; ptID < 5; ptID++)
  {
    c = new TCanvas(Form("c%d", ptID), Form("%.0f < pT < %.0f", ptRange[ptID], ptRange[ptID+1]), 1200, 400);
    c->Divide(3, 1);
  
    TH2* proj[2];
    
    for (Int_t i=0; i<3; i++)
    {
      c->cd(i+1);
      gPad->SetRightMargin(0.2);
      
      if (i < 2)
      {
	h->GetTwoTrackDistance(i)->GetZaxis()->SetRangeUser(ptRange[ptID] + 0.01, ptRange[ptID+1] - 0.01);
	p = (TH2*) h->GetTwoTrackDistance(i)->Project3D(Form("yx_%d_%d", i, ptID));
	
	if (1)
	{
	  // reduce to one quadrant
	  pc = (TH2*) p->Clone(Form("%s_pc", p->GetName()));
	  pc->Reset();
	  pc->Rebin2D(2, 2);
	  for (Int_t x=1; x<=p->GetNbinsX(); x++)
	    for (Int_t y=1; y<=p->GetNbinsY(); y++)
	      pc->Fill(TMath::Abs(p->GetXaxis()->GetBinCenter(x)), TMath::Abs(p->GetYaxis()->GetBinCenter(y)), p->GetBinContent(x, y));
	  pc->GetXaxis()->SetRangeUser(0, 100);
	  pc->GetYaxis()->SetRangeUser(0, 100);
	  p = pc;
	}
	
	if (ptID == 0)
	  p->SetStats(0);
	
	p->DrawCopy("COLZ");
	
	proj[i] = p;
      }
      else
      {
	proj[0]->Divide(proj[1]);
	proj[0]->SetStats(0);

// 	Float_t scale = proj[1]->Integral(0, proj[1]->GetNbinsX()+1, 0, proj[1]->GetNbinsY()+1) / proj[0]->Integral(0, proj[1]->GetNbinsX()+1, 0,  proj[1]->GetNbinsY()+1);
	Float_t scale = proj[0]->Integral(1, proj[1]->GetNbinsX(), 1, proj[1]->GetNbinsY()) / proj[1]->GetNbinsX() / proj[1]->GetNbinsY();
	proj[0]->Scale(1./ scale / 4);
	proj[0]->DrawCopy("COLZ");
      }
    }
  }
}

void SystematicpTResolution(const char* inputYield, Int_t caseId = 18, Int_t triggerId = 1)
{
  //   Study by Jacek comparing TPC only tracks with global tracks for new cuts (crossed rows) (Fwd by Andrew, 07.07.11)
  //   Resolution from tpc only tracks twice as worse than global tracks
  //   Parameterization for tpc-only tracks:
  //     f(pT) = a * pT * sqrt(1+b/(pT^abs(c)))
  //     a = 0.003; b = 2.08; c = 7.07e-7
  Float_t a = 0.003; Float_t b = 2.08; Float_t c = 7.07e-7;
  res = new TF1("res", "[0] * x * sqrt(1+[1]/(x**abs([2])))", 0, 15);
  res->SetParameters(a, b, c);
//   res->Draw(); return;
  
  ReadYields(inputYield);
  
  for (Int_t side = 0; side < 2; side++)
  {
    for (Int_t centrality = 0; centrality < 4; centrality++)
    {
      TGraphErrors** tmp = yields[side][triggerId][centrality]; 
      graph = tmp[caseId];
      
//       graph->DrawClone("A*");
      
      Float_t axisLimits[20];
      for (Int_t i=0; i<graph->GetN(); i++)
      {
	axisLimits[i] = graph->GetX()[i] - graph->GetEX()[i];
	axisLimits[i+1] = graph->GetX()[i] + graph->GetEX()[i];
      }
      
      hist = new TH1F("hist", "", graph->GetN(), axisLimits);
      
      for (Int_t i=0; i<graph->GetN(); i++)
      {
	gaus = new TF1("gaus", "gaus(0)", 0, 15);
	Float_t sigma = graph->GetX()[i] * res->Eval(graph->GetX()[i]);
	Float_t norm = graph->GetY()[i] / TMath::Sqrt(2 * TMath::Pi()) / sigma;
	gaus->SetParameters(norm, graph->GetX()[i], sigma);
// 	gaus->Draw("SAME");
// 	Printf("%f %f", graph->GetY()[i], gaus->Integral(0, 20));
	
	// fill histogram
	for (Int_t j=1; j<=hist->GetNbinsX(); j++)
	  hist->SetBinContent(j, hist->GetBinContent(j) + gaus->Integral(hist->GetBinLowEdge(j), hist->GetXaxis()->GetBinUpEdge(j)));
      }
      
//       hist->Draw("SAME"); return;
      
      for (Int_t i=0; i<graph->GetN(); i++)
	graph->GetY()[i] = hist->GetBinContent(i+1);

//       graph->SetMarkerColor(2); graph->DrawClone("*SAME");
    }
  }
}

void DrawProcessIDPlot(const char* fileName)
{
  if (gFile)
    gFile->Close();
  TFile::Open(fileName);
  list = (TList*) gFile->Get("PWG4_PhiCorrelations/histosPhiCorrelations");
  ((TH1*) list->FindObject("processIDs"))->Draw();
}
 