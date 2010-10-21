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

void CompareBias(const char* mcFile = "PWG4_JetTasksOutput.root", Int_t region, Int_t ueHist)
{
  loadlibs();

  TFile::Open(mcFile);
  list = (TList*) gFile->Get(objectName);
  AliUEHistograms* h = (AliUEHistograms*) list->FindObject("AliUEHistograms");
  h->SetEtaRange(-0.79, 0.79);
  h->SetPtRange(0.51, 100);
  h->SetCombineMinMax(kTRUE);
  
  const char* axis = "z";
  Float_t ptLeadMin = 0;
  Float_t ptLeadMax = -1;
  
  if (ueHist == 2)
  {
    ptLeadMin = 0.51 + 4;
    ptLeadMax = 1.99 + 4;
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

  if (compareUEHist == 2)
  {
    for (Float_t ptLeadMin = 2.51; ptLeadMin < 4; ptLeadMin += 2)
      DrawRatios(TString(Form("UE %d pT %f", compareUEHist, ptLeadMin)), corrected->GetUEHist(compareUEHist), comparison->GetUEHist(compareUEHist), compareStep, compareRegion, ptLeadMin, ptLeadMin + 1.98);      
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
  
  DrawRatios(esd, testSample, compareStep, compareRegion, compareUEHist);
  
  //CompareEventsTracks(esd, testSample, compareStep, compareRegion, compareUEHist);
}

// function to compare only final step for all regions and distributions

void correctData(const char* fileNameCorrections, const char* fileNameESD, Int_t compareStep = -1, Int_t compareRegion = 2, Int_t compareUEHist = 0)
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
  
  DrawRatios(esd, corr, compareStep, compareRegion, compareUEHist);
}

void ITSTPCEfficiency(const char* fileNameData, Int_t itsTPC = 0)
{
  // its = 0; tpc = 1

  // uncertainty from dN/dpT paper
  Double_t pTBins[] =  {0.0, 0.1, 0.15,  0.2,  0.25,  0.3,   0.35,  0.4,   0.45,  0.5,   0.6,   0.7,   0.8,   0.9,   1.0,   1.5,   2.0,   2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 100.0};
  Float_t effITS[] = {0.,  0.,  0.995, 0.98, 0.986, 0.996, 0.999, 0.997, 0.997, 1.008, 1.005, 1.022, 1.021, 1.021, 1.035, 1.04,  1.03  };  // the last three are the same because i don't have entries
  Float_t effTPC[] = {0.,  0,   1.042, 1.026,1.021, 1.018, 1.015, 1.015, 1.012, 1.012, 1.007, 1.0075,1.006, 1.006, 1.004, 1.004, 1.009 }; // the last bins put as if they were the same

  TH1F* effHist = new TH1F("effHist", "effHist", 39, pTBins);
  for (Int_t i=0; i<39; i++)
  {
    Int_t bin = i;
    if (i > 16)
      bin = 16;
    effHist->SetBinContent(i+1, (itsTPC == 0) ? effITS[bin] : effTPC[bin]);
  }

  new TCanvas; effHist->Draw();

  EffectOfModifiedTrackingEfficiency(fileNameData, effHist, (itsTPC == 0) ? "ITS" : "TPC");
} 


void EffectOfModifiedTrackingEfficiency(const char* fileNameData, TH1F* trackingEff, const char* name = "EffectOfModifiedTrackingEfficiency")
{
  // trackingEff should contain the change in tracking efficiency, i.e. between before and after in the eta-pT plane

  loadlibs();
  
  TFile::Open(fileNameData);
  list = (TList*) gFile->Get(objectName);
  AliUEHistograms* corrected = (AliUEHistograms*) list->FindObject("AliUEHistograms");
  corrected->SetEtaRange(-0.79, 0.79);
  corrected->SetPtRange(0.51, 100);
  corrected->SetCombineMinMax(kTRUE);
  
  AliUEHist* ueHist = corrected->GetUEHist(0);
  Int_t region = 2;
  
  // histogram before
  TH1* before = ueHist->GetUEHist(AliUEHist::kCFStepAll, region);

  // copy histogram
  // the CFStepTriggered step is overwritten here and cannot be used for comparison afterwards anymore
  ueHist->CorrectTracks(AliUEHist::kCFStepAll, AliUEHist::kCFStepTriggered, (TH1*) 0, 0, -1);

  // reapply tracking efficiency
  ueHist->CorrectTracks(AliUEHist::kCFStepTriggered, AliUEHist::kCFStepAll, trackingEff, 1);

  // histogram after
  TH1* after = ueHist->GetUEHist(AliUEHist::kCFStepAll, region);
  
  DrawRatio(before, after, name);
}

void ModifyComposition(const char* fileNameData, const char* fileNameCorrections, Int_t id, Bool_t verbose = kFALSE)
{
  loadlibs();
  
  TFile::Open(fileNameData);
  list = (TList*) gFile->Get(objectName);
  AliUEHistograms* corrected = (AliUEHistograms*) list->FindObject("AliUEHistograms");
  corrected->SetEtaRange(-0.79, 0.79);
  corrected->SetPtRange(0.51, 100);
  corrected->SetCombineMinMax(kTRUE);
  
  file2 = TFile::Open(fileNameCorrections);
  list2 = (TList*) file2->Get(objectName);
  AliUEHistograms* corrections = (AliUEHistograms*) list2->FindObject("AliUEHistograms");
  corrections->SetEtaRange(-0.79, 0.79);
  corrections->SetPtRange(0.51, 100);
  corrections->SetCombineMinMax(kTRUE);
  
  ueHist = (AliUEHist*) corrections->GetUEHist(id);
  
  // copy histogram
  // the CFStepTriggered step is overwritten here and cannot be used for comparison afterwards anymore
  ueHist->CorrectTracks(AliUEHist::kCFStepAll, AliUEHist::kCFStepTriggered, (TH1*) 0, 0);
  
  // histogram before
  TH1* before[3];
  for (Int_t region=0; region<3; region++)
    before[region] = ueHist->GetUEHist(AliUEHist::kCFStepAll, region);
  
  defaultEff = ueHist->GetTrackingEfficiency();
  defaultEffpT = ueHist->GetTrackingEfficiency(1);
  defaultContainer = ueHist->GetTrackHistEfficiency();
  
  c = new TCanvas;
  defaultEffpT->Draw("");
  
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
    ueHist->SetTrackHistEfficiency(newContainer);
    
    // ratio
    modifiedEff = ueHist->GetTrackingEfficiency();
    modifiedEff->Divide(modifiedEff, defaultEff);
    //modifiedEff->Draw("COLZ");
    
    c->cd();
    modifiedEffpT = ueHist->GetTrackingEfficiency(1);
    modifiedEffpT->SetLineColor(i+1);
    modifiedEffpT->Draw("SAME");
    
    // apply change in tracking efficiency
    ueHist->CorrectTracks(AliUEHist::kCFStepTriggered, AliUEHist::kCFStepAll, modifiedEff, 0, 1);
  
    for (Int_t region=0; region<3; region++)
    {
      // histogram after
      TH1* after = ueHist->GetUEHist(AliUEHist::kCFStepAll, region);
      
      if (verbose)
        DrawRatio(before[region], (TH1*) after->Clone(), Form("Region %d Composition %d", region, i));
      
      // ratio is flat, extract deviation
      after->Divide(before[region]);
      after->Fit("pol0", "0");
      Printf("Deviation for region %d case %d is %.2f %%", region, i, 100.0 - 100.0 * after->GetFunction("pol0")->GetParameter(0));
    }
    //return;
  }
}    
