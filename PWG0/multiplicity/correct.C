/* $Id$ */

//
// script to correct the multiplicity spectrum + helpers
//

Bool_t is900GeV = 0;
Bool_t is2360TeV = 0;

// 900 GeV, MC
//const Int_t kBinLimits[]   = { 42, 57, 60 };
//const Int_t kTrustLimits[] = { 26, 42, 54 };

// 2.36 TeV
//const Int_t kBinLimits[]   = { 43, 67, 83 };
//const Int_t kTrustLimits[] = { 33, 57, 73 };

// 7 TeV
const Int_t kBinLimits[]   = { 60, 120, 60 };
const Int_t kTrustLimits[] = { 26, 70, 54 };

void SetTPC()
{
  gSystem->Load("libPWG0base");
  AliMultiplicityCorrection::SetQualityRegions(kFALSE);
}

const char* GetMultLabel(Int_t etaR = -1, Bool_t trueM = kTRUE)
{
	if (etaR == -1)
		etaR = etaRange;
		
	TString tmpStr((trueM) ? "True " : "Measured ");
		
	if (etaR == 4)
	{
		tmpStr += "multiplicity (full phase space)";
	}
	else
		tmpStr += Form("multiplicity in |#eta| < %.1f", (etaR+1)* 0.5);
	return Form("%s", tmpStr.Data());
}

void draw(const char* fileName = "multiplicity.root", const char* folder = "Multiplicity")
{
  loadlibs();

  AliMultiplicityCorrection* mult = new AliMultiplicityCorrection(folder, folder);

  TFile::Open(fileName);
  mult->LoadHistograms();
  mult->DrawHistograms();

  TH2* hist = (TH2*) gROOT->FindObject("fCorrelation2_zy");
  canvas = new TCanvas("c1", "c1", 600, 500);
  canvas->SetTopMargin(0.05);
  hist->SetStats(kFALSE);
  hist->Draw("COLZ");
  hist->SetTitle(";true multiplicity in |#eta| < 1.5;measured multiplicity in |#eta| < 1.5");
  hist->GetYaxis()->SetTitleOffset(1.1);
  gPad->SetRightMargin(0.12);
  gPad->SetLogz();

  canvas->SaveAs("responsematrix.eps");
}

void loadlibs()
{
  gSystem->Load("libTree");
  gSystem->Load("libVMC");

  gSystem->Load("libSTEERBase");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libPWG0base");
}

void LoadAndInitialize(void* multVoid, void* esdVoid, void* multTriggerVoid, Int_t histID, Bool_t fullPhaseSpace, Int_t* geneLimits)
{
  AliMultiplicityCorrection* mult = (AliMultiplicityCorrection*) multVoid;
  AliMultiplicityCorrection* esd = (AliMultiplicityCorrection*) esdVoid;
  AliMultiplicityCorrection* multTrigger = (AliMultiplicityCorrection*) multTriggerVoid;
  
  for (Int_t i=0; i<3; i++)
    geneLimits[i] = kBinLimits[i];
  
  // REBINNING
  if (1)
  {
    // 900 GeV
    if (is900GeV)
    {
      if (1)
      {
        Int_t bins05 = 31;
        Double_t newBinsEta05[] = { -0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5,
                  10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 19.5,
                  20.5, 22.5, 24.5, 26.5, 28.5, 30.5, 34.5, 38.5, 42.5, 50.5, 
                  100.5 };
      }

      if (0)
      {
        Int_t bins05 = 29;
        Double_t newBinsEta05[] = { -0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5,
                10.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 19.5,
                20.5, 22.5, 24.5, 26.5, 28.5, 30.5, 34.5, 38.5, 42.5, 50.5, 
                100.5 };
      }
      
      if (0)
      {
        Int_t bins05 = 25;
        Double_t* newBinsEta05 = new Double_t[bins05+1];
      
        //newBinsEta05[0] = -0.5;
        //newBinsEta05[1] = 0.5;
        //newBinsEta05[2] = 1.5;
        //newBinsEta05[3] = 2.5;
        
        for (Int_t i=0; i<bins05+1; i++)
          newBinsEta05[i] = -0.5 + i*2;
          //newBinsEta05[i] = 4.5 + (i-4)*2;
      }

      Int_t bins10 = 54;
      Double_t newBinsEta10[] = { -0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5,
                                  10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 19.5,
                                  20.5, 21.5, 22.5, 23.5, 24.5, 25.5, 26.5, 27.5, 28.5, 29.5,
                                  30.5, 31.5, 32.5, 33.5, 34.5, 35.5, 36.5, 37.5, 38.5, 39.5,
                                  40.5, 42.5, 44.5, 46.5, 48.5, 50.5, 52.5, 54.5, 56.5, 58.5, 
                                  60.5, 65.5, 70.5, 100.5 };
                        
      geneLimits[0] = bins05;
      geneLimits[1] = bins10;
      geneLimits[2] = bins10;
    }
    
    // 2.36 TeV
    if (is2360TeV)
    {
      Int_t bins05 = 36;
      Double_t newBinsEta05[] = { -0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5,
				  10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 19.5,
				  20.5, 22.5, 24.5, 26.5, 28.5, 30.5, 32.5, 34.5, 36.5, 38.5, 
				  40.5, 45.5, 50.5, 55.5, 60.5, 100.5 };
      Int_t bins10 = 64;
      Double_t newBinsEta10[] = { -0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5,
				  10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 19.5,
				  20.5, 21.5, 22.5, 23.5, 24.5, 25.5, 26.5, 27.5, 28.5, 29.5,
				  30.5, 31.5, 32.5, 33.5, 34.5, 35.5, 36.5, 37.5, 38.5, 39.5,
				  40.5, 42.5, 44.5, 46.5, 48.5, 50.5, 52.5, 54.5, 56.5, 58.5, 
				  60.5, 62.5, 64.5, 66.5, 68.5, 70.5, 72.5, 74.5, 76.5, 78.5, 
				  80.5, 85.5, 90.5, 100.5 };

      geneLimits[0] = bins05;
      geneLimits[1] = bins10;
      geneLimits[2] = bins10;
    }
    
    // 7 TeV
    if (!is900GeV && !is2360TeV)
    {
      Int_t bins05 = 36;
      Double_t newBinsEta05[] = { -0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5,
				  10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 19.5,
				  20.5, 22.5, 24.5, 26.5, 28.5, 30.5, 32.5, 34.5, 36.5, 38.5, 
				  40.5, 45.5, 50.5, 55.5, 60.5, 100.5 };
      
      Int_t bins10 = 85;
      Double_t newBinsEta10[] = { 
          -0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5,
				  10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 19.5,
				  20.5, 21.5, 22.5, 23.5, 24.5, 25.5, 26.5, 27.5, 28.5, 29.5,
				  30.5, 31.5, 32.5, 33.5, 34.5, 35.5, 36.5, 37.5, 38.5, 39.5,
				  40.5, 41.5, 42.5, 43.5, 44.5, 45.5, 46.5, 47.5, 48.5, 49.5, 
				  50.5, 51.5, 52.5, 53.5, 54.5, 55.5, 56.5, 57.5, 58.5, 59.5, 
				  60.5, 62.5, 64.5, 66.5, 68.5, 70.5, 72.5, 74.5, 76.5, 78.5, 
				  80.5, 82.5, 84.5, 86.5, 88.5, 90.5, 92.5, 94.5, 96.5, 98.5, 
				  100.5, 105.5, 110.5, 115.5, 120.5 };
      
      geneLimits[0] = bins05;
      geneLimits[1] = bins10;
      geneLimits[2] = bins10;
    }

    esd->RebinGenerated(bins05, newBinsEta05, bins10, newBinsEta10, bins10, newBinsEta10);
    mult->RebinGenerated(bins05, newBinsEta05, bins10, newBinsEta10, bins10, newBinsEta10);
    multTrigger->RebinGenerated(bins05, newBinsEta05, bins10, newBinsEta10, bins10, newBinsEta10);
  }
  
  for (Int_t hID = ((histID == -1) ? 0 : histID); hID <= ((histID == -1) ? 2 : histID); hID++)
  {
    TH2F* hist = esd->GetMultiplicityESD(hID);
  
    mult->SetMultiplicityESD(hID, hist);
    mult->SetTriggeredEvents(hID, esd->GetTriggeredEvents(hID));
    
    // insert trigger efficiency in flat response matrix
    for (Int_t evType = AliMultiplicityCorrection::kTrVtx; evType <= AliMultiplicityCorrection::kNSD; evType++)
      mult->SetMultiplicityMC(hID, evType, multTrigger->GetMultiplicityMC(hID, evType));
    
    mult->SetNoVertexEvents(hID, multTrigger->GetNoVertexEvents(hID));
    mult->FixTriggerEfficiencies(10);
      
    // small hack to get around charge conservation for full phase space ;-)
    if (fullPhaseSpace)
    {
      TH1* corr = mult->GetCorrelation(histID + 4);
  
      for (Int_t i=2; i<=corr->GetNbinsX(); i+=2)
        for (Int_t j=1; j<=corr->GetNbinsY(); ++j)
        {
          corr->SetBinContent(i, j, corr->GetBinContent(i-1, j));
          corr->SetBinError(i, j, corr->GetBinError(i-1, j));
        }
    }
  }
}

void correct(const char* fileNameMC = "multiplicityMC.root", const char* folder = "Multiplicity", const char* fileNameESD = "multiplicityESD.root", Bool_t chi2 = 1, Int_t histID = 1, Bool_t fullPhaseSpace = kFALSE, Float_t beta = -25, Int_t eventType = 2 /* AliMultiplicityCorrection::kTrVtx */, const char* targetFile = "unfolded.root", const char* folderESD = "Multiplicity", Bool_t calcBias = kFALSE)
{
  loadlibs();
  
  Int_t geneLimits[] = { 0, 0, 0 };
  
  AliMultiplicityCorrection* mult = AliMultiplicityCorrection::Open(fileNameMC, folder);
  AliMultiplicityCorrection* esd = AliMultiplicityCorrection::Open(fileNameESD, folderESD);
  AliMultiplicityCorrection* multTrigger = AliMultiplicityCorrection::Open("multiplicityTrigger.root");

  LoadAndInitialize(mult, esd, multTrigger, histID, fullPhaseSpace, geneLimits);
  
  //AliUnfolding::SetDebug(1);

  for (Int_t hID = ((histID == -1) ? 0 : histID); hID <= ((histID == -1) ? 2 : histID); hID++)
  {
    AliUnfolding::SetNbins(kBinLimits[hID], geneLimits[hID]);
  
    TH2F* hist2 = esd->GetMultiplicityMC(hID, eventType);

    if (chi2)
    {
      if (is900GeV)
      {
        //Float_t regParamPol0[] = { 2, 2, 10 };   // TPCITS
        Float_t regParamPol0[] = { 5, 15, 25 }; // SPD
        Float_t regParamPol1[] =  { 0.15, 0.25, 0.5 };
      }
      else if (is2360TeV)
      {
        Float_t regParamPol0[] = { 10, 12, 40 };
        Float_t regParamPol1[] =  { 0.25, 0.25, 2 };
      }
      else
      {
        Float_t regParamPol0[] = { 1, 25, 10 };
        Float_t regParamPol1[] =  { 0.15, 0.5, 0.5 };
        AliUnfolding::SetSkipBinsBegin(3);
      }
      
      Int_t reg = AliUnfolding::kPol0;
      if (beta > 0)
        reg = AliUnfolding::kPol1;
        
      Float_t regParam = TMath::Abs(beta);
      if (histID == -1)
      {
        if (beta > 0)
          regParam = regParamPol1[hID];
        else
          regParam = regParamPol0[hID];
      }
      AliUnfolding::SetChi2Regularization(reg, regParam);
      
      //AliUnfolding::SetChi2Regularization(AliUnfolding::kLog, 1000000);
      //AliUnfolding::SetChi2Regularization(AliUnfolding::kRatio, 10);
      //TVirtualFitter::SetDefaultFitter("Minuit2");
      
      //AliUnfolding::SetCreateOverflowBin(kFALSE);
      //mult->SetRegularizationParameters(AliMultiplicityCorrection::kNone, 0); //mult->SetCreateBigBin(kFALSE);
      //mult->SetRegularizationParameters(AliMultiplicityCorrection::kNone, 0, 125); mult->SetCreateBigBin(kFALSE);
     // AliUnfolding::SetChi2Regularization(AliUnfolding::kEntropy, beta);
      //mult->SetRegularizationParameters(AliMultiplicityCorrection::kLog, 1e5);
      
      if (0)
      {
        // part for checking
        TH1* mcCompare = hist2->ProjectionY("mcmchist", 1, hist2->GetNbinsX());
        mcCompare->Sumw2();
        mult->ApplyMinuitFit(hID, fullPhaseSpace, AliMultiplicityCorrection::kMB, 0, kTRUE, mcCompare);
        mult->SetMultiplicityESDCorrected(hID, (TH1F*) mcCompare);
      }
      else
      {
        // Zero Bin
        Int_t zeroBin = 0;
        if (is900GeV) // from data
        {
          // background subtraction
          Int_t background = 0;
          
          //background = 1091 + 4398; // MB1 for run 104824 - 52 (SPD)
          //background = 1087 + 4398; // MB1 for run 104824 - 52 (TPCITS)
          
          //background = 417 + 1758; // MB1 for run 104867 - 92 (SPD)
          //background = 1162+422;     // MB1 for run 104892 (TPCITS)
          //background = 5830 + 1384;    // MB1 for run 104824,25,45,52,67,92 (TPC runs) (TPCITS)
          
          background = 20 + 0; // V0AND for run 104824 - 52
          //background = 10 + 0; // V0AND for run 104867 - 92
          
          Printf("NOTE: Subtracting %d background events", background);
          gSystem->Sleep(1000);
      
          zeroBin = mult->GetTriggeredEvents(hID)->GetBinContent(1) - background - mult->GetMultiplicityESD(hID)->Integral(0, 2, 1, 1);
        }
        else if (is2360TeV)
        {
          // from MC
          Float_t fractionZeroBin = (multTrigger->GetTriggeredEvents(hID)->GetBinContent(1) - multTrigger->GetMultiplicityESD(hID)->Integral(0, 2, 1, 1)) / multTrigger->GetMultiplicityESD(hID)->Integral(1, 1, 1, mult->GetMultiplicityESD(hID)->GetNbinsY());
          zeroBin = fractionZeroBin * mult->GetMultiplicityESD(hID)->Integral(1, 1, 1, mult->GetMultiplicityESD(hID)->GetNbinsY());
          
          Printf("Zero bin from MC: Estimating %d events with trigger but without vertex", zeroBin);
          gSystem->Sleep(1000);
        }
        else
        {
          AliUnfolding::SetSkip0BinInChi2(kTRUE);
        }
      
        //mult->SetVertexRange(3, 4);
        mult->ApplyMinuitFit(hID, fullPhaseSpace, eventType, zeroBin, kFALSE, 0, calcBias);
      }
    }
    else
    {
      // HACK
      //mult->GetMultiplicityESD(hID)->SetBinContent(1, 1, 0);
      //for (Int_t bin=1; bin<=mult->GetCorrelation(hID)->GetNbinsY(); bin++)
      //  mult->GetCorrelation(hID)->SetBinContent(1, bin, 1, 0);
      AliUnfolding::SetChi2Regularization(AliUnfolding::kNone, 0);
      mult->ApplyBayesianMethod(hID, fullPhaseSpace, eventType, 1, beta, 0, kTRUE);
      //mult->ApplyBayesianMethod(hID, fullPhaseSpace, eventType, 1, 5, 0, kFALSE);
    }
  }
  
  Printf("Writing result to %s", targetFile);
  TFile* file = TFile::Open(targetFile, "RECREATE");
  mult->SaveHistograms("Multiplicity");
  file->Write();
  file->Close();

  if (histID == -1)
    return;

  for (hID = ((histID == -1) ? 0 : histID); hID <= ((histID == -1) ? 2 : histID); hID++)
  {
    TH2F* hist2 = esd->GetMultiplicityMC(hID, eventType);
    TH1* mcCompare = hist2->ProjectionY("mcmchist", 1, hist2->GetNbinsX());
    mult->DrawComparison(Form("%s_%d_%f", (chi2) ? "MinuitChi2" : "Bayesian", hID, beta), hID, fullPhaseSpace, kTRUE, mcCompare, kFALSE, eventType);

    Printf("<n> = %f", mult->GetMultiplicityESDCorrected(hID)->GetMean());
  }
}

void correctAll(Int_t eventType)
{
  const char* suffix = "";
  switch (eventType)
  {
    case 0: suffix = "trvtx"; break;
    case 1: suffix = "mb"; break;
    case 2: suffix = "inel"; break;
    case 3: suffix = "nsd"; break;
  }

  correct("multiplicityMC.root", "Multiplicity", "multiplicityESD.root", 1, -1, 0, -1, eventType, TString(Form("chi2a_%s.root", suffix)));
  correct("multiplicityMC.root", "Multiplicity", "multiplicityESD.root", 1, -1, 0,  1, eventType, TString(Form("chi2b_%s.root", suffix)));
  correct("multiplicityMC.root", "Multiplicity", "multiplicityESD.root", 0, -1, 0, 40, eventType, TString(Form("bayes_%s.root", suffix)));
  
  if (eventType == 3)
    drawAll(1);
  else if (eventType == 2)
    drawAllINEL();
}

void drawAll(Bool_t showUA5 = kFALSE)
{
  const char* files[] = { "chi2a_nsd.root", "chi2b_nsd.root", "bayes_nsd.root" };
  drawAll(files, showUA5);
}

void drawAllINEL()
{
  const char* files[] = { "chi2a_inel.root", "chi2b_inel.root", "bayes_inel.root" };
  drawAll(files);
}

void drawAllMB()
{
  const char* files[] = { "chi2a_mb.root", "chi2b_mb.root", "bayes_mb.root" };
  drawAll(files);
}

void drawAll(const char** files, Bool_t showUA5 = kFALSE, Bool_t normalize = kFALSE)
{
  loadlibs();
  
  Int_t colors[] = { 1, 2, 4 };
  
  c = new TCanvas(Form("c%d", gRandom->Uniform(100)), "c", 1800, 600);
  c->Divide(3, 1);
  
  l = new TLegend(0.6, 0.6, 0.99, 0.99);
  l->SetFillColor(0);
  
  TH1* hist0 = 0;
  
  for (Int_t hID=0; hID<3; hID++)
  {
    c->cd(hID+1)->SetLogy();
    gPad->SetGridx();
    gPad->SetGridy();
    for (Int_t i=0; i<3; i++)
    {
      TFile::Open(files[i]);
      
      hist = (TH1*) gFile->Get(Form("Multiplicity/fMultiplicityESDCorrected%d", hID));

      Float_t average0 = hist->GetMean();
      hist2 = (TH1*) hist->Clone("temp");
      hist2->SetBinContent(1, 0);
      Float_t average1 = hist2->GetMean();
      Printf("%d: <N> = %.2f <N>/(eta) = %.2f | without 0: <N> = %.2f <N>/(eta) = %.2f", hID, average0, average0 / ((hID+1) - 0.4 * (hID / 2)), average1, average1 / ((hID+1) - 0.4 * (hID / 2)));
      
      hist->SetLineColor(colors[i]);
      
      hist->SetStats(0);
      hist->GetXaxis()->SetRangeUser(0, kBinLimits[hID]);
      
      Float_t min = 0.1;
      Float_t max = hist->GetMaximum() * 1.5;
      
      if (normalize)
        min = 1e-6;
      
      hist->GetYaxis()->SetRangeUser(min, max);
      hist->SetTitle(";unfolded multiplicity;events");
      
      if (hID == 0)
      {
        l->AddEntry(hist, files[i]);
      }
      
      if (hist->Integral() <= 0)
        continue;
      
      if (normalize)
        hist->Scale(1.0 / hist->Integral());
      
      AliPWG0Helper::NormalizeToBinWidth(hist);
      
      hist->DrawCopy((i == 0) ? "" : "SAME");
      
      if (!hist0)
        hist0 = (TH1*) hist->Clone("hist0");
    }
      
    if (hist0)
    {
      line = new TLine(kTrustLimits[hID], hist0->GetMinimum(), kTrustLimits[hID], hist0->GetMaximum());
      line->SetLineWidth(3);
      line->Draw();
    }
  }
  
  c->cd(1);
  l->Draw();
  
  if (showUA5)
  {
    TGraphErrors *gre = new TGraphErrors(23);
    gre->SetName("Graph");
    gre->SetTitle("UA5");
    gre->SetFillColor(1);
    gre->SetMarkerStyle(24);
    gre->SetPoint(0,0,0.15);
    gre->SetPointError(0,0.5,0.01);
    gre->SetPoint(1,1,0.171);
    gre->SetPointError(1,0.5,0.007);
    gre->SetPoint(2,2,0.153);
    gre->SetPointError(2,0.5,0.007);
    gre->SetPoint(3,3,0.124);
    gre->SetPointError(3,0.5,0.006);
    gre->SetPoint(4,4,0.099);
    gre->SetPointError(4,0.5,0.005);
    gre->SetPoint(5,5,0.076);
    gre->SetPointError(5,0.5,0.005);
    gre->SetPoint(6,6,0.057);
    gre->SetPointError(6,0.5,0.004);
    gre->SetPoint(7,7,0.043);
    gre->SetPointError(7,0.5,0.004);
    gre->SetPoint(8,8,0.032);
    gre->SetPointError(8,0.5,0.003);
    gre->SetPoint(9,9,0.024);
    gre->SetPointError(9,0.5,0.003);
    gre->SetPoint(10,10,0.018);
    gre->SetPointError(10,0.5,0.002);
    gre->SetPoint(11,11,0.013);
    gre->SetPointError(11,0.5,0.002);
    gre->SetPoint(12,12,0.01);
    gre->SetPointError(12,0.5,0.002);
    gre->SetPoint(13,13,0.007);
    gre->SetPointError(13,0.5,0.001);
    gre->SetPoint(14,14,0.005);
    gre->SetPointError(14,0.5,0.001);
    gre->SetPoint(15,15,0.004);
    gre->SetPointError(15,0.5,0.001);
    gre->SetPoint(16,16,0.0027);
    gre->SetPointError(16,0.5,0.0009);
    gre->SetPoint(17,17,0.0021);
    gre->SetPointError(17,0.5,0.0008);
    gre->SetPoint(18,18,0.0015);
    gre->SetPointError(18,0.5,0.0007);
    gre->SetPoint(19,19,0.0011);
    gre->SetPointError(19,0.5,0.0006);
    gre->SetPoint(20,20,0.0008);
    gre->SetPointError(20,0.5,0.0005);
    gre->SetPoint(21,22,0.00065);
    gre->SetPointError(21,1,0.0003);
    gre->SetPoint(22,27.5,0.00013);
    gre->SetPointError(22,3.5,7.1e-05);
    
    for (Int_t i=0; i<gre->GetN(); i++)
    {
      gre->GetY()[i] *= hist0->Integral();
      gre->GetEY()[i] *= hist0->Integral();
    }
    
    gre->Draw("P");
    
    ratio = (TGraphErrors*) gre->Clone("ratio");
    
    for (Int_t i=0; i<gre->GetN(); i++)
    {
      //Printf("%f %d", hist0->GetBinContent(hist0->FindBin(ratio->GetX()[i])), hist0->FindBin(ratio->GetX()[i]));
      Int_t bin = hist0->FindBin(gre->GetX()[i]);
      if (hist0->GetBinContent(bin) > 0)
      {
        ratio->GetEY()[i] = TMath::Sqrt((ratio->GetEY()[i] / ratio->GetY()[i])**2 + (hist0->GetBinError(bin) / hist0->GetBinContent(bin))**2);
        ratio->GetY()[i] /= hist0->GetBinContent(bin);
        ratio->GetEY()[i] *= ratio->GetY()[i];
      }
    }
    new TCanvas;
    gPad->SetGridx();
    gPad->SetGridy();
    //ratio->Print();
    ratio->Draw("AP");
    ratio->GetYaxis()->SetRangeUser(0.5, 1.5);
  }
  
  c->SaveAs("draw_all.png");
}

TH1* GetChi2Bias(Float_t alpha = -5)
{
  loadlibs();

  AliMultiplicityCorrection::EventType eventType = AliMultiplicityCorrection::kMB;
  Int_t hID = 1;
  
  correct("multiplicityMC.root", "Multiplicity", "multiplicityESD.root", 1, hID, 0, alpha, eventType, "nobias.root", "Multiplicity", kFALSE);

  correct("multiplicityMC.root", "Multiplicity", "multiplicityESD.root", 1, hID, 0, alpha, eventType, "withbias.root", "Multiplicity", kTRUE);

  // without bias calculation
  mult = AliMultiplicityCorrection::Open("nobias.root");
  baseold = (TH1*) mult->GetMultiplicityESDCorrected(hID)->Clone("baseold");
  
  // with bias calculation
  mult = AliMultiplicityCorrection::Open("withbias.root");
  base = (TH1*) mult->GetMultiplicityESDCorrected(hID)->Clone("base");
 
  // relative error plots
  ratio1 = (TH1*) mult->GetMultiplicityESDCorrected(hID)->Clone("hist1");
  ratio1->SetStats(0);
  ratio1->SetTitle(";unfolded multiplicity;relative error");
  ratio1->Reset();
  ratio2 = (TH1*) mult->GetMultiplicityESDCorrected(hID)->Clone("hist1");
  ratio2->Reset();

  for (Int_t t = 0; t<baseold->GetNbinsX(); t++)
  {
    Printf("Bin %d; Content: %f; Chi2 Error: %f; Bias: %f; In percent: %.2f %.2f", t+1, base->GetBinContent(t+1), baseold->GetBinError(t+1), base->GetBinError(t+1), (base->GetBinContent(t+1) > 0) ? 100.0 * baseold->GetBinError(t+1) / base->GetBinContent(t+1) : -1, (base->GetBinContent(t+1) > 0) ? 100.0 * base->GetBinError(t+1) / base->GetBinContent(t+1) : -1);
    if (base->GetBinContent(t+1) <= 0)
      continue;
    ratio1->SetBinContent(t+1, baseold->GetBinError(t+1) / base->GetBinContent(t+1));
    ratio2->SetBinContent(t+1, base->GetBinError(t+1) / base->GetBinContent(t+1));
  }
  
  //new TCanvas; base->Draw(); gPad->SetLogy();

  c = new TCanvas;
  ratio1->GetYaxis()->SetRangeUser(0, 1);
  ratio1->GetXaxis()->SetRangeUser(0, 50);
  ratio1->Draw();
  ratio2->SetLineColor(2);
  ratio2->Draw("SAME");
  
  return base;
}

void CheckBayesianBias()
{
  loadlibs();
  
  const char* fileNameMC = "multiplicityMC.root";
  const char* folder = "Multiplicity";
  const char* fileNameESD = "multiplicityESD.root";
  const char* folderESD = "Multiplicity";

  AliMultiplicityCorrection* mult = AliMultiplicityCorrection::Open(fileNameMC, folder);
  AliMultiplicityCorrection* esd = AliMultiplicityCorrection::Open(fileNameESD, folderESD);
  
  const Int_t kMaxBins = 35;
  AliUnfolding::SetNbins(kMaxBins, kMaxBins);
  //AliUnfolding::SetDebug(1);
  
  Int_t hID = 0;
  
  TH2F* hist = esd->GetMultiplicityESD(hID);
  mult->SetMultiplicityESD(hID, hist);  
  
  // without bias calculation
  mult->ApplyBayesianMethod(hID, kFALSE, 0, 1, 10, 0, 1);
  baseold = (TH1*) mult->GetMultiplicityESDCorrected(hID)->Clone("baseold");
  
  // with bias calculation
  mult->ApplyBayesianMethod(hID, kFALSE, 0, 1, 10, 0, 2);
  base = (TH1*) mult->GetMultiplicityESDCorrected(hID)->Clone("base");
  
  TH1* hist1 = new TH1F("hist1", "", 100, 0, 20);
  TH1* hist2 = new TH1F("hist2", "", 100, 0, 20);
  
  for (Int_t t = 0; t<kMaxBins; t++)
  {
    Printf("Bin %d; Content: %f; Randomization Error: %f; Bias: %f; In percent: %.2f %.2f", t+1, base->GetBinContent(t+1), baseold->GetBinError(t+1), base->GetBinError(t+1), (base->GetBinContent(t+1) > 0) ? 100.0 * baseold->GetBinError(t+1) / base->GetBinContent(t+1) : -1, (base->GetBinContent(t+1) > 0) ? 100.0 * base->GetBinError(t+1) / base->GetBinContent(t+1) : -1);
    
    hist1->Fill((base->GetBinContent(t+1) > 0) ? 100.0 * baseold->GetBinError(t+1) / base->GetBinContent(t+1) : -1);
    hist2->Fill((base->GetBinContent(t+1) > 0) ? 100.0 * base->GetBinError(t+1) / base->GetBinContent(t+1) : -1);
  }
  
  new TCanvas;
  hist1->Draw();
  hist2->SetLineColor(2);
  hist2->Draw("SAME");
}

void DrawUnfoldingLimit(Int_t hID, Float_t min, Float_t max)
{
  line = new TLine(kTrustLimits[hID], min, kTrustLimits[hID], max);
  line->SetLineWidth(2);
  line->Draw();
}

void DataScan(Int_t hID, Bool_t redo = kTRUE)
{
  // makes a set of unfolded spectra and compares
  // don't forget FindUnfoldedLimit in plots.C

  loadlibs();

  // files...
  Bool_t mc = 1;
  const char* fileNameMC = "multiplicityMC.root";
  const char* folder = "Multiplicity";
  const char* fileNameESD = "multiplicityESD.root";
  const char* folderESD = "Multiplicity";
  const Int_t kMaxBins = kBinLimits[hID];
  
  AliMultiplicityCorrection::EventType eventType = AliMultiplicityCorrection::kTrVtx;
  Bool_t evaluteBias = kFALSE;
  
  Int_t referenceCase = 2; // all results are shown as ratio to this case
  
  // chi2 range
  AliUnfolding::RegularizationType regBegin = AliUnfolding::kPol0;
  AliUnfolding::RegularizationType regEnd = AliUnfolding::kPol1;
  Float_t regParamBegin[] = { 0, 1, 0.2, 3 }; 
  Float_t regParamEnd[] = { 0, 11, 0.5, 31 }; 
  Float_t regParamStep[] = { 0, 2, 2, TMath::Sqrt(10) }; 

  // bayesian range
  Int_t iterBegin = 40;
  Int_t iterEnd = 41;
  Int_t iterStep = 20;
  
  TList labels;
  
  AliMultiplicityCorrection* mult = AliMultiplicityCorrection::Open(fileNameMC, folder);
  AliMultiplicityCorrection* esd = AliMultiplicityCorrection::Open(fileNameESD, folderESD);
  
  mult->SetMultiplicityESD(hID, esd->GetMultiplicityESD(hID));
  
  // insert trigger efficiency in flat response matrix
  AliMultiplicityCorrection* multTrigger = AliMultiplicityCorrection::Open("multiplicityTrigger.root");
  for (Int_t evType = AliMultiplicityCorrection::kTrVtx; evType <= AliMultiplicityCorrection::kNSD; evType++)
    mult->SetMultiplicityMC(hID, evType, multTrigger->GetMultiplicityMC(hID, evType));
  
  mult->SetNoVertexEvents(hID, multTrigger->GetNoVertexEvents(hID));
  mult->FixTriggerEfficiencies(10);
  
  Float_t fraction0Generated = multTrigger->GetFraction0Generated(hID);
  
  if (mc)
    mult->SetMultiplicityMC(hID, eventType, esd->GetMultiplicityMC(hID, eventType));
  
  AliUnfolding::SetNbins(kMaxBins, kMaxBins);
  //AliUnfolding::SetDebug(1);
  
  Int_t count = -1;
  
  for (AliUnfolding::RegularizationType reg = regBegin; reg <= regEnd; reg++)
  {
    for (Float_t regParam = regParamBegin[reg]; regParam < regParamEnd[reg]; regParam *= regParamStep[reg])
    {
      count++;
      labels.Add(new TObjString(Form("#chi^{2} Reg %d #beta = %.2f", (Int_t) reg, regParam)));
      
      if (!redo)
        continue;
    
      AliUnfolding::SetChi2Regularization(reg, regParam);
      
      mult->ApplyMinuitFit(hID, kFALSE, eventType, fraction0Generated, kFALSE, 0, evaluteBias);
      
      file = TFile::Open(Form("datascan_%d.root", count), "RECREATE");
      mult->SaveHistograms();
      file->Close();
    }
  }
  
  for (Int_t iter = iterBegin; iter <= iterEnd; iter += iterStep)
  {
    count++;
    labels.Add(new TObjString(Form("Bayesian Iter = %d", iter)));
    
    if (!redo)
      continue;
    
    mult->ApplyBayesianMethod(hID, kFALSE, eventType, 1, iter, 0, kTRUE);
    
    file = TFile::Open(Form("datascan_%d.root", count), "RECREATE");
    mult->SaveHistograms();
    file->Close();
  }
      
  // 1) all distributions
  // 2) ratio to MC
  // 3) ratio to ref point
  // 4) relative error
  // 5) residuals
  c = new TCanvas("c", "c", 1200, 800);
  c->Divide(3, 2);
  
  c->cd(1)->SetLogy();
  
  // get reference point
  mult = AliMultiplicityCorrection::Open(Form("datascan_%d.root", referenceCase));
  if (!mult)
    return;
  refPoint = (TH1*) mult->GetMultiplicityESDCorrected(hID)->Clone("refPoint");
  
  legend = new TLegend(0.5, 0.5, 0.99, 0.99);
  legend->SetFillColor(0);
  
  TH1* residualSum = new TH1F("residualSum", ";;residual squared sum", count + 1, 0.5, count + 1.5);
  
  Int_t allowedCount = count;
  count = 0;
  while (1)
  {
    mult = AliMultiplicityCorrection::Open(Form("datascan_%d.root", count++));
    if (!mult)
      break;
    if (count > allowedCount+1)
      break;
      
    hist = (TH1*) mult->GetMultiplicityESDCorrected(hID);
    hist->SetLineColor((count-1) % 4 + 1);
    hist->SetLineStyle((count-1) / 4 + 1);
    hist->GetXaxis()->SetRangeUser(0, kMaxBins);
    hist->GetYaxis()->SetRangeUser(0.1, hist->GetMaximum() * 1.5);
    hist->SetStats(kFALSE);
    hist->SetTitle("");
    
    legend->AddEntry(hist->Clone(), labels.At(count-1)->GetName());
    
    c->cd(1);
    hist->DrawCopy((count == 1) ? "HIST" : "HISTSAME");
    
    DrawUnfoldingLimit(hID, 0.1, hist->GetMaximum());
    
    if (mc)
    {
      TH2* mcHist2d = mult->GetMultiplicityMC(hID, eventType);
      mcHist = mcHist2d->ProjectionY("mcmchist", 1, mcHist2d->GetNbinsX());
      
      c->cd(1);
      mcHist->SetMarkerStyle(24);
      mcHist->Draw("P SAME");
      
      c->cd(2);
      // calculate ratio using only the error on the mc bin
      ratio = (TH1*) hist->Clone("ratio");
      for (Int_t bin=1; bin<=ratio->GetNbinsX(); bin++)
      {
        if (mcHist->GetBinContent(bin) <= 0)
          continue;
        ratio->SetBinContent(bin, hist->GetBinContent(bin) / mcHist->GetBinContent(bin));
        ratio->SetBinError(bin, mcHist->GetBinError(bin) / mcHist->GetBinContent(bin) * ratio->GetBinContent(bin));
      }
      
      ratio->GetYaxis()->SetRangeUser(0.5, 1.5);
      ratio->GetYaxis()->SetTitle("Ratio unfolded / MC");
      ratio->DrawCopy((count == 1) ? "HIST" : "HISTSAME");
    
      DrawUnfoldingLimit(hID, 0.5, 1.5);
    }
    
    c->cd(3);
    // calculate ratio using no errors for now
    ratio = (TH1*) hist->Clone("ratio");
    for (Int_t bin=1; bin<=ratio->GetNbinsX(); bin++)
    {
      if (refPoint->GetBinContent(bin) <= 0)
        continue;
      ratio->SetBinContent(bin, hist->GetBinContent(bin) / refPoint->GetBinContent(bin));
      ratio->SetBinError(bin, 0);
    }
    
    ratio->GetYaxis()->SetRangeUser(0.5, 1.5);
    ratio->GetYaxis()->SetTitle("Ratio unfolded / unfolded reference case");
    ratio->DrawCopy((count == 1) ? "" : "SAME");
    
    DrawUnfoldingLimit(hID, 0.5, 1.5);
    
    c->cd(4);
    ratio = (TH1*) hist->Clone("ratio");
    for (Int_t bin=1; bin<=ratio->GetNbinsX(); bin++)
    {
      if (ratio->GetBinContent(bin) <= 0)
        continue;
      ratio->SetBinContent(bin, ratio->GetBinError(bin) / ratio->GetBinContent(bin));
      ratio->SetBinError(bin, 0);
    }
    ratio->GetYaxis()->SetRangeUser(0, 1);
    ratio->GetYaxis()->SetTitle("Relative error");
    ratio->DrawCopy((count == 1) ? "" : "SAME");
    
    DrawUnfoldingLimit(hID, 0, 1);
    
    c->cd(5);
    Float_t sum;
    residuals = mult->GetResiduals(hID, AliMultiplicityCorrection::kTrVtx, sum);
    residuals->SetStats(0);
    residuals->GetXaxis()->SetRangeUser(0, kMaxBins);
    residuals->SetStats(kFALSE);
    residuals->SetLineColor((count-1) % 4 + 1);
    residuals->SetLineStyle((count-1) / 4 + 1);
    residuals->GetYaxis()->SetRangeUser(-5, 5);
    residuals->DrawCopy((count == 1) ? "" : "SAME");
    
    residualSum->Fill(count, sum);
    residualSum->GetXaxis()->SetBinLabel(count, labels.At(count-1)->GetName());
  }
  
  c->cd(6);
  residualSum->Draw();
  legend->Draw();
}

void DrawBayesianIterations(const char* fileNameMC = "multiplicityMC.root", const char* fileNameESD = "multiplicityESD.root", Int_t histID = 1, Int_t eventType = 0 /* AliMultiplicityCorrection::kTrVtx */)
{
  loadlibs();
  
  const char* folder = "Multiplicity";

  AliMultiplicityCorrection* mult = new AliMultiplicityCorrection(folder, folder);
  TFile::Open(fileNameMC);
  mult->LoadHistograms();

  AliMultiplicityCorrection* esd = new AliMultiplicityCorrection(folder, folder);
  TFile::Open(fileNameESD);
  esd->LoadHistograms();

  TH2F* hist = esd->GetMultiplicityESD(histID);
  TH2F* hist2 = esd->GetMultiplicityMC(histID, eventType);
  
  mult->SetMultiplicityESD(histID, hist);
  mult->SetMultiplicityMC(histID, eventType, hist2);
  
  TH1* mcCompare = hist2->ProjectionY("mcmchist", 1, hist2->GetNbinsX());
  //mcCompare->Scale(1.0 / mcCompare->Integral());

  TH1* esdHist = (TH1*) hist->ProjectionY("myesd", 1, 1)->Clone("myesd");
  //esdHist->Scale(1.0 / esdHist->Integral());
  
  c = new TCanvas("c", "c", 1200, 600);
  c->Divide(2, 1);
  
  c->cd(1);
  gPad->SetLeftMargin(0.12);
  gPad->SetTopMargin(0.05);
  gPad->SetRightMargin(0.05);
  gPad->SetLogy();
  gPad->SetGridx();
  gPad->SetGridy();
    
  mcCompare->GetXaxis()->SetRangeUser(0, 80);
  mcCompare->SetStats(0);
  mcCompare->SetFillColor(5);
  mcCompare->SetLineColor(5);
  mcCompare->SetTitle(Form(";%s;Entries", GetMultLabel(1)));
  mcCompare->GetYaxis()->SetRangeUser(2, esdHist->GetMaximum() * 2);
  mcCompare->GetYaxis()->SetTitleOffset(1.3);
  mcCompare->Draw("HIST");
  esdHist->SetMarkerStyle(5);
  esdHist->Draw("P HIST SAME");
  
  c->cd(2);
  gPad->SetTopMargin(0.05);
  gPad->SetRightMargin(0.05);
  gPad->SetGridx();
  gPad->SetGridy();
  
  trueMeasuredRatio = (TH1*) mcCompare->Clone("trueMeasuredRatio");
  trueMeasuredRatio->Divide(esdHist);
  trueMeasuredRatio->SetStats(0);
  trueMeasuredRatio->SetTitle(Form(";%s;Ratio", GetMultLabel(1)));
  trueMeasuredRatio->GetYaxis()->SetTitleOffset(1.3);
  trueMeasuredRatio->GetYaxis()->SetRangeUser(0, 2);
  // ROOT displays all values > 2 at 2 which looks weird
  for (Int_t i=1; i<=trueMeasuredRatio->GetNbinsX(); i++)
    if (trueMeasuredRatio->GetBinContent(i) > 2)
      trueMeasuredRatio->SetBinContent(i, 0);
  trueMeasuredRatio->SetMarkerStyle(5);
  trueMeasuredRatio->Draw("P");

  Int_t colors[] = {1, 2, 4, 6};
  
  legend = new TLegend(0.15, 0.13, 0.5, 0.5);
  legend->SetFillColor(0);
  legend->SetTextSize(0.04);
  
  legend->AddEntry(mcCompare, "True", "F");
  legend->AddEntry(esdHist, "Measured", "P");

  legend2 = new TLegend(0.15, 0.13, 0.5, 0.4);
  legend2->SetFillColor(0);
  legend2->SetTextSize(0.04);
    
  legend2->AddEntry(trueMeasuredRatio, "Measured", "P");
  
  Int_t iters[] = {1, 3, 10, -1};
  for (Int_t i = 0; i<4; i++)
  {
    Int_t iter = iters[i];
    mult->ApplyBayesianMethod(histID, kFALSE, eventType, 1, iter, 0, 0);
    corr = mult->GetMultiplicityESDCorrected(histID);
    corr->Scale(1.0 / corr->Integral());
    corr->Scale(mcCompare->Integral());
    corr->GetXaxis()->SetRangeUser(0, 80);
    //corr->SetMarkerStyle(20+iter);
    corr->SetLineColor(colors[i]);
    corr->SetLineStyle(i+1);
    corr->SetLineWidth(2);
    
    c->cd(1);
    corr->DrawCopy("SAME HIST");
    
    c->cd(2);
    clone = (TH1*) corr->Clone("clone");
    clone->Divide(mcCompare, corr);
    clone->GetYaxis()->SetRangeUser(0, 2);
    clone->DrawCopy("SAME HIST");
    
    legend->AddEntry((TH1*) corr->Clone(), (iter > -1) ? Form("%d iteration%s", iter, (iter == 1) ? "" : "s") : "Convergence", "L");
    legend2->AddEntry((TH1*) corr->Clone(), (iter > -1) ? Form("%d iteration%s", iter, (iter == 1) ? "" : "s") : "Convergence", "L");
  }
  
  c->cd(1);
  legend->Draw();
  
  c->cd(2);
  legend2->Draw();
  
  c->SaveAs("bayesian_iterations.eps");
}

void CompareChi2Bayesian(Int_t histID = 2, const char* chi2File = "chi2.root", const char* bayesianFile = "bayesian.root", const char* label1 = "Chi2", const char* label2 = "Bayesian", const char* mcFile = 0, Float_t simpleCorrect = 0)
{
  const char* folder = "Multiplicity";

  loadlibs();

  AliMultiplicityCorrection* chi2 = new AliMultiplicityCorrection(folder, folder);
  TFile::Open(chi2File);
  chi2->LoadHistograms();

  AliMultiplicityCorrection* bayesian = new AliMultiplicityCorrection(folder, folder);
  TFile::Open(bayesianFile);
  bayesian->LoadHistograms();

  histRAW = chi2->GetMultiplicityESD(histID)->ProjectionY("raw", 1, chi2->GetMultiplicityESD(histID)->GetNbinsX());
  histRAW->Scale(1.0 / histRAW->Integral());

  histC = chi2->GetMultiplicityESDCorrected(histID);
  histB = bayesian->GetMultiplicityESDCorrected(histID);

  c = new TCanvas("CompareChi2Bayesian", "CompareChi2Bayesian", 800, 600);
  c->SetRightMargin(0.05);
  c->SetTopMargin(0.05);
  c->SetLogy();
  c->SetGridx();
  c->SetGridy();

  histC->SetTitle(";N;P(N)");
  histC->SetStats(kFALSE);
  histC->GetXaxis()->SetRangeUser(0, 100);

  histC->SetLineColor(1);
  histB->SetLineColor(2);
  histRAW->SetLineColor(3);

  histC->DrawCopy("HISTE");
  histB->DrawCopy("HISTE SAME");
  histRAW->DrawCopy("SAME");

  legend = new TLegend(0.2, 0.2, 0.4, 0.4);
  legend->SetFillColor(0);

  legend->AddEntry(histC, label1);
  legend->AddEntry(histB, label2);
  legend->AddEntry(histRAW, "raw ESD");

  histGraph = 0;
  if (simpleCorrect > 0)
  {
    graph = new TGraph;
    graph->SetMarkerStyle(25);
    graph->SetFillColor(0);
    for (Int_t bin=1; bin<=histRAW->GetNbinsX(); bin++)
      graph->SetPoint(graph->GetN(), histRAW->GetXaxis()->GetBinCenter(bin) * simpleCorrect, histRAW->GetBinContent(bin));

    graph->Draw("PSAME");
    legend->AddEntry(graph, "weighting");

    // now create histogram from graph and normalize
    histGraph = (TH1*) histRAW->Clone();
    histGraph->Reset();
    for (Int_t bin=1; bin<=histGraph->GetNbinsX(); bin++)
    {
      Int_t j=1;
      for (j=1; j<graph->GetN(); j++)
        if (graph->GetX()[j] > histGraph->GetXaxis()->GetBinCenter(bin))
	  break;
      if (j == graph->GetN())
        continue;
      if (histGraph->GetXaxis()->GetBinCenter(bin) - graph->GetX()[j] < graph->GetX()[j-1] - histGraph->GetXaxis()->GetBinCenter(bin))
        j--;
      histGraph->SetBinContent(bin, graph->GetY()[j]);
    }

    Printf("Integral = %f", histGraph->Integral());
    histGraph->Scale(1.0 / histGraph->Integral());

    histGraph->SetLineColor(6);
    histGraph->DrawCopy("SAME");
    legend->AddEntry(histGraph, "weighting normalized");
  }

  if (mcFile)
  {
    AliMultiplicityCorrection* mc = new AliMultiplicityCorrection(folder, folder);
    TFile::Open(mcFile);
    mc->LoadHistograms();

    histMC = mc->GetMultiplicityVtx(histID)->ProjectionY("mc", 1, mc->GetMultiplicityVtx(histID)->GetNbinsX());
    histMC->Sumw2();
    histMC->Scale(1.0 / histMC->Integral());

    histMC->Draw("HISTE SAME");
    histMC->SetLineColor(4);
    legend->AddEntry(histMC, "MC");
  }

  legend->Draw();

  c->SaveAs(Form("%s.png", c->GetName()));
  c->SaveAs(Form("%s.eps", c->GetName()));

  if (!mcFile)
    return;

  // build ratios

  c = new TCanvas("CompareChi2BayesianRatio", "CompareChi2BayesianRatio", 800, 600);
  c->SetRightMargin(0.05);
  c->SetTopMargin(0.05);
  c->SetGridx();
  c->SetGridy();

  for (Int_t bin=1; bin<=histC->GetNbinsX(); bin++)
  {
    if (histMC->GetBinContent(bin) > 0)
    {
      histC->SetBinContent(bin, histC->GetBinContent(bin) / histMC->GetBinContent(bin));
      histB->SetBinContent(bin, histB->GetBinContent(bin) / histMC->GetBinContent(bin));

      /*
      if (histGraph)
      {
        histGraph->SetBinContent(bin, histGraph->GetBinContent(bin) / histMC->GetBinContent(bin));
        histGraph->SetBinError(bin, 0);
      }
      */

      // TODO errors?
      histC->SetBinError(bin, 0);
      histB->SetBinError(bin, 0);
    }
  }

  histC->GetYaxis()->SetRangeUser(0.5, 2);
  histC->GetYaxis()->SetTitle("Unfolded / MC");

  histC->Draw("HIST");
  histB->Draw("HIST SAME");
  /*
  if (histGraph)
    histGraph->Draw("HIST SAME");
  */

  /*
  if (simpleCorrect > 0)
  {
    graph2 = new TGraph;
    graph2->SetMarkerStyle(25);
    graph2->SetFillColor(0);
    for (Int_t i=0; i<graph->GetN(); i++)
    {
      Float_t mcValue = histMC->GetBinContent(histMC->FindBin(graph->GetX()[i]));
      Float_t mcError = histMC->GetBinError(histMC->FindBin(graph->GetX()[i]));
      if (mcValue > 0)
        graph2->SetPoint(graph2->GetN(), graph->GetX()[i], graph->GetY()[i] / mcValue);
    }
    graph2->Draw("PSAME");
  }
  */
  
  c->SaveAs(Form("%s.png", c->GetName()));
  c->SaveAs(Form("%s.eps", c->GetName()));
}


void CompareMC(Int_t histID, Int_t eventType, const char* file1, const char* file2, const char* label1, const char* label2)
{
  const char* folder = "Multiplicity";

  loadlibs();

  AliMultiplicityCorrection* mc1 = new AliMultiplicityCorrection(folder, folder);
  TFile::Open(file1);
  mc1->LoadHistograms();
  histMC1 = mc1->GetMultiplicityMC(histID, eventType)->ProjectionY("mc1", 1, mc1->GetMultiplicityMC(histID, eventType)->GetNbinsX());
  histMC1->Scale(1.0 / histMC1->Integral());

  AliMultiplicityCorrection* mc2 = new AliMultiplicityCorrection(folder, folder);
  TFile::Open(file2);
  mc2->LoadHistograms();
  histMC2 = mc2->GetMultiplicityMC(histID, eventType)->ProjectionY("mc2", 1, mc2->GetMultiplicityMC(histID, eventType)->GetNbinsX());
  histMC2->Scale(1.0 / histMC2->Integral());

  c = new TCanvas("CompareMC", "CompareMC", 800, 600);
  c->SetRightMargin(0.05);
  c->SetTopMargin(0.05);
  c->SetLogy();

  histMC1->SetTitle(";N;P(N)");
  histMC1->SetStats(kFALSE);
  histMC1->GetXaxis()->SetRangeUser(0, 100);

  histMC1->SetLineColor(1);
  histMC2->SetLineColor(2);

  histMC1->Draw();
  histMC2->Draw("SAME");

  legend = new TLegend(0.2, 0.2, 0.4, 0.4);
  legend->SetFillColor(0);

  legend->AddEntry(histMC1, label1);
  legend->AddEntry(histMC2, label2);

  legend->Draw();

  c->SaveAs(Form("%s.gif", c->GetName()));
  c->SaveAs(Form("%s.eps", c->GetName()));
}

void* fit2Step(const char* fileNameMC = "multiplicityMC_2M.root", const char* fileNameESD = "multiplicityMC_1M_3.root", Int_t histID = 3, Bool_t fullPhaseSpace = kFALSE)
{
  gSystem->Load("libPWG0base");

  AliMultiplicityCorrection* mult = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");

  TFile::Open(fileNameMC);
  mult->LoadHistograms("Multiplicity");

  TFile::Open(fileNameESD);
  TH2F* hist = (TH2F*) gFile->Get(Form("Multiplicity/fMultiplicityESD%d", histID));
  TH2F* hist2 = (TH2F*) gFile->Get(Form("Multiplicity/fMultiplicityVtx%d", ((fullPhaseSpace) ? 4 : histID)));
  //hist2 = (TH2F*) gFile->Get(Form("Multiplicity/fMultiplicityINEL%d", histID));

  mult->SetMultiplicityESD(histID, hist);

  // small hack to get around charge conservation for full phase space ;-)
  if (fullPhaseSpace)
  {
    TH1* corr = mult->GetCorrelation(histID + 4);

    for (Int_t i=2; i<=corr->GetNbinsX(); i+=2)
      for (Int_t j=1; j<=corr->GetNbinsY(); ++j)
      {
        corr->SetBinContent(i, j, corr->GetBinContent(i-1, j));
        corr->SetBinError(i, j, corr->GetBinError(i-1, j));
      }
  }

  mult->SetRegularizationParameters(AliMultiplicityCorrection::kPol1, 10000);
  mult->ApplyMinuitFit(histID, fullPhaseSpace, AliMultiplicityCorrection::kTrVtx, kFALSE);
  mult->DrawComparison("MinuitChi2", histID, fullPhaseSpace, kTRUE, hist2->ProjectionY("mymchist"));

  TH1* result = (TH1*) mult->GetMultiplicityESDCorrected((fullPhaseSpace) ? 4 : histID))->Clone("firstresult");

  mult->SetRegularizationParameters(AliMultiplicityCorrection::kEntropy, 100000);
  mult->ApplyMinuitFit(histID, fullPhaseSpace, AliMultiplicityCorrection::kTrVtx, kFALSE, result);
  mult->DrawComparison("MinuitChi2_Step2", histID, fullPhaseSpace, kTRUE, hist2->ProjectionY("mymchist"));

  return mult;
}

const char* GetRegName(Int_t type)
{
  switch (type)
  {
    case AliMultiplicityCorrection::kNone:      return "None"; break;
    case AliMultiplicityCorrection::kPol0:      return "Pol0"; break;
    case AliMultiplicityCorrection::kPol1:      return "Pol1"; break;
    case AliMultiplicityCorrection::kCurvature: return "TotalCurvature"; break;
    case AliMultiplicityCorrection::kEntropy:   return "Reduced cross-entropy"; break;
    case AliMultiplicityCorrection::kLog   :    return "Log"; break;
  }
  return 0;
}

const char* GetEventTypeName(Int_t type)
{
  switch (type)
  {
    case AliMultiplicityCorrection::kTrVtx:   return "trigger, vertex"; break;
    case AliMultiplicityCorrection::kMB:      return "minimum bias"; break;
    case AliMultiplicityCorrection::kINEL:    return "inelastic"; break;
  }
  return 0;
}

void EvaluateBayesianMethodIterationsSmoothing(const char* fileNameMC = "multiplicityMC.root", const char* fileNameESD = "multiplicityESD.root", const char* targetDir = "bayesian", Int_t histID = 1)
{
  gSystem->mkdir(targetDir);

  gSystem->Load("libPWG0base");

  AliMultiplicityCorrection* mult = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");
  TFile::Open(fileNameMC);
  mult->LoadHistograms("Multiplicity");

  AliMultiplicityCorrection* multESD = new AliMultiplicityCorrection("MultiplicityESD", "MultiplicityESD");
  TFile::Open(fileNameESD);
  multESD->LoadHistograms("Multiplicity");
  mult->SetMultiplicityESD(histID, multESD->GetMultiplicityESD(histID));

  Int_t count = 0; // just to order the saved images...

  TFile* graphFile = TFile::Open(Form("%s/EvaluateBayesianMethodIterationsSmoothing.root", targetDir), "RECREATE");

  Int_t colors[] =  {1, 2, 3, 4, 6};
  Int_t markers[] = {20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 3, 4, 5, 6};

  for (AliMultiplicityCorrection::EventType type = AliMultiplicityCorrection::kTrVtx; type <= AliMultiplicityCorrection::kTrVtx; ++type)
  //for (AliMultiplicityCorrection::EventType type = AliMultiplicityCorrection::kTrVtx; type <= AliMultiplicityCorrection::kINEL; ++type)
  {
    TString tmp;
    tmp.Form("EvaluateBayesianMethodIterationsSmoothing_%s", GetEventTypeName(type));

    TCanvas* canvas = new TCanvas(tmp, tmp, 800, 600);

    for (Int_t i = 1; i <= 5; i++)
    {
      Int_t iterArray[5] = {2, 5, 10, 20, -1};
      Int_t iter = iterArray[i-1];

      TGraph* fitResultsMC[3];
      for (Int_t region=0; region<AliMultiplicityCorrection::kQualityRegions; ++region)
      {
        fitResultsMC[region] = new TGraph;
        fitResultsMC[region]->SetTitle(Form("%d iterations - reg. %d", iter, region+1));
        fitResultsMC[region]->GetXaxis()->SetTitle("smoothing parameter #alpha");
        fitResultsMC[region]->GetYaxis()->SetTitle(Form("P_{1} in region %d", region));
        fitResultsMC[region]->SetName(Form("%s_MC_%d", tmp.Data(), i * AliMultiplicityCorrection::kQualityRegions + region - 2));
        fitResultsMC[region]->SetFillColor(0);
        //fitResultsMC[region]->SetMarkerStyle(markers[(i-1) * AliMultiplicityCorrection::kQualityRegions + region]);
        fitResultsMC[region]->SetMarkerStyle(markers[i-1]);
        fitResultsMC[region]->SetLineColor(colors[i-1]);
        fitResultsMC[region]->SetMarkerColor(colors[i-1]);
      }

      TGraph* fitResultsRes = new TGraph;
      fitResultsRes->SetTitle(Form("%d iterations", iter));
      fitResultsRes->SetName(Form("%s_Res_%d", tmp.Data(), i));
      fitResultsRes->GetXaxis()->SetTitle("smoothing parameter");
      fitResultsRes->GetYaxis()->SetTitle("P_{2}");

      fitResultsRes->SetFillColor(0);
      fitResultsRes->SetMarkerStyle(markers[i-1]);
      fitResultsRes->SetMarkerColor(colors[i-1]);
      fitResultsRes->SetLineColor(colors[i-1]);

      for (Float_t weight = 0.0; weight < 1.01; weight += 0.2)
      {
        Float_t chi2MC = 0;
        Float_t residuals = 0;

        mult->ApplyBayesianMethod(histID, kFALSE, type, weight, iter, 0, kFALSE);
        mult->DrawComparison(Form("%s/BayesianIterSmooth_%03d_%d_%d_%f", targetDir, count++, type, iter, weight), histID, kFALSE, kTRUE, multESD->GetMultiplicityMC(histID, type)->ProjectionY());
        mult->GetComparisonResults(&chi2MC, 0, &residuals);

        for (Int_t region=0; region<AliMultiplicityCorrection::kQualityRegions; ++region)
          fitResultsMC[region]->SetPoint(fitResultsMC[region]->GetN(), weight, mult->GetQuality(region));

        fitResultsRes->SetPoint(fitResultsRes->GetN(), weight, residuals);
      }

      graphFile->cd();
      for (Int_t region=0; region<AliMultiplicityCorrection::kQualityRegions; ++region)
        fitResultsMC[region]->Write();

      fitResultsRes->Write();
    }
  }

  graphFile->Close();
}

void EvaluateDrawResult(const char* targetDir, Int_t type = 0, Bool_t plotRes = kTRUE)
{
  gSystem->Load("libPWG0base");

  TString name;
  TFile* graphFile = 0;
  if (type == -1)
  {
    name = "EvaluateChi2Method";
    graphFile = TFile::Open(Form("%s/EvaluateChi2Method.root", targetDir));
  }
  else
  {
    name.Form("EvaluateBayesianMethodIterationsSmoothing_%s", GetEventTypeName(type));
    graphFile = TFile::Open(Form("%s/EvaluateBayesianMethodIterationsSmoothing.root", targetDir));
  }

  TCanvas* canvas = new TCanvas(name, name, 800, 500);
  if (type == -1)
  {
    canvas->SetLogx();
    canvas->SetLogy();
  }
  canvas->SetTopMargin(0.05);
  canvas->SetGridx();
  canvas->SetGridy();

  TLegend* legend = new TLegend(0.8, 0.15, 0.98, 0.98);
  legend->SetFillColor(0);

  Int_t count = 1;

  Float_t xMin = 1e20;
  Float_t xMax = 0;

  Float_t yMin = 1e20;
  Float_t yMax = 0;

  Float_t yMinRegion[3];
  for (Int_t i=0; i<AliMultiplicityCorrection::kQualityRegions; ++i)
    yMinRegion[i] = 1e20;

  TString xaxis, yaxis;

  while (1)
  {
    TGraph* mc = (TGraph*) graphFile->Get(Form("%s_MC_%d", name.Data(), count));
    TGraph* res = (TGraph*) graphFile->Get(Form("%s_Res_%d", name.Data(), count));

    if (!mc)
      break;

    xaxis = mc->GetXaxis()->GetTitle();
    yaxis = mc->GetYaxis()->GetTitle();

    mc->Print();

    if (res)
      res->Print();

    xMin = TMath::Min(xMin, mc->GetXaxis()->GetXmin());
    yMin = TMath::Min(yMin, mc->GetYaxis()->GetXmin());

    xMax = TMath::Max(xMax, mc->GetXaxis()->GetXmax());
    yMax = TMath::Max(yMax, mc->GetYaxis()->GetXmax());

    if (plotRes && res)
    {
      xMin = TMath::Min(xMin, res->GetXaxis()->GetXmin());
      yMin = TMath::Min(yMin, res->GetYaxis()->GetXmin());

      xMax = TMath::Max(xMax, res->GetXaxis()->GetXmax());
      yMax = TMath::Max(yMax, res->GetYaxis()->GetXmax());
    }

    for (Int_t i=0; i<mc->GetN(); ++i)
      yMinRegion[(count-1) % 3] = TMath::Min(yMinRegion[(count-1) % 3], mc->GetY()[i]);

    count++;
  }

  for (Int_t i=0; i<AliMultiplicityCorrection::kQualityRegions; ++i)
    Printf("Minimum for region %d is %f", i, yMinRegion[i]);

  if (type >= 0)
  {
    xaxis = "smoothing parameter";
  }
  else if (type == -1)
  {
    xaxis = "weight parameter";
    xMax *= 5;
  }
  //yaxis = "P_{1} (2 <= t <= 150)";

  printf("%f %f %f %f\n", xMin, xMax, yMin, yMax);

  TGraph* dummy = new TGraph;
  dummy->SetPoint(0, xMin, yMin);
  dummy->SetPoint(1, xMax, yMax);
  dummy->SetTitle(Form(";%s;%s", xaxis.Data(), yaxis.Data()));

  dummy->SetMarkerColor(0);
  dummy->Draw("AP");
  dummy->GetYaxis()->SetMoreLogLabels(1);

  count = 1;

  while (1)
  {
    TGraph* mc = (TGraph*) graphFile->Get(Form("%s_MC_%d", name.Data(), count));
    TGraph* res = (TGraph*) graphFile->Get(Form("%s_Res_%d", name.Data(), count));

    //printf("%s_MC_%d %p %p\n", name.Data(), count, mc, res);

    if (!mc)
      break;

    printf("Loaded %d sucessful.\n", count);

    if (type == -1)
    {
      legend->AddEntry(mc, Form("Eq. (%d) - reg. %d", 10 + (count-1) / 3, 1+ (count-1) % 3));
    }
    else
      legend->AddEntry(mc);

    mc->Draw("SAME PC");

    if (plotRes && res)
    {
      legend->AddEntry(res);
      res->Draw("SAME PC");
    }

    count++;
  }

  legend->Draw();

  canvas->SaveAs(Form("%s/%s.gif", targetDir, canvas->GetName()));
  canvas->SaveAs(Form("%s/%s.eps", targetDir, canvas->GetName()));
}

void EvaluateDrawResultRegions(const char* targetDir, Int_t type = 0)
{
  loadlibs();

  TString name;
  TFile* graphFile = 0;
  if (type == -1)
  {
    name = "EvaluateChi2Method";
    graphFile = TFile::Open(Form("%s/EvaluateChi2Method.root", targetDir));
  }
  else
  {
    name.Form("EvaluateBayesianMethodIterationsSmoothing_%s", GetEventTypeName(type));
    graphFile = TFile::Open(Form("%s/EvaluateBayesianMethodIterationsSmoothing.root", targetDir));
  }

  TCanvas* canvas = new TCanvas(name, name, 1200, 800);
  //canvas->Divide(1, AliMultiplicityCorrection::kQualityRegions, 0, 0);
  canvas->SetTopMargin(0);
  canvas->SetBottomMargin(0);
  canvas->SetRightMargin(0.05);
  canvas->Range(0, 0, 1, 1);

  TPad* pad[4];
  pad[3] = new TPad(Form("%s_pad1", name.Data()), "", 0, 0, 0.5, 0.5); pad[3]->SetTopMargin(0); pad[3]->SetBottomMargin(0.15);
  pad[1] = new TPad(Form("%s_pad2", name.Data()), "", 0, 0.5, 0.5, 1); pad[1]->SetBottomMargin(0); 
  pad[0] = new TPad(Form("%s_pad3", name.Data()), "", 0.5, 0, 1, 0.5); pad[0]->SetTopMargin(0); pad[0]->SetBottomMargin(0.15);
  pad[2] = new TPad(Form("%s_pad4", name.Data()), "", 0.5, 0.5, 1, 1); pad[2]->SetBottomMargin(0); 

  Float_t yMinRegion[3];
  for (Int_t i=0; i<AliMultiplicityCorrection::kQualityRegions; ++i)
    yMinRegion[i] = 1e20;

  for (Int_t region = 0; region <= AliMultiplicityCorrection::kQualityRegions; region++)
  {
    canvas->cd();
    pad[region]->Draw();
    pad[region]->cd();
    pad[region]->SetRightMargin(0.05);
    
    if (type == -1)
    {
      pad[region]->SetLogx();
    }
    pad[region]->SetLogy();
    
    pad[region]->SetGridx();
    pad[region]->SetGridy();

    TLegend* legend = new TLegend(0.5, 0.4, 0.85, 0.85);
    legend->SetFillColor(0);
    //legend->SetNColumns(3);
    legend->SetTextSize(0.06);
    Int_t count = 0;

    Float_t xMin = 1e20;
    Float_t xMax = 0;

    Float_t yMin = 1e20;
    Float_t yMax = 0;

    TString xaxis, yaxis;

    while (1)
    {
      count++;

      if (region > 0)
      {
        TGraph* mc = (TGraph*) graphFile->Get(Form("%s_MC_%d", name.Data(), count));
        if (!mc)
          break;
  
        if (TString(mc->GetTitle()).Contains(Form("reg. %d", region)) == kFALSE)
          continue;
      
        for (Int_t i=0; i<mc->GetN(); ++i)
          yMinRegion[(count-1) % 3] = TMath::Min(yMinRegion[(count-1) % 3], mc->GetY()[i]);
      }
      else
      {
        TGraph* mc = (TGraph*) graphFile->Get(Form("%s_Res_%d", name.Data(), count));
        if (!mc)
          break;
      }

      xaxis = mc->GetXaxis()->GetTitle();
      yaxis = mc->GetYaxis()->GetTitle();

      //mc->Print();

      xMin = TMath::Min(xMin, mc->GetXaxis()->GetXmin());
      yMin = TMath::Min(yMin, TMath::MinElement(mc->GetN(), mc->GetY()));

      xMax = TMath::Max(xMax, mc->GetXaxis()->GetXmax());
      yMax = TMath::Max(yMax, mc->GetYaxis()->GetXmax());
      
      //Printf("%f %f %f %f", xMin, xMax, yMin, yMax);
    }

    if (type >= 0)
    {
      xaxis = "Smoothing parameter #alpha";
    }
    else if (type == -1)
    {
      xaxis = "Weight parameter #beta";
      //xMax *= 5;
    }
    if (region > 0)
    {
      yaxis = Form("Q_{1} in region %d ", region); // (2 <= t <= 150)";
    }
    else
      yaxis = "Q_{2} ";

    printf("%f %f %f %f\n", xMin, xMax, yMin, yMax);

    if (region > 0)
    {
      if (type == -1)
      {
        yMin = 0.534622;
        yMax = 19.228941;
      }
      else
      {
        yMin = 0.759109;
        yMax = 10;
      }
    }
    
    if (type != -1)
      xMax = 1.0;

    TH1* dummy = new TH2F("dummy", "dummy", 100, xMin, xMax, 100, yMin * 0.9, yMax * 1.1);
    dummy->SetTitle(Form(";%s;%s", xaxis.Data(), yaxis.Data()));

    if (region == 0 && type != -1)
      dummy->GetYaxis()->SetMoreLogLabels(1);
    dummy->SetLabelSize(0.06, "xy");
    dummy->SetTitleSize(0.06, "xy");
    dummy->GetXaxis()->SetTitleOffset(1.2);
    dummy->GetYaxis()->SetTitleOffset(0.8);
    dummy->SetStats(0);
    
    dummy->DrawCopy();
   

    count = 0;

    while (1)
    {
      count++;

      if (region > 0)
      {
        TGraph* mc = (TGraph*) graphFile->Get(Form("%s_MC_%d", name.Data(), count));
        if (mc && TString(mc->GetTitle()).Contains(Form("reg. %d", region)) == kFALSE)
          continue;
      }
      else
        TGraph* mc = (TGraph*) graphFile->Get(Form("%s_Res_%d", name.Data(), count));

      if (!mc)
        break;
      //printf("Loaded %d sucessful.\n", count);

      if (type == -1)
      {
        //legend->AddEntry(mc, Form("Eq. (%d) - reg. %d", 10 + (count-1) / 3, 1+ (count-1) % 3));
        //legend->AddEntry(mc, Form("Eq. (7.%d)", 11 + (count-1) / 3));
        const char* names[] = { "Const", "Linear", "Log" };
        legend->AddEntry(mc, names[(count-1) / 3], "PL");
      }
      else
      {
        TString label(mc->GetTitle());
        TString newLabel(label(0, label.Index(" ")));
        //Printf("%s", newLabel.Data());
        if (newLabel.Atoi() == -1)
        {
          newLabel = "Convergence";
        }
        else
          newLabel = label(0, label.Index("iterations") + 10);
          
        legend->AddEntry(mc, newLabel, "PL");
      }

      mc->Draw("SAME PC");

    }

    if (region == 2)
      legend->Draw();
  }

  for (Int_t i=0; i<AliMultiplicityCorrection::kQualityRegions; ++i)
    Printf("Minimum for region %d is %f", i, yMinRegion[i]);

  canvas->Modified();
  canvas->SaveAs(Form("%s/%s.gif", targetDir, canvas->GetName()));
  canvas->SaveAs(Form("%s/%s.eps", targetDir, canvas->GetName()));
}

void EvaluateChi2Method(const char* fileNameMC = "multiplicityMC.root", const char* fileNameESD = "multiplicityESD.root", const char* targetDir = "chi2compare", Int_t histID = 1)
{
  gSystem->mkdir(targetDir);

  gSystem->Load("libPWG0base");

  AliMultiplicityCorrection* mult = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");

  TFile::Open(fileNameMC);
  mult->LoadHistograms("Multiplicity");

  TFile::Open(fileNameESD);
  TH2F* hist = (TH2F*) gFile->Get(Form("Multiplicity/fMultiplicityESD%d", histID));
  TH2F* hist2 = (TH2F*) gFile->Get(Form("Multiplicity/fMultiplicityVtx%d", histID));

  mult->SetMultiplicityESD(histID, hist);

  Int_t count = 0; // just to order the saved images...
  Int_t colors[] = {1, 2, 4, 3};
  Int_t markers[] = {20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 3};

  TGraph* fitResultsRes = 0;

  TFile* graphFile = TFile::Open(Form("%s/EvaluateChi2Method.root", targetDir), "RECREATE");

  //for (AliMultiplicityCorrection::RegularizationType type = AliMultiplicityCorrection::kPol0; type <= AliMultiplicityCorrection::kLog; ++type)
  //for (AliMultiplicityCorrection::RegularizationType type = AliMultiplicityCorrection::kEntropy; type <= AliMultiplicityCorrection::kEntropy; ++type)
  for (AliMultiplicityCorrection::RegularizationType type = AliMultiplicityCorrection::kPol1; type <= AliMultiplicityCorrection::kPol1; ++type)
  {
    TGraph* fitResultsMC[3];
    for (Int_t region=0; region<AliMultiplicityCorrection::kQualityRegions; ++region)
    {
      fitResultsMC[region] = new TGraph;
      fitResultsMC[region]->SetTitle(Form("Eq. (%d) - reg. %d", type+9, region+1));
      fitResultsMC[region]->GetXaxis()->SetTitle("weight parameter #alpha");
      fitResultsMC[region]->GetYaxis()->SetTitle(Form("P_{1} in region %d", region));
      fitResultsMC[region]->SetName(Form("EvaluateChi2Method_MC_%d", type * AliMultiplicityCorrection::kQualityRegions + region - 2));
      fitResultsMC[region]->SetFillColor(0);
      //fitResultsMC[region]->SetMarkerStyle(markers[(type-1) * AliMultiplicityCorrection::kQualityRegions + region]);
      //fitResultsMC[region]->SetLineColor(colors[region]);
      fitResultsMC[region]->SetMarkerStyle(markers[type-1]);
      fitResultsMC[region]->SetMarkerColor(colors[type-1]);
      fitResultsMC[region]->SetLineColor(colors[type-1]);
    }

    fitResultsRes = new TGraph;
    fitResultsRes->SetTitle(Form("%s residual chi2", GetRegName(type)));
    fitResultsRes->SetName(Form("EvaluateChi2Method_Res_%d", type));
    fitResultsRes->GetXaxis()->SetTitle("Weight Parameter");

    fitResultsRes->SetFillColor(0);
    fitResultsRes->SetMarkerStyle(markers[type-1]);
    fitResultsRes->SetMarkerColor(colors[type-1]);
    fitResultsRes->SetLineColor(colors[type-1]);

    for (Int_t i=0; i<15; ++i)
    {
      Float_t weight = TMath::Power(TMath::Sqrt(10), i+1);
      //Float_t weight = TMath::Power(10, i+2);

      //if (type == AliMultiplicityCorrection::kEntropy) weight = 1e4 * (i+1) * 1.5;

      Float_t chi2MC = 0;
      Float_t residuals = 0;
      Float_t chi2Limit = 0;

      TString runName;
      runName.Form("MinuitChi2_%02d_%d_%f", count++, type, weight);

      mult->SetRegularizationParameters(type, weight);
      Int_t status = mult->ApplyMinuitFit(histID, kFALSE, AliMultiplicityCorrection::kTrVtx);
      mult->DrawComparison(Form("%s/%s", targetDir, runName.Data()), histID, kFALSE, kTRUE, hist2->ProjectionY("mymc", 1, hist2->GetNbinsX()));
      if (status != 0)
      {
        printf("MINUIT did not succeed. Skipping...\n");
        continue;
      }

      mult->GetComparisonResults(&chi2MC, 0, &residuals);
      TH1* result = mult->GetMultiplicityESDCorrected(histID);
      result->SetName(runName);
      result->Write();

      for (Int_t region=0; region<AliMultiplicityCorrection::kQualityRegions; ++region)
        fitResultsMC[region]->SetPoint(fitResultsMC[region]->GetN(), weight, mult->GetQuality(region));

      fitResultsRes->SetPoint(fitResultsRes->GetN(), weight, residuals);
    }

    graphFile->cd();
    for (Int_t region=0; region<AliMultiplicityCorrection::kQualityRegions; ++region)
      fitResultsMC[region]->Write();
    fitResultsRes->Write();
  }

  graphFile->Close();
}

void EvaluateChi2MethodAll()
{
  EvaluateChi2Method("multiplicityMC_3M.root", "multiplicityMC_3M.root", "eval-3M-3M");
  EvaluateChi2Method("multiplicityMC_2M.root", "multiplicityMC_1M_3.root", "eval-2M-1M");
  EvaluateChi2Method("multiplicityMC_3M.root", "multiplicityMC_3M_NBD.root", "eval-3M-NBD");
  EvaluateChi2Method("multiplicityMC_2M_smoothed.root", "multiplicityMC_1M_3.root", "eval-2MS-1M");
  EvaluateChi2Method("multiplicityMC_2M_smoothed.root", "multiplicityMC_3M_NBD.root", "eval-2MS-NBD");
}

void EvaluateBayesianMethodAll()
{
  EvaluateBayesianMethod("multiplicityMC_3M.root", "multiplicityMC_3M.root", "eval-3M-3M");
  EvaluateBayesianMethod("multiplicityMC_2M.root", "multiplicityMC_1M_3.root", "eval-2M-1M");
  EvaluateBayesianMethod("multiplicityMC_3M.root", "multiplicityMC_3M_NBD.root", "eval-3M-NBD");
  EvaluateBayesianMethod("multiplicityMC_2M_smoothed.root", "multiplicityMC_1M_3.root", "eval-2MS-1M");
  EvaluateBayesianMethod("multiplicityMC_2M_smoothed.root", "multiplicityMC_3M_NBD.root", "eval-2MS-NBD");
}

void CompareMethods(const char* fileNameMC = "multiplicityMC.root", const char* fileNameESD = "multiplicityMC.root", const char* targetDir, Int_t histID = 3)
{
  gSystem->mkdir(targetDir);

  gSystem->Load("libPWG0base");

  AliMultiplicityCorrection* mult = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");

  TFile::Open(fileNameMC);
  mult->LoadHistograms("Multiplicity");

  TFile::Open(fileNameESD);
  AliMultiplicityCorrection* multESD = new AliMultiplicityCorrection("MultiplicityESD", "MultiplicityESD");
  multESD->LoadHistograms("Multiplicity");

  mult->SetMultiplicityESD(histID, multESD->GetMultiplicityESD(histID));

  TCanvas* canvas = new TCanvas("CompareMethods", "CompareMethods", 1200, 1200);
  canvas->Divide(3, 3);

  Int_t count = 0;

  for (AliMultiplicityCorrection::EventType type = AliMultiplicityCorrection::kTrVtx; type <= AliMultiplicityCorrection::kTrVtx; ++type)
  {
    TH1* mc = multESD->GetMultiplicityMC(histID, type)->ProjectionY("mymc");
    mc->Sumw2();

    mult->SetRegularizationParameters(AliMultiplicityCorrection::kPol1, 10000);
    mult->ApplyMinuitFit(histID, kFALSE, type);
    mult->DrawComparison(Form("%s/CompareMethods_%d_MinuitChi2", targetDir, type), histID, kFALSE, kTRUE, mc);
    TH1* chi2Result = (TH1*) mult->GetMultiplicityESDCorrected(histID)->Clone("chi2Result");

    mult->ApplyBayesianMethod(histID, kFALSE, type, 0.1);
    mult->DrawComparison(Form("%s/CompareMethods_%d_Bayesian", targetDir, type), histID, kFALSE, kTRUE, mc);
    TH1* bayesResult = (TH1*) mult->GetMultiplicityESDCorrected(histID)->Clone("bayesResult");

    mc->GetXaxis()->SetRangeUser(0, 150);
    chi2Result->GetXaxis()->SetRangeUser(0, 150);

/*    // skip errors for now
    for (Int_t i=1; i<=chi2Result->GetNbinsX(); ++i)
    {
      chi2Result->SetBinError(i, 0);
      bayesResult->SetBinError(i, 0);
    }*/

    canvas->cd(++count);
    mc->SetFillColor(kYellow);
    mc->DrawCopy();
    chi2Result->SetLineColor(kRed);
    chi2Result->DrawCopy("SAME");
    bayesResult->SetLineColor(kBlue);
    bayesResult->DrawCopy("SAME");
    gPad->SetLogy();

    canvas->cd(count + 3);
    chi2ResultRatio = (TH1*) chi2Result->Clone("chi2ResultRatio");
    bayesResultRatio = (TH1*) bayesResult->Clone("bayesResultRatio");
    chi2ResultRatio->Divide(chi2Result, mc, 1, 1, "");
    bayesResultRatio->Divide(bayesResult, mc, 1, 1, "");

    chi2ResultRatio->GetYaxis()->SetRangeUser(0.5, 1.5);

    chi2ResultRatio->DrawCopy("HIST");
    bayesResultRatio->DrawCopy("SAME HIST");

    canvas->cd(count + 6);
    chi2Result->Divide(chi2Result, bayesResult, 1, 1, "");
    chi2Result->GetYaxis()->SetRangeUser(0.5, 1.5);
    chi2Result->DrawCopy("HIST");
  }

  canvas->SaveAs(Form("%s/%s.gif", targetDir, canvas->GetName()));
  canvas->SaveAs(Form("%s/%s.C", targetDir, canvas->GetName()));
}

void StatisticsPlot(const char* fileNameMC = "multiplicityMC_2M.root", Int_t histID = 3)
{
  gSystem->Load("libPWG0base");

  AliMultiplicityCorrection* mult = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");

  TFile::Open(fileNameMC);
  mult->LoadHistograms("Multiplicity");

  const char* files[] = { "multiplicityMC_100k_1.root", "multiplicityMC_200k.root", "multiplicityMC_400k.root", "multiplicityMC_600k.root", "multiplicityMC_800k.root" };

  TGraph* fitResultsChi2[3];
  TGraph* fitResultsBayes[3];

  for (Int_t region=0; region<AliMultiplicityCorrection::kQualityRegions; ++region)
  {
    fitResultsChi2[region] = new TGraph;
    fitResultsChi2[region]->SetTitle(";Nevents;Chi2");
    fitResultsChi2[region]->SetName(Form("fitResultsChi2_%d", region));
    fitResultsChi2[region]->SetMarkerStyle(20+region);

    fitResultsBayes[region] = new TGraph;
    fitResultsBayes[region]->SetTitle(";Nevents;Chi2");
    fitResultsBayes[region]->SetName(Form("fitResultsBayes_%d", region));
    fitResultsBayes[region]->SetMarkerStyle(20+region);
    fitResultsBayes[region]->SetMarkerColor(2);
  }

  TGraph* fitResultsChi2Limit = new TGraph;  fitResultsChi2Limit->SetTitle(";Nevents;Multiplicity reach");
  TGraph* fitResultsBayesLimit = new TGraph; fitResultsBayesLimit->SetTitle(";Nevents;Multiplicity reach");
  TGraph* fitResultsChi2Res = new TGraph;       fitResultsChi2Res->SetTitle(";Nevents;Chi2");
  TGraph* fitResultsBayesRes = new TGraph;      fitResultsBayesRes->SetTitle(";Nevents;Chi2");

  fitResultsChi2Limit->SetName("fitResultsChi2Limit");
  fitResultsBayesLimit->SetName("fitResultsBayesLimit");
  fitResultsChi2Res->SetName("fitResultsChi2Res");
  fitResultsBayesRes->SetName("fitResultsBayesRes");

  TCanvas* canvas = new TCanvas("StatisticsPlot", "StatisticsPlot", 1200, 600);
  canvas->Divide(5, 2);

  Float_t min = 1e10;
  Float_t max = 0;

  TFile* file = TFile::Open("StatisticsPlot.root", "RECREATE");
  file->Close();

  for (Int_t i=0; i<5; ++i)
  {
    TFile::Open(files[i]);
    AliMultiplicityCorrection* multESD = new AliMultiplicityCorrection("MultiplicityESD", "MultiplicityESD");
    multESD->LoadHistograms("Multiplicity");

    mult->SetMultiplicityESD(histID, multESD->GetMultiplicityESD(histID));
    Int_t nEntries = multESD->GetMultiplicityESD(histID)->GetEntries();
    TH1* mc = multESD->GetMultiplicityMC(histID, AliMultiplicityCorrection::kTrVtx)->ProjectionY(Form("mc_%d", i));

    mult->SetRegularizationParameters(AliMultiplicityCorrection::kPol1, 10000);
    mult->ApplyMinuitFit(histID, kFALSE, AliMultiplicityCorrection::kTrVtx);
    mult->DrawComparison(Form("StatisticsPlot_%d_MinuitChi2", i), histID, kFALSE, kTRUE, mc);

    Int_t chi2MCLimit = 0;
    Float_t chi2Residuals = 0;
    mult->GetComparisonResults(0, &chi2MCLimit, &chi2Residuals);
    for (Int_t region=0; region<AliMultiplicityCorrection::kQualityRegions; ++region)
    {
      fitResultsChi2[region]->SetPoint(fitResultsChi2[region]->GetN(), nEntries, mult->GetQuality(region));
      min = TMath::Min(min, mult->GetQuality(region));
      max = TMath::Max(max, mult->GetQuality(region));
    }
    fitResultsChi2Limit->SetPoint(fitResultsChi2Limit->GetN(), nEntries, chi2MCLimit);
    fitResultsChi2Res->SetPoint(fitResultsChi2Res->GetN(), nEntries, chi2Residuals);

    TH1* chi2Result = (TH1*) mult->GetMultiplicityESDCorrected(histID)->Clone(Form("chi2Result_%d", i));

    mult->ApplyBayesianMethod(histID, kFALSE, AliMultiplicityCorrection::kTrVtx, 1, 100, 0, kFALSE);
    mult->DrawComparison(Form("StatisticsPlot_%d_Bayesian", i), histID, kFALSE, kTRUE, mc);
    TH1* bayesResult = (TH1*) mult->GetMultiplicityESDCorrected(histID)->Clone(Form("bayesResult_%d", i));
    mult->GetComparisonResults(0, &chi2MCLimit, &chi2Residuals);
    for (Int_t region=0; region<AliMultiplicityCorrection::kQualityRegions; ++region)
    {
      fitResultsBayes[region]->SetPoint(fitResultsBayes[region]->GetN(), nEntries, mult->GetQuality(region));
      min = TMath::Min(min, mult->GetQuality(region));
      max = TMath::Max(max, mult->GetQuality(region));
    }
    fitResultsBayesLimit->SetPoint(fitResultsBayesLimit->GetN(), nEntries, chi2MCLimit);
    fitResultsBayesRes->SetPoint(fitResultsBayesRes->GetN(), nEntries, chi2Residuals);

    mc->GetXaxis()->SetRangeUser(0, 150);
    chi2Result->GetXaxis()->SetRangeUser(0, 150);

    // skip errors for now
    for (Int_t j=0; j<=chi2Result->GetNbinsX(); ++j)
    {
      chi2Result->SetBinError(j, 0);
      bayesResult->SetBinError(j, 0);
    }

    canvas->cd(i+1);
    mc->SetFillColor(kYellow);
    mc->DrawCopy();
    chi2Result->SetLineColor(kRed);
    chi2Result->DrawCopy("SAME");
    bayesResult->SetLineColor(kBlue);
    bayesResult->DrawCopy("SAME");
    gPad->SetLogy();

    canvas->cd(i+6);
    chi2Result->Divide(chi2Result, mc, 1, 1, "B");
    bayesResult->Divide(bayesResult, mc, 1, 1, "B");

    // skip errors for now
    for (Int_t j=0; j<=chi2Result->GetNbinsX(); ++j)
    {
      chi2Result->SetBinError(j, 0);
      bayesResult->SetBinError(j, 0);
    }

    chi2Result->SetTitle("Ratios;Npart;unfolded measured/MC");
    chi2Result->GetYaxis()->SetRangeUser(0.5, 1.5);

    chi2Result->DrawCopy("");
    bayesResult->DrawCopy("SAME");

    TFile* file = TFile::Open("StatisticsPlot.root", "UPDATE");
    mc->Write();
    chi2Result->Write();
    bayesResult->Write();
    file->Close();
  }

  canvas->SaveAs(Form("%s.gif", canvas->GetName()));
  canvas->SaveAs(Form("%s.C", canvas->GetName()));

  TCanvas* canvas2 = new TCanvas("StatisticsPlot2", "StatisticsPlot2", 800, 400);
  canvas2->Divide(2, 1);

  canvas2->cd(1);

  for (Int_t region=0; region<AliMultiplicityCorrection::kQualityRegions; ++region)
  {
    fitResultsChi2[region]->GetYaxis()->SetRangeUser(0.5 * min, 1.5 * max);
    fitResultsChi2[region]->Draw(((region == 0) ? "AP" : "P SAME"));

    fitResultsBayes[region]->Draw("P SAME");
  }

  gPad->SetLogy();

  canvas2->cd(2);
  fitResultsChi2Limit->SetMarkerStyle(20);
  fitResultsChi2Limit->GetYaxis()->SetRangeUser(0.9 * TMath::Min(fitResultsChi2Limit->GetYaxis()->GetXmin(), fitResultsBayesLimit->GetYaxis()->GetXmin()), 1.1 * TMath::Max(fitResultsChi2Limit->GetYaxis()->GetXmax(), fitResultsBayesLimit->GetYaxis()->GetXmax()));
  fitResultsChi2Limit->Draw("AP");

  fitResultsBayesLimit->SetMarkerStyle(3);
  fitResultsBayesLimit->SetMarkerColor(2);
  fitResultsBayesLimit->Draw("P SAME");

  canvas2->SaveAs(Form("%s.gif", canvas2->GetName()));
  canvas2->SaveAs(Form("%s.C", canvas2->GetName()));

  TFile* file = TFile::Open("StatisticsPlot.root", "UPDATE");

  for (Int_t region=0; region<AliMultiplicityCorrection::kQualityRegions; ++region)
  {
    fitResultsChi2[region]->Write();
    fitResultsBayes[region]->Write();
  }
  fitResultsChi2Limit->Write();
  fitResultsBayesLimit->Write();
  fitResultsChi2Res->Write();
  fitResultsBayesRes->Write();
  file->Close();
}

void StartingConditions(const char* fileNameMC = "multiplicityMC.root", Int_t histID = 1)
{
  loadlibs();

  AliMultiplicityCorrection* mult = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");

  TFile::Open(fileNameMC);
  mult->LoadHistograms("Multiplicity");

  const char* files[] = { "multiplicityESD.root", "multiplicityESD_100k_1.root", "multiplicityESD_100k_2.root", "multiplicityESD_100k_1.root", "multiplicityESD_100k_2.root" }

  // this one we try to unfold
  TFile::Open(files[0]);
  AliMultiplicityCorrection* multESD = new AliMultiplicityCorrection("MultiplicityESD", "MultiplicityESD");
  multESD->LoadHistograms("Multiplicity");
  mult->SetMultiplicityESD(histID, multESD->GetMultiplicityESD(histID));
  TH1* mc = multESD->GetMultiplicityMC(histID, AliMultiplicityCorrection::kTrVtx)->ProjectionY("mc", 1, 1);

  TGraph* fitResultsChi2 = new TGraph;
  fitResultsChi2->SetTitle(";Input Dist ID;Chi2");
  TGraph* fitResultsBayes = new TGraph;
  fitResultsBayes->SetTitle(";Input Dist ID;Chi2");
  TGraph* fitResultsChi2Limit = new TGraph;
  fitResultsChi2Limit->SetTitle(";Input Dist ID;Multiplicity reach");
  TGraph* fitResultsBayesLimit = new TGraph;
  fitResultsBayesLimit->SetTitle(";Input Dist ID;Multiplicity reach");

  TCanvas* canvas = new TCanvas("StartingConditions", "StartingConditions", 1200, 600);
  canvas->Divide(8, 2);

  TCanvas* canvas3 = new TCanvas("StartingConditions3", "StartingConditions3", 1000, 400);
  canvas3->Divide(2, 1);

  Float_t min = 1e10;
  Float_t max = 0;

  TH1* firstChi = 0;
  TH1* firstBayesian = 0;
  TH1* startCond = multESD->GetMultiplicityESD(histID)->ProjectionY("startCond");

  TLegend* legend = new TLegend(0.7, 0.7, 1, 1);

  TFile* file = TFile::Open("StartingConditions.root", "RECREATE");
  mc->Write();
  file->Close();

  for (Int_t i=0; i<8; ++i)
  {
    if (i == 0)
    {
      startCond = (TH1*) mc->Clone("startCond2");
    }
    else if (i < 6)
    {
      TFile::Open(files[i-1]);
      AliMultiplicityCorrection* multESD2 = new AliMultiplicityCorrection("MultiplicityESD2", "MultiplicityESD2");
      multESD2->LoadHistograms("Multiplicity");
      startCond = multESD2->GetMultiplicityESD(histID)->ProjectionY("startCond");
    }
    else if (i == 6)
    {
      func = new TF1("nbd", "exp(log([0]) + TMath::LnGamma([2]+x) - TMath::LnGamma([2]) - TMath::LnGamma(x+1) + log([1] / ([1]+[2])) * x + log(1.0 + [1]/[2]) * -[2])", 0, 50);
      
      func->SetParNames("scaling", "averagen", "k");
      func->SetParLimits(0, 1e-3, 1e10);
      func->SetParLimits(1, 0.001, 1000);
      func->SetParLimits(2, 0.001, 1000);

      func->SetParameters(1, 10, 2);
      for (Int_t j=2; j<=startCond->GetNbinsX(); j++)
        startCond->SetBinContent(j, func->Eval(j-1));
    }
    else if (i == 7)
    {
      for (Int_t j=1; j<=startCond->GetNbinsX(); j++)
        startCond->SetBinContent(j, 1);
    }

    mult->SetRegularizationParameters(AliMultiplicityCorrection::kPol0, 1e5);
    mult->ApplyMinuitFit(histID, kFALSE, AliMultiplicityCorrection::kTrVtx, kFALSE, startCond);
    //mult->DrawComparison(Form("StartingConditions_%d_MinuitChi2", i), histID, kFALSE, kTRUE, mc);

    Float_t chi2MC = 0;
    Int_t chi2MCLimit = 0;
    mult->GetComparisonResults(&chi2MC, &chi2MCLimit, 0);
    fitResultsChi2->SetPoint(fitResultsChi2->GetN(), i, chi2MC);
    fitResultsChi2Limit->SetPoint(fitResultsChi2Limit->GetN(), i, chi2MCLimit);
    min = TMath::Min(min, chi2MC);
    max = TMath::Max(max, chi2MC);

    TH1* chi2Result = (TH1*) mult->GetMultiplicityESDCorrected(histID)->Clone(Form("chi2Result_%d", i));
    if (!firstChi)
      firstChi = (TH1*) chi2Result->Clone("firstChi");

    mult->ApplyBayesianMethod(histID, kFALSE, AliMultiplicityCorrection::kTrVtx, 1, 10, startCond);
    //mult->DrawComparison(Form("StartingConditions_%d_Bayesian", i), histID, kFALSE, kTRUE, mc);
    TH1* bayesResult = (TH1*) mult->GetMultiplicityESDCorrected(histID)->Clone(Form("bayesResult_%d", i));
    if (!firstBayesian)
      firstBayesian = (TH1*) bayesResult->Clone("firstBayesian");

    mult->GetComparisonResults(&chi2MC, &chi2MCLimit, 0);
    fitResultsBayes->SetPoint(fitResultsBayes->GetN(), i, chi2MC);
    fitResultsBayesLimit->SetPoint(fitResultsBayesLimit->GetN(), i, chi2MCLimit);

    TFile* file = TFile::Open("StartingConditions.root", "UPDATE");
    chi2Result->Write();
    bayesResult->Write();
    file->Close();

    min = TMath::Min(min, chi2MC);
    max = TMath::Max(max, chi2MC);
    mc->GetXaxis()->SetRangeUser(0, 150);
    chi2Result->GetXaxis()->SetRangeUser(0, 150);

    // skip errors for now
    for (Int_t j=0; j<=chi2Result->GetNbinsX(); ++j)
    {
      chi2Result->SetBinError(j, 0);
      bayesResult->SetBinError(j, 0);
    }

    canvas3->cd(1);
    TH1* tmp = (TH1*) chi2Result->Clone("tmp");
    tmp->SetTitle("Difference to best initial conditions;Npart;Ratio");
    tmp->Divide(firstChi);
    tmp->GetYaxis()->SetRangeUser(0.5, 1.5);
    tmp->GetXaxis()->SetRangeUser(0, 200);
    tmp->SetLineColor(i+1);
    legend->AddEntry(tmp, Form("%d", i));
    tmp->DrawCopy((i > 0) ? "SAME HIST" : "HIST");

    canvas3->cd(2);
    tmp = (TH1*) bayesResult->Clone("tmp");
    tmp->SetTitle("Difference to best initial conditions;Npart;Ratio");
    tmp->Divide(firstBayesian);
    tmp->SetLineColor(i+1);
    tmp->GetYaxis()->SetRangeUser(0.5, 1.5);
    tmp->GetXaxis()->SetRangeUser(0, 200);
    tmp->DrawCopy((i > 0) ? "SAME HIST" : "HIST");

    canvas->cd(i+1);
    mc->SetFillColor(kYellow);
    mc->DrawCopy();
    chi2Result->SetLineColor(kRed);
    chi2Result->DrawCopy("SAME");
    bayesResult->SetLineColor(kBlue);
    bayesResult->DrawCopy("SAME");
    gPad->SetLogy();

    canvas->cd(i+9);
    chi2Result->Divide(chi2Result, mc, 1, 1, "B");
    bayesResult->Divide(bayesResult, mc, 1, 1, "B");

    // skip errors for now
    for (Int_t j=0; j<=chi2Result->GetNbinsX(); ++j)
    {
      chi2Result->SetBinError(j, 0);
      bayesResult->SetBinError(j, 0);
    }

    chi2Result->SetTitle("Ratios;Npart;unfolded measured/MC");
    chi2Result->GetYaxis()->SetRangeUser(0.5, 1.5);

    chi2Result->DrawCopy("");
    bayesResult->DrawCopy("SAME");
  }

  canvas3->cd(1);
  legend->Draw();

  canvas->SaveAs(Form("%s.gif", canvas->GetName()));

  TCanvas* canvas2 = new TCanvas("StartingConditions2", "StartingConditions2", 800, 400);
  canvas2->Divide(2, 1);

  canvas2->cd(1);
  fitResultsChi2->SetMarkerStyle(20);
  fitResultsChi2->GetYaxis()->SetRangeUser(0.5 * min, 1.5 * max);
  fitResultsChi2->Draw("AP");

  fitResultsBayes->SetMarkerStyle(3);
  fitResultsBayes->SetMarkerColor(2);
  fitResultsBayes->Draw("P SAME");

  gPad->SetLogy();

  canvas2->cd(2);
  fitResultsChi2Limit->SetMarkerStyle(20);
  fitResultsChi2Limit->GetYaxis()->SetRangeUser(0.9 * TMath::Min(fitResultsChi2Limit->GetYaxis()->GetXmin(), fitResultsBayesLimit->GetYaxis()->GetXmin()), 1.1 * TMath::Max(fitResultsChi2Limit->GetYaxis()->GetXmax(), fitResultsBayesLimit->GetYaxis()->GetXmax()));
  fitResultsChi2Limit->Draw("AP");

  fitResultsBayesLimit->SetMarkerStyle(3);
  fitResultsBayesLimit->SetMarkerColor(2);
  fitResultsBayesLimit->Draw("P SAME");

  canvas2->SaveAs(Form("%s.gif", canvas2->GetName()));
  canvas3->SaveAs(Form("%s.gif", canvas3->GetName()));
}

void Merge(Int_t n, const char** files, const char* output)
{
  // const char* files[] = { "multiplicityMC_100k_1.root",  "multiplicityMC_100k_2.root",  "multiplicityMC_100k_3.root", "multiplicityMC_100k_4.root",  "multiplicityMC_100k_5.root",  "multiplicityMC_100k_6.root",  "multiplicityMC_100k_7.root",  "multiplicityMC_100k_8.root" };


  gSystem->Load("libPWG0base");

  AliMultiplicityCorrection** data = new AliMultiplicityCorrection*[n];
  TList list;
  for (Int_t i=0; i<n; ++i)
  {
    TString name("Multiplicity");
    if (i > 0)
      name.Form("Multiplicity%d", i);

    TFile::Open(files[i]);
    Printf("Loading from file %s", files[i]);
    data[i] = new AliMultiplicityCorrection(name, name);
    data[i]->LoadHistograms("Multiplicity");
    if (i > 0)
      list.Add(data[i]);
  }

  data[0]->Merge(&list);

  //data[0]->DrawHistograms();

  TFile::Open(output, "RECREATE");
  data[0]->SaveHistograms();
  gFile->Close();
}

void testMethod(Int_t caseNo, const char* fileName = "multiplicityMC.root")
{
  loadlibs();

  AliMultiplicityCorrection* mult = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");

  TFile::Open(fileName);
  mult->LoadHistograms("Multiplicity");

  TF1* func = 0;

  if (caseNo >= 4)
  {
    TF1* func = new TF1("nbd", "exp(log([0]) + TMath::LnGamma([2]+x) - TMath::LnGamma([2]) - TMath::LnGamma(x+1) + log([1] / ([1]+[2])) * x + log(1.0 + [1]/[2]) * -[2])", 0, 200);
    //func = new TF1("nbd", "[0] * TMath::Binomial([2]+TMath::Nint(x)-1, [2]-1) * pow([1] / ([1]+[2]), TMath::Nint(x)) * pow(1 + [1]/[2], -[2])", 0, 200);
    func->SetParNames("scaling", "averagen", "k");
  }

  switch (caseNo)
  {
    case 0: func = new TF1("flat", "1000"); break;
    case 1: func = new TF1("flat", "501-x"); break;
    case 2: func = new TF1("flat", "1000 * 1/(x+1)"); break;
    case 3: func = new TF1("flat", "1000 * TMath::Landau(x, 10, 5)"); break;
    case 4: func->SetParameters(2e5, 15, 2); break;
    case 5: func->SetParameters(1, 13, 7); break;
    case 6: func->SetParameters(1e7, 30, 4); break;
    case 7: func->SetParameters(1e7, 30, 2); break; // ***
    case 8: func = new TF1("testlaszlo", "10*1000*x*exp(-0.1*x)"); break;

    default: return;
  }

  new TCanvas;
  func->Draw();

  mult->SetGenMeasFromFunc(func, 1);

  TFile::Open("out.root", "RECREATE");
  mult->SaveHistograms();

  new TCanvas; mult->GetMultiplicityESD(1)->ProjectionY()->DrawCopy();
  new TCanvas; mult->GetMultiplicityVtx(1)->ProjectionY("mc", 1, mult->GetMultiplicityVtx(1)->GetNbinsX())->DrawCopy();

  //mult->ApplyBayesianMethod(2, kFALSE);
  //mult->ApplyMinuitFit(2, kFALSE);
  //mult->ApplyGaussianMethod(2, kFALSE);
  //mult->ApplyLaszloMethod(2, kFALSE, AliMultiplicityCorrection::kTrVtx);
}

void smoothCorrelationMap(const char* fileName = "multiplicity.root", Int_t corrMatrix = 1)
{
  gSystem->Load("libPWG0base");

  AliMultiplicityCorrection* mult = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");

  TFile::Open(fileName);
  mult->LoadHistograms("Multiplicity");

  // empty under/overflow bins in x, otherwise Project3D takes them into account
  TH3* corr = mult->GetCorrelation(corrMatrix);
  for (Int_t j=1; j<=corr->GetYaxis()->GetNbins(); ++j)
  {
    for (Int_t k=1; k<=corr->GetZaxis()->GetNbins(); ++k)
    {
      corr->SetBinContent(0, j, k, 0);
      corr->SetBinContent(corr->GetXaxis()->GetNbins()+1, j, k, 0);
    }
  }

  TH2* proj = (TH2*) corr->Project3D("zy");

  // normalize correction for given nPart
  for (Int_t i=1; i<=proj->GetNbinsX(); ++i)
  {
    Double_t sum = proj->Integral(i, i, 1, proj->GetNbinsY());
    if (sum <= 0)
      continue;

    for (Int_t j=1; j<=proj->GetNbinsY(); ++j)
    {
      // npart sum to 1
      proj->SetBinContent(i, j, proj->GetBinContent(i, j) / sum);
      proj->SetBinError(i, j, proj->GetBinError(i, j) / sum);
    }
  }

  new TCanvas;
  proj->DrawCopy("COLZ");

  TH1* scaling = proj->ProjectionY("scaling", 1, 1);
  scaling->Reset();
  scaling->SetMarkerStyle(3);
  //scaling->GetXaxis()->SetRangeUser(0, 50);
  TH1* mean = (TH1F*) scaling->Clone("mean");
  TH1* width = (TH1F*) scaling->Clone("width");

  TF1* lognormal = new TF1("lognormal", "[0]*exp(-(log(x)-[1])^2/(2*[2]^2))/(x*[2]*TMath::Sqrt(2*TMath::Pi()))", 0.01, 500);
  lognormal->SetParNames("scaling", "mean", "sigma");
  lognormal->SetParameters(1, 1, 1);
  lognormal->SetParLimits(0, 1, 1);
  lognormal->SetParLimits(1, 0, 100);
  lognormal->SetParLimits(2, 1e-3, 1);

  TF1* nbd = new TF1("nbd", "[0] * TMath::Binomial([2]+TMath::Nint(x)-1, [2]-1) * pow([1] / ([1]+[2]), TMath::Nint(x)) * pow(1 + [1]/[2], -[2])", 0, 50);
  nbd->SetParNames("scaling", "averagen", "k");
  nbd->SetParameters(1, 13, 5);
  nbd->SetParLimits(0, 1, 1);
  nbd->SetParLimits(1, 1, 100);
  nbd->SetParLimits(2, 1, 1e8);

  TF1* poisson = new TF1("poisson", "[0] * exp(-(x+[2])) * (x+[2])**[1] / TMath::Factorial([1])", 0.01, 50);
  poisson->SetParNames("scaling", "k", "deltax");
  poisson->SetParameters(1, 1, 0);
  poisson->SetParLimits(0, 0, 10);
  poisson->SetParLimits(1, 0.01, 100);
  poisson->SetParLimits(2, 0, 10);

  TF1* mygaus = new TF1("mygaus", "[0] * exp(-(x-[1])**2 / 2 / [2] - [3] * log(x + [4]) / [5])", 0.01, 50);
  mygaus->SetParNames("scaling", "mean", "width", "scale2log", "logmean", "logwidth");
  mygaus->SetParameters(1, 0, 1, 1, 0, 1);
  mygaus->SetParLimits(2, 1e-5, 10);
  mygaus->SetParLimits(4, 1, 1);
  mygaus->SetParLimits(5, 1e-5, 10);

  //TF1* sqrt = new TF1("sqrt", "[0] + [1] * sqrt((x + [3]) * [2])", 0, 50);
  TF1* sqrt = new TF1("sqrt", "[0] + (x + [1])**[2]", 0, 50);
  sqrt->SetParNames("ydelta", "exp", "xdelta");
  sqrt->SetParameters(0, 0, 1);
  sqrt->SetParLimits(1, 0, 10);

  const char* fitWith = "gaus";

  for (Int_t i=1; i<=80; ++i)
  {
    printf("Fitting %d...\n", i);

    TH1* hist = proj->ProjectionY(Form("proj%d", i), i, i, "e");

    //hist->GetXaxis()->SetRangeUser(0, 50);
    //lognormal->SetParameter(0, hist->GetMaximum());
    hist->Fit(fitWith, "0 M", "");

    TF1* func = hist->GetListOfFunctions()->FindObject(fitWith);

    if (0 && (i % 5 == 0))
    {
      pad = new TCanvas;
      hist->Draw();
      func->Clone()->Draw("SAME");
      pad->SetLogy();
    }

    scaling->SetBinContent(i, func->GetParameter(0));
    mean->SetBinContent(i, func->GetParameter(1));
    width->SetBinContent(i, func->GetParameter(2));
    
    scaling->SetBinError(i, func->GetParError(0));
    mean->SetBinError(i, func->GetParError(1));
    width->SetBinError(i, func->GetParError(2));
  }

  TF1* log = new TF1("log", "[0] + [1] * log([2] * x)", 0.01, 500);
  log->SetParameters(0, 1, 1);
  log->SetParLimits(1, 0, 100);
  log->SetParLimits(2, 1e-3, 10);

  TF1* over = new TF1("over", "[0] + [1] / (x+[2])", 0.01, 500);
  over->SetParameters(0, 1, 0);
  //over->SetParLimits(0, 0, 100);
  over->SetParLimits(1, 1e-3, 10);
  over->SetParLimits(2, 0, 100);

  c1 = new TCanvas("fitparams", "fitparams", 1200, 400);
  c1->Divide(3, 1);

  c1->cd(1);
  scaling->Draw("P");

  //TF1* scalingFit = new TF1("mypol0", "[0]");
  TF1* scalingFit = over;
  scaling->Fit(scalingFit, "", "", 3, 140);
  scalingFit->SetRange(0, 200);
  scalingFit->Draw("SAME");
  
  c1->cd(2);
  mean->Draw("P");

  //TF1* meanFit = log;
  TF1* meanFit = new TF1("mypol1", "[0]+[1]*x");
  mean->Fit(meanFit, "", "", 3, 140);
  meanFit->SetRange(0, 200);
  meanFit->Draw("SAME");

  c1->cd(3);
  width->Draw("P");

  //TF1* widthFit = over;
  TF1* widthFit = new TF1("mypol", "[0]+[1]*TMath::Sqrt([2]*x)");
  widthFit->SetParLimits(2, 1e-5, 1e5);
  width->Fit(widthFit, "", "", 5, 140);
  widthFit->SetRange(0, 200);
  widthFit->Draw("SAME");
  
  // build new correction matrix
  TH2* new = (TH2*) proj->Clone("new");
  new->Reset();
  Float_t x, y;
  for (Int_t i=1; i<=new->GetXaxis()->GetNbins(); i+=1)
  {
    TF1* func = (TF1*) gROOT->FindObject(fitWith);
    x = new->GetXaxis()->GetBinCenter(i);
    //if (i == 1)
    //  x = 0.1;
    x++;
    func->SetParameters(scalingFit->Eval(x), meanFit->Eval(x), widthFit->Eval(x));
    printf("%f %f %f %f\n", x, scalingFit->Eval(x), meanFit->Eval(x), widthFit->Eval(x));

    for (Int_t j=1; j<=new->GetYaxis()->GetNbins(); j+=1)
    {
      if (i < 11)
      {
        // leave bins 1..20 untouched
        new->SetBinContent(i, j, corr->Integral(1, corr->GetNbinsX(), i, i, j, j));
      }
      else
      {
        y = new->GetYaxis()->GetBinCenter(j);
        if (j == 1)
          y = 0.1;
        if (func->Eval(y) > 1e-4)
          new->SetBinContent(i, j, func->Eval(y));
      }
    }
  }

  // fill 0 multiplicity bins, this cannot be done with the function because it does not accept 0
  // we take the values from the old response matrix
  //for (Int_t i=1; i<=new->GetXaxis()->GetNbins(); i+=1)
  //  new->SetBinContent(i, 1, proj->GetBinContent(i, 1));

  //for (Int_t j=1; j<=new->GetYaxis()->GetNbins(); j+=1)
  //  new->SetBinContent(1, j, proj->GetBinContent(1, j));

  // normalize correction for given nPart
  for (Int_t i=1; i<=new->GetNbinsX(); ++i)
  {
    Double_t sum = new->Integral(i, i, 1, proj->GetNbinsY());
    if (sum <= 0)
      continue;

    for (Int_t j=1; j<=new->GetNbinsY(); ++j)
    {
      // npart sum to 1
      new->SetBinContent(i, j, new->GetBinContent(i, j) / sum);
      new->SetBinError(i, j, new->GetBinError(i, j) / sum);
    }
  }

  new TCanvas;
  new->Draw("COLZ");

  TH2* diff = (TH2*) new->Clone("diff");
  diff->Add(proj, -1);

  new TCanvas;
  diff->Draw("COLZ");
  diff->SetMinimum(-0.05);
  diff->SetMaximum(0.05);

  corr->Reset();

  for (Int_t i=1; i<=new->GetNbinsX(); ++i)
    for (Int_t j=1; j<=new->GetNbinsY(); ++j)
      corr->SetBinContent(corr->GetXaxis()->GetNbins() / 2 + 1, i, j, new->GetBinContent(i, j));

  new TCanvas;
  corr->Project3D("zy")->Draw("COLZ");

  TFile::Open("out.root", "RECREATE");
  mult->SaveHistograms();

  TH1* proj1 = proj->ProjectionY("proj1", 36, 36);
  TH1* proj2 = new->ProjectionY("proj2", 36, 36);
  proj2->SetLineColor(2);

  new TCanvas;
  proj1->Draw();
  proj2->Draw("SAME");
}

void buildCorrelationMap(const char* fileName = "multiplicityMC_2M.root", Int_t corrMatrix = 3)
{
  gSystem->Load("libPWG0base");

  AliMultiplicityCorrection* mult = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");

  TFile::Open(fileName);
  mult->LoadHistograms("Multiplicity");

  TH3F* new = mult->GetCorrelation(corrMatrix);
  new->Reset();

  TF1* func = new TF1("func", "gaus(0)");

  Int_t vtxBin = new->GetNbinsX() / 2;
  if (vtxBin == 0)
    vtxBin = 1;

  Float_t sigma = 2;
  for (Int_t i=1; i<=new->GetYaxis()->GetNbins(); i+=1)
  {
    Float_t x = new->GetYaxis()->GetBinCenter(i);
    func->SetParameters(1, x * 0.8, sigma);
    //func->SetParameters(1, x, sigma);

    for (Int_t j=1; j<=new->GetZaxis()->GetNbins(); j+=1)
    {
      Float_t y = new->GetYaxis()->GetBinCenter(j);

      // cut at 1 sigma
      if (TMath::Abs(y-x*0.8) < sigma)
        new->SetBinContent(vtxBin, i, j, func->Eval(y));

      // test only bin 40 has smearing
      //if (x != 40)
      //  new->SetBinContent(vtxBin, i, j, (i == j));
    }
  }

  new TCanvas;
  new->Project3D("zy")->DrawCopy("COLZ");

  TFile* file = TFile::Open("out.root", "RECREATE");
  mult->SetCorrelation(corrMatrix, new);
  mult->SaveHistograms();
  file->Close();
}

void GetCrossSections(const char* fileName)
{
  gSystem->Load("libPWG0base");

  AliMultiplicityCorrection* mult = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");

  TFile::Open(fileName);
  mult->LoadHistograms("Multiplicity");

  TH1* xSection2 = mult->GetCorrelation(3)->Project3D("y")->Clone("xSection2");
  xSection2->Sumw2();
  xSection2->Scale(1.0 / xSection2->Integral());

  TH1* xSection15 = mult->GetCorrelation(2)->Project3D("y")->Clone("xSection15");
  xSection15->Sumw2();
  xSection15->Scale(1.0 / xSection15->Integral());

  TFile::Open("crosssection.root", "RECREATE");
  xSection2->Write();
  xSection15->Write();
  gFile->Close();
}

void AnalyzeSpeciesTree(const char* fileName)
{
  //
  // prints statistics about fParticleSpecies
  //

  gSystem->Load("libPWG0base");

  TFile::Open(fileName);
  TNtuple* fParticleSpecies = (TNtuple*) gFile->Get("fParticleSpecies");

  const Int_t nFields = 8;
  Long_t totals[8];
  for (Int_t i=0; i<nFields; i++)
    totals[i] = 0;

  for (Int_t i=0; i<fParticleSpecies->GetEntries(); i++)
  {
    fParticleSpecies->GetEvent(i);

    Float_t* f = fParticleSpecies->GetArgs();

    for (Int_t j=0; j<nFields; j++)
      totals[j] += f[j+1];
  }

  for (Int_t i=0; i<nFields; i++)
    Printf("%d --> %ld", i, totals[i]);
}

void BuildResponseFromTree(const char* fileName, const char* target)
{
  //
  // builds several response matrices with different particle ratios (systematic study)
  // particle-species study
  // 
  // WARNING doesn't work uncompiled, see test.C

  loadlibs();

  TFile::Open(fileName);
  TNtuple* fParticleSpecies = (TNtuple*) gFile->Get("fParticleSpecies");

  TFile* file = TFile::Open(target, "RECREATE");
  file->Close();

  Int_t tracks = 0; // control variables
  Int_t noLabel = 0;
  Int_t secondaries = 0;
  Int_t doubleCount = 0;

  for (Int_t num = 0; num < 7; num++)
  {
    Printf("%d", num);
  
    TString name;
    name.Form("Multiplicity_%d", num);
    AliMultiplicityCorrection* fMultiplicity = new AliMultiplicityCorrection(name, name);

    Float_t ratio[4]; // pi, K, p, other
    for (Int_t i = 0; i < 4; i++)
      ratio[i] = 1;

    if (num == 1)
      ratio[1] = 0.7;
    if (num == 2)
      ratio[1] = 1.3;
    
    if (num == 3)
      ratio[2] = 0.7;
    if (num == 4)
      ratio[2] = 1.3;
    
    if (num == 5)
      ratio[3] = 0.7;
    if (num == 6)
      ratio[3] = 1.3;
      
    for (Int_t i=0; i<fParticleSpecies->GetEntries(); i++)
    {
      fParticleSpecies->GetEvent(i);
      

      Float_t* f = fParticleSpecies->GetArgs();

      Float_t gene = 0;
      Float_t meas = 0;

      for (Int_t j = 0; j < 4; j++)
      {
        if (ratio[j] == 1)
        {
          gene += f[j+1];
        }
        else
        {
          for (Int_t k=0; k<f[j+1]; k++)
          {
            if (ratio[j] < 1 && gRandom->Uniform() < ratio[j])
            {
              gene += 1;
            }
            else if (ratio[j] > 1)
            {
              gene += 1;
              if (gRandom->Uniform() < ratio[j] - 1)
                gene += 1;
            }
          }
        }
          
        if (ratio[j] == 1)
        {
          meas += f[j+1+4];
        }
        else
        {
          for (Int_t k=0; k<f[j+1+4]; k++)
          {
            if (ratio[j] < 1 && gRandom->Uniform() < ratio[j])
            {
              meas += 1;
            }
            else if (ratio[j] > 1)
            {
              meas += 1;
              if (gRandom->Uniform() < ratio[j] - 1)
                meas += 1;
            }
          }
        }
        
        tracks += f[j+1+4];
      }

      // add the ones w/o label
      tracks += f[9];
      noLabel += f[9];

      // secondaries are already part of meas!
      secondaries += f[10];

      // double counted are already part of meas!
      doubleCount += f[11];

      // ones w/o label are added without weighting to allow comparison to default analysis. however this is only valid when their fraction is low enough!
      meas += f[9];

      //Printf("%.f %.f %.f %.f %.f", f[5], f[6], f[7], f[8], f[9]);

      fMultiplicity->FillCorrection(f[0], gene, gene, gene, gene, meas, meas, meas);
      // HACK all as kND = 1
      fMultiplicity->FillGenerated(f[0], kTRUE, kTRUE, 1, gene, gene, gene, gene);
      fMultiplicity->FillMeasured(f[0], meas, meas, meas);
    }

    //fMultiplicity->DrawHistograms();

    TFile* file = TFile::Open(target, "UPDATE");
    fMultiplicity->SaveHistograms();
    file->Close();

    if (num == 0)
    {
      Printf("%d total tracks, %d w/o label = %.2f %%, %d double counted = %.2f %%, secondaries = %.2f %%", tracks, noLabel, 100.0 * noLabel / tracks, doubleCount, 100.0 * doubleCount / tracks, 100.0 * secondaries / tracks);
      if ((Float_t) noLabel / tracks > 0.02)
        Printf("WARNING: More than 2%% of tracks without label, this might bias the study!");
    }
  }
}

void ParticleRatioStudy()
{
  loadlibs();

  for (Int_t num = 0; num < 7; num++)
  {
    TString target;
    target.Form("chi2a_inel_species_%d.root", num);
    
    correct("compositions.root", Form("Multiplicity_%d", num), "multiplicityESD.root", 1, -1, 0, -1, 2, target);
  }
}

void MergeCrossSection(Int_t xsectionID, const char* output = "multiplicityMC_xsection.root")
{
  const char* files[] = { "multiplicitySD.root", "multiplicityDD.root", "multiplicityND.root" };

  loadlibs();

  TFile::Open(output, "RECREATE");
  gFile->Close();

  AliMultiplicityCorrection* data[3];
  TList list;

  Float_t ref_SD = -1;
  Float_t ref_DD = -1;
  Float_t ref_ND = -1;

  Float_t error_SD = -1;
  Float_t error_DD = -1;
  Float_t error_ND = -1;

  gROOT->ProcessLine(gSystem->ExpandPathName(".L $ALICE_ROOT/PWG0/dNdEta/drawSystematics.C"));
  GetRelativeFractions(xsectionID, ref_SD, ref_DD, ref_ND, error_SD, error_DD, error_ND);
  
  for (Int_t i=0; i<3; ++i)
  {
    TString name;
    name.Form("Multiplicity");
    if (i > 0)
      name.Form("Multiplicity_%d", i);

    TFile::Open(files[i]);
    data[i] = new AliMultiplicityCorrection(name, name);
    data[i]->LoadHistograms("Multiplicity");
  }
  
  // TODO is the under/overflow bin scaled as well? --> seems like, verify anyway!

  // calculating relative
  Float_t sd = data[0]->GetMultiplicityINEL(0)->Integral(0, data[0]->GetMultiplicityINEL(0)->GetNbinsX()+1);
  Float_t dd = data[1]->GetMultiplicityINEL(0)->Integral(0, data[1]->GetMultiplicityINEL(0)->GetNbinsX()+1);
  Float_t nd = data[2]->GetMultiplicityINEL(0)->Integral(0, data[2]->GetMultiplicityINEL(0)->GetNbinsX()+1);
  Float_t total = nd + dd + sd;
  
  nd /= total;
  sd /= total;
  dd /= total;
  
  Printf("Ratios in the correction map are: ND=%f, DD=%f, SD=%f", nd, dd, sd);
  
  Float_t ratio[3];
  ratio[0] = ref_SD / sd;
  ratio[1] = ref_DD / dd;
  ratio[2] = ref_ND / nd;
  
  Printf("SD=%.2f, DD=%.2f, ND=%.2f",ratio[0], ratio[1], ratio[2]);
  
  for (Int_t i=0; i<3; ++i)
  {
    // modify x-section
    for (Int_t j=0; j<AliMultiplicityCorrection::kMCHists; j++)
    {
      data[i]->GetMultiplicityVtx(j)->Scale(ratio[i]);
      data[i]->GetMultiplicityMB(j)->Scale(ratio[i]);
      data[i]->GetMultiplicityINEL(j)->Scale(ratio[i]);
      data[i]->GetMultiplicityNSD(j)->Scale(ratio[i]);
    }

    for (Int_t j=0; j<AliMultiplicityCorrection::kESDHists; j++)
    {
      data[i]->GetMultiplicityESD(j)->Scale(ratio[i]);
      data[i]->GetTriggeredEvents(j)->Scale(ratio[i]);
      data[i]->GetNoVertexEvents(j)->Scale(ratio[i]);
    }
    
    for (Int_t j=0; j<AliMultiplicityCorrection::kCorrHists; j++)
      data[i]->GetCorrelation(j)->Scale(ratio[i]);

    if (i > 0)
      list.Add(data[i]);
  }

  printf("Case %d, %s: Entries in response matrix: SD: %.2f DD: %.2f ND: %.2f", xsectionID, data[0]->GetName(), data[0]->GetCorrelation(0)->Integral(), data[1]->GetCorrelation(0)->Integral(), data[2]->GetCorrelation(0)->Integral());

  data[0]->Merge(&list);

  Printf(" Total: %.2f", data[0]->GetCorrelation(0)->Integral());

  TFile::Open(output, "RECREATE");
  data[0]->SaveHistograms();
  gFile->Close();

  list.Clear();

  for (Int_t i=0; i<3; ++i)
    delete data[i];
}

void Rebin(const char* fileName = "multiplicityMC_3M.root", Int_t corrMatrix = 3)
{
  // rebins MC axis of correlation map, MC and histogram for corrected (for evaluation of effect of regularization)
  // rebin does not exist for 3D hists, so we convert to 2D and then back to 3D (loosing the vertex information)

  Printf("WARNING: Vertex information is lost in this process. Use result only for evaluation of errors.");

  gSystem->Load("libPWG0base");

  AliMultiplicityCorrection* mult = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");

  TFile::Open(fileName);
  mult->LoadHistograms("Multiplicity");

  // rebin correlation
  TH3* old = mult->GetCorrelation(corrMatrix);

  // empty under/overflow bins in x, otherwise Project3D takes them into account
  for (Int_t y=1; y<=old->GetYaxis()->GetNbins(); ++y)
  {
    for (Int_t z=1; z<=old->GetZaxis()->GetNbins(); ++z)
    {
      old->SetBinContent(0, y, z, 0);
      old->SetBinContent(old->GetXaxis()->GetNbins()+1, y, z, 0);
    }
  }

  TH2* response = (TH2*) old->Project3D("zy");
  response->RebinX(2);

  TH3F* new = new TH3F(old->GetName(), old->GetTitle(),
    old->GetXaxis()->GetNbins(), old->GetXaxis()->GetBinLowEdge(1), old->GetXaxis()->GetBinUpEdge(old->GetXaxis()->GetNbins()),
    old->GetYaxis()->GetNbins() / 2, old->GetYaxis()->GetBinLowEdge(1), old->GetYaxis()->GetBinUpEdge(old->GetYaxis()->GetNbins()),
    old->GetZaxis()->GetNbins(), old->GetZaxis()->GetBinLowEdge(1), old->GetZaxis()->GetBinUpEdge(old->GetZaxis()->GetNbins()));
  new->Reset();

  Int_t vtxBin = new->GetNbinsX() / 2;
  if (vtxBin == 0)
    vtxBin = 1;

  for (Int_t i=1; i<=new->GetYaxis()->GetNbins(); i+=1)
    for (Int_t j=1; j<=new->GetZaxis()->GetNbins(); j+=1)
      new->SetBinContent(vtxBin, i, j, response->GetBinContent(i, j));

  // rebin MC + hist for corrected
  for (AliMultiplicityCorrection::EventType eventType = AliMultiplicityCorrection::kTrVtx; eventType <= AliMultiplicityCorrection::kINEL; eventType++)
    mult->GetMultiplicityMC(corrMatrix, eventType)->RebinY(2);

  mult->GetMultiplicityESDCorrected(corrMatrix)->Rebin(2);

  // recreate measured from correlation matrix to get rid of vertex shift effect
  TH2* newMeasured = (TH2*) old->Project3D("zx");
  TH2* esd = mult->GetMultiplicityESD(corrMatrix);
  esd->Reset();

  // transfer from TH2D to TH2F
  for (Int_t i=0; i<=new->GetXaxis()->GetNbins()+1; i+=1)
    for (Int_t j=0; j<=new->GetYaxis()->GetNbins()+1; j+=1)
      esd->SetBinContent(i, j, newMeasured->GetBinContent(i, j));

  new TCanvas;
  new->Project3D("zy")->DrawCopy("COLZ");

  TFile* file = TFile::Open("out.root", "RECREATE");
  mult->SetCorrelation(corrMatrix, new);
  mult->SaveHistograms();
  file->Close();
}

void EvaluateRegularizationEffect(Int_t step, const char* fileNameRebinned = "multiplicityMC_3M_rebinned.root", const char* fileNameNormal = "multiplicityMC_3M.root", Int_t histID = 3)
{
  // due to some static members in AliMultiplicityCorrection, the session has to be restarted after changing the number of parameters, to be fixed
  // that is why this is done in 2 steps

  gSystem->Load("libPWG0base");

  Bool_t fullPhaseSpace = kFALSE;

  if (step == 1)
  {
    // first step: unfold without regularization and rebinned histogram ("N=M")
    AliMultiplicityCorrection* mult = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");
    TFile::Open(fileNameRebinned);
    mult->LoadHistograms();

    mult->SetRegularizationParameters(AliMultiplicityCorrection::kNone, 0, 125);
    mult->SetCreateBigBin(kFALSE);

    mult->ApplyMinuitFit(histID, fullPhaseSpace, AliMultiplicityCorrection::kTrVtx, kFALSE);
    mult->DrawComparison("MinuitChi2", histID, fullPhaseSpace, kTRUE, mult->GetMultiplicityVtx(histID)->ProjectionY("mymchist"));

    TFile* file = TFile::Open("EvaluateRegularizationEffect1.root", "RECREATE");
    mult->SaveHistograms();
    file->Close();
  }
  else if (step == 2)
  {
    // second step: unfold with regularization and normal histogram
    AliMultiplicityCorrection* mult2 = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");
    TFile::Open(fileNameNormal);
    mult2->LoadHistograms();

    mult2->SetRegularizationParameters(AliMultiplicityCorrection::kPol1, 1e4);
    mult2->SetCreateBigBin(kTRUE);
    mult2->ApplyMinuitFit(histID, fullPhaseSpace, AliMultiplicityCorrection::kTrVtx, kFALSE);
    mult2->DrawComparison("MinuitChi2", histID, fullPhaseSpace, kTRUE, mult2->GetMultiplicityVtx(histID)->ProjectionY("mymchist"));

    TH1* result2 = mult2->GetMultiplicityESDCorrected(histID);

    AliMultiplicityCorrection* mult = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");
    TFile* file = TFile::Open("EvaluateRegularizationEffect1.root");
    mult->LoadHistograms();

    TH1* result1 = mult->GetMultiplicityESDCorrected(histID);

    // compare results
    TCanvas* canvas = new TCanvas("EvaluateRegularizationEffect", "EvaluateRegularizationEffect", 1000, 800);
    canvas->Divide(2, 2);

    canvas->cd(1);
    result1->SetLineColor(1);
    result1->DrawCopy();
    result2->SetLineColor(2);
    result2->DrawCopy("SAME");
    gPad->SetLogy();

    result2->Rebin(2);
    result1->Scale(1.0 / result1->Integral());
    result2->Scale(1.0 / result2->Integral());

    canvas->cd(2);
    result1->DrawCopy();
    result2->DrawCopy("SAME");
    gPad->SetLogy();

    TH1* diff = (TH1*) result1->Clone("diff");
    diff->Add(result2, -1);

    canvas->cd(3);
    diff->DrawCopy("HIST");

    canvas->cd(4);
    diff->Divide(result1);
    diff->GetYaxis()->SetRangeUser(-0.3, 0.3);
    diff->DrawCopy("HIST");

    Double_t chi2 = 0;
    for (Int_t i=1; i<=diff->GetNbinsX(); i++)
      chi2 += diff->GetBinContent(i) * diff->GetBinContent(i);

    Printf("Chi2 is %e", chi2);

    canvas->SaveAs(Form("%s.eps", canvas->GetName()));
  }
}

void MergeDistributions()
{
  loadlibs();

  const char* dir1 = "run000104824-52_pass4";
  const char* dir2 = "run000104867_90_92_pass4";
  
  for (Int_t evType = 0; evType < 2; evType++)
  {
    Printf("%d", evType);
    
    const char* evTypeStr = ((evType == 0) ? "inel" : "nsd");

    const char* id = "chi2a";
    //const char* id = "bayes";
    
    TString suffix;
    suffix.Form("/all/spd/%s_", id);
    if (evType == 1)
      suffix.Form("/v0and/spd/%s_", id);
  
    TString file1, file2;
    file1.Form("%s%s%%s.root", dir1, suffix.Data());
    file2.Form("%s%s%%s.root", dir2, suffix.Data());
    
    if (1)
    {
      const char* files[] = { Form(file1.Data(), evTypeStr), Form(file2.Data(), evTypeStr) };
      Merge(2, files, Form("merged/%s_%s.root", id, evTypeStr));
    }
    else
    {
      AliMultiplicityCorrection* mult1 = AliMultiplicityCorrection::Open(Form(file1.Data(), evTypeStr));
      AliMultiplicityCorrection* mult2 = AliMultiplicityCorrection::Open(Form(file2.Data(), evTypeStr));
      
      AliMultiplicityCorrection* target = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");
      
      for (Int_t etaRange = 0; etaRange < 3; etaRange++)
      {
        hist1 = mult1->GetMultiplicityESDCorrected(etaRange);
        hist2 = mult2->GetMultiplicityESDCorrected(etaRange);
        targetHist = target->GetMultiplicityESDCorrected(etaRange);
        
        //hist1->Scale(1.0 / hist1->Integral());
        //hist2->Scale(1.0 / hist2->Integral());
        
        //targetHist->SetBinContent(1, hist1->GetBinContent(1) * 0.5 + hist2->GetBinContent(1) * 0.5);
        targetHist->SetBinContent(1, hist1->GetBinContent(1) + hist2->GetBinContent(1));
        for (Int_t i=1; i<=hist1->GetNbinsX(); i++)
        {
          if (hist1->GetBinError(i) > 0 && hist2->GetBinError(i) > 0)
          {
            Float_t w1 = 1.0 / hist1->GetBinError(i) / hist1->GetBinError(i);
            Float_t w2 = 1.0 / hist2->GetBinError(i) / hist2->GetBinError(i);
            
            Float_t average = (hist1->GetBinContent(i) * w1 + hist2->GetBinContent(i) * w2) / (w1 + w2);
            
            //targetHist->SetBinContent(i, average);
            //targetHist->SetBinError(i, TMath::Max(mult1->GetBinError(i), mult2->GetBinError(i)));
            
            targetHist->SetBinContent(i, hist1->GetBinContent(i) + hist2->GetBinContent(i));
            targetHist->SetBinError(i, hist1->GetBinError(i) + hist2->GetBinError(i));
          }
        }
      }
      
      file = TFile::Open(Form("merged/%s_%s.root", id, evTypeStr), "RECREATE");
      target->SaveHistograms();
      file->Close();
    }

    const char* files2[] = { Form(file1.Data(), evTypeStr), Form(file2.Data(), evTypeStr), Form("merged/%s_%s.root", id, evTypeStr) };
    drawAll(files2, (evType == 1), kTRUE);
  }
}

void CompareDistributions(Int_t evType)
{
  loadlibs();
  gROOT->ProcessLine(gSystem->ExpandPathName(".L $ALICE_ROOT/PWG0/dNdEta/drawPlots.C"));

  const char* dir1 = "run000104824-52_pass4";
  const char* dir2 = "run000104867_90_92_pass4";
  
  const char* evTypeStr = (evType == 0) ? "inel" : "nsd";

  const char* suffix = "/all/spd/chi2a_";
  if (evType == 1)
    suffix = "/v0and/spd/chi2a_";
  
  TString file1, file2;
  file1.Form("%s%s%s.root", dir1, suffix, evTypeStr);
  file2.Form("%s%s%s.root", dir2, suffix, evTypeStr);
    
  for (Int_t hID = 0; hID < 3; hID++)
    CompareQualityHists(file1, file2, Form("Multiplicity/fMultiplicityESDCorrected%d", hID), 1, 1);
}

void StatisticalUncertainty(Int_t methodType, Int_t etaRange, Bool_t mc = kFALSE)
{
  loadlibs();
  
  AliMultiplicityCorrection* mult = AliMultiplicityCorrection::Open("multiplicityMC.root");
  AliMultiplicityCorrection* esd = AliMultiplicityCorrection::Open("multiplicityESD.root");
  AliMultiplicityCorrection* multTrigger = AliMultiplicityCorrection::Open("multiplicityTrigger.root");

  TH1* mcHist = esd->GetMultiplicityMB(etaRange)->ProjectionY("mymc", 1, 1);

  Int_t geneLimits[] = { 0, 0, 0 };

  LoadAndInitialize(mult, esd, multTrigger, etaRange, kFALSE, geneLimits);

  AliUnfolding::SetNbins(kBinLimits[etaRange], geneLimits[etaRange]);
  AliUnfolding::SetChi2Regularization(AliUnfolding::kPol0, 5);
  AliUnfolding::SetBayesianParameters(1, 10);

  // background subtraction
  Int_t background = 0;
  
  //background = 1091 + 4398; // MB1 for run 104824 - 52
  background = 417 + 1758; // MB1 for run 104867 - 92
  
  //background = 20 + 0; // V0AND for run 104824 - 52
  //background = 10 + 0; // V0AND for run 104867 - 92
  
  Printf("NOTE: Subtracting %d background events", background);

  Int_t zeroBin = mult->GetTriggeredEvents(etaRange)->GetBinContent(1) - background - mult->GetMultiplicityESD(etaRange)->Integral(0, 2, 1, 1);
  
  TH1* errorMeasured = (TH1*) mult->StatisticalUncertainty((AliUnfolding::MethodType) methodType, etaRange, kFALSE, AliMultiplicityCorrection::kMB, zeroBin, kTRUE, kFALSE, ((mc) ? mcHist : 0))->Clone("errorMeasured");
  
  new TCanvas; errorMeasured->Draw();
  new TCanvas; 

  AliPWG0Helper::NormalizeToBinWidth(mult->GetMultiplicityESDCorrected(etaRange));
  mult->GetMultiplicityESDCorrected(etaRange)->Draw();
  
  if (mc)
  {
    AliPWG0Helper::NormalizeToBinWidth(mcHist);
    mcHist->Scale(mult->GetMultiplicityESDCorrected(etaRange)->Integral() / mcHist->Integral()); 
    mcHist->SetLineColor(2);
    mcHist->Draw("SAME");
  }
  gPad->SetLogy();
  
  TH1* errorResponse = (TH1*) mult->StatisticalUncertainty((AliUnfolding::MethodType) methodType, etaRange, kFALSE, AliMultiplicityCorrection::kMB, zeroBin, kFALSE, kTRUE, ((mc) ? mcHist : 0))->Clone("errorResponse");

  TH1* errorBoth = (TH1*) mult->StatisticalUncertainty((AliUnfolding::MethodType) methodType, etaRange, kFALSE, AliMultiplicityCorrection::kMB, zeroBin, kTRUE, kTRUE, ((mc) ? mcHist : 0))->Clone("errorBoth");

  if (mc)
  {
    TH1* result = mult->GetMultiplicityESDCorrected(etaRange);
    DrawResultRatio(mcHist, result, "StatisticalUncertainty2.eps");
  }

  TFile* file = new TFile(Form("StatisticalUncertaintySPD%s.root", (methodType == 0) ? "Chi2" : "Bayesian"), "RECREATE");
  errorResponse->Write();
  errorMeasured->Write();
  errorBoth->Write();
  file->Close();
}

void ChangeResponseMatrixEfficiency(Bool_t reduceEff = kTRUE, Float_t factor = 0.01625)
{
  loadlibs();
  
  AliMultiplicityCorrection* mult = AliMultiplicityCorrection::Open("multiplicity.root");
  
  for (Int_t etaR = 0; etaR < 3; etaR++)
  {
    TH3* corr = mult->GetCorrelation(etaR);
    TH3* corrOld = (TH3*) corr->Clone("corrOld");
    
    corr->Reset();
    
    Int_t total = 0;
    Int_t change = 0;
    for (Int_t x=0; x<=corr->GetNbinsX()+1; x++)
    {
      for (Int_t y=0; y<=corr->GetNbinsY()+1; y++)
      {
        Printf("%d", y);
        for (Int_t z=0; z<=corr->GetNbinsZ()+1; z++)
        {
          Float_t tracklets = corr->GetZaxis()->GetBinCenter(z);
          
          for (Int_t n=0; n<corrOld->GetBinContent(x, y, z); n++)
          {
            total += tracklets;
            Float_t trackletsNew = tracklets;
            
            for (Int_t i=0; i<tracklets; i++)
            {
              if (gRandom->Uniform() < factor)
              {
                trackletsNew += (reduceEff) ? -1 : 1;
                change += (reduceEff) ? -1 : 1;
              }
            }
                
            corr->Fill(corr->GetXaxis()->GetBinCenter(x), corr->GetYaxis()->GetBinCenter(y), trackletsNew);
          }
        }
      }
    }
    
    Printf("%d: Changed %d out of %d total tracklets", etaR, change, total);
  }

  file = TFile::Open("out.root", "RECREATE");
  mult->SaveHistograms();
  file->Close();
}

