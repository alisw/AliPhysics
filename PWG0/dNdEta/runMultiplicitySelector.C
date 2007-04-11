/* $Id$ */

//
// script to run the AliMultiplicityESDSelector
//

#include "../CreateESDChain.C"
#include "../PWG0Helper.C"

void runMultiplicitySelector(Char_t* data, Int_t nRuns=20, Int_t offset=0, Bool_t aMC = kFALSE, Bool_t aDebug = kFALSE, Bool_t aProof = kFALSE, const char* option = "")
{
  if (aProof)
    connectProof("jgrosseo@lxb6046");

  TString libraries("libEG;libGeom;libESD;libPWG0base");
  TString packages("PWG0base");

  if (aMC != kFALSE)
  {
    libraries += ";libVMC;libMinuit;libSTEER;libPWG0dep;libEVGEN;libFASTSIM;libmicrocern;libpdf;libpythia6;libEGPythia6;libAliPythia6";
    packages += ";PWG0dep";
  }

  if (!prepareQuery(libraries, packages, kTRUE))
    return;

  gROOT->ProcessLine(".L CreateCuts.C");
  gROOT->ProcessLine(".L drawPlots.C");

  // selection of esd tracks
  AliESDtrackCuts* esdTrackCuts = CreateTrackCuts();
  if (!esdTrackCuts)
  {
    printf("ERROR: esdTrackCuts could not be created\n");
    return;
  }

  TList inputList;
  inputList.Add(esdTrackCuts);

  TChain* chain = CreateESDChain(data, nRuns, offset, kFALSE, kFALSE);

  TString selectorName = ((aMC == kFALSE) ? "AliMultiplicityESDSelector" : "AliMultiplicityMCSelector");
  AliLog::SetClassDebugLevel(selectorName, AliLog::kInfo);

  selectorName += ".cxx+";

  if (aDebug != kFALSE)
    selectorName += "g";

  Int_t result = executeQuery(chain, &inputList, selectorName, option);

  if (result != 0)
  {
    printf("ERROR: Executing process failed with %d.\n", result);
    return;
  }
}

void draw(const char* fileName = "multiplicityMC.root")
{
  gSystem->Load("libPWG0base");

  AliMultiplicityCorrection* mult = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");

  TFile::Open(fileName);
  mult->LoadHistograms("Multiplicity");

  mult->DrawHistograms();
}

void* fit(const char* fileName = "multiplicityMC.root", Int_t hist = 2, Int_t eventType = 0)
{
  gSystem->Load("libPWG0base");

  AliMultiplicityCorrection* mult = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");

  TFile::Open(fileName);
  mult->LoadHistograms("Multiplicity");

  //mult->ApplyLaszloMethod(hist, kFALSE, AliMultiplicityCorrection::kTrVtx);

  //return;


  mult->ApplyBayesianMethod(hist, kFALSE, AliMultiplicityCorrection::kTrVtx);
  mult->DrawComparison("Bayesian", hist, kFALSE, kTRUE, mult->GetMultiplicityMC(hist, AliMultiplicityCorrection::kTrVtx)->ProjectionY());

  return;

  TStopwatch timer;

  timer.Start();

  mult->SetRegularizationParameters(AliMultiplicityCorrection::kTest, 1);
  mult->ApplyMinuitFit(hist, kFALSE, AliMultiplicityCorrection::kTrVtx);
  mult->DrawComparison("MinuitChi2", hist, kFALSE, kFALSE, mult->GetMultiplicityMC(hist, AliMultiplicityCorrection::kTrVtx)->ProjectionY());

  timer.Stop();
  timer.Print();

  return 0;

  //mult->ApplyGaussianMethod(hist, kFALSE);

  mult->SetRegularizationParameters(AliMultiplicityCorrection::kNone, 0);
  mult->ApplyNBDFit(hist, kFALSE);
  mult->DrawComparison("NBDChi2Fit", hist, kFALSE, kTRUE, mult->GetMultiplicityMC(hist, eventType)->ProjectionY());

  return mult;
}

void* fitOther(const char* fileNameMC = "multiplicityMC_3M.root", const char* fileNameESD = "multiplicityMC_3M.root", Int_t histID = 2)
{
  gSystem->Load("libPWG0base");

  AliMultiplicityCorrection* mult = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");

  TFile::Open(fileNameMC);
  mult->LoadHistograms("Multiplicity");

  TFile::Open(fileNameESD);
  TH2F* hist = (TH2F*) gFile->Get(Form("Multiplicity/fMultiplicityESD%d", histID));
  TH2F* hist2 = (TH2F*) gFile->Get(Form("Multiplicity/fMultiplicityVtx%d", histID));
  //hist2 = (TH2F*) gFile->Get(Form("Multiplicity/fMultiplicityMB%d", histID));

  mult->SetMultiplicityESD(histID, hist);

  mult->SetMultiplicityVtx(histID, hist2);
  mult->ApplyLaszloMethod(histID, kFALSE, AliMultiplicityCorrection::kTrVtx);
  return;

  mult->SetRegularizationParameters(AliMultiplicityCorrection::kTest, 1.1);
  mult->ApplyMinuitFit(histID, kFALSE, AliMultiplicityCorrection::kTrVtx);
  mult->DrawComparison("MinuitChi2", histID, kFALSE, kTRUE, hist2->ProjectionY());

  return;

  //mult->ApplyGaussianMethod(histID, kFALSE);

  for (Float_t f=0.1; f<=0.11; f+=0.05)
  {
    mult->ApplyBayesianMethod(histID, kFALSE, AliMultiplicityCorrection::kTrVtx, f);
    mult->DrawComparison(Form("Bayesian_%f", f), histID, kFALSE, kTRUE, hist2->ProjectionY());
  }

  //mult->SetRegularizationParameters(AliMultiplicityCorrection::kEntropy, 1e7);
  //mult->ApplyMinuitFit(histID, kFALSE);
  //mult->DrawComparison("MinuitChi2", histID, kFALSE, kTRUE, hist2->ProjectionY());


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
    case AliMultiplicityCorrection::kTest   :   return "Test"; break;
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

void EvaluateBayesianMethod(const char* fileNameMC = "multiplicityMC.root", const char* fileNameESD = "multiplicityMC.root", const char* targetDir, Int_t histID = 2)
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

  TCanvas* canvas = new TCanvas("EvaluateBayesianMethod", "EvaluateBayesianMethod", 800, 600);
  TLegend* legend = new TLegend(0.2, 0.7, 0.58, 0.98);
  legend->SetFillColor(0);

  Float_t min = 1e10;
  Float_t max = 0;

  TGraph* first = 0;
  Int_t count = 0; // just to order the saved images...

  for (AliMultiplicityCorrection::EventType type = AliMultiplicityCorrection::kTrVtx; type <= AliMultiplicityCorrection::kINEL; ++type)
  {
    TGraph* fitResultsMC = new TGraph;
    fitResultsMC->SetTitle(";Weight Parameter");
    TGraph* fitResultsRes = new TGraph;
    fitResultsRes->SetTitle(";Weight Parameter");

    fitResultsMC->SetFillColor(0);
    fitResultsRes->SetFillColor(0);
    fitResultsMC->SetMarkerStyle(20+type);
    fitResultsRes->SetMarkerStyle(24+type);
    fitResultsRes->SetMarkerColor(kRed);
    fitResultsRes->SetLineColor(kRed);

    legend->AddEntry(fitResultsMC, Form("%s MC chi2", GetEventTypeName(type)));
    legend->AddEntry(fitResultsRes, Form("%s residual chi2", GetEventTypeName(type)));

    for (Float_t weight = 0; weight < 0.301; weight += 0.02)
    {
      Float_t chi2MC = 0;
      Float_t residuals = 0;

      mult->ApplyBayesianMethod(histID, kFALSE, type, weight);
      mult->DrawComparison(Form("%s/Bayesian_%02d_%d_%f", targetDir, count++, type, weight), histID, kFALSE, kTRUE, multESD->GetMultiplicityMC(histID, type)->ProjectionY());
      mult->GetComparisonResults(&chi2MC, 0, &residuals);

      fitResultsMC->SetPoint(fitResultsMC->GetN(), weight, chi2MC);
      fitResultsRes->SetPoint(fitResultsRes->GetN(), weight, residuals);

      min = TMath::Min(min, TMath::Min(chi2MC, residuals));
      max = TMath::Max(max, TMath::Max(chi2MC, residuals));
    }

    fitResultsMC->Print();
    fitResultsRes->Print();

    canvas->cd();
    fitResultsMC->Draw(Form("%s CP", (first == 0) ? "A" : "SAME"));
    fitResultsRes->Draw("SAME CP");

    if (first == 0)
      first = fitResultsMC;
  }

  gPad->SetLogy();
  printf("min = %f, max = %f\n", min, max);
  if (min <= 0)
    min = 1e-5;
  first->GetYaxis()->SetRangeUser(min * 0.5, max * 1.5);

  legend->Draw();

  canvas->SaveAs(Form("%s/%s.gif", targetDir, canvas->GetName()));
  canvas->SaveAs(Form("%s/%s.C", targetDir, canvas->GetName()));
}

void EvaluateBayesianMethodIterations(const char* fileNameMC = "multiplicityMC.root", const char* fileNameESD = "multiplicityMC.root", const char* targetDir, Int_t histID = 2)
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

  TCanvas* canvas = new TCanvas("EvaluateBayesianMethodIterations", "EvaluateBayesianMethodIterations", 800, 600);
  TLegend* legend = new TLegend(0.2, 0.7, 0.58, 0.98);
  legend->SetFillColor(0);

  Float_t min = 1e10;
  Float_t max = 0;

  TGraph* first = 0;
  Int_t count = 0; // just to order the saved images...

  for (AliMultiplicityCorrection::EventType type = AliMultiplicityCorrection::kTrVtx; type <= AliMultiplicityCorrection::kINEL; ++type)
  {
    TGraph* fitResultsMC = new TGraph;
    fitResultsMC->SetTitle(";Iterations");
    TGraph* fitResultsRes = new TGraph;
    fitResultsRes->SetTitle(";Iterations");

    fitResultsMC->SetFillColor(0);
    fitResultsRes->SetFillColor(0);
    fitResultsMC->SetMarkerStyle(20+type);
    fitResultsRes->SetMarkerStyle(24+type);
    fitResultsRes->SetMarkerColor(kRed);
    fitResultsRes->SetLineColor(kRed);

    legend->AddEntry(fitResultsMC, Form("%s MC chi2", GetEventTypeName(type)));
    legend->AddEntry(fitResultsRes, Form("%s residual chi2", GetEventTypeName(type)));

    for (Int_t iter = 5; iter <= 50; iter += 5)
    {
      Float_t chi2MC = 0;
      Float_t residuals = 0;

      mult->ApplyBayesianMethod(histID, kFALSE, type, 0.1, iter);
      mult->DrawComparison(Form("%s/BayesianIter_%02d_%d_%d", targetDir, count++, type, iter), histID, kFALSE, kTRUE, multESD->GetMultiplicityMC(histID, type)->ProjectionY());
      mult->GetComparisonResults(&chi2MC, 0, &residuals);

      fitResultsMC->SetPoint(fitResultsMC->GetN(), iter, chi2MC);
      fitResultsRes->SetPoint(fitResultsRes->GetN(), iter, residuals);

      min = TMath::Min(min, TMath::Min(chi2MC, residuals));
      max = TMath::Max(max, TMath::Max(chi2MC, residuals));
    }

    fitResultsMC->Print();
    fitResultsRes->Print();

    canvas->cd();
    fitResultsMC->Draw(Form("%s CP", (first == 0) ? "A" : "SAME"));
    fitResultsRes->Draw("SAME CP");

    if (first == 0)
      first = fitResultsMC;
  }

  gPad->SetLogy();
  printf("min = %f, max = %f\n", min, max);
  if (min <= 0)
    min = 1e-5;
  first->GetYaxis()->SetRangeUser(min * 0.5, max * 1.5);

  legend->Draw();

  canvas->SaveAs(Form("%s/%s.gif", targetDir, canvas->GetName()));
  canvas->SaveAs(Form("%s/%s.C", targetDir, canvas->GetName()));
}

void EvaluateChi2Method(const char* fileNameMC = "multiplicityMC.root", const char* fileNameESD = "multiplicityMC.root", const char* targetDir, Int_t histID = 2)
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

  TCanvas* canvas = new TCanvas("EvaluateChi2Method", "EvaluateChi2Method", 800, 600);
  TLegend* legend = new TLegend(0.2, 0.7, 0.58, 0.98);
  legend->SetFillColor(0);

  Float_t min = 1e10;
  Float_t max = 0;

  TGraph* first = 0;
  Int_t count = 0; // just to order the saved images...

  Bool_t firstLoop = kTRUE;

  for (AliMultiplicityCorrection::RegularizationType type = AliMultiplicityCorrection::kPol0; type <= AliMultiplicityCorrection::kEntropy; ++type)
  //for (AliMultiplicityCorrection::RegularizationType type = AliMultiplicityCorrection::kPol1; type <= AliMultiplicityCorrection::kPol1; ++type)
  {
    TGraph* fitResultsMC = new TGraph;
    fitResultsMC->SetTitle(";Weight Parameter");
    TGraph* fitResultsRes = new TGraph;
    fitResultsRes->SetTitle(";Weight Parameter");

    fitResultsMC->SetFillColor(0);
    fitResultsRes->SetFillColor(0);
    fitResultsMC->SetMarkerStyle(19+type);
    fitResultsRes->SetMarkerStyle(23+type);
    fitResultsRes->SetMarkerColor(kRed);
    fitResultsRes->SetLineColor(kRed);

    legend->AddEntry(fitResultsMC, Form("%s MC chi2", GetRegName(type)));
    legend->AddEntry(fitResultsRes, Form("%s residual chi2", GetRegName(type)));

    if (first == 0)
      first = fitResultsMC;

    for (Float_t weight = 1e-4; weight < 2e4; weight *= 10)
    //for (Float_t weight = 0.1; weight < 10; weight *= TMath::Sqrt(TMath::Sqrt(10)))
    {
      Float_t chi2MC = 0;
      Float_t residuals = 0;

      mult->SetRegularizationParameters(type, weight);
      mult->ApplyMinuitFit(histID, kFALSE, AliMultiplicityCorrection::kTrVtx);
      mult->DrawComparison(Form("%s/MinuitChi2_%02d_%d_%f", targetDir, count++, type, weight), histID, kFALSE, kTRUE, hist2->ProjectionY());
      mult->GetComparisonResults(&chi2MC, 0, &residuals);

      fitResultsMC->SetPoint(fitResultsMC->GetN(), weight, chi2MC);
      fitResultsRes->SetPoint(fitResultsRes->GetN(), weight, residuals);

      min = TMath::Min(min, TMath::Min(chi2MC, residuals));
      max = TMath::Max(max, TMath::Max(chi2MC, residuals));
    }

    fitResultsMC->Print();
    fitResultsRes->Print();

    canvas->cd();
    fitResultsMC->Draw(Form("%s CP", (firstLoop) ? "A" : "SAME"));
    fitResultsRes->Draw("SAME CP");

    firstLoop = kFALSE;
  }

  gPad->SetLogx();
  gPad->SetLogy();
  printf("min = %f, max = %f\n", min, max);
  if (min <= 0)
    min = 1e-5;
  first->GetYaxis()->SetRangeUser(min * 0.5, max * 1.5);

  legend->Draw();

  canvas->SaveAs(Form("%s/%s.gif", targetDir, canvas->GetName()));
  canvas->SaveAs(Form("%s/%s.C", targetDir, canvas->GetName()));
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

void CompareMethods(const char* fileNameMC = "multiplicityMC.root", const char* fileNameESD = "multiplicityMC.root", const char* targetDir, Int_t histID = 2)
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

  TCanvas* canvas = new TCanvas("CompareMethods", "CompareMethods", 1200, 800);
  canvas->Divide(3, 2);

  Int_t count = 0;

  for (AliMultiplicityCorrection::EventType type = AliMultiplicityCorrection::kTrVtx; type <= AliMultiplicityCorrection::kINEL; ++type)
  {
    TH1* mc = multESD->GetMultiplicityMC(histID, type)->ProjectionY();

    mult->SetRegularizationParameters(AliMultiplicityCorrection::kPol1, 0.3);
    mult->ApplyMinuitFit(histID, kFALSE, type);
    mult->DrawComparison(Form("%s/CompareMethods_%d_MinuitChi2", targetDir, type), histID, kFALSE, kTRUE, mc);
    TH1* chi2Result = (TH1*) mult->GetMultiplicityESDCorrected(histID)->Clone("chi2Result");

    mult->ApplyBayesianMethod(histID, kFALSE, type, 0.1);
    mult->DrawComparison(Form("%s/CompareMethods_%d_Bayesian", targetDir, type), histID, kFALSE, kTRUE, mc);
    TH1* bayesResult = (TH1*) mult->GetMultiplicityESDCorrected(histID)->Clone("bayesResult");

    mc->GetXaxis()->SetRangeUser(0, 150);
    chi2Result->GetXaxis()->SetRangeUser(0, 150);

    // skip errors for now
    for (Int_t i=1; i<=chi2Result->GetNbinsX(); ++i)
    {
      chi2Result->SetBinError(i, 0);
      bayesResult->SetBinError(i, 0);
    }

    canvas->cd(++count);
    mc->SetFillColor(kYellow);
    mc->DrawCopy();
    chi2Result->SetLineColor(kRed);
    chi2Result->DrawCopy("SAME");
    bayesResult->SetLineColor(kBlue);
    bayesResult->DrawCopy("SAME");
    gPad->SetLogy();

    canvas->cd(count + 3);
    chi2Result->Divide(chi2Result, mc, 1, 1, "B");
    bayesResult->Divide(bayesResult, mc, 1, 1, "B");

    // skip errors for now
    for (Int_t i=1; i<=chi2Result->GetNbinsX(); ++i)
    {
      chi2Result->SetBinError(i, 0);
      bayesResult->SetBinError(i, 0);
    }

    chi2Result->GetYaxis()->SetRangeUser(0.5, 1.5);

    chi2Result->DrawCopy("");
    bayesResult->DrawCopy("SAME");
  }

  canvas->SaveAs(Form("%s/%s.gif", targetDir, canvas->GetName()));
  canvas->SaveAs(Form("%s/%s.C", targetDir, canvas->GetName()));
}

void StatisticsPlot(const char* fileNameMC = "multiplicityMC.root", Int_t histID = 2)
{
  gSystem->Load("libPWG0base");

  AliMultiplicityCorrection* mult = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");

  TFile::Open(fileNameMC);
  mult->LoadHistograms("Multiplicity");

  const char* files[] = { "multiplicityMC_0.5M.root", "multiplicityMC_0.75M.root", "multiplicityMC_1M_3.root", "multiplicityMC_1.25M.root", "multiplicityMC_1.5M.root" };

  TGraph* fitResultsChi2 = new TGraph;
  fitResultsChi2->SetTitle(";Nevents;Chi2");
  TGraph* fitResultsBayes = new TGraph;
  fitResultsBayes->SetTitle(";Nevents;Chi2");
  TGraph* fitResultsChi2Limit = new TGraph;
  fitResultsChi2Limit->SetTitle(";Nevents;Multiplicity reach");
  TGraph* fitResultsBayesLimit = new TGraph;
  fitResultsBayesLimit->SetTitle(";Nevents;Multiplicity reach");

  TCanvas* canvas = new TCanvas("StatisticsPlot", "StatisticsPlot", 1200, 600);
  canvas->Divide(5, 2);

  Float_t min = 1e10;
  Float_t max = 0;

  for (Int_t i=0; i<5; ++i)
  {
    TFile::Open(files[i]);
    AliMultiplicityCorrection* multESD = new AliMultiplicityCorrection("MultiplicityESD", "MultiplicityESD");
    multESD->LoadHistograms("Multiplicity");

    mult->SetMultiplicityESD(histID, multESD->GetMultiplicityESD(histID));
    Int_t nEntries = multESD->GetMultiplicityESD(histID)->GetEntries();
    TH1* mc = multESD->GetMultiplicityMC(histID, AliMultiplicityCorrection::kTrVtx)->ProjectionY();

    mult->SetRegularizationParameters(AliMultiplicityCorrection::kPol1, 0.3);
    mult->ApplyMinuitFit(histID, kFALSE, AliMultiplicityCorrection::kTrVtx);
    mult->DrawComparison(Form("StatisticsPlot_%d_MinuitChi2", i), histID, kFALSE, kTRUE, mc);

    Float_t chi2MC = 0;
    Int_t chi2MCLimit = 0;
    mult->GetComparisonResults(&chi2MC, &chi2MCLimit, 0);
    fitResultsChi2->SetPoint(fitResultsChi2->GetN(), nEntries, chi2MC);
    fitResultsChi2Limit->SetPoint(fitResultsChi2Limit->GetN(), nEntries, chi2MCLimit);
    min = TMath::Min(min, chi2MC);
    max = TMath::Max(max, chi2MC);

    TH1* chi2Result = (TH1*) mult->GetMultiplicityESDCorrected(histID)->Clone("chi2Result");

    mult->ApplyBayesianMethod(histID, kFALSE, AliMultiplicityCorrection::kTrVtx, 0.1);
    mult->DrawComparison(Form("StatisticsPlot_%d_Bayesian", i), histID, kFALSE, kTRUE, mc);
    TH1* bayesResult = (TH1*) mult->GetMultiplicityESDCorrected(histID)->Clone("bayesResult");
    mult->GetComparisonResults(&chi2MC, &chi2MCLimit, 0);
    fitResultsBayes->SetPoint(fitResultsBayes->GetN(), nEntries, chi2MC);
    fitResultsBayesLimit->SetPoint(fitResultsBayesLimit->GetN(), nEntries, chi2MCLimit);

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
  }

  canvas->SaveAs(Form("%s.gif", canvas->GetName()));
  canvas->SaveAs(Form("%s.C", canvas->GetName()));

  TCanvas* canvas2 = new TCanvas("StatisticsPlot2", "StatisticsPlot2", 800, 400);
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
  canvas2->SaveAs(Form("%s.C", canvas2->GetName()));
}

void StartingConditions(const char* fileNameMC = "multiplicityMC.root", Int_t histID = 2)
{
  gSystem->Load("libPWG0base");

  AliMultiplicityCorrection* mult = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");

  TFile::Open(fileNameMC);
  mult->LoadHistograms("Multiplicity");

  const char* files[] = { "multiplicityMC_0.5M.root", "multiplicityMC_0.75M.root", "multiplicityMC_1M_3.root", "multiplicityMC_1.25M.root", "multiplicityMC_1.5M.root" };

  // this one we try to unfold
  TFile::Open(files[0]);
  AliMultiplicityCorrection* multESD = new AliMultiplicityCorrection("MultiplicityESD", "MultiplicityESD");
  multESD->LoadHistograms("Multiplicity");
  mult->SetMultiplicityESD(histID, multESD->GetMultiplicityESD(histID));
  TH1* mc = multESD->GetMultiplicityMC(histID, AliMultiplicityCorrection::kTrVtx)->ProjectionY();

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
      func = new TF1("nbd", "[0] * TMath::Binomial([2]+TMath::Nint(x)-1, [2]-1) * pow([1] / ([1]+[2]), TMath::Nint(x)) * pow(1 + [1]/[2], -[2])", 0, 50);
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

    mult->SetRegularizationParameters(AliMultiplicityCorrection::kPol1, 0.3);
    mult->ApplyMinuitFit(histID, kFALSE, AliMultiplicityCorrection::kTrVtx, kFALSE, startCond);
    mult->DrawComparison(Form("StartingConditions_%d_MinuitChi2", i), histID, kFALSE, kTRUE, mc);

    Float_t chi2MC = 0;
    Int_t chi2MCLimit = 0;
    mult->GetComparisonResults(&chi2MC, &chi2MCLimit, 0);
    fitResultsChi2->SetPoint(fitResultsChi2->GetN(), i, chi2MC);
    fitResultsChi2Limit->SetPoint(fitResultsChi2Limit->GetN(), i, chi2MCLimit);
    min = TMath::Min(min, chi2MC);
    max = TMath::Max(max, chi2MC);

    TH1* chi2Result = (TH1*) mult->GetMultiplicityESDCorrected(histID)->Clone("chi2Result");
    if (!firstChi)
      firstChi = (TH1*) chi2Result->Clone("firstChi");

    mult->ApplyBayesianMethod(histID, kFALSE, AliMultiplicityCorrection::kTrVtx, 0.1, 30, startCond);
    mult->DrawComparison(Form("StartingConditions_%d_Bayesian", i), histID, kFALSE, kTRUE, mc);
    TH1* bayesResult = (TH1*) mult->GetMultiplicityESDCorrected(histID)->Clone("bayesResult");
    if (!firstBayesian)
      firstBayesian = (TH1*) bayesResult->Clone("firstBayesian");

    mult->GetComparisonResults(&chi2MC, &chi2MCLimit, 0);
    fitResultsBayes->SetPoint(fitResultsBayes->GetN(), i, chi2MC);
    fitResultsBayesLimit->SetPoint(fitResultsBayesLimit->GetN(), i, chi2MCLimit);

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
  gSystem->Load("libPWG0base");

  AliMultiplicityCorrection** data = new AliMultiplicityCorrection*[n];
  TList list;
  for (Int_t i=0; i<n; ++i)
  {
    TString name("Multiplicity");
    if (i > 0)
      name.Form("Multiplicity%d", i);

    TFile::Open(files[i]);
    data[i] = new AliMultiplicityCorrection(name, name);
    data[i]->LoadHistograms("Multiplicity");
    if (i > 0)
      list.Add(data[i]);
  }

  data[0]->Merge(&list);

  data[0]->DrawHistograms();

  TFile::Open(output, "RECREATE");
  data[0]->SaveHistograms();
  gFile->Close();
}

void testMethod(Int_t caseNo, const char* fileName = "multiplicityMC.root")
{
  gSystem->Load("libPWG0base");

  AliMultiplicityCorrection* mult = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");

  TFile::Open(fileName);
  mult->LoadHistograms("Multiplicity");

  TF1* func = 0;

  if (caseNo >= 4)
  {
    func = new TF1("nbd", "[0] * TMath::Binomial([2]+TMath::Nint(x)-1, [2]-1) * pow([1] / ([1]+[2]), TMath::Nint(x)) * pow(1 + [1]/[2], -[2])", 0, 50);
    func->SetParNames("scaling", "averagen", "k");
  }

  switch (caseNo)
  {
    case 0: func = new TF1("flat", "1"); break;
    case 1: func = new TF1("flat", "501-x"); break;
    case 2: func = new TF1("flat", "1000 * 1/(x+1)"); break;
    case 3: func = new TF1("flat", "1000 * TMath::Landau(x, 10, 5)"); break;
    case 4: func->SetParameters(1e7, 10, 2); break;
    case 5: func->SetParameters(1e7, 20, 3); break;
    case 6: func->SetParameters(1e7, 30, 4); break;
    case 7: func->SetParameters(1e7, 70, 2); break;
    case 8: func = new TF1("testlaszlo", "10*1000*x*exp(-0.1*x)"); break;

    default: return;
  }

  mult->SetGenMeasFromFunc(func, 2);

  TFile::Open("out.root", "RECREATE");
  mult->SaveHistograms();

  //mult->ApplyBayesianMethod(2, kFALSE);
  //mult->ApplyMinuitFit(2, kFALSE);
  //mult->ApplyGaussianMethod(2, kFALSE);
  mult->ApplyLaszloMethod(2, kFALSE, AliMultiplicityCorrection::kTrVtx);
}

void smoothCorrelationMap(const char* fileName = "multiplicityMC.root", Int_t corrMatrix = 2)
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
  proj->Draw("COLZ");

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

  for (Int_t i=1; i<=150; ++i)
  {
    printf("Fitting %d...\n", i);

    TH1* hist = proj->ProjectionY(Form("proj%d", i), i, i, "e");
    //hist->GetXaxis()->SetRangeUser(0, 50);
    //lognormal->SetParameter(0, hist->GetMaximum());
    hist->Fit(fitWith, "0 M", "");

    TF1* func = hist->GetListOfFunctions()->FindObject(fitWith);

    if (((i-1) % 15 == 0) || ((i % 5 == 0) && i < 30))
    {
      new TCanvas;
      hist->Draw();
      func->Clone()->Draw("SAME");
      gPad->SetLogy();
    }

    scaling->Fill(i, func->GetParameter(0));
    mean->Fill(i, func->GetParameter(1));
    width->Fill(i, func->GetParameter(2));
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
  scaling->Fit(scalingFit, "", "", 3, 100);

  c1->cd(2);
  mean->Draw("P");

  //TF1* meanFit = log;
  TF1* meanFit = new TF1("mypol1", "[0]+[1]*x");
  mean->Fit(meanFit, "", "", 3, 100);

  c1->cd(3);
  width->Draw("P");

  //TF1* widthFit = over;
  TF1* widthFit = new TF1("mypol2", "[0]+[1]*x+[2]*x*x");
  width->Fit(widthFit, "", "", 5, 100);

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
      if (i < 21)
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
      corr->SetBinContent(corr->GetXaxis()->GetNbins() / 2, i, j, new->GetBinContent(i, j));

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
