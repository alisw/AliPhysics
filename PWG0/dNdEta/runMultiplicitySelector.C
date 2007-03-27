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

  gProof->SetParallel(0);
  gProof->SetLogLevel(4);
  gProof->SetParallel(9999);

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

  mult->ApplyBayesianMethod(hist, kFALSE, eventType);
  mult->DrawComparison("Bayesian", hist, kFALSE, kTRUE, mult->GetMultiplicityMC(hist, eventType)->ProjectionY());

  //mult->ApplyMinuitFit(hist, kFALSE);
  //mult->DrawComparison("MinuitChi2", hist, kFALSE, kTRUE, mult->GetMultiplicityMC(hist, AliMultiplicityCorrection::kTrVtx)->ProjectionY());

  //mult->ApplyGaussianMethod(hist, kFALSE);

  return mult;
}

void* fitOther(const char* fileNameMC = "multiplicityMC.root", const char* fileNameESD = "multiplicityESD.root", Int_t histID = 2)
{
  gSystem->Load("libPWG0base");

  AliMultiplicityCorrection* mult = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");

  TFile::Open(fileNameMC);
  mult->LoadHistograms("Multiplicity");

  TFile::Open(fileNameESD);
  TH2F* hist = (TH2F*) gFile->Get(Form("Multiplicity/fMultiplicityESD%d", histID));
  TH2F* hist2 = (TH2F*) gFile->Get(Form("Multiplicity/fMultiplicityVtx%d", histID));

  mult->SetMultiplicityESD(histID, hist);

  //mult->ApplyGaussianMethod(histID, kFALSE);
  //for (Float_t f=0; f<0.1; f+=0.01)
  //mult->ApplyBayesianMethod(histID, kFALSE, AliMultiplicityCorrection::kTrVtx);
  mult->ApplyMinuitFit(histID, kFALSE);
  mult->DrawComparison("MinuitChi2", hist, kFALSE, kTRUE, hist2->ProjectionY());

  //mult->ApplyLaszloMethod(histID, kFALSE, AliMultiplicityCorrection::kTrVtx);

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
  }
  return 0;
}

void EvaluateChi2Method(const char* fileNameMC = "multiplicityMC.root", const char* fileNameESD = "multiplicityMC.root", Int_t histID = 2)
{
  gSystem->Load("libPWG0base");

  AliMultiplicityCorrection* mult = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");

  TFile::Open(fileNameMC);
  mult->LoadHistograms("Multiplicity");

  TFile::Open(fileNameESD);
  TH2F* hist = (TH2F*) gFile->Get(Form("Multiplicity/fMultiplicityESD%d", histID));
  TH2F* hist2 = (TH2F*) gFile->Get(Form("Multiplicity/fMultiplicityVtx%d", histID));

  mult->SetMultiplicityESD(histID, hist);

  TCanvas* canvas = new TCanvas("EvaluateChi2Method", "EvaluateChi2Method", 800, 600);
  TLegend* legend = new TLegend(0.6, 0.1, 0.98, 0.4);
  legend->SetFillColor(0);

  Float_t min = 1e10;
  Float_t max = 0;

  TGraph* first = 0;

  for (AliMultiplicityCorrection::RegularizationType type = AliMultiplicityCorrection::kPol0; type <= AliMultiplicityCorrection::kEntropy; ++type)
  {
    TGraph* fitResultsMC = new TGraph;
    fitResultsMC->SetTitle(";Weight Parameter");
    TGraph* fitResultsRes = new TGraph;
    fitResultsRes->SetTitle(";Weight Parameter");

    fitResultsMC->SetFillColor(0);
    fitResultsRes->SetFillColor(0);
    fitResultsMC->SetMarkerStyle(19+type);
    fitResultsRes->SetMarkerStyle(19+type);
    fitResultsRes->SetMarkerColor(kRed);
    fitResultsRes->SetLineColor(kRed);

    legend->AddEntry(fitResultsMC, Form("%s MC chi2", GetRegName(type)));
    legend->AddEntry(fitResultsRes, Form("%s residual chi2", GetRegName(type)));

    if (first == 0)
      first = fitResultsMC;

    for (Float_t weight = 1e-2; weight < 2e4; weight *= 1e2)
    {
      Float_t chi2MC = 0;
      Float_t residuals = 0;

      mult->SetRegularizationParameters(type, weight);
      mult->ApplyMinuitFit(histID, kFALSE);
      mult->DrawComparison(Form("MinuitChi2_%d_%f", type, weight), histID, kFALSE, kTRUE, hist2->ProjectionY());
      mult->GetComparisonResults(chi2MC, residuals);

      fitResultsMC->SetPoint(fitResultsMC->GetN(), weight, chi2MC);
      fitResultsRes->SetPoint(fitResultsRes->GetN(), weight, residuals);

      min = TMath::Min(min, TMath::Min(chi2MC, residuals));
      max = TMath::Max(max, TMath::Max(chi2MC, residuals));
    }

    fitResultsMC->Print();
    fitResultsRes->Print();

    canvas->cd();
    fitResultsMC->Draw(Form("%s CP", (type == AliMultiplicityCorrection::kPol0) ? "A" : "SAME"));
    fitResultsRes->Draw("SAME CP");
  }

  gPad->SetLogx();
  gPad->SetLogy();
  printf("min = %f, max = %f\n", min, max);
  if (min <= 0)
    min = 1e-5;
  first->GetYaxis()->SetRangeUser(min * 0.5, max * 1.5);

  legend->Draw();

  canvas->SaveAs(Form("%s.gif", canvas->GetName()));
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
    case 4: func->SetParameters(1000, 10, 2); break;
    case 5: func->SetParameters(1000, 20, 3); break;
    case 6: func->SetParameters(1000, 30, 4); break;
    case 7: func->SetParameters(1000, 40, 5); break;

    default: return;
  }

  mult->SetGenMeasFromFunc(func, 2);

  mult->ApplyBayesianMethod(2, kFALSE);
  mult->ApplyMinuitFit(2, kFALSE);
  //mult->ApplyGaussianMethod(2, kFALSE);
}
