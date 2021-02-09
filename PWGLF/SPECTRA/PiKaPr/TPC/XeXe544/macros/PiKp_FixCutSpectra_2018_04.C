//
// raw spectranormalized with (1 / Nev)*  dN/(dy*dpt)
// applied efficiency correction
// positive and negative particles both
// pT range cut applied
// multigaussian fit applied
// feed down correction applied
//
#include "TArrayD.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TH2D.h"
#include "TH3.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TStyle.h"

void Cosmetics(TH1D* hist);
void CutPtRange(TH1D* inputHist, const Char_t* particleType = "Pion");
TH1D* GetRawSpectrum(TH2D* pTvsNsigma, const Char_t* particleType = "Pion",Float_t binlow = 0, Float_t binUp = 5);
TH1D* GetSpectrumInCentralityBinMultGaus(TH3D* histPartCent, const Char_t* particleType = "Pion", Float_t binlow = 0.5, Float_t binUp = 4.5);
TH1D* GetSpectrumInCentralityBinFixCut(TH3D* histPartCent, const Char_t* particleType = "Pion", Float_t binlow = 0.5, Float_t binUp = 4.5);
void NormalizeSpectrum(TH1D* spectrum, Float_t dy, Float_t numberOfEvents);
TH1D* GetFraction(const Char_t* particleType = "Pion", Float_t binlow = 0, Float_t binUp = 5);
TH1D* GetEfficiency(const Char_t* particleType = "Muon");
void CalcProjPerCent(const Char_t* particleType, TH1D* efficiency);
//void CalcProjPerCent(const Char_t * particleType , TH1D * efficiency); //for kaons

void PiKp_FixCutSpectra_2018_04()
{

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  //
  TFile* inFile = TFile::Open("./Output_2018_03_17data/MergeData.root"); // open analysis results folder
  TList* list = (TList*)inFile->Get("chist");                            // make my output container to store results
  //
  // QA plot for centrality
  //
  TH1F* fHistCent = (TH1F*)list->FindObject("fHistCent");
  TCanvas* canvCentrality = new TCanvas("canvCentrality", "Centrality");
  fHistCent->DrawCopy("EP");
  inFile->Close();
  //
  // smooth primary fraction with a fit
  //
  // calculate efficiencies
  //
  TH1D * effPartPion = GetEfficiency("Pion");
  TH1D * effPartAntiPion = GetEfficiency("AntiPion");
  TH1D * effPartKaon = GetEfficiency("Kaon");
  TH1D * effPartAntiKaon = GetEfficiency("AntiKaon");
  TH1D* effPartProt = GetEfficiency("Prot");
  TH1D* effPartAntiProt = GetEfficiency("AntiProt");
  //
  // Generate and correct Raw Spectra
  //
  CalcProjPerCent("Pion", effPartPion);
  CalcProjPerCent("AntiPion", effPartAntiPion);
  CalcProjPerCent("Kaon", effPartKaon);
  CalcProjPerCent("AntiKaon", effPartAntiKaon);
  CalcProjPerCent("Proton", effPartProt);
  CalcProjPerCent("AntiProton", effPartAntiProt);
}
//******************************************************************************************************************

void CalcProjPerCent(const Char_t * particleType, TH1D * efficiency)
{

  TFile* inFile = TFile::Open("./Output_2018_03_17data/MergeData.root"); // open analysis results folder
  TList* list = (TList*)inFile->Get("chist");                            // make my output container to store results
  //
  // Plot centrality
  //
  TH1F* fHistCent = (TH1F*)list->FindObject("fHistCent");
  //
  // Plot 3D raw spectra histograms
  //
  //TCanvas* canv3D = new TCanvas(Form("canv%s3D", particleType),Form("canv%s3D", particleType));
  TH3D* fPartCent = (TH3D*)list->FindObject(Form("f%sDeDxCent", particleType));
  //fPartCent->GetXaxis()->SetTitle("pT(GeV/c)");
  //fPartCent->DrawCopy();
  fPartCent->Sumw2();

    TH1D* frac_0to5 = GetFraction(particleType, 0, 5);
    TH1D* frac_5to10 = GetFraction(particleType, 5, 10);
    TH1D* frac_10to20 = GetFraction(particleType, 10, 20);
    TH1D* frac_20to30 = GetFraction(particleType, 20, 30);
    TH1D* frac_30to40 = GetFraction(particleType, 30, 40);
    TH1D* frac_40to50 = GetFraction(particleType, 40, 50);
    TH1D* frac_50to60 = GetFraction(particleType, 50, 60);
    TH1D* frac_60to70 = GetFraction(particleType, 60, 70);
    TH1D* frac_70to90 = GetFraction(particleType, 70, 90);
    

    TString strType(particleType);
    if (strType.EqualTo("AntiPion")) {
      frac_0to5 = GetFraction("Pion", 0, 5);
      frac_5to10 = GetFraction("Pion", 5, 10);
      frac_10to20 = GetFraction("Pion", 10, 20);
      frac_20to30 = GetFraction("Pion", 20, 30);
      frac_30to40 = GetFraction("Pion", 30, 40);
      frac_40to50 = GetFraction("Pion", 40, 50);
      frac_50to60 = GetFraction("Pion", 50, 60);
      frac_60to70 = GetFraction("Pion", 60, 70);
      frac_70to90 = GetFraction("Pion", 70, 90);
    }
        if (strType.Contains("Kaon")) {
      frac_0to5 = 0x0;
      frac_5to10 = 0x0;
      frac_10to20 = 0x0;
      frac_20to30 = 0x0;
      frac_30to40 = 0x0;
      frac_40to50 = 0x0;
      frac_50to60 = 0x0;
      frac_60to70 = 0x0;
      frac_70to90 = 0x0;
    }

    TH1D * hist0to5Part = GetSpectrumInCentralityBinFixCut(fPartCent, particleType, 0.5, 4.5);
//TH1D* hist0to5Part = GetSpectrumInCentralityBinMultGaus(fPartCent, particleType, 0.5, 4.5);
    NormalizeSpectrum(hist0to5Part, 1.0, fHistCent->Integral(fHistCent->GetXaxis()->FindBin(0.5), fHistCent->GetXaxis()->FindBin(4.5)));
    hist0to5Part->Divide(efficiency);
    hist0to5Part->Multiply(frac_0to5);
    CutPtRange(hist0to5Part, particleType);
    //hist0to5Part->Scale(1e9);
    hist0to5Part->SetLineColor(2);
    //hist0to5Part->Print("all");
    //

    TH1D * hist5to10Part = GetSpectrumInCentralityBinFixCut(fPartCent, particleType, 5.5, 9.5);
//TH1D* hist5to10Part = GetSpectrumInCentralityBinMultGaus(fPartCent, particleType, 5.5, 9.5);
    NormalizeSpectrum(hist5to10Part, 1.0, fHistCent->Integral(fHistCent->GetXaxis()->FindBin(5.5), fHistCent->GetXaxis()->FindBin(9.5)));
    hist5to10Part->Divide(efficiency);
    hist5to10Part->Multiply(frac_5to10);
    CutPtRange(hist5to10Part, particleType);
    //hist5to10Part->Scale(1e8);
    hist5to10Part->SetLineColor(3);

    //
    TH1D * hist10to20Part = GetSpectrumInCentralityBinFixCut(fPartCent, particleType, 10.5, 19.5);
//TH1D* hist10to20Part = GetSpectrumInCentralityBinMultGaus(fPartCent, particleType, 10.5, 19.5);
    NormalizeSpectrum(hist10to20Part, 1.0, fHistCent->Integral(fHistCent->GetXaxis()->FindBin(10.5), fHistCent->GetXaxis()->FindBin(19.5)));
    hist10to20Part->Divide(efficiency);
    hist10to20Part->Multiply(frac_10to20);
    CutPtRange(hist10to20Part, particleType);
    //hist10to20Part->Scale(1e7);
    hist10to20Part->SetLineColor(4);

    //
    TH1D * hist20to30Part = GetSpectrumInCentralityBinFixCut(fPartCent, particleType, 20.5, 29.5);
//TH1D* hist20to30Part = GetSpectrumInCentralityBinMultGaus(fPartCent, particleType, 20.5, 29.5);
    NormalizeSpectrum(hist20to30Part, 1.0, fHistCent->Integral(fHistCent->GetXaxis()->FindBin(20.5), fHistCent->GetXaxis()->FindBin(29.5)));
    hist20to30Part->Divide(efficiency);
    hist20to30Part->Multiply(frac_20to30);
    CutPtRange(hist20to30Part, particleType);
    //hist20to30Part->Scale(1e6);
    hist20to30Part->SetLineColor(5);

    //
    TH1D * hist30to40Part = GetSpectrumInCentralityBinFixCut(fPartCent, particleType,  30.5, 39.5);
//TH1D* hist30to40Part = GetSpectrumInCentralityBinMultGaus(fPartCent, particleType, 30.5, 39.5);
    NormalizeSpectrum(hist30to40Part, 1.0, fHistCent->Integral(fHistCent->GetXaxis()->FindBin(30.5), fHistCent->GetXaxis()->FindBin(39.5)));
    hist30to40Part->Divide(efficiency);
    hist30to40Part->Multiply(frac_30to40);
    CutPtRange(hist30to40Part, particleType);
    //hist30to40Part->Scale(1e5);
    hist30to40Part->SetLineColor(6);
    //
    TH1D * hist40to50Part = GetSpectrumInCentralityBinFixCut(fPartCent, particleType, 40.5, 49.5);
//TH1D* hist40to50Part = GetSpectrumInCentralityBinMultGaus(fPartCent, particleType, 40.5, 49.5);
    NormalizeSpectrum(hist40to50Part, 1.0, fHistCent->Integral(fHistCent->GetXaxis()->FindBin(40.5), fHistCent->GetXaxis()->FindBin(49.5)));
    hist40to50Part->Divide(efficiency);
    hist40to50Part->Multiply(frac_40to50);
    CutPtRange(hist40to50Part, particleType);
    //hist40to50Part->Scale(1e4);
    hist40to50Part->SetLineColor(7);

    //
    TH1D * hist50to60Part = GetSpectrumInCentralityBinFixCut(fPartCent, particleType, 50.5, 59.5);
//TH1D* hist50to60Part = GetSpectrumInCentralityBinMultGaus(fPartCent, particleType, 50.5, 59.5);
    NormalizeSpectrum(hist50to60Part, 1.0, fHistCent->Integral(fHistCent->GetXaxis()->FindBin(50.5), fHistCent->GetXaxis()->FindBin(59.5)));
    hist50to60Part->Divide(efficiency);
    hist50to60Part->Multiply(frac_50to60);
    CutPtRange(hist50to60Part, particleType);
    //hist50to60Part->Scale(1e3);
    hist50to60Part->SetLineColor(8);

    //
    TH1D * hist60to70Part = GetSpectrumInCentralityBinFixCut(fPartCent, particleType, 60.5, 69.5);
//TH1D* hist60to70Part = GetSpectrumInCentralityBinMultGaus(fPartCent, particleType, 60.5, 69.5);
    NormalizeSpectrum(hist60to70Part, 1.0, fHistCent->Integral(fHistCent->GetXaxis()->FindBin(60.5), fHistCent->GetXaxis()->FindBin(69.5)));
    hist60to70Part->Divide(efficiency);
    hist60to70Part->Multiply(frac_60to70);
    CutPtRange(hist60to70Part, particleType);
    //hist60to70Part->Scale(1e2);
    hist60to70Part->SetLineColor(9);

    //
    TH1D * hist70to90Part = GetSpectrumInCentralityBinFixCut(fPartCent, particleType, 70.5, 89.5);
    //TH1D* hist70to90Part = GetSpectrumInCentralityBinMultGaus(fPartCent, particleType, 70.5, 89.5);
    NormalizeSpectrum(hist70to90Part, 1.0, fHistCent->Integral(fHistCent->GetXaxis()->FindBin(70.5), fHistCent->GetXaxis()->FindBin(89.5)));
    hist70to90Part->Divide(efficiency);
    hist70to90Part->Multiply(frac_70to90);
    CutPtRange(hist70to90Part, particleType);
    //hist70to90Part->Scale(10);
    hist70to90Part->SetLineColor(46);
    //
    //Plot the raw spectra
    //
    TCanvas* canvPtDist = new TCanvas(Form("canvPtDist%s", particleType), Form("canvPtDist%s", particleType));
    canvPtDist->SetTitle(Form("%sCentralityClasses", particleType));
    gPad->SetLogy();

    hist0to5Part->GetYaxis()->SetTitle("1/N_{ev} #frac{d^{2}N}{dp_{T}d#it{y}}(GeV^{-1}/c)");
    hist0to5Part->GetXaxis()->SetTitle("p_{T}(GeV/c)");
    hist0to5Part->GetYaxis()->SetRangeUser(1e-2, 1000);
    hist0to5Part->GetXaxis()->SetRangeUser(0, 1.0);

    hist0to5Part->SetMarkerColor(2);
    hist5to10Part->SetMarkerColor(3);
    hist10to20Part->SetMarkerColor(4);
    hist20to30Part->SetMarkerColor(5);
    hist30to40Part->SetMarkerColor(6);
    hist40to50Part->SetMarkerColor(7);
    hist50to60Part->SetMarkerColor(8);
    hist60to70Part->SetMarkerColor(9);
    hist70to90Part->SetMarkerColor(46);

    Cosmetics(hist0to5Part);
    Cosmetics(hist5to10Part);
    Cosmetics(hist10to20Part);
    Cosmetics(hist20to30Part);
    Cosmetics(hist30to40Part);
    Cosmetics(hist40to50Part);
    Cosmetics(hist50to60Part);
    Cosmetics(hist60to70Part);
    Cosmetics(hist60to70Part);
    Cosmetics(hist70to90Part);

    hist0to5Part->DrawCopy("Ep");
    hist5to10Part->DrawCopy("Ep,SAME");
    hist10to20Part->DrawCopy("Ep,SAME");
    hist20to30Part->DrawCopy("Ep,SAME");
    hist30to40Part->DrawCopy("Ep,SAME");
    hist40to50Part->DrawCopy("Ep,SAME");
    hist50to60Part->DrawCopy("Ep,SAME");
    hist60to70Part->DrawCopy("Ep,SAME");
    hist70to90Part->DrawCopy("Ep,SAME");

    //Int_t Char=0;
    TLegend* legend = new TLegend(0.75, 0.6, 0.9, 0.9, "Xe-Xe,5.44TeV");
    legend->SetHeader(Form("p_{T} dist for %s in centrality classes", particleType));
     legend->SetTextSize(0.02);
    legend->AddEntry(hist0to5Part, "0-5%", "lp");
    legend->AddEntry(hist5to10Part, "5-10%", "lp");
    legend->AddEntry(hist10to20Part, "10-20%", "lp");
    legend->AddEntry(hist20to30Part, "20-30%", "lp");
    legend->AddEntry(hist30to40Part, "30-40%", "lp");
    legend->AddEntry(hist40to50Part, "40-50%", "lp");
    legend->AddEntry(hist50to60Part, "50-60%", "lp");
    legend->AddEntry(hist60to70Part, "60-70%", "lp");
    legend->AddEntry(hist70to90Part, "70-90%", "lp");

    legend->Draw();
    //canvPtDist->SaveAs(Form("FinalSpectra_forQM/%s.root", canvPtDist->GetName()));
    //canvPtDist->SaveAs(Form("FinalSpectra_forQM/%s.png", canvPtDist->GetName()));
    canvPtDist->SaveAs(Form("Fix_Cut_Spectra/%s.root", canvPtDist->GetName()));
    canvPtDist->SaveAs(Form("Fix_Cut_Spectra/%s.png", canvPtDist->GetName()));
}

TH1D* GetEfficiency(const Char_t* particleType)
{
  //
  // get efficiency for a certain particle type
  //
  TFile* inFileMC = TFile::Open("./Output_2018_03_17_MC/MergeOutput_MC.root"); // open analysis results folder
  TList* listMC = (TList*)inFileMC->Get("chist");
  //
  TH1F* histPartGen = (TH1F*)listMC->FindObject(Form("fHist%sGen", particleType));
  TH1F* histPartReco = (TH1F*)listMC->FindObject(Form("fHist%sReco", particleType));
  histPartGen->Sumw2();
  histPartReco->Sumw2();
  //
  TH1D* effPart = (TH1D*)histPartReco->Clone(Form("eff%s", particleType));
  effPart->GetXaxis()->SetTitle("p_{T}(GeV/c)");
  effPart->GetYaxis()->SetTitle("All Reco. tracks/All Gen.tracks");
  effPart->GetXaxis()->SetRangeUser(0, 10);
  effPart->Divide(histPartReco, histPartGen, 1.0, 1.0, "B");
  TCanvas* canvPartEff = new TCanvas(Form("canv%sEff", particleType), Form("canv%sEff", particleType));
  effPart->DrawCopy("EP");
  //canvPartEff->SaveAs(Form("%s.eps", effPart->GetName()));

  effPart->SetDirectory(0);
  delete listMC;
  inFileMC->Close();
  return effPart;
}

TH1D* GetFraction(const Char_t* particleType, Float_t binlow, Float_t binUp)
{
  TVirtualPad* currpad = gPad;
  Printf("%f %f", binlow, binUp);
  TFile newfile(Form("FeedDownPlots/%s_%.0f_%.0f_fraction.root", particleType, binlow, binUp), "READ");
  Printf("%s", newfile.GetName());
  newfile.ls();
  TCanvas* c = (TCanvas*)newfile.Get(Form("%s_%.0f_%.0f_fraction", particleType, binlow, binUp));
  Printf("%s", c->GetName());
  TList* l = c->GetListOfPrimitives();
  l->ls();
  TH1D* hFrac = (TH1D*)l->At(2)->Clone(Form("%s_%s", l->At(2)->GetName(), particleType));
  //TH1D * hFrac = (TH1D*)l->At(1)->Clone(Form("%s", particleType));
  hFrac->SetDirectory(0);
  newfile.Close();
  currpad->cd();
  return hFrac;
}

void NormalizeSpectrum(TH1D* spectrum, Float_t dy, Float_t numberOfEvents)
{
  //
  // make (1 / Nev)*  dN/(dy*dpt)
  //
  spectrum->Scale(1. / dy);
  spectrum->Scale(1. / numberOfEvents);
  //
  Int_t nBins = spectrum->GetNbinsX();
  for (Int_t i = 0; i < nBins + 1; i++) {
    Double_t Content = spectrum->GetBinContent(i);
    Double_t error = spectrum->GetBinError(i);
    if (spectrum->GetBinWidth(i) == 0)
      continue;
    spectrum->SetBinContent(i, Content / spectrum->GetBinWidth(i));
    spectrum->SetBinError(i, error / spectrum->GetBinWidth(i));
  }
}

TH1D* GetSpectrumInCentralityBinFixCut(TH3D* histPartCent, const Char_t* particleType, Float_t binlow, Float_t binUp)
{

  TH1D* hist = histPartCent->ProjectionX(Form("hist%ito%i_%s", TMath::Nint(binlow - 0.5), TMath::Nint(binUp + 0.5), particleType), histPartCent->GetYaxis()->FindBin(-3.5), histPartCent->GetYaxis()->FindBin(3.5), histPartCent->GetZaxis()->FindBin(binlow), histPartCent->GetZaxis()->FindBin(binUp));
  return hist;
}

TH1D* GetSpectrumInCentralityBinMultGaus(TH3D* histPartCent, const Char_t* particleType, Float_t binlow, Float_t binUp)
{

  histPartCent->GetZaxis()->SetRangeUser(binlow, binUp);
  TH2D* pTvsNsigma = (TH2D*)histPartCent->Project3D("yx");
  pTvsNsigma->SetNameTitle(Form("hist%ito%i_pTvsNsigma_%s", TMath::Nint(binlow - 0.5), TMath::Nint(binUp + 0.5), particleType), Form("hist%ito%i_pTvsNsigma_%s", TMath::Nint(binlow - 0.5), TMath::Nint(binUp + 0.5), particleType));
  TH1D* hist = GetRawSpectrum(pTvsNsigma,particleType, binlow, binUp);

  return hist;
}

TH1D* GetRawSpectrum(TH2D* pTvsNsigma, const Char_t* particleType, Float_t binlow, Float_t binUp)
{
  //
  // get the raw yield from pT vs Nisgma for a given centrality class
  //
  //TCanvas * canvQAmultGausFit = new TCanvas(Form("canvQAmultGausFit_%s",pTvsNsigma->GetName()),Form("canvQAmultGausFit_%s",pTvsNsigma->GetName()));
  //pTvsNsigma->Draw("colz");
  //canvQAmultGausFit->SaveAs(Form("canvQAmultGausFit_%s.root",pTvsNsigma->GetName()));
  TH1D* rawSpectrum = pTvsNsigma->ProjectionX(Form("rawSpectrum_%s", pTvsNsigma->GetName()));
  rawSpectrum->Reset();
  //
  TCanvas* canvQAmultGausFitSingleBin = new TCanvas(Form("canvQAmultGausFitSingleBin_%s", pTvsNsigma->GetName()), Form("canvQAmultGausFitSingleBin_%s", pTvsNsigma->GetName()));
  canvQAmultGausFitSingleBin->Print(Form("Gaus_Fit/canvQAmultGausFitSingleBin_%s.pdf(", pTvsNsigma->GetName()), "pdf");
  for (Int_t iBin = 0; iBin < pTvsNsigma->GetXaxis()->GetNbins(); iBin++) {
    Double_t pT = pTvsNsigma->GetXaxis()->GetBinCenter(iBin);
    if (pT < 0.1 || pT > 0.8)
      continue;

    TH1D* multGausHist = pTvsNsigma->ProjectionY(Form("multGaus_%s_%f", pTvsNsigma->GetName(), pT), iBin, iBin);

    canvQAmultGausFitSingleBin->cd();
    multGausHist->DrawCopy("Ep");

    TF1* funcGausFit = new TF1("funcGausFit", "(x <= ([3] + [1])) * TMath::Abs([0]) * TMath::Gaus(x, [1], [2]) + (x > ([3] + [1])) * [0] * TMath::Gaus([3] + [1], [1], [2]) * TMath::Exp(-([3]) * (x - [3] - [1]) / ([2] * [2]))    + (x <= ([7] + [5])) * TMath::Abs([4]) * TMath::Gaus(x, [5], [6]) + (x > ([7] + [5])) * [4] * TMath::Gaus([7] + [5], [5], [6]) * TMath::Exp(-([7]) * (x - [7] - [5]) / ([6] * [6]))  ", -8.2, 8.2); //0:norm, 1:mean, 2:sigma, 3:tail
    
    TString t[3] = { "Pion.", "Kaon.", "Prot." };
    TString name = multGausHist->GetName();
    if (name.Contains("Pion")) {
      funcGausFit->SetParameters(24924.038834, -0.172245, 1.139421, 1.394517, 1550.601596, 6.411475, 1.646372, 5); //Pion
      if (pT < 0.2) funcGausFit->SetParLimits(5, -10.0, -2.0); } // Pion

    if (name.Contains("Kaon")) {
     funcGausFit->SetParameters(3638.901406, -0.231609, -1.333264, 2.337703, 53302.352113, -6.889051, 1.141718, 2.522051); // Kaon
     if ((binlow >= 5) && (binUp <10)) funcGausFit->SetParameters(3177.248170,-0.056969,-1.144069,1.342167,39845.347373,-5.941882,1.021963,1.751);
     if ((binlow >= 10) && (binUp <20)) funcGausFit->SetParameters(4885.054574,-0.090713,-1.119390,1.326326,59672.801550,-6.026855,1.034954,1.908);
     if ((binlow >= 20) && (binUp <30)) funcGausFit->SetParameters(3391.500582,-0.077089,-1.108721,1.371476,40180.088961,-6.070771,1.048102,2.087);
     if ((binlow >= 30) && (binUp <40)) funcGausFit->SetParameters(2172.455562,-0.185090,-1.224220,2.519541,28726.819913,-7.242756,1.217320,4.654);
     if ((binlow >= 40) && (binUp <50)) funcGausFit->SetParameters(1438.102784,0.005925,-1.128071,1.515512,16395.436530,-6.043770,1.069208,4.4235);
     if ((binlow >= 50) && (binUp <60)) funcGausFit->SetParameters( 895.042847,0.058939,-1.055676,1.334239,9715.939687,-6.021401,1.061027,2.211260 );
     if ((binlow >= 60) && (binUp <70)) funcGausFit->SetParameters(474.762605,0.006047,-1.168819,2.608116,5810.827642,-7.203814,1.234378,3.629573);
     else if ((binlow >= 70) && (binUp <90)) funcGausFit->SetParameters(94.650060,0.089907,-1.006088,1.479073,877203.314232,-18.939989,0.903836,0.7324);
     if (pT < 0.2) funcGausFit->SetParLimits(5, -80.0, -40.0); 
      }// Kaon

    if (name.Contains("_AntiProt")) {
      funcGausFit->SetParameters(1000, -0.223973, -1.248492, 2.014887, 23738.374711, -9.443739, 2.016047, 3.0);  // Proton
      if ((binlow >= 30) && (binUp <40)) funcGausFit->SetParameters(644.175855,-0.056694,-1.052726,1.980896,8728.640906,-9.761502,1.688535,3.03758);
      if ((binlow >= 40) && (binUp <50)) {
        funcGausFit->SetParameters(543.132403,-0.107918,2.014887,2.107382,6985.531383,-9.441000,2.011205);
        funcGausFit->SetRange(-5,5);
      }
      if ((binlow >= 50) && (binUp <60)) funcGausFit->SetParameters(355.220795,-0.056542,1.095427,2.106012,4188.763432,-9.509362,2.031393,5.00000);
      if ((binlow >= 60) && (binUp <70)) funcGausFit->SetParameters(215.734719,0.045289,-0.962056,1.587409,25962.450783,-12.249531,1.643087,1.9079);
      else if ((binlow >= 70) && (binUp <90)) funcGausFit->SetParameters(110.200770,0.000506,-1.045747,1.917268,371.896561,-8.318614,0.659736,0.561846);
      //else if ((binlow >= 80) && (binUp <90)) funcGausFit->SetParameters(50.462281,0.005254,-1.077229,2.279959,484.952001,-9.322721,2.077086,5.000000);
      if (pT < 0.4) funcGausFit->SetParLimits(5, -80.0, -40.0); 
    } // Anti-Proton

    if (name.Contains("_Prot")) {
      funcGausFit->SetParameters(1000, -0.223973, -1.248492, 2.014887, 23738.374711, -9.443739, 2.016047, 3.0);  // Proton
      if ((binlow >= 30) && (binUp <40)) funcGausFit->SetParameters(644.175855,-0.056694,-1.052726,1.980896,8728.640906,-9.761502,1.688535,3.03758);
      if ((binlow >= 40) && (binUp <50)) {
        funcGausFit->SetParameters(543.132403,-0.107918,2.014887,2.107382,6985.531383,-9.441000,2.011205);
        //funcGausFit->SetRange(-5,5);
      }
      if ((binlow >= 50) && (binUp <60)) funcGausFit->SetParameters(355.220795,-0.056542,1.095427,2.106012,4188.763432,-9.509362,2.031393,5.00000);
      if ((binlow >= 60) && (binUp <70)) funcGausFit->SetParameters(215.734719,0.045289,-0.962056,1.587409,25962.450783,-12.249531,1.643087,1.9079);
      else if ((binlow >= 70) && (binUp <90)) funcGausFit->SetParameters(110.200770,0.000506,-1.045747,1.917268,371.896561,-8.318614,0.659736,0.561846);
      //else if ((binlow >= 80) && (binUp <90)) funcGausFit->SetParameters(50.462281,0.005254,-1.077229,2.279959,484.952001,-9.322721,2.077086,5.000000);
      if (pT < 0.4) funcGausFit->SetParLimits(5, -80.0, -40.0); 
    } // Proton


    funcGausFit->SetParLimits(1, -0.3, 0.3);
    funcGausFit->SetParLimits(3, 0, 3.0);
    funcGausFit->SetParLimits(7, 0, 5.0);
    multGausHist->Fit(funcGausFit, "QNRL");
    funcGausFit->Draw("SAME");
    Printf(Form("Parameters: %4.3f, %4.3f,%4.3f, %4.3f,%4.3f, %4.3f,%4.3f, %4.3f,%4.3f", pT, funcGausFit->GetParameter(0), funcGausFit->GetParameter(1), funcGausFit->GetParameter(2), funcGausFit->GetParameter(3), funcGausFit->GetParameter(4), funcGausFit->GetParameter(5), funcGausFit->GetParameter(6), funcGausFit->GetParameter(7)));
    //
    TF1* funcSignal = new TF1("funcSignal", "(x <= ([3] + [1])) * [0] * TMath::Gaus(x, [1], [2]) + (x > ([3] + [1])) * [0] * TMath::Gaus([3] + [1], [1], [2]) * TMath::Exp(-([3]) * (x - [3] - [1]) / ([2] * [2]))", -10.0, 10.0);
    funcSignal->SetParameters(funcGausFit->GetParameter(0), funcGausFit->GetParameter(1), funcGausFit->GetParameter(2), funcGausFit->GetParameter(3));
    for (int i = 0; i < 4; i++)
      funcSignal->SetParError(i, funcGausFit->GetParError(i));
    funcSignal->SetLineColor(kBlue);
    funcSignal->Draw("SAME");
    Double_t yieldSignal = funcSignal->Integral(-10.0, 10.0) / multGausHist->GetBinWidth(10); // random
   //Double_t yieldSignalError = funcSignal->IntegralError(-10.0, 10.0) / multGausHist->GetBinWidth(10);

    TF1* funcBackGround = new TF1("funcBackGround", "(x <= ([3] + [1])) * TMath::Abs([0]) * TMath::Gaus(x, [1], [2]) + (x > ([3] + [1])) * TMath::Abs([0]) * TMath::Gaus([3] + [1], [1], [2]) * TMath::Exp(-([3]) * (x - [3] - [1]) / ([2] * [2])) ", -10.0, 10.0);
    funcBackGround->SetParameters(funcGausFit->GetParameter(4), funcGausFit->GetParameter(5), funcGausFit->GetParameter(6), funcGausFit->GetParameter(7));
    for (int i = 0; i < 4; i++)
      funcBackGround->SetParError(i + 4, funcGausFit->GetParError(i + 4));
    funcBackGround->SetLineColor(kOrange);
    funcBackGround->Draw("SAME");

    TLatex* ptlabel = new TLatex(0.7, 0.8, Form("pT = %f", pT));
    ptlabel->SetNDC();
    ptlabel->SetTextSize(0.04);
    ptlabel->SetTextFont(42);
    ptlabel->SetLineWidth(2);
    ptlabel->Draw();

    TLatex* params = new TLatex(0.01, 0.7, Form("%4.3f,%4.3f,%4.3f,%4.3f,%4.3f,%4.3f,%4.3f,%4.3f", funcGausFit->GetParameter(0), funcGausFit->GetParameter(1), funcGausFit->GetParameter(2), funcGausFit->GetParameter(3), funcGausFit->GetParameter(4), funcGausFit->GetParameter(5), funcGausFit->GetParameter(6), funcGausFit->GetParameter(7)));
    params->SetNDC();
    params->SetTextSize(0.04);
    params->SetTextFont(42);
    params->SetLineWidth(2);
    params->Draw();

    //canvQAmultGausFitSingleBin->SaveAs(Form("canvQAmultGausFitSingleBin_%s_%f.png",pTvsNsigma->GetName(),pT));
    Double_t yieldBck = funcBackGround->Integral(-4.0, 4.0) / multGausHist->GetBinWidth(4); // random
    //Double_t yieldBckError = funcBackGround->IntegralError(-4.0, 4.0) / multGausHist->GetBinWidth(4);
    //
    Double_t yieldError=0;
    Double_t yield = multGausHist->IntegralAndError(multGausHist->GetXaxis()->FindBin(-4.0), multGausHist->GetXaxis()->FindBin(+4.0),yieldError) ;
    yield -= yieldBck;
    //yieldError= TMath::Sqrt(yieldError*yieldError+yieldBckError*yieldBckError);
    yieldError= TMath::Sqrt(yieldError + 2*yieldBck);
    //
    rawSpectrum->SetBinContent(iBin, yield);
    rawSpectrum->SetBinError(iBin, yieldError);
    canvQAmultGausFitSingleBin->SetLogy();
    canvQAmultGausFitSingleBin->Print(Form("Gaus_Fit/canvQAmultGausFitSingleBin_%s.pdf", pTvsNsigma->GetName()), "pdf");

    //canvQAmultGausFitSingleBin->SaveAs(Form("Gaus_Fit/canvQAmultGausFitSingleBin_%s_%f.pdf", pTvsNsigma->GetName(), pT), "pdf");
  }
  canvQAmultGausFitSingleBin->DrawFrame(0, 0, 1, 1, "EMPTY;EMPTY;EMPTY");
  canvQAmultGausFitSingleBin->Print(Form("Gaus_Fit/canvQAmultGausFitSingleBin_%s.pdf)", pTvsNsigma->GetName()), "pdf");
  //delete canvQAmultGausFitSingleBin;
  //
  return rawSpectrum;
  //return hist;
}



void CutPtRange(TH1D* inputHist, const Char_t* particleType)
{
  //
  // remove parts were PID does not work
  //
  TH1D* oldHist = (TH1D*)inputHist->Clone("oldHist");
  inputHist->Reset();
  TString particleString(particleType);
  Double_t pTmin = 0.0;
  Double_t pTmax = 10.0;
  if (particleString.Contains("Pion")) {
    pTmin = 0.1;
    pTmax = 0.75;
  }
  if (particleString.Contains("Kaon")) {
    pTmin = 0.2;
    pTmax = 0.55;
  }
  if (particleString.Contains("Proton")) {
    pTmin = 0.4;
    pTmax = 0.85;
  }
  //
  for (Int_t iBin = 0; iBin < inputHist->GetXaxis()->GetNbins(); iBin++) {
    Double_t pT = inputHist->GetXaxis()->GetBinCenter(iBin);
    if (pT > pTmax || pT < pTmin)
      continue;
    inputHist->SetBinContent(iBin, oldHist->GetBinContent(iBin));
    inputHist->SetBinError(iBin, oldHist->GetBinError(iBin));
  }
  delete oldHist;
}

void Cosmetics(TH1D* hist)
{

  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(0.4);
  // hist->GetXaxis()->SetRangeUser(0,1.5);
}
