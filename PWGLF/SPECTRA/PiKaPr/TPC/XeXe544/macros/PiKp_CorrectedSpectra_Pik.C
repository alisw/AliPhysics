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
Double_t TrackingPtGeantFlukaCorrectionKaMinus(Double_t pTmc= 0);
//void CalcProjPerCent(const Char_t * particleType , TH1D * efficiency); //for kaons

void PiKp_CorrectedSpectra_Pik()
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
  // calculate efficiencies
  //
  TH1D * effPartPion = GetEfficiency("Pion");
  TH1D * effPartAntiPion = GetEfficiency("AntiPion");
  TH1D * effPartKaon = GetEfficiency("Kaon");
  TH1D* effPartAntiKaon = GetEfficiency("AntiKaon");
   for (Int_t iBin = 0; iBin < effPartAntiKaon->GetXaxis()->GetNbins(); iBin++) {
    Double_t pTmc = effPartAntiKaon->GetXaxis()->GetBinCenter(iBin);
    Double_t yvalue = effPartAntiKaon->GetBinContent(iBin);
    Double_t corrFactor = TrackingPtGeantFlukaCorrectionKaMinus(pTmc);
    Double_t newYvalue = yvalue/corrFactor;
    effPartAntiKaon->SetBinContent(iBin,newYvalue);
         } 
  //TH1D* effPartProt = GetEfficiency("Prot");
  //TH1D* effPartAntiProt = GetEfficiency("AntiProt");
    //
    // Generate and correct Raw Spectra
    //
    CalcProjPerCent("Pion", effPartPion);
    CalcProjPerCent("AntiPion", effPartAntiPion);
    CalcProjPerCent("Kaon", effPartKaon);
    CalcProjPerCent("AntiKaon", effPartAntiKaon);
    //CalcProjPerCent("Proton", effPartProt);
    //CalcProjPerCent("AntiProton", effPartAntiProt);
}

Double_t TrackingPtGeantFlukaCorrectionKaMinus(Double_t pTmc)
{
  return TMath::Min((0.972865 + 0.0117093 * pTmc), 1.);
}

void CalcProjPerCent(const Char_t* particleType, TH1D* efficiency)
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
  TH3D* fPartCent = (TH3D*)list->FindObject(Form("f%sDeDxCent", particleType));
  fPartCent->Sumw2();

  const Int_t nCents = 9;
  Double_t centMin[nCents] = { 0.5, 5.5, 10.5, 20.5, 30.5, 40.5, 50.5, 60.5, 70.5 };
  Double_t centMax[nCents] = { 4.5, 9.5, 19.5, 29.5, 39.5, 49.5, 59.5, 69.5, 89.5 };
  const TString col[10] = { "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a" }; //http://colorbrewer2.org/#type=qualitative&scheme=Paired&n=10
  TH1D* histPart[nCents] = { 0x0 };
  TH1D* frac[nCents] = { 0x0 };
  TString strType(particleType);
  if (strType.EqualTo("AntiPion"))
    strType = "Pion";

  for (Int_t icent = 0; icent < nCents && !strType.Contains("Kaon"); icent++) {
    frac[icent] = GetFraction(strType, centMin[icent] - .5, centMax[icent] + .5);
  }

  for (Int_t icent = 0; icent < nCents; icent++) {
      //  histPart[icent] = GetSpectrumInCentralityBinFixCut(fPartCent,
    histPart[icent] = GetSpectrumInCentralityBinMultGaus(fPartCent,
        particleType,
        centMin[icent],
        centMax[icent]);
    NormalizeSpectrum(histPart[icent],
        1.0,
        fHistCent->Integral(
            fHistCent->GetXaxis()->FindBin(centMin[icent]),
            fHistCent->GetXaxis()->FindBin(centMax[icent])));
    histPart[icent]->Divide(efficiency);
    if (frac[icent])
      histPart[icent]->Multiply(frac[icent]);
    CutPtRange(histPart[icent], particleType);
    //histPart[icent]->Scale(1e9);
    histPart[icent]->SetLineColor(TColor::GetColor(col[icent]));
    histPart[icent]->SetMarkerColor(TColor::GetColor(col[icent]));
    Cosmetics(histPart[icent]);
    //histPart[icent]->Print("all");
  }
  //
  //Plot the raw spectra
  //
  TCanvas* canvPtDist = new TCanvas(Form("canvPtDist%s", particleType), Form("canvPtDist%s", particleType));
  canvPtDist->SetTitle(Form("%sCentralityClasses", particleType));
  gPad->SetLogy();

  TLegend* legend = new TLegend(0.75, 0.6, 0.9, 0.9, "Xe-Xe,5.44TeV");
  legend->SetHeader(Form("centrality classes"));
  legend->SetTextSize(0.02);
  for (Int_t icent = 0; icent < nCents; icent++) {
    histPart[icent]->GetYaxis()->SetTitle("1/N_{ev} #frac{d^{2}N}{dp_{T}d#it{y}}(GeV^{-1}/c)");
    histPart[icent]->GetXaxis()->SetTitle("p_{T}(GeV/c)");
    histPart[icent]->GetYaxis()->SetRangeUser(1, 1000);
    histPart[icent]->GetXaxis()->SetRangeUser(0.1, 1.0);
    histPart[icent]->DrawCopy(icent == 0 ? "Ep" : "Ep,SAME");
    legend->AddEntry(histPart[icent], Form("%.0f-%.0f%%", centMin[icent] - .5, centMax[icent] + .5), "lp");
  }

  legend->Draw();
  canvPtDist->SaveAs(Form("FinalSpectra_forQM/%s.root", canvPtDist->GetName()));
  canvPtDist->SaveAs(Form("FinalSpectra_forQM/%s.png", canvPtDist->GetName()));
  //canvPtDist->SaveAs(Form("Fix_Cut_Spectra/%s.root", canvPtDist->GetName()));
  //canvPtDist->SaveAs(Form("Fix_Cut_Spectra/%s.png", canvPtDist->GetName()));
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
  TH1D* hist = GetRawSpectrum(pTvsNsigma, particleType, binlow, binUp);

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
      if (pT < 0.2)
        funcGausFit->SetParLimits(5, -10.0, -2.0);
    } // Pion

    if (name.Contains("Kaon")) {
     funcGausFit->SetParameters(3638.901406, -0.231609, -1.333264, 2.337703, 53302.352113, -6.889051, 1.141718, 2.522051); // Kaon
     if ((binlow >= 5) && (binUp <10)) funcGausFit->SetParameters(3177.248170,-0.056969,-1.144069,1.342167,39845.347373,-5.941882,1.021963,1.751);
     if ((binlow >= 10) && (binUp <20)) funcGausFit->SetParameters(4885.054574,-0.090713,-1.119390,1.326326,59672.801550,-6.026855,1.034954,1.908);
     if ((binlow >= 20) && (binUp <30)) funcGausFit->SetParameters(3391.500582,-0.077089,-1.108721,1.371476,40180.088961,-6.070771,1.048102,2.087);
     if ((binlow >= 30) && (binUp <40)) funcGausFit->SetParameters(2172.455562,-0.185090,-1.224220,2.519541,28726.819913,-7.242756,1.217320,4.654);
     if ((binlow >= 40) && (binUp <50)) funcGausFit->SetParameters(1438.102784,0.005925,-1.128071,1.515512,16395.436530,-6.043770,1.069208,4.4235);
     if ((binlow >= 50) && (binUp <60)) funcGausFit->SetParameters( 895.042847,0.058939,-1.055676,1.334239,9715.939687,-6.021401,1.061027,2.211260 );
     if ((binlow >= 60) && (binUp <70)) funcGausFit->SetParameters(474.762605,0.006047,-1.168819,2.608116,5810.827642,-7.203814,1.234378,3.629573);
     else if ((binlow >= 70) && (binUp <90)) funcGausFit->SetParameters(474.762605,0.006047,-1.168819,2.608116,5810.827642,-7.203814,1.234378,3.629573);
     if (pT < 0.2) funcGausFit->SetParLimits(5, -80.0, -40.0); 
               //  funcGausFit->SetRange(-7,8);
      }

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
      if (pT < 0.4) funcGausFit->SetParLimits(5, -80.0, -40.0); 
    } 

    if (name.Contains("_Prot")) {
      funcGausFit->SetParameters(1000, -0.223973, -1.248492, 2.014887, 23738.374711, -9.443739, 2.016047, 3.0);  // Proton
      //if ((binlow >= 20) && (binUp <30)) funcGausFit->SetParameters(958.159,-0.107,-1.067,1.814,13185.099,-9.920,1.423,2.573);
      if ((binlow >= 30) && (binUp <40)) funcGausFit->SetParameters(644.175855,-0.056694,-1.052726,1.980896,8728.640906,-9.761502,1.688535,3.03758);
      if ((binlow >= 40) && (binUp <50)) {
        funcGausFit->SetParameters( 543.132403,-0.107918,2.014887,2.107382,6985.531383,-9.441000,2.011205);
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
    Printf("Parameters: %4.3f, %4.3f,%4.3f, %4.3f,%4.3f, %4.3f,%4.3f, %4.3f,%4.3f", pT, funcGausFit->GetParameter(0), funcGausFit->GetParameter(1), funcGausFit->GetParameter(2), funcGausFit->GetParameter(3), funcGausFit->GetParameter(4), funcGausFit->GetParameter(5), funcGausFit->GetParameter(6), funcGausFit->GetParameter(7));
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
