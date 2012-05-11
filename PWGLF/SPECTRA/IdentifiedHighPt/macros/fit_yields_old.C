
#include <TFile.h>
#include <TH1.h>
#include <TProfile.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TGraphAsymmErrors.h>
#include <TList.h>
#include <TMath.h>
#include <TSystem.h>
#include <TLatex.h>
#include <TF1.h>
#include <TF2.h>

#include "AliHighPtDeDxData.h"

#include "my_tools.C"
#include "my_functions.C"

#include <iostream>
#include <fstream>
#include <string>

using namespace std;

/*
  Ideas to improve:
  =================

  Use the real mean p -> Might give some effect
  Effect: Push the <dE/dx> down (let'a see if it is still needed in the fits)

  Use the real <nCl> 
  Effect: Increase sigma for p<15. For p>15 seems to saturate.
  
  To use:
  =======
  root is enough
  .L libMyDeDxAnalysis.so 
  .L my_tools.C+
  .L my_functions.C+
  .L fit_yields.C+


  mkdir results
  mkdir results/fits
  mkdir results/eta04
  mkdir results/eta48
  mkdir results/etaneg
  mkdir results/etapos
  mkdir results/neg
  mkdir results/pos
  mkdir results/comparison
  mkdir results/kaon
  mkdir results/lambda


  *********************
       ALL ETAS
  *********************
  FitYields("data/lhc10h_aod_all.root", 3.0, 20.0, 0, kTRUE, 0, 1, kTRUE, "lhc10h_aod_all.root", 0, "results/fits") 

  FitYields("data/lhc10h_aod_all.root", 3.0, 20.0, -1, kTRUE, 0, 1, kTRUE, "lhc10h_aod_all_neg.root", 0, "results/neg") 

  FitYields("data/lhc10h_aod_all.root", 3.0, 20.0, +1, kTRUE, 0, 1, kTRUE, "lhc10h_aod_all_pos.root", 0, "results/pos") 

  Compare("fit_yields_results/lhc10h_aod_all.root", "fit_yields_results/lhc10h_aod_all_neg.root", "fit_yields_results/lhc10h_aod_all_pos.root", "q<0", "q>0", "results/comparison/neg_vs_pos")

  CompareYields("data/lhc10h_aod_all.root", 0, 3.0, 20.0, -1, +1, 0, 0, kFALSE, 0, 1, kFALSE)



  *********************
     ETA LOW & HIGH
  *********************

  FitYields("data/lhc10h_aod_all.root", 3.0, 20.0, 0, kTRUE, 0, 1, kTRUE, "lhc10h_aod_all_eta_low.root", "etaabs04", "results/eta04") 

  FitYields("data/lhc10h_aod_all.root", 3.0, 20.0, 0, kTRUE, 0, 1, kTRUE, "lhc10h_aod_all_eta_high.root", "etaabs48", "results/eta48") 

  Compare("fit_yields_results/lhc10h_aod_all.root", "fit_yields_results/lhc10h_aod_all_eta_low.root", "fit_yields_results/lhc10h_aod_all_eta_high.root", "|#eta|<0.4", "0.4<|#eta|<0.8", "results/comparison/low_vs_high")

  CompareYields("data/lhc10h_aod_all.root", 0, 3.0, 20.0, 0, 0, "etaabs04", "etaabs48")


  *********************
     ETA NEG & POS
  *********************
  FitYields("data/lhc10h_aod_all.root", 3.0, 20.0, 0, kTRUE, 0, 1, kTRUE, "lhc10h_aod_all_eta_neg.root", "eta-80", "results/etaneg") 

  FitYields("data/lhc10h_aod_all.root", 3.0, 20.0, 0, kTRUE, 0, 1, kTRUE, "lhc10h_aod_all_eta_pos.root", "eta08", "results/etapos") 

  Compare("fit_yields_results/lhc10h_aod_all.root", "fit_yields_results/lhc10h_aod_all_eta_neg.root", "fit_yields_results/lhc10h_aod_all_eta_pos.root", "#eta < 0", "#eta > 0", "results/comparison/etaneg_vs_etapos")

  CompareYields("data/lhc10h_aod_all.root", 0, 3.0, 20.0, 0, 0, "eta-80", "eta08")


  *********************
     V0 ETA LOW & HIGH
  *********************
  FitYieldsV0("data/lhc10h_aod_all.root", "data/lhc10h_v0_all_loose.root", "kaon", 3.0, 20.0, 0, kTRUE, 1, 0, kTRUE, "etaabs04") 
  FitYieldsV0("data/lhc10h_aod_all.root", "data/lhc10h_v0_all_loose.root", "kaon", 3.0, 20.0, 0, kTRUE, 1, 0, kTRUE, "etaabs48") 
  FitYieldsV0("data/lhc10h_aod_all.root", "data/lhc10h_v0_all_loose.root", "lambda", 3.0, 20.0, 0, kTRUE, 1, 0, kTRUE, "etaabs04") 
  FitYieldsV0("data/lhc10h_aod_all.root", "data/lhc10h_v0_all_loose.root", "lambda", 3.0, 20.0, 0, kTRUE, 1, 0, kTRUE, "etaabs48") 

  *********************

  *********************
     ETA FINE
  *********************
  FitYields("data/lhc10h_aod_all.root", 3.0, 20.0, 0, kTRUE, 0, 1, kTRUE, "lhc10h_aod_all_eta_86.root", "eta-8-6", "results/eta86") 
  FitYields("data/lhc10h_aod_all.root", 3.0, 20.0, 0, kTRUE, 0, 1, kTRUE, "lhc10h_aod_all_eta_64.root", "eta-6-4", "results/eta64") 
  FitYields("data/lhc10h_aod_all.root", 3.0, 20.0, 0, kTRUE, 0, 1, kTRUE, "lhc10h_aod_all_eta_42.root", "eta-4-2", "results/eta42") 
  FitYields("data/lhc10h_aod_all.root", 3.0, 20.0, 0, kTRUE, 0, 1, kTRUE, "lhc10h_aod_all_eta_20.root", "eta-20", "results/eta20") 
  FitYields("data/lhc10h_aod_all.root", 3.0, 20.0, 0, kTRUE, 0, 1, kTRUE, "lhc10h_aod_all_eta_02.root", "eta02", "results/eta02") 
  FitYields("data/lhc10h_aod_all.root", 3.0, 20.0, 0, kTRUE, 0, 1, kTRUE, "lhc10h_aod_all_eta_24.root", "eta24", "results/eta24") 
  FitYields("data/lhc10h_aod_all.root", 3.0, 20.0, 0, kTRUE, 0, 1, kTRUE, "lhc10h_aod_all_eta_46.root", "eta46", "results/eta46") 
  FitYields("data/lhc10h_aod_all.root", 3.0, 20.0, 0, kTRUE, 0, 1, kTRUE, "lhc10h_aod_all_eta_68.root", "eta68", "results/eta68") 

  Compare(3.1)
  Compare(5.1)
  Compare(7.1)
  Compare(9.1)

 FitYields("data/lhc10h_mc_aod_all.root", 3.0, 20.0, 0, kTRUE)   

  FitYields("data/lhc10h_aod_test.root", 3.0, 20.0, 0, kTRUE)

  FitYields("data/lhc10h_aod_all.root", 3.0, 20.0, 0, kTRUE)

  FitYieldsV0("data/lhc10h_aod_all.root", "data/lhc10h_v0_all_loose.root", "kaon", 3.0, 20.0, 0, kTRUE)
  FitYieldsV0("data/lhc10h_aod_all.root", "data/lhc10h_v0_all_loose.root", "lambda", 3.0, 20.0, 0, kTRUE)



  MakeRatios("../ptspectraPP/fit_yields_results/lhc10h_aod_all_neg.root", "../ptspectraPP/fit_yields_results/lhc10h_aod_all_pos.root", "../ptspectraPP/results/comparison/")

  MakeRatios("../ptspectraPbPb_60_80/fit_yields_results/lhc10h_aod_all_neg.root", "../ptspectraPbPb_60_80/fit_yields_results/lhc10h_aod_all_pos.root", "../ptspectraPbPb_60_80/results/comparison/")

  MakeRatios("../ptspectraPbPb_40_60/fit_yields_results/lhc10h_aod_all_neg.root", "../ptspectraPbPb_40_60/fit_yields_results/lhc10h_aod_all_pos.root", "../ptspectraPbPb_40_60/results/comparison/")

  MakeRatios("../ptspectraPbPb_20_40/fit_yields_results/lhc10h_aod_all_neg.root", "../ptspectraPbPb_20_40/fit_yields_results/lhc10h_aod_all_pos.root", "../ptspectraPbPb_20_40/results/comparison/")

  MakeRatios("../ptspectraPbPb_10_20/fit_yields_results/lhc10h_aod_all_neg.root", "../ptspectraPbPb_10_20/fit_yields_results/lhc10h_aod_all_pos.root", "../ptspectraPbPb_10_20/results/comparison/")

  MakeRatios("../ptspectraPbPb_0_5/fit_yields_results/lhc10h_aod_all_neg.root", "../ptspectraPbPb_0_5/fit_yields_results/lhc10h_aod_all_pos.root", "../ptspectraPbPb_0_5/results/comparison/")

  MakeRatios("fit_yields_results/lhc10h_aod_all_neg.root", "fit_yields_results/lhc10h_aod_all_pos.root", "results/comparison/")


 */

void MakeRatios(const Char_t* file1Name, const Char_t* file2Name, 
		Bool_t drawFractionRatios,
		Bool_t drawYieldRatios);


//___________________________________________________________________________________________
void MakeRatios(const Char_t* fileNeg, const Char_t* filePos, 
		const Char_t* outdirname)
{
  /*
    For yields we assume that file 1 is negative charge and file
    2 is positive charge.
   */
  
  TFile* file1 = FindFileFresh(fileNeg);
  if(!file1)
    return;
  
  TH1D* hPionRatio1   = (TH1D*)file1->Get("hPionYield");
  SetHistError(hPionRatio1, 0.0);
  // TH1D* hKaonRatio1   = (TH1D*)file1->Get("hKaonRatio");
  // TH1D* hProtonRatio1 = (TH1D*)file1->Get("hProtonRatio");


  TFile* file2 = FindFileFresh(filePos);
  if(!file2)
    return;
  
  TH1D* hPionRatio2   = (TH1D*)file2->Get("hPionYield");
  // TH1D* hKaonRatio2   = (TH1D*)file2->Get("hKaonRatio");
  // TH1D* hProtonRatio2 = (TH1D*)file2->Get("hProtonRatio");

  
  TH1D* hPionAssymetry = (TH1D*)hPionRatio1->Clone("hPionAssymetry");
  hPionAssymetry->Divide(hPionRatio2);
  TF1 f1("f1", "-1.0", 0.0, 50.0);
  hPionAssymetry->Add(&f1); // correct for muons and electrons
  CutHistogram(hPionAssymetry, 3.0, 20.0);
  hPionAssymetry->GetYaxis()->SetRangeUser(-0.25, 0.25);
  hPionAssymetry->SetTitle("Syst. error evaluation: Assymetry; p_{T} [GeV/c]; pion assymmetry: (neg-pos)/pos");
  
  
  hPionRatio1->Divide(hPionRatio2);
  hPionRatio1->GetYaxis()->SetRangeUser(0.75, 1.25);
  hPionRatio1->SetTitle("Syst. error evaluation: Ratio; p_{T} [GeV/c]; pion fraction ratio: neg/pos");


  TCanvas* cPionFractionRatio = new TCanvas("cPionFractionRatio", "pion fraction ratio", 400, 300);
  cPionFractionRatio->Clear();
  cPionFractionRatio->SetGridy();
  cPionFractionRatio->cd();
  hPionRatio1->Draw();
  if(strstr(fileNeg, "PP"))
    gROOT->ProcessLine(".x ../ptspectraPP/drawText.C");
  else if(strstr(fileNeg, "0_5"))
    gROOT->ProcessLine(".x ../ptspectraPbPb_0_5/drawText.C");
  else if(strstr(fileNeg, "5_10"))
    gROOT->ProcessLine(".x ../ptspectraPbPb_5_10/drawText.C");
  else if(strstr(fileNeg, "10_20"))
    gROOT->ProcessLine(".x ../ptspectraPbPb_10_20/drawText.C");
  else if(strstr(fileNeg, "20_40"))
    gROOT->ProcessLine(".x ../ptspectraPbPb_20_40/drawText.C");
  else if(strstr(fileNeg, "40_60"))
    gROOT->ProcessLine(".x ../ptspectraPbPb_40_60/drawText.C");
  else if(strstr(fileNeg, "60_80"))
    gROOT->ProcessLine(".x ../ptspectraPbPb_60_80/drawText.C");

  cPionFractionRatio->SaveAs(Form("%s/pion_frac_ratio_neg_over_pos.gif", outdirname));
  cPionFractionRatio->SaveAs(Form("%s/pion_frac_ratio_neg_over_pos.pdf", outdirname));


  // TCanvas* cPionAssymetry = new TCanvas("cPionAssymetry", "pion fraction assymetry", 400, 300);
  // cPionAssymetry->Clear();
  // cPionAssymetry->SetGridy();
  // cPionAssymetry->cd();
  // hPionAssymetry->Draw();
  // if(strstr(fileNeg, "PP"))
  //   gROOT->ProcessLine(".x ../ptspectraPP/drawText.C");
  // else if(strstr(fileNeg, "0_5"))
  //   gROOT->ProcessLine(".x ../ptspectraPbPb_0_5/drawText.C");
  // else if(strstr(fileNeg, "5_10"))
  //   gROOT->ProcessLine(".x ../ptspectraPbPb_5_10/drawText.C");
  // else if(strstr(fileNeg, "10_20"))
  //   gROOT->ProcessLine(".x ../ptspectraPbPb_10_20/drawText.C");
  // else if(strstr(fileNeg, "20_40"))
  //   gROOT->ProcessLine(".x ../ptspectraPbPb_20_40/drawText.C");
  // else if(strstr(fileNeg, "40_60"))
  //   gROOT->ProcessLine(".x ../ptspectraPbPb_40_60/drawText.C");
  // else if(strstr(fileNeg, "60_80"))
  //   gROOT->ProcessLine(".x ../ptspectraPbPb_60_80/drawText.C");
  // cPionAssymetry->SaveAs(Form("%s/pion_assymetry_between_neg_and_pos.gif", outdirname));
  // cPionAssymetry->SaveAs(Form("%s/pion_assymetry_between_neg_and_pos.pdf", outdirname));  


}

//____________________________________________________________________________
void FitYields(const Char_t* dataFileName,
	       Double_t ptStart, Double_t ptStop,
	       Int_t charge,
	       Bool_t performFit = kFALSE,
	       Int_t run    = 0,
	       Int_t filter = 1,
	       Bool_t usePhiCut = kTRUE,
	       const Char_t* outFileName=0,
	       const Char_t* endName=0,
	       const Char_t* dirName="debugfits")
{
  gStyle->SetOptStat(0);

  
  TFile* dataFile = FindFileFresh(dataFileName);
  if(!dataFile)
    return;
  AliHighPtDeDxData* data = (AliHighPtDeDxData*)GetObject(dataFile, filter, usePhiCut, run, "filter", endName);
  data->Print();

  gSystem->Exec(Form("mv %s/* old/", dirName));
  if(data->IsMc())
    gSystem->Exec("mv debugfitsmc/* olddebugfitsmc/");


  TH2D* hDeltaPiVsPt = data->GetHistDeltaPiVsPt(0, charge);
  hDeDxVsP = hDeltaPiVsPt; // for the 2d fit to pick up the right bin

  TH2D* hDeltaPiVsPtPiGen = data->GetHistDeltaPiVsPt(1, charge);
  TH2D* hDeltaPiVsPtKGen  = data->GetHistDeltaPiVsPt(2, charge);
  TH2D* hDeltaPiVsPtPGen  = data->GetHistDeltaPiVsPt(3, charge);

  TH2D* hDeltaPiVsPtPi = 0;
  TH2D* hDeltaPiVsPtK  = 0;
  TH2D* hDeltaPiVsPtP  = 0;

  if(data->IsMc()) {

    hDeltaPiVsPtPi = data->GetHistDeltaPiVsPtMc(1, charge);
    hDeltaPiVsPtK  = data->GetHistDeltaPiVsPtMc(2, charge);
    hDeltaPiVsPtP  = data->GetHistDeltaPiVsPtMc(3, charge);
  }



  TProfile* hPiGenProfile = hDeltaPiVsPtPiGen->ProfileX();
  hPiGenProfile->SetMarkerStyle(29);
  TProfile* hKGenProfile = hDeltaPiVsPtKGen->ProfileX();
  hKGenProfile->SetMarkerStyle(29);
  TProfile* hPGenProfile = hDeltaPiVsPtPGen->ProfileX();
  hPGenProfile->SetMarkerStyle(29);

  TCanvas* cDeltaPiVsPt = new TCanvas("cDeltaPiVsPt", "dE/dx vs p", 400, 300);
  cDeltaPiVsPt->Clear();
  cDeltaPiVsPt->cd();
  cDeltaPiVsPt->SetLogz();
  hDeltaPiVsPt->Draw("COLZ");
  hPiGenProfile->Draw("SAME P");
  hKGenProfile->Draw("SAME P");
  hPGenProfile->Draw("SAME P");
  gROOT->ProcessLine(".x drawText.C");
  cDeltaPiVsPt->SaveAs(Form("%s/deltapi_vs_pt.gif", dirName));
  cDeltaPiVsPt->SaveAs(Form("%s/deltapi_vs_pt.pdf", dirName));

  TCanvas* cDeltaPiVsPtLogX = new TCanvas("cDeltaPiVsPtLogX", "dE/dx vs p", 400, 300);
  cDeltaPiVsPtLogX->Clear();
  cDeltaPiVsPtLogX->cd();
  cDeltaPiVsPtLogX->SetLogz();
  cDeltaPiVsPtLogX->SetLogx();
  hDeltaPiVsPt->Draw("COLZ");
  hPiGenProfile->Draw("SAME P");
  hKGenProfile->Draw("SAME P");
  hPGenProfile->Draw("SAME P");
  gROOT->ProcessLine(".x drawText.C");
  cDeltaPiVsPtLogX->SaveAs(Form("%s/deltapi_vs_pt_logx.gif", dirName));
  cDeltaPiVsPtLogX->SaveAs(Form("%s/deltapi_vs_pt_logx.pdf", dirName));

  // Root is a bit stupid with finidng bins so we have to add and subtract a
  // little to be sure we get the right bin as we typically put edges as
  // limits
  const Int_t binStart = hDeltaPiVsPt->GetXaxis()->FindBin(ptStart+0.01);
  ptStart = hDeltaPiVsPt->GetXaxis()->GetBinLowEdge(binStart);
  const Int_t binStop  = hDeltaPiVsPt->GetXaxis()->FindBin(ptStop-0.01);
  ptStop = hDeltaPiVsPt->GetXaxis()->GetBinUpEdge(binStop);
  //  const Int_t nBins    = binStop - binStart + 1;

  cout << "Doing fits from pTlow = " << ptStart << " (bin: " << binStart
       << ") to pThigh = " << ptStop << " (bin: " << binStop << ")" << endl;
  

  //cross check
  TCanvas* cFits = new TCanvas("cFits", "Fit comparison to data", 1200, 800);
  cFits->Clear();
  cFits->Divide(7, 4);

  // TF1* pionGen = new TF1("pionGen", "[0]/sqrt(6.2832*[2]*[2])*exp(-(x-[1])*(x-[1])/2.0/[2]/[2])", -30, 20);
  // pionGen->SetLineWidth(2);
  // pionGen->SetLineColor(kRed);

  TF1* pion = new TF1("pion", "gausn", -30, 20);
  pion->SetLineWidth(2);
  pion->SetLineColor(kRed);
  TF1* kaon = new TF1("kaon", "gausn", -30, 20);
  kaon->SetLineWidth(2);
  kaon->SetLineColor(kGreen);
  TF1* proton = new TF1("proton", "gausn", -30, 20);
  proton->SetLineWidth(2);
  proton->SetLineColor(kBlue);
  TF1* total = new TF1("total", "gausn(0)+gausn(3)+gausn(6)", -30, 20);
  total->SetLineColor(kBlack);
  total->SetLineWidth(2);
  total->SetLineStyle(2);

  TLegend* legend = new TLegend(0.11, 0.6, 0.35, 0.85);    
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->AddEntry(total, "3-Gauss fit", "L");
  legend->AddEntry(pion, "#pi", "L");
  legend->AddEntry(kaon, "K", "L");
  legend->AddEntry(proton, "p", "L");

  TCanvas* cSingleFit = new TCanvas("cSingleFit", "single fit");
  //  cSingleFit->SetLogy();

  TH1D* hPionRatio =(TH1D*)hDeltaPiVsPt->ProjectionX("hPionRatio", 1, 1);
  hPionRatio->SetTitle("particle fractions; p_{T} [GeV/c]; particle fraction");
  hPionRatio->Reset();
  TH1D* hKaonRatio   = (TH1D*)hPionRatio->Clone("hKaonRatio");
  TH1D* hProtonRatio = (TH1D*)hPionRatio->Clone("hProtonRatio");
  TH1D* hPionRatioMc = (TH1D*)hPionRatio->Clone("hPionRatioMc");
  TH1D* hKaonRatioMc = (TH1D*)hPionRatio->Clone("hKaonRatioMc");
  TH1D* hProtonRatioMc = (TH1D*)hPionRatio->Clone("hProtonRatioMc");

  TH1D* hPionYield =(TH1D*)hDeltaPiVsPt->ProjectionX("hPionYield", 1, 1);
  hPionYield->SetTitle("particle fractions; p_{T} [GeV/c]; particle fraction");
  hPionYield->Reset();
  TH1D* hKaonYield   = (TH1D*)hPionYield->Clone("hKaonYield");
  TH1D* hProtonYield = (TH1D*)hPionYield->Clone("hProtonYield");
  TH1D* hPionYieldMc = (TH1D*)hPionYield->Clone("hPionYieldMc");
  TH1D* hKaonYieldMc = (TH1D*)hPionYield->Clone("hKaonYieldMc");
  TH1D* hProtonYieldMc = (TH1D*)hPionYield->Clone("hProtonYieldMc");
  
  for(Int_t bin = binStart; bin <= binStop; bin++){
    
    cout << "Making projection for bin: " << bin << endl;
    
    const Int_t j = bin-binStart;
    
    TH1D* hDeltaPiVsPtProj =(TH1D*)hDeltaPiVsPt->ProjectionY(Form("hDeltaPiVsPtProj%d", bin), bin, bin);
    hDeltaPiVsPtProj->GetXaxis()->SetRangeUser(-25, 20);
    hDeltaPiVsPtProj->SetTitle(Form("%.1f<p_{T}<%.1f GeV/c", 
  				hDeltaPiVsPt->GetXaxis()->GetBinLowEdge(bin),
  				hDeltaPiVsPt->GetXaxis()->GetBinUpEdge(bin)));

    Double_t all =  hDeltaPiVsPtProj->GetEntries();

    TH1D* hDeltaPiVsPtPiGenProj =(TH1D*)hDeltaPiVsPtPiGen->ProjectionY(Form("hDeltaPiVsPtPiGenProj%d", bin), bin, bin);
    TH1D* hDeltaPiVsPtKGenProj =(TH1D*)hDeltaPiVsPtKGen->ProjectionY(Form("hDeltaPiVsPtKGenProj%d", bin), bin, bin);
    TH1D* hDeltaPiVsPtPGenProj =(TH1D*)hDeltaPiVsPtPGen->ProjectionY(Form("hDeltaPiVsPtPGenProj%d", bin), bin, bin);
    
    Double_t gausParams[9] = { 
      0.6*all,
      hDeltaPiVsPtPiGenProj->GetMean(), 
      hDeltaPiVsPtPiGenProj->GetRMS(), 
      0.2*all,
      hDeltaPiVsPtKGenProj->GetMean(), 
      hDeltaPiVsPtKGenProj->GetRMS(), 
      0.2*all,
      hDeltaPiVsPtPGenProj->GetMean(), 
      hDeltaPiVsPtPGenProj->GetRMS(), 
    };

    cFits->cd();
    cFits->cd(j + 1);

    total->SetParameters(gausParams);
    for(Int_t i = 0; i < 9; i++) {

      if((i%3) > 0)
	total->FixParameter(i, gausParams[i]);
      // else {
      // 	if(bin>56)
      // 	  total->SetParLimits(i, 0, 10*all);
      // }
    }
    
    if(performFit) {

      hDeltaPiVsPtProj->Fit(total, "0L");

    } 

    hDeltaPiVsPtProj->DrawCopy();
    total->DrawCopy("same");    
    
    Double_t parametersOut[9];
    total->GetParameters(parametersOut);
    const Double_t* parameterErrorsOut = total->GetParErrors();

    for(Int_t i = 0; i < 9; i++) 
      cout << parametersOut[i] << ", ";
    cout << endl;


    pion->SetParameters(&parametersOut[0]);
    pion->DrawCopy("same");
    hPionRatio->SetBinContent(bin, parametersOut[0]/all);
    hPionRatio->SetBinError(bin, parameterErrorsOut[0]/all);
    hPionYield->SetBinContent(bin, parametersOut[0]);
    hPionYield->SetBinError(bin, parameterErrorsOut[0]);

    kaon->SetParameters(&parametersOut[3]);
    kaon->DrawCopy("same");
    hKaonRatio->SetBinContent(bin, parametersOut[3]/all);
    hKaonRatio->SetBinError(bin, parameterErrorsOut[3]/all);
    hKaonYield->SetBinContent(bin, parametersOut[3]);
    hKaonYield->SetBinError(bin, parameterErrorsOut[3]);
    
    proton->SetParameters(&parametersOut[6]);
    proton->DrawCopy("same");
    hProtonRatio->SetBinContent(bin, parametersOut[6]/all);
    hProtonRatio->SetBinError(bin, parameterErrorsOut[6]/all);
    hProtonYield->SetBinContent(bin, parametersOut[6]);
    hProtonYield->SetBinError(bin, parameterErrorsOut[6]);


    cSingleFit->cd();
    cSingleFit->Clear();
    //    cSingleFit->SetLogy();
    hDeltaPiVsPtProj->Draw();
    pion->DrawCopy("same");
    kaon->DrawCopy("same");
    proton->DrawCopy("same");
    total->DrawCopy("same");
    
    if(bin==42 || bin==49) {

      TFile* fileOut = new TFile(Form("%s/ptspectrum_bin%d.root", dirName, bin), "RECREATE");

      TH1D* hist = (TH1D*)hDeltaPiVsPtProj->Clone("hist");

      TF1* sumFit    = (TF1*)total->Clone("sumFit");
      TF1* pionFit   = (TF1*)pion->Clone("pionFit");
      TF1* kaonFit   = (TF1*)kaon->Clone("kaonFit");
      TF1* protonFit = (TF1*)proton->Clone("protonFit");

      hist->Write();

      sumFit->Write();
      pionFit->Write();
      kaonFit->Write();
      protonFit->Write();

      fileOut->Close();
    }

    gROOT->ProcessLine(".x drawText.C");
    cSingleFit->SaveAs(Form("%s/ptspectrum_bin%d.gif", dirName, bin));
    cSingleFit->SaveAs(Form("%s/ptspectrum_bin%d.pdf", dirName, bin));
    //    legend->Draw();

    if(data->IsMc()) {

      cSingleFit->cd();
      cSingleFit->Clear();
      TH1D* hDeltaPiVsPtPiProj =(TH1D*)hDeltaPiVsPtPi->ProjectionY(Form("hDeltaPiVsPtPiProj%d", bin), bin, bin);
      hDeltaPiVsPtPiProj->SetMarkerStyle(20);
      hDeltaPiVsPtPiProj->SetMarkerColor(2);
      hDeltaPiVsPtPiProj->SetTitle(Form("%.1f<p_{T}<%.1f GeV/c", 
					hDeltaPiVsPt->GetXaxis()->GetBinLowEdge(bin),
					hDeltaPiVsPt->GetXaxis()->GetBinUpEdge(bin)));
      hDeltaPiVsPtPiProj->Draw("P");
      Double_t integral = hDeltaPiVsPtPiProj->Integral();
      hPionRatioMc->SetBinContent(bin, integral/all);
      hPionRatioMc->SetBinError(bin, TMath::Sqrt(integral)/all);
      hPionYieldMc->SetBinContent(bin, integral);
      hPionYieldMc->SetBinError(bin, TMath::Sqrt(integral));
      TH1D* hDeltaPiVsPtKProj =(TH1D*)hDeltaPiVsPtK->ProjectionY(Form("hDeltaPiVsPtKProj%d", bin), bin, bin);
      hDeltaPiVsPtKProj->SetMarkerStyle(21);
      hDeltaPiVsPtKProj->SetMarkerColor(3);
      hDeltaPiVsPtKProj->Draw("SAME P");
      integral = hDeltaPiVsPtKProj->Integral();
      hKaonRatioMc->SetBinContent(bin, integral/all);
      hKaonRatioMc->SetBinError(bin, TMath::Sqrt(integral)/all);
      hKaonYieldMc->SetBinContent(bin, integral);
      hKaonYieldMc->SetBinError(bin, TMath::Sqrt(integral));
      TH1D* hDeltaPiVsPtPProj =(TH1D*)hDeltaPiVsPtP->ProjectionY(Form("hDeltaPiVsPtPProj%d", bin), bin, bin);
      hDeltaPiVsPtPProj->SetMarkerStyle(22);
      hDeltaPiVsPtPProj->SetMarkerColor(4);
      hDeltaPiVsPtPProj->Draw("SAME P");
      integral = hDeltaPiVsPtPProj->Integral();
      hProtonRatioMc->SetBinContent(bin, integral/all);
      hProtonRatioMc->SetBinError(bin, TMath::Sqrt(integral)/all);
      hProtonYieldMc->SetBinContent(bin, integral);
      hProtonYieldMc->SetBinError(bin, TMath::Sqrt(integral));

      pion->DrawCopy("same");
      kaon->DrawCopy("same");
      proton->DrawCopy("same");
      cSingleFit->SaveAs(Form("debugfitsmc/ptspectrum_bin%d.gif", bin));
    }
  }

  TCanvas* cRatio = new TCanvas("cRatio", "ratios/all vs p", 600, 400);
  cRatio->Clear();
  hPionRatio->SetMarkerStyle(20);
  hPionRatio->SetMarkerColor(2);
  hPionRatio->GetXaxis()->SetRangeUser(0.0, ptStop-0.001);
  hPionRatio->GetYaxis()->SetRangeUser(0.0, 1.0);
  hPionRatio->DrawCopy("P E");
  hPionYield->SetMarkerStyle(20);
  hPionYield->SetMarkerColor(2);
  hPionYield->GetXaxis()->SetRangeUser(0.0, ptStop-0.001);

  hKaonRatio->SetMarkerStyle(20);
  hKaonRatio->SetMarkerColor(3);
  hKaonRatio->GetXaxis()->SetRangeUser(0.0, ptStop-0.001);
  hKaonRatio->DrawCopy("SAME P E");
  hKaonYield->SetMarkerStyle(20);
  hKaonYield->SetMarkerColor(3);
  hKaonYield->GetXaxis()->SetRangeUser(0.0, ptStop-0.001);

  hProtonRatio->SetMarkerStyle(20);
  hProtonRatio->SetMarkerColor(4);
  hProtonRatio->GetXaxis()->SetRangeUser(0.0, ptStop-0.001);
  hProtonRatio->DrawCopy("SAME P E");
  hProtonYield->SetMarkerStyle(20);
  hProtonYield->SetMarkerColor(4);
  hProtonYield->GetXaxis()->SetRangeUser(0.0, ptStop-0.001);

  gROOT->ProcessLine(".x drawText.C");
  cRatio->SaveAs(Form("%s/particle_ratios.gif", dirName));
  cRatio->SaveAs(Form("%s/particle_ratios.pdf", dirName));

  if(data->IsMc()) {
    
    hPionRatioMc->SetMarkerStyle(24);
    hPionRatioMc->SetMarkerColor(2);
    hPionRatioMc->GetXaxis()->SetRangeUser(0.0, ptStop-0.001);
    hPionYieldMc->GetXaxis()->SetRangeUser(0.0, ptStop-0.001);
    hPionYieldMc->SetMarkerStyle(24);
    hPionYieldMc->SetMarkerColor(2);
    hPionRatioMc->DrawCopy("SAME P");
    
    hKaonRatioMc->SetMarkerStyle(24);
    hKaonRatioMc->SetMarkerColor(3);
    hKaonRatioMc->GetXaxis()->SetRangeUser(0.0, ptStop-0.001);
    hKaonYieldMc->GetXaxis()->SetRangeUser(0.0, ptStop-0.001);
    hKaonYieldMc->SetMarkerStyle(24);
    hKaonYieldMc->SetMarkerColor(3);
    hKaonRatioMc->DrawCopy("SAME P");
    
    hProtonRatioMc->SetMarkerStyle(24);
    hProtonRatioMc->SetMarkerColor(4);
    hProtonRatioMc->GetXaxis()->SetRangeUser(0.0, ptStop-0.001);
    hProtonYieldMc->GetXaxis()->SetRangeUser(0.0, ptStop-0.001);
    hProtonYieldMc->SetMarkerStyle(24);
    hProtonYieldMc->SetMarkerColor(4);
    hProtonRatioMc->DrawCopy("SAME P");

    cRatio->SaveAs("debugfitsmc/particle_ratios.gif");
  }

  if(outFileName) {

    TFile* fileOut = new TFile(Form("fit_yields_results/%s", outFileName), "RECREATE");

    // Histograms for normalization
    data->GetHistVtxStatus()->Write();
    TH1D* hPhiCutEff = data->GetHistNclVsPhiVsPtAfter()->ProjectionX("hPhiCutEff");
    TH1D* hHelp = data->GetHistNclVsPhiVsPtBefore()->ProjectionX("hHelp");
    hPhiCutEff->Divide(hPhiCutEff, hHelp, 1.0, 1.0, "B");
    hPhiCutEff->Write();
    delete hHelp;
    delete hPhiCutEff;

    hPionRatio->Write();
    hKaonRatio->Write();
    hProtonRatio->Write();

    hPionYield->Write();
    hKaonYield->Write();
    hProtonYield->Write();

    if(data->IsMc()) {
      hPionRatioMc->Write();
      hKaonRatioMc->Write();
      hProtonRatioMc->Write();

      hPionYieldMc->Write();
      hKaonYieldMc->Write();
      hProtonYieldMc->Write();			       
    }

    fileOut->Close();
  }
}

//____________________________________________________________________________
void FitYieldsV0(const Char_t* dataFileName,
		 const Char_t* v0FileName,
		 const Char_t* v0Name,
		 Double_t ptStart, Double_t ptStop,
		 Int_t charge,
		 Bool_t performFit = kFALSE,
		 Int_t filter = 1,
		 Int_t run    = 0,
		 Bool_t usePhiCut = kTRUE,
		 const Char_t* endName=0)
{
  gStyle->SetOptStat(0);
  
  TFile* dataFile = FindFileFresh(dataFileName);
  if(!dataFile)
    return;
  AliHighPtDeDxData* data = (AliHighPtDeDxData*)GetObject(dataFile, filter, usePhiCut, run, "filter", endName);
  data->Print();

  TFile* v0File = FindFileFresh(v0FileName);
  if(!v0File)
    return;
  AliHighPtDeDxData* v0 = (AliHighPtDeDxData*)GetObject(v0File, 0, usePhiCut, run, v0Name, endName);
  v0->Print();

  TH2D* hDeltaPiVsPt = v0->GetHistDeltaPiVsPt(0, charge);
  hDeDxVsP = hDeltaPiVsPt; // for the 2d fit to pick up the right bin

  TH2D* hDeltaPiVsPtPiGen = data->GetHistDeltaPiVsPt(1, charge);
  TH2D* hDeltaPiVsPtKGen  = data->GetHistDeltaPiVsPt(2, charge);
  TH2D* hDeltaPiVsPtPGen  = data->GetHistDeltaPiVsPt(3, charge);

  TProfile* hPiGenProfile = hDeltaPiVsPtPiGen->ProfileX();
  hPiGenProfile->SetMarkerStyle(29);
  TProfile* hKGenProfile = hDeltaPiVsPtKGen->ProfileX();
  hKGenProfile->SetMarkerStyle(29);
  TProfile* hPGenProfile = hDeltaPiVsPtPGen->ProfileX();
  hPGenProfile->SetMarkerStyle(29);

  TCanvas* cDeltaPiVsPt = new TCanvas("cDeltaPiVsPt", "dE/dx vs p", 400, 300);
  cDeltaPiVsPt->Clear();
  cDeltaPiVsPt->cd();
  cDeltaPiVsPt->SetLogz();
  hDeltaPiVsPt->Draw("COLZ");
  hPiGenProfile->Draw("SAME P");
  hKGenProfile->Draw("SAME P");
  hPGenProfile->Draw("SAME P");

  TCanvas* cDeltaPiVsPtLogX = new TCanvas("cDeltaPiVsPtLogX", "dE/dx vs p", 400, 300);
  cDeltaPiVsPtLogX->Clear();
  cDeltaPiVsPtLogX->cd();
  cDeltaPiVsPtLogX->SetLogz();
  cDeltaPiVsPtLogX->SetLogx();
  hDeltaPiVsPt->Draw("COLZ");
  hPiGenProfile->Draw("SAME P");
  hKGenProfile->Draw("SAME P");
  hPGenProfile->Draw("SAME P");

  // Root is a bit stupid with finidng bins so we have to add and subtract a
  // little to be sure we get the right bin as we typically put edges as
  // limits
  const Int_t binStart = hDeltaPiVsPt->GetXaxis()->FindBin(ptStart+0.01);
  ptStart = hDeltaPiVsPt->GetXaxis()->GetBinLowEdge(binStart);
  const Int_t binStop  = hDeltaPiVsPt->GetXaxis()->FindBin(ptStop-0.01);
  ptStop = hDeltaPiVsPt->GetXaxis()->GetBinUpEdge(binStop);
  //  const Int_t nBins    = binStop - binStart + 1;

  cout << "Doing fits from pTlow = " << ptStart << " (bin: " << binStart
       << ") to pThigh = " << ptStop << " (bin: " << binStop << ")" << endl;
  

  //cross check
  TCanvas* cFits = new TCanvas("cFits", "Fit comparison to data", 1200, 800);
  cFits->Clear();
  cFits->Divide(7, 4);

  TF1* pion = new TF1("pion", "[0]/sqrt(6.2832*[2]*[2])*exp(-(x-[1])*(x-[1])/2.0/[2]/[2])", -30, 20);
    //  TF1* pion = new TF1("pion", "gausn", -30, 20);
  pion->SetLineWidth(2);
  pion->SetLineColor(kRed);
  // TF1* kaon = new TF1("kaon", "gausn", -30, 20);
  // kaon->SetLineWidth(2);
  // kaon->SetLineColor(kGreen);
  TF1* proton = new TF1("proton", "[0]/sqrt(6.2832*[2]*[2])*exp(-(x-[1])*(x-[1])/2.0/[2]/[2])", -30, 20);
  //  TF1* proton = new TF1("proton", "gausn", -30, 20);
  proton->SetLineWidth(2);
  proton->SetLineColor(kBlue);

  TLegend* legend = new TLegend(0.11, 0.6, 0.35, 0.85);    
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->AddEntry(pion, "#pi", "L");
  //  legend->AddEntry(kaon, "K", "L");
  legend->AddEntry(proton, "p", "L");

  if(strstr(v0Name, "lambda")) {
    gSystem->Exec("mv debugfitslambda/* olddebugfitslambda/");
  } else {
    gSystem->Exec("mv debugfitskaon/* olddebugfitskaon/");
  }
  TCanvas* cSingleFit = new TCanvas("cSingleFit", "single fit");

  for(Int_t bin = binStart; bin <= binStop; bin++){
    
    cout << "Making projection for bin: " << bin << endl;
    
    const Int_t j = bin-binStart;
    
    TH1D* hDeltaPiVsPtProj =(TH1D*)hDeltaPiVsPt->ProjectionY(Form("hDeltaPiVsPtProj%d", bin), bin, bin);
    //    hDeltaPiVsPtProj->GetXaxis()->SetRangeUser(40, 85);
    hDeltaPiVsPtProj->SetTitle(Form("%.1f<p_{T}<%.1f GeV/c", 
  				hDeltaPiVsPt->GetXaxis()->GetBinLowEdge(bin),
  				hDeltaPiVsPt->GetXaxis()->GetBinUpEdge(bin)));

    Double_t all =  hDeltaPiVsPtProj->GetEntries();

    TH1D* hDeltaPiVsPtPiGenProj =(TH1D*)hDeltaPiVsPtPiGen->ProjectionY(Form("hDeltaPiVsPtPiGenProj%d", bin), bin, bin);
    //    TH1D* hDeltaPiVsPtKGenProj =(TH1D*)hDeltaPiVsPtKGen->ProjectionY(Form("hDeltaPiVsPtKGenProj%d", bin), bin, bin);
    TH1D* hDeltaPiVsPtPGenProj =(TH1D*)hDeltaPiVsPtPGen->ProjectionY(Form("hDeltaPiVsPtPGenProj%d", bin), bin, bin);
    
    cFits->cd();
    cFits->cd(j + 1);

    hDeltaPiVsPtProj->Draw();

    if(strstr(v0Name, "lambda")) {
    
      Double_t mean  = hDeltaPiVsPtPGenProj->GetMean();
      Double_t sigma = hDeltaPiVsPtPGenProj->GetRMS();

      cout << "Bin: " << bin 
	   <<   " (" << hDeltaPiVsPtProj->GetTitle()				
	   << "), sigma: " << sigma << endl;
      proton->SetParameters(all, mean, sigma);
      proton->FixParameter(1, mean);
      proton->FixParameter(2, sigma);

      if(performFit) {
	
      	hDeltaPiVsPtProj->Fit(proton, "L0", "", mean-sigma, mean+sigma);
	
      	cout << "Bin: " << bin 
      	     <<   " (" << hDeltaPiVsPtProj->GetTitle()				
      	     << "), sigma: " << proton->GetParameter(2) << endl;
      } 

      proton->DrawCopy("same");

      cSingleFit->cd();
      cSingleFit->Clear();
      hDeltaPiVsPtProj->Draw();
      proton->DrawCopy("same");
    } else {

      Double_t mean  = hDeltaPiVsPtPiGenProj->GetMean();
      Double_t sigma = hDeltaPiVsPtPiGenProj->GetRMS();
      cout << "Bin: " << bin 
	   <<   " (" << hDeltaPiVsPtProj->GetTitle()				
	   << "), sigma: " << sigma << endl;
      pion->SetParameters(all, mean, sigma);
      pion->FixParameter(1, mean);
      pion->FixParameter(2, sigma);
      
      if(performFit) {
	
      	hDeltaPiVsPtProj->Fit(pion, "L0", "", mean-sigma, mean+sigma);

      	cout << "Bin: " << bin 
      	     <<   " (" << hDeltaPiVsPtProj->GetTitle()				
      	     << "), sigma: " << pion->GetParameter(2) << endl;	
      } 
      pion->DrawCopy("same");

      cSingleFit->cd();
      cSingleFit->Clear();
      hDeltaPiVsPtProj->Draw();
      pion->DrawCopy("same");
    }
    
    //    cSingleFit->SetLogy();
    
    gROOT->ProcessLine(".x drawText.C");

    if(strstr(v0Name, "lambda")) {

      if(!endName) {
	cSingleFit->SaveAs(Form("results/lambda/ptspectrum_bin%d.gif", bin));
	cSingleFit->SaveAs(Form("results/lambda/ptspectrum_bin%d.pdf", bin));
      } else {
	cSingleFit->SaveAs(Form("results/lambda%s/ptspectrum_bin%d.gif", endName, bin));
	cSingleFit->SaveAs(Form("results/lambda%s/ptspectrum_bin%d.pdf", endName, bin));
      }
    } else {
      if(!endName) {
	cSingleFit->SaveAs(Form("results/kaon/ptspectrum_bin%d.gif", bin));
	cSingleFit->SaveAs(Form("results/kaon/ptspectrum_bin%d.pdf", bin));
      }else {
	cSingleFit->SaveAs(Form("results/kaon%s/ptspectrum_bin%d.gif", endName, bin));
	cSingleFit->SaveAs(Form("results/kaon%s/ptspectrum_bin%d.pdf", endName, bin));
      }
    }
    //    legend->Draw();
  }

}

//____________________________________________________________________________
void CompareYields(const Char_t* dataFileName1,
		   const Char_t* dataFileName2,
		   Double_t ptStart, Double_t ptStop,
		   Int_t charge1,
		   Int_t charge2,
		   const Char_t* endName1=0,
		   const Char_t* endName2=0,
		   Bool_t performFit = kFALSE,
		   Int_t run    = 0,
		   Int_t filter = 1,
		   Bool_t usePhiCut = kTRUE)
{
  gStyle->SetOptStat(0);

  
  TFile* dataFile1 = FindFileFresh(dataFileName1);
  if(!dataFile1)
    return;
  AliHighPtDeDxData* data1 = (AliHighPtDeDxData*)GetObject(dataFile1, filter, usePhiCut, run, "filter", endName1);
  data1->Print();

  gSystem->Exec("mv debugfits/* olddebugfits/");
  gSystem->Exec("mv debugfitsratio/* olddebugfitsratio/");
  // if(data1->IsMc())
  //   gSystem->Exec("mv debugfitsmc/* olddebugfitsmc/");


  TH2D* hDeltaPiVsPt1 = data1->GetHistDeltaPiVsPt(0, charge1);
  hDeDxVsP = hDeltaPiVsPt1; // for the 2d fit to pick up the right bin

  TH2D* hDeltaPiVsPtPiGen1 = data1->GetHistDeltaPiVsPt(1, charge1);
  TH2D* hDeltaPiVsPtKGen1  = data1->GetHistDeltaPiVsPt(2, charge1);
  TH2D* hDeltaPiVsPtPGen1  = data1->GetHistDeltaPiVsPt(3, charge1);

  // TH2D* hDeltaPiVsPtPi1 = 0;
  // TH2D* hDeltaPiVsPtK1  = 0;
  // TH2D* hDeltaPiVsPtP1  = 0;

  // if(data1->IsMc()) {

  //   hDeltaPiVsPtPi1 = data1->GetHistDeltaPiVsPtPiMc();
  //   hDeltaPiVsPtK1  = data1->GetHistDeltaPiVsPtKMc();
  //   hDeltaPiVsPtP1  = data1->GetHistDeltaPiVsPtPMc();
  // }

  AliHighPtDeDxData* data2 = data1;
  if(dataFileName2) {
    TFile* dataFile2 = FindFileFresh(dataFileName2);
    if(!dataFile2)
      return;
    data2 = (AliHighPtDeDxData*)GetObject(dataFile2, filter, usePhiCut, run, "filter", endName2);
    data2->Print();
  } else if (endName2) {

    data2 = (AliHighPtDeDxData*)GetObject(dataFile1, filter, usePhiCut, run, "filter", endName2);
  }

  TH2D* hDeltaPiVsPt2 = data2->GetHistDeltaPiVsPt(0, charge2);
  hDeDxVsP = hDeltaPiVsPt2; // for the 2d fit to pick up the right bin

  // TH2D* hDeltaPiVsPtPiGen2 = data2->GetHistDeltaPiVsPt(1, charge2);
  // TH2D* hDeltaPiVsPtKGen2  = data2->GetHistDeltaPiVsPt(2, charge2);
  // TH2D* hDeltaPiVsPtPGen2  = data2->GetHistDeltaPiVsPt(3, charge2);

  // TH2D* hDeltaPiVsPtPi2 = 0;
  // TH2D* hDeltaPiVsPtK2  = 0;
  // TH2D* hDeltaPiVsPtP2  = 0;

  // if(data2->IsMc()) {

  //   hDeltaPiVsPtPi2 = data2->GetHistDeltaPiVsPtPiMc();
  //   hDeltaPiVsPtK2  = data2->GetHistDeltaPiVsPtKMc();
  //   hDeltaPiVsPtP2  = data2->GetHistDeltaPiVsPtPMc();
  // }


  // Root is a bit stupid with finidng bins so we have to add and subtract a
  // little to be sure we get the right bin as we typically put edges as
  // limits
  const Int_t binStart = hDeltaPiVsPt1->GetXaxis()->FindBin(ptStart+0.01);
  ptStart = hDeltaPiVsPt1->GetXaxis()->GetBinLowEdge(binStart);
  const Int_t binStop  = hDeltaPiVsPt1->GetXaxis()->FindBin(ptStop-0.01);
  ptStop = hDeltaPiVsPt1->GetXaxis()->GetBinUpEdge(binStop);
  //  const Int_t nBins    = binStop - binStart + 1;

  cout << "Doing fits from pTlow = " << ptStart << " (bin: " << binStart
       << ") to pThigh = " << ptStop << " (bin: " << binStop << ")" << endl;
  

  //cross check
  TCanvas* cFits = new TCanvas("cFits", "Fit comparison to data", 1200, 800);
  cFits->Clear();
  cFits->Divide(7, 4);

  TF1* pion = new TF1("pion", "gausn", -30, 20);
  pion->SetLineWidth(2);
  pion->SetLineColor(kRed);
  TF1* kaon = new TF1("kaon", "gausn", -30, 20);
  kaon->SetLineWidth(2);
  kaon->SetLineColor(kGreen);
  TF1* proton = new TF1("proton", "gausn", -30, 20);
  proton->SetLineWidth(2);
  proton->SetLineColor(kBlue);
  TF1* total = new TF1("total", "gausn(0)+gausn(3)+gausn(6)", -30, 20);
  total->SetLineColor(kBlack);
  total->SetLineWidth(2);
  total->SetLineStyle(2);

  TLegend* legend = new TLegend(0.11, 0.6, 0.35, 0.85);    
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->AddEntry(total, "3-Gauss fit", "L");
  legend->AddEntry(pion, "#pi", "L");
  legend->AddEntry(kaon, "K", "L");
  legend->AddEntry(proton, "p", "L");

  TCanvas* cSingleFit = new TCanvas("cSingleFit", "single fit");

  for(Int_t bin = binStart; bin <= binStop; bin++){
    
    cout << "Making projection for bin: " << bin << endl;
    
    const Int_t j = bin-binStart;
    
    TH1D* hDeltaPiVsPtProj1 =(TH1D*)hDeltaPiVsPt1->ProjectionY(Form("hDeltaPiVsPtProj1%d", bin), bin, bin);
    hDeltaPiVsPtProj1->GetXaxis()->SetRangeUser(-25, 20);
    hDeltaPiVsPtProj1->SetTitle(Form("%.1f<p_{T}<%.1f GeV/c", 
  				hDeltaPiVsPt1->GetXaxis()->GetBinLowEdge(bin),
  				hDeltaPiVsPt1->GetXaxis()->GetBinUpEdge(bin)));

    TH1D* hDeltaPiVsPtProj2 =(TH1D*)hDeltaPiVsPt2->ProjectionY(Form("hDeltaPiVsPtProj2%d", bin), bin, bin);

    Double_t all1 =  hDeltaPiVsPtProj1->GetEntries();

    TH1D* hDeltaPiVsPtPiGenProj1 =(TH1D*)hDeltaPiVsPtPiGen1->ProjectionY(Form("hDeltaPiVsPtPiGenProj1%d", bin), bin, bin);
    TH1D* hDeltaPiVsPtKGenProj1 =(TH1D*)hDeltaPiVsPtKGen1->ProjectionY(Form("hDeltaPiVsPtKGenProj1%d", bin), bin, bin);
    TH1D* hDeltaPiVsPtPGenProj1 =(TH1D*)hDeltaPiVsPtPGen1->ProjectionY(Form("hDeltaPiVsPtPGenProj1%d", bin), bin, bin);
    
    Double_t gausParams[9] = { 
      0.6*all1,
      hDeltaPiVsPtPiGenProj1->GetMean(), 
      hDeltaPiVsPtPiGenProj1->GetRMS(), 
      0.2*all1,
      hDeltaPiVsPtKGenProj1->GetMean(), 
      hDeltaPiVsPtKGenProj1->GetRMS(), 
      0.2*all1,
      hDeltaPiVsPtPGenProj1->GetMean(), 
      hDeltaPiVsPtPGenProj1->GetRMS(), 
    };

    cFits->cd();
    cFits->cd(j + 1);

    total->SetParameters(gausParams);
    for(Int_t i = 0; i < 9; i++) {

      if((i%3) > 0)
	total->FixParameter(i, gausParams[i]);
    }
    
    if(performFit) {

      hDeltaPiVsPtProj1->Fit(total, "0L");

    } 

    hDeltaPiVsPtProj1->SetLineColor(4);
    hDeltaPiVsPtProj2->SetLineColor(2);
    hDeltaPiVsPtProj1->DrawCopy();
    hDeltaPiVsPtProj2->DrawCopy("SAME");
    if(performFit)
      total->DrawCopy("same");    
    
    Double_t parametersOut[9];
    total->GetParameters(parametersOut);
    //    const Double_t* parameterErrorsOut = total->GetParErrors();

    for(Int_t i = 0; i < 9; i++) 
      cout << parametersOut[i] << ", ";
    cout << endl;


    if(performFit) {
      pion->SetParameters(&parametersOut[0]);
      pion->DrawCopy("same");
      
      kaon->SetParameters(&parametersOut[3]);
      kaon->DrawCopy("same");
      
      proton->SetParameters(&parametersOut[6]);
      proton->DrawCopy("same");
    }

    cSingleFit->cd();
    cSingleFit->Clear();
    //    cSingleFit->SetLogy();
    hDeltaPiVsPtProj1->DrawCopy();
    hDeltaPiVsPtProj2->DrawCopy("SAME");
    if(performFit) {
      pion->DrawCopy("same");
      kaon->DrawCopy("same");
      proton->DrawCopy("same");
      total->DrawCopy("same");
    }

    cSingleFit->SaveAs(Form("debugfits/ptspectrum_bin%d.gif", bin));

    cSingleFit->cd();
    cSingleFit->Clear();
    //    cSingleFit->SetLogy();
    hDeltaPiVsPtProj1->Divide(hDeltaPiVsPtProj2);
    hDeltaPiVsPtProj1->GetYaxis()->SetRangeUser(0.8, 1.2);
    hDeltaPiVsPtProj1->DrawCopy();

    cSingleFit->SaveAs(Form("debugfitsratio/ptspectrum_bin%d.gif", bin));
    //    legend->Draw();

    // if(data->IsMc()) {

    //   cSingleFit->cd();
    //   cSingleFit->Clear();
    //   TH1D* hDeltaPiVsPtPiProj =(TH1D*)hDeltaPiVsPtPi->ProjectionY(Form("hDeltaPiVsPtPiProj%d", bin), bin, bin);
    //   hDeltaPiVsPtPiProj->SetMarkerStyle(20);
    //   hDeltaPiVsPtPiProj->SetMarkerColor(2);
    //   hDeltaPiVsPtPiProj->SetTitle(Form("%.1f<p_{T}<%.1f GeV/c", 
    // 					hDeltaPiVsPt->GetXaxis()->GetBinLowEdge(bin),
    // 					hDeltaPiVsPt->GetXaxis()->GetBinUpEdge(bin)));
    //   hDeltaPiVsPtPiProj->Draw("P");
    //   hPionRatioMc->SetBinContent(bin, hDeltaPiVsPtPiProj->Integral()/all);
    //   hPionYieldMc->SetBinContent(bin, hDeltaPiVsPtPiProj->Integral());
    //   TH1D* hDeltaPiVsPtKProj =(TH1D*)hDeltaPiVsPtK->ProjectionY(Form("hDeltaPiVsPtKProj%d", bin), bin, bin);
    //   hDeltaPiVsPtKProj->SetMarkerStyle(21);
    //   hDeltaPiVsPtKProj->SetMarkerColor(3);
    //   hDeltaPiVsPtKProj->Draw("SAME P");
    //   hKaonRatioMc->SetBinContent(bin, hDeltaPiVsPtKProj->Integral()/all);
    //   hKaonYieldMc->SetBinContent(bin, hDeltaPiVsPtKProj->Integral());
    //   TH1D* hDeltaPiVsPtPProj =(TH1D*)hDeltaPiVsPtP->ProjectionY(Form("hDeltaPiVsPtPProj%d", bin), bin, bin);
    //   hDeltaPiVsPtPProj->SetMarkerStyle(22);
    //   hDeltaPiVsPtPProj->SetMarkerColor(4);
    //   hDeltaPiVsPtPProj->Draw("SAME P");
    //   hProtonRatioMc->SetBinContent(bin, hDeltaPiVsPtPProj->Integral()/all);
    //   hProtonYieldMc->SetBinContent(bin, hDeltaPiVsPtPProj->Integral());

    //   pion->DrawCopy("same");
    //   kaon->DrawCopy("same");
    //   proton->DrawCopy("same");
    //   cSingleFit->SaveAs(Form("debugfitsmc/ptspectrum_bin%d.gif", bin));
    // }
  }
}

//___________________________________________________________________________________________
void MakeRatios(const Char_t* file1Name, const Char_t* file2Name, 
		Bool_t drawFractionRatios,
		Bool_t drawYieldRatios)
{
  /*
    For yields we assume that file 1 is negative charge and file
    2 is positive charge.
   */
  
  TFile* file1 = FindFileFresh(file1Name);
  if(!file1)
    return;
  
  TH1D* hPionRatio1   = (TH1D*)file1->Get("hPionRatio");
  TH1D* hKaonRatio1   = (TH1D*)file1->Get("hKaonRatio");
  TH1D* hProtonRatio1 = (TH1D*)file1->Get("hProtonRatio");

  TH1D* hPionYield1   = (TH1D*)file1->Get("hPionYield");
  TH1D* hKaonYield1   = (TH1D*)file1->Get("hKaonYield");
  TH1D* hProtonYield1 = (TH1D*)file1->Get("hProtonYield");

  TH1D* hPionRatioMc1   = (TH1D*)file1->Get("hPionRatioMc");
  TH1D* hKaonRatioMc1   = (TH1D*)file1->Get("hKaonRatioMc");
  TH1D* hProtonRatioMc1 = (TH1D*)file1->Get("hProtonRatioMc");

  TH1D* hPionYieldMc1   = (TH1D*)file1->Get("hPionYieldMc");
  TH1D* hKaonYieldMc1   = (TH1D*)file1->Get("hKaonYieldMc");
  TH1D* hProtonYieldMc1 = (TH1D*)file1->Get("hProtonYieldMc");

  TFile* file2 = FindFileFresh(file2Name);
  if(!file2)
    return;
  
  TH1D* hPionRatio2   = (TH1D*)file2->Get("hPionRatio");
  TH1D* hKaonRatio2   = (TH1D*)file2->Get("hKaonRatio");
  TH1D* hProtonRatio2 = (TH1D*)file2->Get("hProtonRatio");

  TH1D* hPionYield2   = (TH1D*)file2->Get("hPionYield");
  TH1D* hKaonYield2   = (TH1D*)file2->Get("hKaonYield");
  TH1D* hProtonYield2 = (TH1D*)file2->Get("hProtonYield");

  TH1D* hPionRatioMc2   = (TH1D*)file2->Get("hPionRatioMc");
  TH1D* hKaonRatioMc2   = (TH1D*)file2->Get("hKaonRatioMc");
  TH1D* hProtonRatioMc2 = (TH1D*)file2->Get("hProtonRatioMc");

  TH1D* hPionYieldMc2   = (TH1D*)file2->Get("hPionYieldMc");
  TH1D* hKaonYieldMc2   = (TH1D*)file2->Get("hKaonYieldMc");
  TH1D* hProtonYieldMc2 = (TH1D*)file2->Get("hProtonYieldMc");
  
  hPionRatio1->Divide(hPionRatio2);
  hPionRatio1->GetYaxis()->SetRangeUser(0.8, 1.2);
  hPionRatio1->SetTitle("Xcheck: pion fraction ratios; p_{T} [GeV/c]; pion fraction ratio");

  hKaonRatio1->Divide(hKaonRatio2);
  hKaonRatio1->GetYaxis()->SetRangeUser(0.8, 1.2);
  hKaonRatio1->SetTitle("Xcheck: kaon fraction ratios; p_{T} [GeV/c]; kaon fraction ratio");

  hProtonRatio1->Divide(hProtonRatio2);
  hProtonRatio1->GetYaxis()->SetRangeUser(0.8, 1.2);
  hProtonRatio1->SetTitle("Xcheck: proton fraction ratios; p_{T} [GeV/c]; proton fraction ratio");

  hPionYield1->Divide(hPionYield2);
  hPionYield1->GetYaxis()->SetRangeUser(0.8, 1.2);
  hPionYield1->SetTitle("#pi^{-}/#pi^{+} vs p_{T}; p_{T} [GeV/c]; #pi^{-}/#pi^{+}");

  hKaonYield1->Divide(hKaonYield2);
  hKaonYield1->GetYaxis()->SetRangeUser(0.8, 1.2);
  hKaonYield1->SetTitle("K^{-}/K^{+} vs p_{T}; p_{T} [GeV/c]; K^{-}/K^{+}");

  hProtonYield1->Divide(hProtonYield2);
  hProtonYield1->GetYaxis()->SetRangeUser(0.8, 1.2);
  hProtonYield1->SetTitle("#bar{p}/p vs p_{T}; p_{T} [GeV/c]; #bar{p}/p");

  if(hPionRatioMc1) {
    hPionRatioMc1->Divide(hPionRatioMc2);
    hPionRatioMc1->GetYaxis()->SetRangeUser(0.8, 1.2);
    hPionRatioMc1->SetTitle("Xcheck: pion fraction ratios (MC TRUTH); p_{T} [GeV/c]; pion fraction ratio");
    
    hKaonRatioMc1->Divide(hKaonRatioMc2);
    hKaonRatioMc1->GetYaxis()->SetRangeUser(0.8, 1.2);
    hKaonRatioMc1->SetTitle("Xcheck: kaon fraction ratios (MC TRUTH); p_{T} [GeV/c]; kaon fraction ratio");
    
    hProtonRatioMc1->Divide(hProtonRatioMc2);
    hProtonRatioMc1->GetYaxis()->SetRangeUser(0.8, 1.2);
    hProtonRatioMc1->SetTitle("Xcheck: proton fraction ratios (MC TRUTH); p_{T} [GeV/c]; proton fraction ratio");
    
    hPionYieldMc1->Divide(hPionYieldMc2);
    hPionYieldMc1->GetYaxis()->SetRangeUser(0.8, 1.2);
    hPionYieldMc1->SetTitle("#pi^{-}/#pi^{+} vs p_{T} (MC TRUTH); p_{T} [GeV/c]; #pi^{-}/#pi^{+}");
    
    hKaonYieldMc1->Divide(hKaonYieldMc2);
    hKaonYieldMc1->GetYaxis()->SetRangeUser(0.8, 1.2);
    hKaonYieldMc1->SetTitle("K^{-}/K^{+} vs p_{T} (MC TRUTH); p_{T} [GeV/c]; K^{-}/K^{+}");
    
    hProtonYieldMc1->Divide(hProtonYieldMc2);
    hProtonYieldMc1->GetYaxis()->SetRangeUser(0.8, 1.2);
    hProtonYieldMc1->SetTitle("#bar{p}/p vs p_{T} (MC TRUTH); p_{T} [GeV/c]; #bar{p}/p");
  }

  if(drawFractionRatios) {
    TCanvas* cPionFractionRatio = new TCanvas("cPionFractionRatio", "pion fraction ratio", 400, 300);
    cPionFractionRatio->Clear();
    cPionFractionRatio->SetGridy();
    cPionFractionRatio->cd();
    hPionRatio1->Draw();
    cPionFractionRatio->SaveAs("pion_frac_ratio.gif");

    TCanvas* cKaonFractionRatio = new TCanvas("cKaonFractionRatio", "kaon fraction ratio", 400, 300);
    cKaonFractionRatio->Clear();
    cKaonFractionRatio->SetGridy();
    cKaonFractionRatio->cd();
    hKaonRatio1->Draw();
    cKaonFractionRatio->SaveAs("kaon_frac_ratio.gif");

    TCanvas* cProtonFractionRatio = new TCanvas("cProtonFractionRatio", "proton fraction ratio", 400, 300);
    cProtonFractionRatio->Clear();
    cProtonFractionRatio->SetGridy();
    cProtonFractionRatio->cd();
    hProtonRatio1->Draw();
    cProtonFractionRatio->SaveAs("proton_frac_ratio.gif");

    if(hPionRatioMc1) {
      TCanvas* cPionFractionRatioMc = new TCanvas("cPionFractionRatioMc", "pion fraction ratio", 400, 300);
      cPionFractionRatioMc->Clear();
      cPionFractionRatioMc->SetGridy();
      cPionFractionRatioMc->cd();
      hPionRatioMc1->Draw();
      cPionFractionRatioMc->SaveAs("pion_frac_ratio_mc.gif");
      
      TCanvas* cKaonFractionRatioMc = new TCanvas("cKaonFractionRatioMc", "kaon fraction ratio", 400, 300);
      cKaonFractionRatioMc->Clear();
      cKaonFractionRatioMc->SetGridy();
      cKaonFractionRatioMc->cd();
      hKaonRatioMc1->Draw();
      cKaonFractionRatioMc->SaveAs("kaon_frac_ratio_mc.gif");
      
      TCanvas* cProtonFractionRatioMc = new TCanvas("cProtonFractionRatioMc", "proton fraction ratio", 400, 300);
      cProtonFractionRatioMc->Clear();
      cProtonFractionRatioMc->SetGridy();
      cProtonFractionRatioMc->cd();
      hProtonRatioMc1->Draw();
      cProtonFractionRatioMc->SaveAs("proton_frac_ratio_mc.gif");
    }
  }
  
  if(drawYieldRatios) {
    TCanvas* cPionRatio = new TCanvas("cPionRatio", "pion yields ratio", 400, 300);
    cPionRatio->Clear();
    cPionRatio->cd();
    cPionRatio->SetGridy();
    hPionYield1->Draw();
    cPionRatio->SaveAs("pion_ratio.gif");

    TCanvas* cKaonRatio = new TCanvas("cKaonRatio", "kaon yields ratio", 400, 300);
    cKaonRatio->Clear();
    cKaonRatio->cd();
    cKaonRatio->SetGridy();
    hKaonYield1->Draw();
    cKaonRatio->SaveAs("kaon_ratio.gif");

    TCanvas* cProtonRatio = new TCanvas("cProtonRatio", "proton yields ratio", 400, 300);
    cProtonRatio->Clear();
    cProtonRatio->cd();
    cProtonRatio->SetGridy();
    hProtonYield1->Draw();
    cProtonRatio->SaveAs("proton_ratio.gif");

    if(hPionRatioMc1) {
      
      TCanvas* cPionRatioMc = new TCanvas("cPionRatioMc", "pion yields ratio", 400, 300);
      cPionRatioMc->Clear();
      cPionRatioMc->cd();
      cPionRatioMc->SetGridy();
      hPionYieldMc1->Draw();
      cPionRatioMc->SaveAs("pion_ratio_mc.gif");
      
      TCanvas* cKaonRatioMc = new TCanvas("cKaonRatioMc", "kaon yields ratio", 400, 300);
      cKaonRatioMc->Clear();
      cKaonRatioMc->cd();
      cKaonRatioMc->SetGridy();
      hKaonYieldMc1->Draw();
      cKaonRatioMc->SaveAs("kaon_ratio_mc.gif");
      
      TCanvas* cProtonRatioMc = new TCanvas("cProtonRatioMc", "proton yields ratio", 400, 300);
      cProtonRatioMc->Clear();
      cProtonRatioMc->cd();
      cProtonRatioMc->SetGridy();
      hProtonYieldMc1->Draw();
      cProtonRatioMc->SaveAs("proton_ratio_mc.gif");
    }
  }
}

//___________________________________________________________________________________________
void Compare(const Char_t* file1Name, const Char_t* file2Name, const Char_t* file3Name, 
	     const Char_t* legend2, const Char_t* legend3, const Char_t* outfilename)
{
  /*
    filename1 is the default
   */
  
  TLegend* legend = new TLegend(0.11, 0.68, 0.35, 0.88);    
  legend->SetBorderSize(0);
  legend->SetFillColor(0);

  TFile* file1 = FindFileFresh(file1Name);
  if(!file1)
    return;
  
  TH1D* hPionRatio1   = (TH1D*)file1->Get("hPionRatio");
  hPionRatio1->SetMarkerStyle(29);
  TH1D* hPionRatio1Clone = (TH1D*)hPionRatio1->Clone("hPionsRatio1Clone");
  hPionRatio1Clone->SetMarkerColor(1);
  legend->AddEntry(hPionRatio1Clone, "default", "P");
  TH1D* hKaonRatio1   = (TH1D*)file1->Get("hKaonRatio");
  hKaonRatio1->SetMarkerStyle(29);
  TH1D* hProtonRatio1 = (TH1D*)file1->Get("hProtonRatio");
  hProtonRatio1->SetMarkerStyle(29);

  TFile* file2 = FindFileFresh(file2Name);
  if(!file2)
    return;
  
  TH1D* hPionRatio2   = (TH1D*)file2->Get("hPionRatio");
  hPionRatio2->SetMarkerStyle(20);
  TH1D* hPionRatio2Clone = (TH1D*)hPionRatio2->Clone("hPionsRatio2Clone");
  hPionRatio2Clone->SetMarkerColor(1);
  legend->AddEntry(hPionRatio2Clone, legend2, "P");
  TH1D* hKaonRatio2   = (TH1D*)file2->Get("hKaonRatio");
  hKaonRatio2->SetMarkerStyle(20);
  TH1D* hProtonRatio2 = (TH1D*)file2->Get("hProtonRatio");
  hProtonRatio2->SetMarkerStyle(20);

  TFile* file3 = FindFileFresh(file3Name);
  if(!file3)
    return;
  
  TH1D* hPionRatio3   = (TH1D*)file3->Get("hPionRatio");
  hPionRatio3->SetMarkerStyle(24);
  TH1D* hPionRatio3Clone = (TH1D*)hPionRatio3->Clone("hPionsRatio3Clone");
  hPionRatio3Clone->SetMarkerColor(1);
  legend->AddEntry(hPionRatio3Clone, legend3, "P");
  TH1D* hKaonRatio3   = (TH1D*)file3->Get("hKaonRatio");
  hKaonRatio3->SetMarkerStyle(24);
  TH1D* hProtonRatio3 = (TH1D*)file3->Get("hProtonRatio");
  hProtonRatio3->SetMarkerStyle(24);

  TCanvas* cRatios = new TCanvas("cRatios", "pion fraction ratio", 400, 300);
  cRatios->Clear();
  cRatios->SetGridy();
  cRatios->cd();
  hPionRatio1->Draw("EP");
  hKaonRatio1->Draw("SAME EP");
  hProtonRatio1->Draw("SAME EP");
  hPionRatio2->Draw("HIST SAME P");
  hKaonRatio2->Draw("HIST SAME P");
  hProtonRatio2->Draw("HIST SAME P");
  hPionRatio3->Draw("HIST SAME P");
  hKaonRatio3->Draw("HIST SAME P");
  hProtonRatio3->Draw("HIST SAME P");
  legend->Draw();
  gROOT->ProcessLine(".x drawText.C");
  cRatios->SaveAs(Form("%s.gif", outfilename));
  cRatios->SaveAs(Form("%s.pdf", outfilename));
  
}

//___________________________________________________________________________________________
void Compare(Double_t x)
{
  /*
    filename1 is the default
   */
  gStyle->SetOptStat(0);

  TCanvas* cRatios = new TCanvas("cRatios", "pion fraction ratio", 400, 300);

  const Int_t nEtaBins = 8;
  Double_t etaLimits[nEtaBins+1] = {-0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8};

  TGraphErrors* graphPionFractions = new TGraphErrors(nEtaBins);
  graphPionFractions->SetMarkerStyle(29);
  graphPionFractions->SetMarkerColor(2);

  TGraphErrors* graphKaonFractions = new TGraphErrors(nEtaBins);
  graphKaonFractions->SetMarkerStyle(29);
  graphKaonFractions->SetMarkerColor(3);

  TGraphErrors* graphProtonFractions = new TGraphErrors(nEtaBins);
  graphProtonFractions->SetMarkerStyle(29);
  graphProtonFractions->SetMarkerColor(4);

  TH1F* hHelp = new TH1F("hHelp", "particle fractions vs #eta; #eta; particle fraction",
			 nEtaBins, etaLimits[0], etaLimits[nEtaBins]);
  hHelp->SetMinimum(0.0);
  hHelp->SetMaximum(1.0);
  hHelp->SetDirectory(0);

  for(Int_t i = 0; i < nEtaBins; i++) {

    Double_t etaCenter = (etaLimits[i] + etaLimits[i+1])/2.0;
    Double_t etaWidth  = TMath::Abs(etaCenter - etaLimits[i]);

    Int_t etaLow  = Int_t(TMath::Abs(etaLimits[i]*10.0));
    Int_t etaHigh = Int_t(TMath::Abs(etaLimits[i+1]*10.0));

    cout << etaCenter << ", " << etaWidth << ", " << etaLow << ", " << etaHigh << endl;

    TFile* file = FindFileFresh(Form("fit_yields_results/lhc10h_aod_all_eta_%d%d.root", 
				     etaLow, etaHigh));
    if(!file)
      return;
    
    TH1D* hPionRatio   = (TH1D*)file->Get("hPionRatio");

    Int_t bin = hPionRatio->FindBin(x);
    if(i==0)
      hHelp->SetTitle(Form("%s (%.1f<p_{T}<%.1f GeV/c)", 
			   hHelp->GetTitle(),
			   hPionRatio->GetXaxis()->GetBinLowEdge(bin),
			   hPionRatio->GetXaxis()->GetBinUpEdge(bin)));

    
    graphPionFractions->SetPoint(i, etaCenter, hPionRatio->GetBinContent(bin)); 
    graphPionFractions->SetPointError(i, etaWidth, hPionRatio->GetBinError(bin)); 
    
    TH1D* hKaonRatio   = (TH1D*)file->Get("hKaonRatio");
    graphKaonFractions->SetPoint(i, etaCenter, hKaonRatio->GetBinContent(bin)); 
    graphKaonFractions->SetPointError(i, etaWidth, hKaonRatio->GetBinError(bin)); 

    TH1D* hProtonRatio = (TH1D*)file->Get("hProtonRatio");
    graphProtonFractions->SetPoint(i, etaCenter, hProtonRatio->GetBinContent(bin)); 
    graphProtonFractions->SetPointError(i, etaWidth, hProtonRatio->GetBinError(bin));     
  }
  

  cRatios->Clear();
  cRatios->cd();
  hHelp->DrawCopy();
  graphPionFractions->Draw("P");
  graphKaonFractions->Draw("P");
  graphProtonFractions->Draw("P");
  gROOT->ProcessLine(".x drawText.C");
  cRatios->SaveAs(Form("results/comparison/ratios_vs_eta_%.1f.gif", x));
  cRatios->SaveAs(Form("results/comparison/ratios_vs_eta_%.1f.pdf", x));
  
  delete hHelp;
}
