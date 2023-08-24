#include <iostream>
#include <TH3D.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TVirtualFitter.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TMath.h>
#include "TLatex.h"

void DCA_plots_pT_bins(){
  const char *part_type[6] = {"Proton", "ProtonMinus", "Pion", "PionMinus", "Kaon", "KaonMinus"};
  const char *part_type_char[6] = {"p", "#bar{p}", "#pi^{+}", "#pi^{-}", "K^{+}", "K^{-}"};
  int type = 2;

  double pT_min = 0.3;
  double pT_max = 2.5;

  ofstream myfile;
  myfile.open(Form("output_pT_%.2f_%.2f.txt", pT_min, pT_max));

  double MC_XY_sum = 0.0;
  double MC_XY_prim = 0.0;
  double MC_XY_sW = 0.0;
  double MC_XY_sM = 0.0;

  double MC_Z_sum = 0.0;
  double MC_Z_prim = 0.0;
  double MC_Z_sW = 0.0;
  double MC_Z_sM = 0.0;

  
  TFile *ifile_mc = new TFile("Mergeds.dir.root");

  for (int ii = 0; ii < 6; ii++){
    type = ii;

    TH2D *primXY, *secMatXY, *secWeakXY, *fakeXY, *primZ, *secMatZ, *secWeakZ, *fakeZ;
    TH1D *primXYproj, *secMatXYproj, *secWeakXYproj, *fakeXYproj, *primZproj, *secMatZproj, *secWeakZproj, *fakeZproj;

    primXY = (TH2D*) ifile_mc->Get(Form("hPrimVsDCAXYM0%s", part_type[type]));
    secMatXY = (TH2D*) ifile_mc->Get(Form("hSecMatVsDCAXYM0%s", part_type[type]));
    secWeakXY = (TH2D*) ifile_mc->Get(Form("hSecWeakVsDCAXYM0%s", part_type[type]));
    fakeXY = (TH2D*) ifile_mc->Get(Form("hFakeVsDCAXYM0%s", part_type[type]));

    primZ = (TH2D*) ifile_mc->Get(Form("hPrimVsDCAZM0%s", part_type[type]));
    secMatZ = (TH2D*) ifile_mc->Get(Form("hSecMatVsDCAZM0%s", part_type[type]));
    secWeakZ = (TH2D*) ifile_mc->Get(Form("hSecWeakVsDCAZM0%s", part_type[type]));
    fakeZ = (TH2D*) ifile_mc->Get(Form("hFakeVsDCAZM0%s", part_type[type]));

    primXY->GetYaxis()->SetRangeUser(pT_min, pT_max);
    secMatXY->GetYaxis()->SetRangeUser(pT_min, pT_max);
    secWeakXY->GetYaxis()->SetRangeUser(pT_min, pT_max);
    primZ->GetYaxis()->SetRangeUser(pT_min, pT_max);
    secMatZ->GetYaxis()->SetRangeUser(pT_min, pT_max);
    secWeakZ->GetYaxis()->SetRangeUser(pT_min, pT_max);

    MC_XY_prim = primXY->Integral();
    MC_XY_sW = secWeakXY->Integral();
    MC_XY_sM = secMatXY->Integral();
    MC_XY_sum = MC_XY_prim + MC_XY_sW + MC_XY_sM;
    myfile << part_type[ii] << endl << endl;
    myfile << "DCA_XY: " << endl;
    myfile << "Prim. fraction: " << MC_XY_prim/MC_XY_sum << endl;
    myfile << "Sec. weak fraction: " << MC_XY_sW/MC_XY_sum << endl;
    myfile << "Sec. mat. fraction: " << MC_XY_sM/MC_XY_sum << endl << endl;

    MC_Z_prim = primZ->Integral();
    MC_Z_sW = secWeakZ->Integral();
    MC_Z_sM = secMatZ->Integral();
    MC_Z_sum = MC_Z_prim + MC_Z_sW + MC_Z_sM;
    myfile << "DCA_Z: " << endl;
    myfile << "Prim. fraction: " << MC_Z_prim/MC_Z_sum << endl;
    myfile << "Sec. weak fraction: " << MC_Z_sW/MC_Z_sum << endl;
    myfile << "Sec. mat. fraction: " << MC_Z_sM/MC_Z_sum << endl << endl;


    primXYproj = primXY->ProjectionX();
    secMatXYproj = secMatXY->ProjectionX();
    secWeakXYproj = secWeakXY->ProjectionX();
    fakeXYproj = fakeXY->ProjectionX();

    primZproj = primZ->ProjectionX();
    secMatZproj = secMatZ->ProjectionX();
    secWeakZproj = secWeakZ->ProjectionX();
    fakeZproj = fakeZ->ProjectionX();

  
    TCanvas *canvXY =  new TCanvas(Form("DCA_XY_%d",1),Form("DCA_XY_%d",1),1200,600);
    canvXY->Divide(2,1);
    canvXY->cd(1);
    TPad* cXY = new TPad(Form("GraphsXY_%d",1),Form("GraphsXY_%d",1),0.01,0.05,0.95,0.95);
    gStyle->SetOptStat(000);
    gPad->SetLogy();

    int bin = primXYproj->GetNbinsX();
    int mirror_bin = 1;
    for(int oo = 0; oo < primXYproj->GetNbinsX()/2 - 1; oo++){
      primXYproj->SetBinContent(bin-1, primXYproj->GetBinContent(mirror_bin));
      primXYproj->SetBinError(bin-1, primXYproj->GetBinError(mirror_bin)); 
      secMatXYproj->SetBinContent(bin-1, secMatXYproj->GetBinContent(mirror_bin));
      secMatXYproj->SetBinError(bin-1, secMatXYproj->GetBinError(mirror_bin));
      secWeakXYproj->SetBinContent(bin-1, secWeakXYproj->GetBinContent(mirror_bin));
      secWeakXYproj->SetBinError(bin-1, secWeakXYproj->GetBinError(mirror_bin));
      fakeXYproj->SetBinContent(bin-1, fakeXYproj->GetBinContent(mirror_bin));
      bin--;
      mirror_bin++;
    }

    bin = primZproj->GetNbinsX();
    mirror_bin = 1;
    for(int oo = 0; oo < primZproj->GetNbinsX()/2 - 1; oo++){
      primZproj->SetBinContent(bin-1, primZproj->GetBinContent(mirror_bin));
      primZproj->SetBinError(bin-1, primZproj->GetBinError(mirror_bin)); 
      secMatZproj->SetBinContent(bin-1, secMatZproj->GetBinContent(mirror_bin));
      secMatZproj->SetBinError(bin-1, secMatZproj->GetBinError(mirror_bin));
      secWeakZproj->SetBinContent(bin-1, secWeakZproj->GetBinContent(mirror_bin));
      secWeakZproj->SetBinError(bin-1, secWeakZproj->GetBinError(mirror_bin));
      fakeZproj->SetBinContent(bin-1, fakeZproj->GetBinContent(mirror_bin));
      bin--;
      mirror_bin++;
    }

  
    primXYproj->SetTitle(Form("%s DCA_{XY} distribution", part_type_char[type]));
    primXYproj->GetXaxis()->SetTitle("DCA_{XY} (cm)");
    primXYproj->GetYaxis()->SetTitle("Entries");
    primXYproj->SetLineColor(kBlack);
    secMatXYproj->SetLineColor(kGreen+2);
    secWeakXYproj->SetLineColor(kRed);
    fakeXYproj->SetLineColor(kBlue);
    // primXYproj->Rebin(5);
    // secMatXYproj->Rebin(5);
    // secWeakXYproj->Rebin(5);

    primXYproj->GetXaxis()->SetRangeUser(-3.0, 3.0);
    if(ii == 0) primXYproj->GetYaxis()->SetRangeUser(10, 500000);
    else if(ii == 1) primXYproj->GetYaxis()->SetRangeUser(1, 500000);
    else if(ii == 2 || ii == 3) primXYproj->GetYaxis()->SetRangeUser(100, 50000000);
    else primXYproj->GetYaxis()->SetRangeUser(1, 10000000);

    
    primXYproj->Draw("hist");
    secMatXYproj->Draw("same hist");
    secWeakXYproj->Draw("same hist");
    //fakeXYproj->Draw("same hist");

    auto legend = new TLegend(0.55,0.67,0.88,0.88);
    legend->SetHeader("#sqrt{s_{NN}} = 5.02 TeV, Pb-Pb","C"); // option "C" allows to center the header
    legend->AddEntry(primXYproj,"Primaries","l");
    legend->AddEntry(secMatXYproj, "Sec. from material","l");
    legend->AddEntry(secWeakXYproj,"Sec. form weak decays","l");
    //legend->AddEntry(fakeXYproj, "Fake", "l");
    legend->SetBorderSize(0);
    legend->Draw();

    canvXY->cd(2);
    TPad* cZ = new TPad(Form("GraphsZ_%d",1),Form("GraphsZ_%d",1),0.01,0.05,0.95,0.95);
    gStyle->SetOptStat(000);
    gPad->SetLogy();
  
    primZproj->SetTitle(Form("%s DCA_{Z} distribution", part_type_char[type]));
    primZproj->GetXaxis()->SetTitle("DCA_{Z} (cm)");
    primZproj->GetYaxis()->SetTitle("Entries");
    primZproj->SetLineColor(kBlack);
    secMatZproj->SetLineColor(kGreen+2);
    secWeakZproj->SetLineColor(kRed);
    fakeZproj->SetLineColor(kBlue);
    
    primZproj->GetXaxis()->SetRangeUser(-3.0, 3.0);
    if(ii == 0) primZproj->GetYaxis()->SetRangeUser(10, 500000);
    else if(ii == 1) primZproj->GetYaxis()->SetRangeUser(1, 500000);
    else if(ii == 2 || ii == 3) primZproj->GetYaxis()->SetRangeUser(100, 50000000);
    else primZproj->GetYaxis()->SetRangeUser(1, 10000000);

    primZproj->Draw("hist");
    secMatZproj->Draw("same hist");
    secWeakZproj->Draw("same hist");

    canvXY->SaveAs(Form("DCA_plots_%s_pT_%.2f_%.2f.pdf", part_type[type], pT_min, pT_max));
  }
}
