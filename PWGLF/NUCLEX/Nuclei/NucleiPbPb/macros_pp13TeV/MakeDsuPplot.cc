#include "src/Common.h"
#include "src/Utils.h"
using namespace utils;
#include <TFile.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TColor.h>
#include <TLegendEntry.h>

void MakeDsuPplot(){
  gStyle->SetOptStat(0);

  TFile file_coalescence(Form("%sd_to_p_for_Max.root",kBaseOutputDir.data()));
  TCanvas* cvCoal = (TCanvas*)file_coalescence.Get("canvas");
  TGraph* grCoal = (TGraph*)cvCoal->GetPrimitive("two_body_coal");
  TGraph* grCoalUp = (TGraph*)cvCoal->GetPrimitive("two_body_coal_uncert2");
  TGraph* grCoalDown = (TGraph*)cvCoal->GetPrimitive("two_body_coal_uncert1");
  // grCoal->SetLineColor(TColor::GetColor("#8EA604"));
  grCoal->SetLineColor(kMagenta+2);
  grCoal->SetLineWidth(2);
  grCoal->SetFillStyle(0);
  grCoal->SetFillColor(0);
  grCoalUp->SetLineColor(kMagenta+2);
  grCoalUp->SetLineWidth(1);
  grCoalUp->SetFillStyle(0);
  grCoalUp->SetFillColor(0);
  grCoalDown->SetLineColor(kMagenta+2);
  grCoalDown->SetLineWidth(1);
  grCoalDown->SetFillStyle(0);
  grCoalDown->SetFillColor(0);
  const int n = grCoal->GetN();
  TGraph *grCoalShade = new TGraph(2*n);
  for(int i=0;i<n;i++) {
    grCoalShade->SetPoint(i,grCoal->GetX()[i],grCoalUp->GetY()[i]);
    grCoalShade->SetPoint(n+i,grCoal->GetX()[n-i-1],grCoalDown->GetY()[n-i-1]);
  }
  grCoalShade->SetFillStyle(3002);
  grCoalShade->SetFillColor(kMagenta+2);

  //grCoal->SetLineStyle();

  TFile* input_file = new TFile(Form("%sdoverpPaperProp.root",kBaseOutputDir.data()));
  TFile* prediction_file = new TFile(Form("%scsm_file.root",kBaseOutputDir.data()));
  TGraph* upper = (TGraph*)prediction_file->Get("upper");
  upper->SetLineStyle(kDashed);
  TGraph* lower = (TGraph*)prediction_file->Get("lower");
  TGraph* shade = (TGraph*)prediction_file->Get("shade");
  TCanvas* cv_input = (TCanvas*)input_file->Get("c1_n19");
  //cv_input->GetListOfPrimitives()->ls();
  TPad* pad = (TPad*)cv_input->GetPrimitive("c1_n19_1");
  //Previous results
  double dndeta_pp7TeV[5] = {3.295484, 7.539474, 10.76316, 13.5, 17.47021};
  double dndeta_pp7TeV_err[5] = {0.1345161, 0.2289474, 0.3, 0.4, 0.5202128};
  double dndeta_pp7TeV_err_null[5] = {0., 0., 0., 0., 0.};
  double dndeta_pp7TeV_dop[5] = {0.001283449, 0.001701064, 0.001970563, 0.001947286, 0.002206761};
  double dndeta_pp7TeV_dop_err_stat[5] = {1.515455e-05, 2.026621e-05, 2.621776e-05, 3.531506e-05, 3.161472e-05};
  double dndeta_pp7TeV_dop_err_syst[5] = {0.000215048, 0.0002194244, 0.0002575069, 0.0002606413, 0.0002978742};
  //
  TGraphErrors *pp7TeVStat = new TGraphErrors(5, dndeta_pp7TeV, dndeta_pp7TeV_dop, dndeta_pp7TeV_err_null, dndeta_pp7TeV_dop_err_stat);
  Requires(pp7TeVStat, "pp7TeVStat");
  pp7TeVStat->SetName("pp7TeVStat");
  pp7TeVStat->SetMarkerColor(TColor::GetColor("#006600"));
  pp7TeVStat->SetMarkerStyle(21);
  pp7TeVStat->SetMarkerSize(1.4);
  //
  TGraphErrors *pp7TeVSyst = new TGraphErrors(5, dndeta_pp7TeV, dndeta_pp7TeV_dop, dndeta_pp7TeV_err, dndeta_pp7TeV_dop_err_syst);
  Requires(pp7TeVSyst, "pp7TeVSyst");
  pp7TeVSyst->SetName("pp7TeVSyst");
  pp7TeVSyst->SetMarkerColor(TColor::GetColor("#006600"));
  pp7TeVSyst->SetLineColor(TColor::GetColor("#006600"));
  pp7TeVSyst->SetMarkerStyle(21);
  pp7TeVSyst->SetMarkerSize(1.4);
  pp7TeVSyst->SetFillStyle(0);
  pp7TeVSyst->SetFillColor(0);
  //
  TGraphErrors* PbPb5TeVStat = (TGraphErrors*)pad->GetPrimitive("Graph6");
  Requires(PbPb5TeVStat,"PbPb5TeVStat");
  PbPb5TeVStat->SetName("PbPb5TeVStat");
  TGraphErrors* PbPb5TeVSyst = (TGraphErrors*)pad->GetPrimitive("Graph5");
  Requires(PbPb5TeVSyst,"PbPb5TeVSyst");
  PbPb5TeVSyst->SetName("PbPb5TeVSyst");
  //
  double dndeta_pPb5TeV[5] = {40.6, 30.5, 23.3, 16.1, 7.1};
  double dndeta_pPb5TeV_err[5] = {0.9, 0.7, 0.5, 0.4, 0.2};
  double dndeta_pPb5TeV_err_null[5] = {0., 0., 0., 0., 0.};
  double dndeta_pPb5TeV_dop[5] = {0.00276558, 0.00264059, 0.00234453, 0.00207092, 0.00150940};
  double dndeta_pPb5TeV_dop_err_stat[5] = {2.87700e-05, 3.12207e-05, 2.45079e-05, 1.77377e-05, 2.68787e-05};
  double dndeta_pPb5TeV_dop_err_syst[5] = {0.000407579, 0.000403995, 0.000359305, 0.000337693, 0.000285351};
  //
  TGraphErrors *pPb5TeVStat = new TGraphErrors(5, dndeta_pPb5TeV, dndeta_pPb5TeV_dop, dndeta_pPb5TeV_err_null, dndeta_pPb5TeV_dop_err_stat);
  Requires(pPb5TeVStat,"pPb5TeVStat");
  pPb5TeVStat->SetName("pPb5TeVStat");
  pPb5TeVStat->SetMarkerColor(600);
  pPb5TeVStat->SetMarkerStyle(33);
  pPb5TeVStat->SetMarkerSize(1.4);
  //
  TGraphErrors *pPb5TeVSyst = new TGraphErrors(5, dndeta_pPb5TeV, dndeta_pPb5TeV_dop, dndeta_pPb5TeV_err, dndeta_pPb5TeV_dop_err_syst);
  Requires(pPb5TeVSyst,"pPb5TeVSyst");
  pPb5TeVSyst->SetName("pPb5TeVSyst");
  pPb5TeVSyst->SetMarkerColor(600);
  pPb5TeVSyst->SetLineColor(600);
  pPb5TeVSyst->SetMarkerStyle(33);
  pPb5TeVSyst->SetMarkerSize(1.4);
  pPb5TeVSyst->SetFillStyle(0);
  pPb5TeVSyst->SetFillColor(0);
  //
  TGraphErrors* PbPb2d7TeVStat = (TGraphErrors*)pad->GetPrimitive("ALICE_d");
  Requires(PbPb2d7TeVStat,"PbPb2d7TeVStat");
  PbPb2d7TeVStat->SetName("PbPb2d7TeVStat");
  TGraphAsymmErrors* PbPb2d7TeVSyst = (TGraphAsymmErrors*)pad->GetPrimitive("Graph4");
  Requires(PbPb2d7TeVSyst,"PbPb2d7TeVSyst");
  PbPb2d7TeVSyst->SetName("PbPb2d7TeVSyst");
  //
  TGraphAsymmErrors* pp0d9TeVSyst = new TGraphAsymmErrors(1);
  Requires(pp0d9TeVSyst,"pp0d9TeVSyst");
  pp0d9TeVSyst->SetName("pp0d9TeVSyst");
  pp0d9TeVSyst->SetMarkerColor(kYellow+1);
  pp0d9TeVSyst->SetLineColor(kYellow+1);
  pp0d9TeVSyst->SetFillStyle(0);
  pp0d9TeVSyst->SetMarkerStyle(29);
  pp0d9TeVSyst->SetMarkerSize(1.8);
  pp0d9TeVSyst->SetPoint(0,3.81,0.00138);
  pp0d9TeVSyst->SetPointError(0,0.07,0.07,0.00015,0.00015);
  TGraphAsymmErrors* pp0d9TeVStat = new TGraphAsymmErrors(1);
  Requires(pp0d9TeVStat,"pp0d9TeVStat");
  pp0d9TeVStat->SetName("pp0d9TeVStat");
  pp0d9TeVStat->SetMarkerColor(kYellow+1);
  pp0d9TeVStat->SetLineColor(kYellow+1);
  pp0d9TeVStat->SetFillStyle(0);
  pp0d9TeVStat->SetMarkerStyle(29);
  pp0d9TeVStat->SetMarkerSize(1.8);
  pp0d9TeVStat->SetPoint(0,3.81,0.00138);
  pp0d9TeVStat->SetPointError(0,0,0,0.00011,0.00011);
  //
  TGraphAsymmErrors* pp2d7TeVSyst = new TGraphAsymmErrors(1);
  pp2d7TeVSyst->SetName("pp2d7TeVSyst");
  pp2d7TeVSyst->SetMarkerColor(kOrange-7);
  pp2d7TeVSyst->SetLineColor(kOrange-7);
  pp2d7TeVSyst->SetFillStyle(0);
  pp2d7TeVSyst->SetMarkerStyle(45);
  pp2d7TeVSyst->SetMarkerSize(1.8);
  pp2d7TeVSyst->SetPoint(0,4.88,0.001482);
  pp2d7TeVSyst->SetPointError(0,0.09,0.13,0.00016,0.00016);
  TGraphAsymmErrors* pp2d7TeVStat = new TGraphAsymmErrors(1);
  pp2d7TeVStat->SetName("pp2d7TeVStat");
  pp2d7TeVStat->SetMarkerColor(kOrange-7);
  pp2d7TeVStat->SetLineColor(kOrange-7);
  pp2d7TeVStat->SetFillStyle(0);
  pp2d7TeVStat->SetMarkerStyle(45);
  pp2d7TeVStat->SetMarkerSize(1.8);
  pp2d7TeVStat->SetPoint(0,4.88,0.001482);
  pp2d7TeVStat->SetPointError(0,0,0,4.7e-05,4.7e-05);
  //
  //
  TGraphAsymmErrors* MBpp7TeVSyst = new TGraphAsymmErrors(1);
  MBpp7TeVSyst->SetName("MBpp7TeVSyst");
  MBpp7TeVSyst->SetMarkerColor(kGreen+3);
  MBpp7TeVSyst->SetLineColor(kGreen+3);
  MBpp7TeVSyst->SetFillStyle(0);
  MBpp7TeVSyst->SetMarkerStyle(31);
  MBpp7TeVSyst->SetMarkerSize(1.8);
  MBpp7TeVSyst->SetPoint(0,6.01,0.001626);
  MBpp7TeVSyst->SetPointError(0,0.12,0.2,0.00017,0.00017);
  TGraphAsymmErrors* MBpp7TeVStat = new TGraphAsymmErrors(1);
  MBpp7TeVStat->SetName("MBpp7TeVStat");
  MBpp7TeVStat->SetMarkerColor(kGreen+3);
  MBpp7TeVStat->SetLineColor(kGreen+3);
  MBpp7TeVStat->SetFillStyle(0);
  MBpp7TeVStat->SetMarkerStyle(31);
  MBpp7TeVStat->SetMarkerSize(1.8);
  MBpp7TeVStat->SetPoint(0,6.01,0.001626);
  MBpp7TeVStat->SetPointError(0,0,0,1.3e-05,1.3e-05);
  //

  //my work
  TFile* my_file = new TFile(Form("%sfinal_doverp.root",kBaseOutputDir.data()));
  TGraphErrors* pp13TeVStat = (TGraphErrors*)my_file->Get("ratio_gr_stat");
  Requires(pp13TeVStat,"pp13TeVStat");
  pp13TeVStat->SetName("pp13TeVStat");
  pp13TeVStat->SetMarkerStyle(47);
  pp13TeVStat->SetMarkerSize(1.5);
  TGraphErrors* pp13TeVSyst = (TGraphErrors*)my_file->Get("ratio_gr_syst");
  Requires(pp13TeVSyst,"pp13TeVSyst");
  pp13TeVSyst->SetName("pp13TeVSyst");
  pp13TeVSyst->SetMarkerStyle(47);
  pp13TeVSyst->SetMarkerSize(1.5);
  TGraphErrors* pp13TeVSystCorr = (TGraphErrors*)my_file->Get("ratio_gr_syst_corr");
  Requires(pp13TeVSystCorr,"pp13TeVSystCorr");
  pp13TeVSystCorr->SetName("pp13TeVSystCorr");
  pp13TeVSystCorr->SetMarkerStyle(47);
  pp13TeVSystCorr->SetMarkerSize(1.5);
  pp13TeVSystCorr->SetFillStyle(1001);

  TCanvas* cv_output = new TCanvas("cDoverPmy","cDoverPmy",800,600);
  cv_output->cd();
  cv_output->SetLeftMargin(0.14);
  cv_output->SetBottomMargin(0.15);
  cv_output->SetLogx();
  TH1F* hFrame = new TH1F("hFrame",";#LTd#it{N}_{ch} / d#it{#eta}_{lab}#GT_{|#it{#eta}_{lab}|<0.5}; 2d / (p + #bar{p})",9000,1,9001);
  hFrame->GetYaxis()->SetRangeUser(0.,0.006);
  hFrame->GetYaxis()->SetTitleOffset(1.3);
  hFrame->GetXaxis()->SetTitleOffset(1.2);
  hFrame->Draw();
  upper->Draw("same");
  lower->Draw("same");
  shade->Draw("samef");
  grCoal->Draw("L");
  grCoalShade->Draw("f");
  PbPb2d7TeVStat->Draw("samepz");
  PbPb2d7TeVSyst->Draw("e2same");
  pPb5TeVStat->Draw("samepz");
  pPb5TeVSyst->Draw("e2same");
  //PbPb5TeVStat->Draw("samepz");
  //PbPb5TeVSyst->Draw("e2same");
  // pp0d9TeVStat->Draw("samepz");
  // pp0d9TeVSyst->Draw("e2same");
  // pp2d7TeVStat->Draw("samepz");
  // pp2d7TeVSyst->Draw("e2same");
  // MBpp7TeVStat->Draw("samepz");
  // MBpp7TeVSyst->Draw("e2same");
  pp7TeVStat->Draw("samepz");
  pp7TeVSyst->Draw("e2same");
  pp13TeVSyst->Draw("e2same");
  pp13TeVSystCorr->Draw("e2same");
  pp13TeVStat->Draw("samepz");

  TLatex text;
  text.SetTextFont(63);
  text.SetTextSize(22);
  //text.DrawText(247.573,0.00529739,"ALICE Preliminary");

  TLegend* legPrel = new TLegend(0.18,0.56,0.305764,0.88,NULL,"brNDC");
  legPrel->SetBorderSize(0);
  legPrel->SetTextSize(0.03);
  legPrel->SetLineColor(1);
  legPrel->SetLineStyle(1);
  legPrel->SetLineWidth(1);
  legPrel->SetFillColor(0);
  legPrel->SetFillStyle(0);
  legPrel->SetHeader("#bf{ALICE}");
  legPrel->AddEntry(pPb5TeVSyst,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV", "pf");// (arXiv:1906.03136)","pf");
  legPrel->AddEntry((TObject*)nullptr,"V0A Multiplicity Classes (Pb-side)","");
  //legPrel->AddEntry(PbPb5TeVSyst,"Pb-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV","pf");
  legPrel->AddEntry(PbPb2d7TeVSyst,"Pb-Pb, #sqrt{#it{s}_{NN}} = 2.76 TeV", "pf");// (PRC 93 (2015) 024917)","pf");
  legPrel->AddEntry(pp7TeVSyst,"pp, #sqrt{#it{s}} = 7 TeV", "pf");// (PLB (2019) 50-63)","pf");
  legPrel->AddEntry(pp13TeVSyst,"pp, #sqrt{#it{s}} = 13 TeV","pf");
  legPrel->AddEntry((TObject*)nullptr,"V0M Multiplicity Classes","");
  legPrel->Draw();

  TLegend* legPub = new TLegend(0.47,0.22,0.55,0.42,NULL,"brNDC");
  legPub->SetBorderSize(0);
  legPub->SetTextSize(0.03);
  legPub->SetLineColor(1);
  legPub->SetLineStyle(1);
  legPub->SetLineWidth(1);
  legPub->SetFillColor(0);
  legPub->SetFillStyle(0);
  legPub->AddEntry("","Thermal-FIST CSM (PLB 785 (2018) 171-174)","");
  legPub->AddEntry(upper, "T_{ch} = 155 MeV, #it{V}_{c} = 3 d#it{V}/d#it{y}");
  legPub->AddEntry(lower, "T_{ch} = 155 MeV, #it{V}_{c} = d#it{V}/d#it{y}");
  legPub->AddEntry(grCoal, "Coalescence (PLB 792 (2019) 132-137)","l");
  // legPub->AddEntry(pp0d9TeVSyst,"pp, #sqrt{#it{s}} = 900 GeV, d/p (PRC 97 (2018) 024615)","pf");
  // legPub->AddEntry(pp2d7TeVSyst,"pp, #sqrt{#it{s}} = 2.76 TeV, d/p (PRC 97 (2018) 024615)","pf");
  // legPub->AddEntry(MBpp7TeVSyst,"pp, #sqrt{#it{s}} = 7 TeV, d/p (PRC 97 (2018) 024615)","pf");
  legPub->Draw();

  TFile output_file(Form("%scustomdoverp.root",kBaseOutputDir.data()),"recreate");
  cv_output->Print(Form("%sdoverp.pdf",kFiguresFolder.data()));
  cv_output->Write();
  hFrame->Write();
  pp7TeVStat->Write();
  pp7TeVSyst->Write();
  PbPb5TeVStat->Write();
  PbPb5TeVSyst->Write();
  PbPb2d7TeVStat->Write();
  PbPb2d7TeVSyst->Write();
  pPb5TeVStat->Write();
  pPb5TeVSyst->Write();
  pp13TeVStat->Write();
  pp13TeVSyst->Write();
  pp0d9TeVStat->Write();
  pp0d9TeVSyst->Write();
  pp2d7TeVStat->Write();
  pp2d7TeVSyst->Write();
  MBpp7TeVStat->Write();
  MBpp7TeVSyst->Write();
}
