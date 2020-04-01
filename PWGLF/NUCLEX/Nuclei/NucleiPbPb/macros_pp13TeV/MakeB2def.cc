#include "src/Common.h"
#include "src/Utils.h"
using namespace utils;
#include <TFile.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TStyle.h>
#include <TLatex.h>

void MakeB2def(){

  gStyle->SetOptStat(0);
  //
  TFile* input_file = new TFile(Form("%sB2vsMult_w13TeV.root",kBaseOutputDir.data()));
  TCanvas* cv_input = (TCanvas*)input_file->Get("cb2vsdNdeta");
  //Previous results
  double dndeta_pp7TeV[5] = {3.295484, 7.539474, 10.76316, 13.5, 17.47021};
  double dndeta_pp7TeV_err[5] = {0.1345161, 0.2289474, 0.3, 0.4, 0.5202128};
  double dndeta_pp7TeV_err_null[5] = {0., 0., 0., 0., 0.};
  double dndeta_pp7TeV_b2[5] = {0.01252897, 0.01074527, 0.00991707, 0.009443706, 0.009336069};
  double dndeta_pp7TeV_b2_err_stat[5] = {0.0005042418, 0.0003745772, 0.0003789578, 0.0004422433, 0.000362355};
  double dndeta_pp7TeV_b2_err_syst[5] = {0.001805432, 0.001365899, 0.00139078, 0.001299496, 0.001227568};
  //
  TGraphErrors *pp7TeVStat = new TGraphErrors(5, dndeta_pp7TeV, dndeta_pp7TeV_b2, dndeta_pp7TeV_err_null, dndeta_pp7TeV_b2_err_stat);
  Requires(pp7TeVStat,"pp7TeVStat");
  pp7TeVStat->SetName("pp7TeVStat");
	pp7TeVStat->SetFillStyle(0);
  pp7TeVStat->SetMarkerColor(TColor::GetColor("#006600"));
  pp7TeVStat->SetMarkerStyle(21);
  pp7TeVStat->SetMarkerSize(1.4);
  //
  TGraphErrors *pp7TeVSyst = new TGraphErrors(5, dndeta_pp7TeV, dndeta_pp7TeV_b2, dndeta_pp7TeV_err, dndeta_pp7TeV_b2_err_syst);
  pp7TeVSyst->SetName("pp7TeVSyst");
  pp7TeVSyst->SetMarkerColor(TColor::GetColor("#006600"));
  pp7TeVSyst->SetLineColor(TColor::GetColor("#006600"));
  pp7TeVSyst->SetMarkerStyle(21);
  pp7TeVSyst->SetMarkerSize(1.4);
  pp7TeVSyst->SetFillStyle(0);
  pp7TeVSyst->SetFillColor(0);
  //
  double dndeta_pPb5TeV[5] = {40.6, 30.5, 23.3, 16.1, 7.1};
  double dndeta_pPb5TeV_err[5] = {0.9, 0.7, 0.5, 0.4, 0.2};
  double dndeta_pPb5TeV_err_null[5] = {0., 0., 0., 0., 0.};
  double dndeta_pPb5TeV_b2[5] = {0.007025, 0.007429, 0.008189, 0.008699, 0.012181};
  double dndeta_pPb5TeV_b2_err_stat[5] = {0.000223, 0.000285, 0.000264, 0.000352, 0.000305};
  double dndeta_pPb5TeV_b2_err_syst[5] = {0.000718, 0.000985, 0.001091, 0.001123, 0.001399};
  //
  TGraphErrors *pPb5TeVStat = new TGraphErrors(5, dndeta_pPb5TeV, dndeta_pPb5TeV_b2, dndeta_pPb5TeV_err_null, dndeta_pPb5TeV_b2_err_stat);
  Requires(pPb5TeVStat,"pPb5TeVStat");
  pPb5TeVStat->SetName("pPb5TeVStat");
	pPb5TeVStat->SetFillStyle(0);
  pPb5TeVStat->SetMarkerColor(600);
  pPb5TeVStat->SetMarkerStyle(33);
  pPb5TeVStat->SetMarkerSize(1.4);
  //
  TGraphErrors *pPb5TeVSyst = new TGraphErrors(5, dndeta_pPb5TeV, dndeta_pPb5TeV_b2, dndeta_pPb5TeV_err, dndeta_pPb5TeV_b2_err_syst);
  Requires(pPb5TeVSyst,"pPb5TeVSyst");
  pPb5TeVSyst->SetName("pPb5TeVSyst");
	pPb5TeVSyst->SetMarkerColor(600);
  pPb5TeVSyst->SetLineColor(600);
  pPb5TeVSyst->SetMarkerStyle(33);
  pPb5TeVSyst->SetMarkerSize(1.4);
  pPb5TeVSyst->SetFillStyle(0);
  pPb5TeVSyst->SetFillColor(0);
  //
  //
  //   TGraphErrors* PbPb5TeVStat = (TGraphErrors*)cv_input->GetPrimitive("B2_PbPb15_pToA=0.750");
  //   Requires(PbPb5TeVStat,"PbPb5TeVStat");
  //   PbPb5TeVStat->SetName("PbPb5TeVStat");
  //   TGraphErrors* PbPb5TeVSyst = (TGraphErrors*)cv_input->GetPrimitive("B2_PbPb15_pToA=0.750_sys");
  //   Requires(PbPb5TeVSyst,"PbPb5TeVSyst");
  //   PbPb5TeVSyst->SetName("PbPb5TeVSyst");
  TGraphErrors* PbPb2d7TeVStat = (TGraphErrors*)cv_input->GetPrimitive("B2_PbPb10_pToA=0.750");
  Requires(PbPb2d7TeVStat,"PbPb2d7TeVStat");
  PbPb2d7TeVStat->SetName("PbPb2d7TeVStat");
	PbPb2d7TeVStat->SetFillStyle(0);
  TGraphAsymmErrors* PbPb2d7TeVSyst = (TGraphAsymmErrors*)cv_input->GetPrimitive("B2_PbPb10_pToA=0.750_sys");
  Requires(PbPb2d7TeVSyst,"PbPb2d7TeVSyst");
  PbPb2d7TeVSyst->SetName("PbPb2d7TeVSyst");
	PbPb2d7TeVSyst->SetFillStyle(0);
  //
  TGraphErrors* hB2_coalescenceParam0 = (TGraphErrors*)cv_input->GetPrimitive("hB2_coalescenceParam0");
  Requires(hB2_coalescenceParam0,"hB2_coalescenceParam0");
  hB2_coalescenceParam0->SetName("hB2_coalescenceParam0");
  //
  TGraphErrors* hB2_coalescenceParam1 = (TGraphErrors*)cv_input->GetPrimitive("hB2_coalescenceParam1");
  Requires(hB2_coalescenceParam1,"hB2_coalescenceParam1");
  hB2_coalescenceParam1->SetName("hB2_coalescenceParam1");
  //
  //myWork
  TFile* my_file = new TFile(Form("%sB2.root",kBaseOutputDir.data()));
  TGraphErrors* pp13TeVStat = (TGraphErrors*)my_file->Get("grB2atPtFixedStat_A_75");
  Requires(pp13TeVStat,"pp13TeVStat");
  pp13TeVStat->SetName("pp13TeVStat");
  pp13TeVStat->SetMarkerStyle(47);
  pp13TeVStat->SetMarkerSize(1.5);
  TGraphErrors* pp13TeVSyst = (TGraphErrors*)my_file->Get("grB2atPtFixedSyst_A_75");
  Requires(pp13TeVSyst,"pp13TeVSyst");
  pp13TeVSyst->SetName("pp13TeVSyst");
  pp13TeVSyst->SetMarkerStyle(47);
  pp13TeVSyst->SetMarkerSize(1.5);
  // MB
  TGraphErrors* pp13TeVStatMB = (TGraphErrors*)my_file->Get("grB2atPtFixedStatMB_A_75");
  Requires(pp13TeVStatMB,"pp13TeVStatMB");
  pp13TeVStatMB->SetName("pp13TeVStatMB");
  pp13TeVStatMB->SetMarkerStyle(47);
  pp13TeVStatMB->SetMarkerSize(1.5);
  TGraphErrors* pp13TeVSystMB = (TGraphErrors*)my_file->Get("grB2atPtFixedSystMB_A_75");
  Requires(pp13TeVSystMB,"pp13TeVSystMB");
  pp13TeVSystMB->SetName("pp13TeVSystMB");
  pp13TeVSystMB->SetMarkerStyle(47);
  pp13TeVSystMB->SetMarkerSize(1.5);



  TCanvas* cv_output = new TCanvas("cB2atPtFixed75_A","cB2atPtFixed75_A",800,600);
  cv_output->SetLogx();
  cv_output->SetLogy();
  cv_output->SetLeftMargin(0.14);
  cv_output->SetBottomMargin(0.14);
  TH1F* hFrame = new TH1F("hFrame",";#LTd#it{N}_{ch} / d#it{#eta}_{lab}#GT_{|#it{#eta}_{lab}| < 0.5}; #it{B}_{2} (GeV^{2}/#it{c}^{3})",3000,1,3001);
  hFrame->GetYaxis()->SetRangeUser(5.e-5,3e-2);
  hFrame->GetYaxis()->SetTitleOffset(1.1);
  hFrame->GetXaxis()->SetTitleOffset(1.2);
  hFrame->Draw();
  //
  hB2_coalescenceParam0->Draw("same");
	hB2_coalescenceParam1->Draw("same");
  //
  PbPb2d7TeVStat->Draw("samepz");
  PbPb2d7TeVSyst->Draw("e2same");
  pPb5TeVStat->Draw("samepz");
  pPb5TeVSyst->Draw("e2same");
  //PbPb5TeVStat->Draw("samepz");
  //PbPb5TeVSyst->Draw("e2same");
  pp7TeVStat->Draw("samepz");
  pp7TeVSyst->Draw("e2same");
  pp13TeVStat->Draw("samepz");
  pp13TeVSyst->Draw("e2same");
  // pp13TeVStatMB->Draw("samepz");
  // pp13TeVSystMB->Draw("e2same");

  TLegend* legPrel = new TLegend(0.192982,0.394783,0.302005,0.68,NULL,"brNDC");
  legPrel->SetBorderSize(0);
  legPrel->SetTextSize(0.03);
  legPrel->SetLineColor(1);
  legPrel->SetLineStyle(1);
  legPrel->SetLineWidth(1);
  legPrel->SetFillColor(0);
  legPrel->SetFillStyle(0);
  legPrel->SetHeader("#bf{ALICE}");
	//
  legPrel->AddEntry(pp13TeVSyst,"#bar{d}, pp, #sqrt{#it{s}} = 13 TeV","pf");
  //legPrel->AddEntry(pp13TeVSystMB,"#bar{d} (INEL), pp, #sqrt{#it{s}} = 13 TeV","pf");
  legPrel->AddEntry(pp7TeVSyst,"d+#bar{d}, pp, #sqrt{#it{s}} = 7 TeV", "pf");// (PLB 794 (2019) 50-63)","pf");
  //legPrel->AddEntry((TObject*)nullptr,"V0M Multiplicity Classes","");
  legPrel->AddEntry(pPb5TeVSyst,"d+#bar{d}, p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV", "pf");// (arXiv:1906.03136)","pf");
  //legPrel->AddEntry((TObject*)nullptr,"V0A Multiplicity Classes (Pb-side)","");
  //legPrel->AddEntry(PbPb5TeVSyst,"d, Pb-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV","pf");
  legPrel->AddEntry(PbPb2d7TeVSyst,"d, Pb-Pb, #sqrt{#it{s}_{NN}} = 2.76 TeV", "pf");// (PRC 93 (2015) 024917)","pf");
	TLegendEntry *header = (TLegendEntry*)legPrel->GetListOfPrimitives()->First();
	header->SetTextSize(0.04);
  legPrel->Draw();

	TLegend* legModel = new TLegend(0.196742,0.206957,0.305764,0.366957,NULL,"brNDC");
  legModel->SetBorderSize(0);
  legModel->SetTextSize(0.03);
  legModel->SetLineColor(1);
  legModel->SetLineStyle(1);
  legModel->SetLineWidth(1);
  legModel->SetFillColor(0);
  legModel->SetFillStyle(0);  
  legModel->SetHeader("#it{B}_{2} coalesc. #it{r}(d) = 3.2 fm (PRC 99 (2019) 054905)");
  legModel->AddEntry(hB2_coalescenceParam0,"Param. A (fit to HBT radii)","l");
  legModel->AddEntry(hB2_coalescenceParam1,"Param. B (constrained to ALICE Pb--Pb #it{B}_{2})","l");
  legModel->Draw();

  TLatex ptLatex;
  ptLatex.SetTextAlign(12);
  ptLatex.SetTextSize(0.04);
  ptLatex.DrawLatexNDC(0.65,0.78,"#bf{#it{p}_{T}/#it{A} = 0.75 GeV/#it{c}}");

  TLatex text;
  text.SetTextFont(63);
  text.SetTextSize(25);
  //text.DrawTextNDC(0.635,0.8,"ALICE Preliminary");

  TFile output_file(Form("%scustomB2PtFixed_def.root",kBaseOutputDir.data()),"recreate");
  cv_output->Write();
  cv_output->SaveAs(Form("%sB2PtFixed_def.pdf",kFiguresFolder.data()));
  hFrame->Write();
  pp7TeVStat->Write();
  pp7TeVSyst->Write();
//   PbPb5TeVStat->Write();
//   PbPb5TeVSyst->Write();
  PbPb2d7TeVStat->Write();
  PbPb2d7TeVSyst->Write();
  pPb5TeVStat->Write();
  pPb5TeVSyst->Write();
  pp13TeVStat->Write();
  pp13TeVSyst->Write();
  pp13TeVStatMB->Write();
  pp13TeVSystMB->Write();
}
