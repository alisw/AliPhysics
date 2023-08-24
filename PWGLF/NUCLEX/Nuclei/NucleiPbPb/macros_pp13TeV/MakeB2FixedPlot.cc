#include "src/Common.h"
#include "src/Utils.h"
using namespace utils;
#include <TFile.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLatex.h>

const char* file_name[3]={"B2vsCollsAll_pToA_75.root","B2vsCollsAll_pToA_55.root","B2vsCollsAll_pToA_105.root"};

const char* kFixedPtLabelMyFile[3] = {"75","55","105"};
const char* kFixedPtLabelOtherFile[3] = {"0.750","0.550","1.050"};
const float kPtLatex[3] = {0.75, 0.55, 1.05};
const char* kLabelPart[2] = {"d","#bar{d}"};

void MakeB2FixedPlot(int iFile=0, int iS=1){

  gStyle->SetOptStat(0);
  TFile* input_file = new TFile(Form("%s%s",kBaseOutputDir.data(),file_name[iFile]));
  TCanvas* cv_input = (TCanvas*)input_file->Get("c1_n5");
  //cv_input->GetListOfPrimitives()->ls();
  TPad* pad = (TPad*)cv_input->GetPrimitive("c1_n5_1");
  //Previous results
  TGraphErrors* pp7TeVStat = (TGraphErrors*)pad->GetPrimitive(Form("B2_pToA=%s",kFixedPtLabelOtherFile[iFile]));
  Requires(pp7TeVStat,"pp7TeVStat");
  pp7TeVStat->SetName("pp7TeVStat");
  TGraphErrors* pp7TeVSyst = (TGraphErrors*)pad->GetPrimitive(Form("B2_pToA=%s_sys",kFixedPtLabelOtherFile[iFile]));
  Requires(pp7TeVSyst,"pp7TeVSyst");
  pp7TeVStat->SetName("pp7TeVSyst");
  //
  TGraphErrors* pPb5TeVStat = (TGraphErrors*)pad->GetPrimitive(Form("B2_pPb_pToA=%s",kFixedPtLabelOtherFile[iFile]));
  Requires(pPb5TeVStat,"pPb5TeVStat");
  pPb5TeVStat->SetName("pPb5TeVStat");
  TGraphAsymmErrors* pPb5TeVSyst = (TGraphAsymmErrors*)pad->GetPrimitive(Form("B2_pPb_pToA=%s_sys",kFixedPtLabelOtherFile[iFile]));
  Requires(pPb5TeVSyst,"pPb5TeVSyst");
  pPb5TeVSyst->SetName("pPb5TeVSyst");
  //
  //
  TGraphErrors* PbPb5TeVStat = (TGraphErrors*)pad->GetPrimitive(Form("B2_PbPb15_pToA=%s",kFixedPtLabelOtherFile[iFile]));
  Requires(PbPb5TeVStat,"PbPb5TeVStat");
  PbPb5TeVStat->SetName("PbPb5TeVStat");
  TGraphErrors* PbPb5TeVSyst = (TGraphErrors*)pad->GetPrimitive(Form("B2_PbPb15_pToA=%s_sys",kFixedPtLabelOtherFile[iFile]));
  Requires(PbPb5TeVSyst,"PbPb5TeVSyst");
  PbPb5TeVSyst->SetName("PbPb5TeVSyst");
  //
  TGraphErrors* PbPb2d7TeVStat = (TGraphErrors*)pad->GetPrimitive(Form("B2_PbPb10_pToA=%s",kFixedPtLabelOtherFile[iFile]));
  Requires(PbPb2d7TeVStat,"PbPb2d7TeVStat");
  PbPb2d7TeVStat->SetName("PbPb2d7TeVStat");
  TGraphAsymmErrors* PbPb2d7TeVSyst = (TGraphAsymmErrors*)pad->GetPrimitive(Form("B2_PbPb10_pToA=%s_sys",kFixedPtLabelOtherFile[iFile]));
  Requires(PbPb2d7TeVSyst,"PbPb2d7TeVSyst");
  PbPb2d7TeVSyst->SetName("PbPb2d7TeVSyst");
  //myWork
  TFile* my_file = new TFile(Form("%sB2.root",kBaseOutputDir.data()));
  TGraphErrors* pp13TeVStat = (TGraphErrors*)my_file->Get(Form("grB2atPtFixedStat_%c_%s",kLetter[iS],kFixedPtLabelMyFile[iFile]));
  Requires(pp13TeVStat,"pp13TeVStat");
  pp13TeVStat->SetName("pp13TeVStat");
  pp13TeVStat->SetMarkerStyle(47);
  pp13TeVStat->SetMarkerSize(1.5);
  TGraphErrors* pp13TeVSyst = (TGraphErrors*)my_file->Get(Form("grB2atPtFixedSyst_%c_%s",kLetter[iS],kFixedPtLabelMyFile[iFile]));
  Requires(pp13TeVSyst,"pp13TeVSyst");
  pp13TeVSyst->SetName("pp13TeVSyst");
  pp13TeVSyst->SetMarkerStyle(47);
  pp13TeVSyst->SetMarkerSize(1.5);

  TCanvas* cv_output = new TCanvas(Form("cB2atPtFixed%s_%c",kFixedPtLabelMyFile[iFile],kLetter[iS]),Form("cB2atPt%s_%c",kFixedPtLabelMyFile[iFile],kLetter[iS]),800,600);
  cv_output->SetLogx();
  cv_output->SetLogy();
  cv_output->SetLeftMargin(0.14);
  cv_output->SetBottomMargin(0.14);
  TH1F* hFrame = new TH1F("hFrame",";#LTd#it{N}_{ch} / d#it{#eta}_{lab}#GT_{|#it{#eta}_{lab}| < 0.5}; #it{B}_{2} (GeV^{2}/#it{c}^{3})",3000,1,3001);
  hFrame->GetYaxis()->SetRangeUser(1.e-4,3e-2);
  hFrame->GetYaxis()->SetTitleOffset(1.1);
  hFrame->GetXaxis()->SetTitleOffset(1.2);
  hFrame->Draw();
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

  TLegend* legPrel = new TLegend(0.180451,0.206957,0.289474,0.657391,NULL,"brNDC");
  legPrel->SetBorderSize(0);
  legPrel->SetTextSize(0.04);
  legPrel->SetLineColor(1);
  legPrel->SetLineStyle(1);
  legPrel->SetLineWidth(1);
  legPrel->SetFillColor(0);
  legPrel->SetFillStyle(0);
  //legPrel->SetHeader("#bf{ALICE Preliminary}");

  legPrel->AddEntry(pp13TeVSyst,Form("%s, pp, #sqrt{#it{s}} = 13 TeV",kLabelPart[iS]),"pf");
  legPrel->AddEntry(pp7TeVSyst,"d+#bar{d}, pp, #sqrt{#it{s}} = 7 TeV","pf");
  legPrel->AddEntry((TObject*)nullptr,"V0M Multiplicity Classes","");
  legPrel->AddEntry(pPb5TeVSyst,"d+#bar{d}, p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV","pf");
  legPrel->AddEntry((TObject*)nullptr,"V0A Multiplicity Classes (Pb-side)","");
  //legPrel->AddEntry(PbPb5TeVSyst,"d, Pb-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV","pf");
  legPrel->AddEntry(PbPb2d7TeVSyst,"d, Pb-Pb, #sqrt{#it{s}_{NN}} = 2.76 TeV (PRC 93 (2015) 024917)","pf");
  legPrel->Draw();

  TLatex ptLatex;
  ptLatex.SetTextAlign(12);
  ptLatex.SetTextSize(0.04);
  ptLatex.DrawLatexNDC(0.65,0.75,Form("#bf{#it{p}_{T}/#it{A} = %.2f GeV/#it{c}}",kPtLatex[iFile]));

  TLatex text;
  text.SetTextFont(63);
  text.SetTextSize(25);
  //text.DrawTextNDC(0.635,0.8,"ALICE Preliminary");

  TFile output_file(Form("%scustomB2PtFixed%s_%c.root",kBaseOutputDir.data(),kFixedPtLabelMyFile[iFile],kLetter[iS]),"recreate");
  cv_output->Write();
  cv_output->SaveAs(Form("%sB2PtFixed%s_%c.pdf",kFiguresFolder.data(),kFixedPtLabelMyFile[iFile],kLetter[iS]));
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
}
