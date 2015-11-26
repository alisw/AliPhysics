#include "syst_mu-tkl.C"
#include "TLatex.h"
#include "TPad.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TString.h"
#include "TLine.h"


Int_t nc = 3;
//Float_t dphi[5] = {1,2,3,4,5};
//Float_t xerr[6] = {0.25,0.25,0.25,0.25,0.5,0.5};
const char* kSystemEnergy = "p-Pb #sqrt{s_{NN}} = 5.02 TeV";
Float_t fontSize = 0.05;
const Bool_t gSystfromfit=1;

Double_t binst_centers[100];
Double_t binst_centersTrk[100];
Double_t binst_centersTkl[100];
Double_t binst_widths[100];

TH1F* draw_frame(Double_t xmax, Double_t ymax){
  TH1F* fr = gPad->DrawFrame(0,0,xmax,ymax);
  fr->SetTitle(";#it{p}_{T} (GeV/#it{c});v_{2}{2PC,sub}");
  fr->GetXaxis()->SetTitleOffset(1.3);
  fr->GetYaxis()->SetTitleOffset(1.5);
  fr->GetYaxis()->SetLabelSize(fontSize);
  fr->GetXaxis()->SetLabelSize(fontSize);
  fr->GetXaxis()->SetTitleSize(fontSize);
  fr->GetYaxis()->SetTitleSize(fontSize);
  return fr;
}

void draw_muon_factorization(TGraphErrors** g, TString period, TString centMethod){
  TCanvas* c2 = new TCanvas(Form("cm_%s",period.Data()),Form("muon factorization %s",period.Data()),800,600);
  TH1F* fr2 = draw_frame(4,0.15);
  TLegend* l = new TLegend(0.4,0.7,0.98,0.98);
  l->SetFillColor(0);
  l->SetBorderSize(0);
  for (Int_t i=0;i<5;i++){
    g[i]->SetLineWidth(2);
    g[i]->DrawClone("psame");
    l->AddEntry(g[i],Form("LHC13%s, %s, #Delta#varphi = %.0f mrad",period.Data(),centMethod.Data(),bins_tkl[i]));
  }
  l->Draw("same");
  gPad->Print(Form("eps/v2_muons_vs_dphi_%s_%s.eps",period.Data(),centMethod.Data()));
  gPad->Print(Form("png/v2_muons_vs_dphi_%s_%s.png",period.Data(),centMethod.Data()));
}

void draw_tracklet_factorization(TGraphErrors** g, TString period, TString centMethod){
  TCanvas* c2 = new TCanvas(Form("ct_%s",period.Data()),Form("tracklet factorization %s",period.Data()),800,600);
  TH1F* fr4 = draw_frame(6,0.2);
  fr4->SetTitle(";d#varphi_{cut},trig (mrad);v_{2}{2PC,sub}");
  fr4->GetXaxis()->SetTitleOffset(1.3);
  fr4->GetYaxis()->SetTitleOffset(1.4);
  TLegend* l = new TLegend(0.2,0.7,0.98,0.98);
  l->SetFillColor(0);
  l->SetBorderSize(0);
  for (Int_t i=0;i<5;i++)   {
    g[i]->SetLineWidth(2);
    g[i]->DrawClone("psame");
    l->AddEntry(g[i],Form("LHC13%s, %s, d #varphi_{cut},assoc = %.0f mrad",period.Data(),centMethod.Data(),bins_tkl[i]));
  }
  l->Draw("same");
  gPad->Print(Form("eps/v2_tracklets_vs_dphi_%s_%s.eps",period.Data(),centMethod.Data()));
  gPad->Print(Form("png/v2_tracklets_vs_dphi_%s_%s.png",period.Data(),centMethod.Data()));
}

void compare_multiplicity_estimators(TGraphErrors* g1, TGraphErrors* g2, TGraphErrors* g3, TGraphErrors* g4, TGraphErrors* g5, TString period, Bool_t isTracklets=0, Bool_t isRatio=0){
  TCanvas* c1 = new TCanvas(Form("%s_%i",period.Data(),isTracklets),Form("%s_%i",period.Data(),isTracklets),800,600);
  TH1F* fr = draw_frame(isTracklets?6:4,isRatio?2.5:0.09);//0.15);
  if (isTracklets) fr->SetTitle(";#Delta#varphi (mrad);v_{2}{2PC,sub}");
  g1->SetMarkerColor(kBlack);
  g2->SetMarkerColor(kRed);
  g3->SetMarkerColor(kGreen);
  g4->SetMarkerColor(kMagenta);
  g5->SetMarkerColor(kBlue);
  g1->SetMarkerStyle(21);
  g2->SetMarkerStyle(22);
  g3->SetMarkerStyle(23);
  g4->SetMarkerStyle(24);
  g5->SetMarkerStyle(25);
  g1->SetLineColor(kBlack);
  g2->SetLineColor(kRed);
  g3->SetLineColor(kGreen);
  g4->SetLineColor(kMagenta);
  g5->SetLineColor(kBlue);
  g1->SetLineWidth(2);
  g2->SetLineWidth(2);
  g3->SetLineWidth(2);
  g4->SetLineWidth(2);
  g5->SetLineWidth(2);
  g4->RemovePoint(1);
  g4->Fit("pol0");
  g1->DrawClone("plsame");
  //  g2->DrawClone("plsame");
  //  g3->DrawClone("plsame");
  g4->Draw("plsame");
  g5->DrawClone("plsame");
  TLegend* l = new TLegend(0.5,0.15,0.98,0.4);
  l->SetFillColor(0);
  l->SetBorderSize(0);
  l->AddEntry(g1,Form("LHC13%s, %s",period.Data(),centMethods[3].Data()));
  //  l->AddEntry(g2,Form("LHC13%s, %s",period.Data(),centMethods[4].Data()));
  //  l->AddEntry(g3,Form("LHC13%s, %s",period.Data(),centMethods[5].Data()));
  l->AddEntry(g4,Form("LHC13%s, %s",period.Data(),centMethods[6].Data()));
  l->AddEntry(g5,Form("LHC13%s, %s",period.Data(),centMethods[7].Data()));
  l->Draw("same");
  gPad->Print(Form("eps/v2_%s_%s_DiffEstimators.eps",isTracklets ? "tkl_tkl" : "mu_tkl",period.Data()));
  gPad->Print(Form("png/v2_%s_%s_DiffEstimators.png",isTracklets ? "tkl_tkl" : "mu_tkl",period.Data()));
}

void compare_ampt(TGraphErrors* g1, TGraphErrors* g2, TGraphErrors* g3, TGraphErrors* g4, TGraphErrors* g5, TString period,Bool_t isRatio=0){
  TCanvas* c1 = new TCanvas(Form("c_ampt_%s_",period.Data()),Form("c_ampt_%s",period.Data()),800,600);
  TH1F* fr = draw_frame(4,isRatio?2.5:0.15);
  g1->SetMarkerColor(kBlack);
  g2->SetMarkerColor(kRed);
  g3->SetMarkerColor(kGreen);
  g4->SetMarkerColor(kMagenta);
  g5->SetMarkerColor(kBlue);
  g1->SetMarkerStyle(21);
  g2->SetMarkerStyle(22);
  g3->SetMarkerStyle(23);
  g4->SetMarkerStyle(24);
  g5->SetMarkerStyle(25);
  g1->SetLineColor(kBlack);
  g2->SetLineColor(kRed);
  g3->SetLineColor(kGreen);
  g4->SetLineColor(kMagenta);
  g5->SetLineColor(kBlue);
  g1->SetLineWidth(2);
  g2->SetLineWidth(2);
  g3->SetLineWidth(2);
  g4->SetLineWidth(2);
  g5->SetLineWidth(2);
  g1->DrawClone("plsame");
  g2->DrawClone("plsame");
  g3->DrawClone("plsame");
  g4->Draw("plsame");
  g5->DrawClone("plsame");
  TLegend* l = new TLegend(0.5,0.15,0.98,0.4);
  l->SetFillColor(0);
  l->SetBorderSize(0);
  l->AddEntry(g1,Form("%s, %s",period.Data(),"all"));
  l->AddEntry(g2,Form("%s, %s",period.Data(),"#pi"));
  l->AddEntry(g3,Form("%s, %s",period.Data(),"K"));
  l->AddEntry(g4,Form("%s, %s",period.Data(),"#mu"));
  l->AddEntry(g5,Form("%s, %s",period.Data(),"reco"));
  l->Draw("same");
  gPad->Print(Form("eps/v2_mu_tkl_%s_ampt.eps",period.Data()));
  gPad->Print(Form("png/v2_mu_tkl_%s_ampt.png",period.Data()));
}

void compare_eta(TGraphErrors* g0,TGraphErrors* g0s,TGraphErrors* g1,TGraphErrors* g1s, TGraphErrors* g2,TGraphErrors* g2s, TGraphErrors* g3,TGraphErrors* g3s, TString period,Bool_t isRatio=0){
  TCanvas* c1 = new TCanvas(Form("c_eta_%s_",period.Data()),Form("c_eta_%s",period.Data()),800,600);
  TH1F* fr = draw_frame(4,isRatio?2.5:0.15);
  g0->SetMarkerColor(kBlue);
  g1->SetMarkerColor(kBlack);
  g2->SetMarkerColor(kRed);
  g3->SetMarkerColor(kGreen);
  g0->SetLineColor(kBlack);
  g1->SetLineColor(kBlue);
  g2->SetLineColor(kRed);
  g3->SetLineColor(kGreen);
  g0->SetMarkerStyle(20);
  g1->SetMarkerStyle(21);
  g2->SetMarkerStyle(22);
  g3->SetMarkerStyle(23);
  g0->SetLineWidth(2);
  g1->SetLineWidth(2);
  g2->SetLineWidth(2);
  g3->SetLineWidth(2);
  g0s->SetMarkerColor(kBlue);
  g1s->SetMarkerColor(kBlack);
  g2s->SetMarkerColor(kRed);
  g3s->SetMarkerColor(kGreen);
  g0s->SetLineColor(kBlue);
  g1s->SetLineColor(kBlack);
  g2s->SetLineColor(kRed);
  g3s->SetLineColor(kGreen);
  g0s->SetMarkerStyle(20);
  g1s->SetMarkerStyle(21);
  g2s->SetMarkerStyle(22);
  g3s->SetMarkerStyle(23);
  g0s->SetLineWidth(2);
  g1s->SetLineWidth(2);
  g2s->SetLineWidth(2);
  g3s->SetLineWidth(2);
  g0s->SetFillColor(kBlue-8);
  g1s->SetFillColor(15);
  g2s->SetFillColor(kRed-8);
  g3s->SetFillColor(kGreen-8);
  g0s->SetFillStyle(1001);
  g1s->SetFillStyle(1001);
  g2s->SetFillStyle(1001);
  g3s->SetFillStyle(1001);
  g0s->DrawClone("2 lsame");
  g1s->DrawClone("2 same");
  g2s->DrawClone("2 same");
  g3s->DrawClone("2 same");
  g0->DrawClone("psame");
  g1->DrawClone("psame");
  g2->DrawClone("psame");
  g3->DrawClone("psame");
  TLegend* l = new TLegend(0.15,0.7,0.4,0.98);
  l->SetFillColor(0);
  l->SetBorderSize(0);
  l->AddEntry(g0s,Form("%s, %s",period.Data(),"-4. < #eta_{#mu} < -2.5"));
  l->AddEntry(g1s,Form("%s, %s",period.Data(),"-3. < #eta_{#mu} < -2.5"));
  l->AddEntry(g2s,Form("%s, %s",period.Data(),"-3.5 < #eta_{#mu} <-3."));
  l->AddEntry(g3s,Form("%s, %s",period.Data(),"-4. < #eta_{#mu} < -3.5"));
  l->Draw("same");
  gPad->Print(Form("eps/v2_mu_tkl_%s_eta.eps",period.Data()));
  gPad->Print(Form("png/v2_mu_tkl_%s_eta.png",period.Data()));
}


void compare_trk_tkl(TGraphErrors* g1, TGraphErrors* g2,TGraphErrors* g1syst, TGraphErrors* g2syst){
  TCanvas* c1 = new TCanvas("trk_tkl","trk_tkl",800,600);
  TH1F* fr = draw_frame(4,0.13);
  g1->SetMarkerStyle(kFullSquare);
  g2->SetMarkerStyle(kOpenCircle);
  g2->SetMarkerColor(kBlue);
  g1->SetMarkerColor(kRed);
  g2->SetLineColor(kBlue);
  g1->SetLineColor(kRed);
  g1->SetLineWidth(2);
  g2->SetLineWidth(2);
  g1syst->SetMarkerStyle(kFullSquare);
  g2syst->SetMarkerStyle(kOpenCircle);
  g2syst->SetMarkerColor(kBlue);
  g1syst->SetMarkerColor(kRed);
  g2syst->SetLineColor(kBlue);
  g1syst->SetLineColor(kRed);
  g1syst->SetLineWidth(2);
  g2syst->SetLineWidth(2);
  g1syst->SetFillColor(kRed-8);
  g2syst->SetFillColor(kBlue-8);
  g2syst->SetFillStyle(1001);
  g1syst->SetFillStyle(1001);
  g2syst->DrawClone("2 same");
  g1syst->DrawClone("2 same");
  g1->DrawClone("psame");
  g2->DrawClone("psame");
  TLegend* l = new TLegend(0.6,0.80,0.98,0.98);
  l->SetFillColor(0);
  l->SetBorderSize(0);
  l->AddEntry(g1syst,"muon-tracks");
  l->AddEntry(g2syst,"muon-tracklets");
  l->Draw("same");

  TLatex* latex = new TLatex();
  latex->SetTextSize(0.04);
  latex->SetTextFont(42);
  latex->SetTextAlign(21);
  latex->SetNDC();
  latex->DrawLatex(0.35,0.93,"ALICE");
  latex->DrawLatex(0.35,0.86,kSystemEnergy);
  latex->DrawLatex(0.35,0.79,"V0S mult: (0-20%) - (60-100%)");
  TFile* fout = new TFile("trk_tkl_final.root","recreate");
  g1->Write("g_trk");
  g2->Write("g_tkl");
  g1syst->Write("gsyst_trk");
  g2syst->Write("gsyst_tkl");
  fout->Close();

  gPad->Print("eps/trk_vs_tkl.eps");
  gPad->Print("png/trk_vs_tkl.png");
}

void compare_trk_bug(TGraphErrors* g1, TGraphErrors* g2,TGraphErrors* g1syst, TGraphErrors* g2syst){
  TCanvas* c1 = new TCanvas("trk_bug","trk_bug",800,600);
  TH1F* fr = draw_frame(4,0.13);
  g1->SetMarkerStyle(kFullSquare);
  g2->SetMarkerStyle(kOpenCircle);
  g2->SetMarkerColor(kBlue);
  g1->SetMarkerColor(kRed);
  g2->SetLineColor(kBlue);
  g1->SetLineColor(kRed);
  g1->SetLineWidth(2);
  g2->SetLineWidth(2);
  g1syst->SetMarkerStyle(kFullSquare);
  g2syst->SetMarkerStyle(kOpenCircle);
  g2syst->SetMarkerColor(kBlue);
  g1syst->SetMarkerColor(kRed);
  g2syst->SetLineColor(kBlue);
  g1syst->SetLineColor(kRed);
  g1syst->SetLineWidth(2);
  g2syst->SetLineWidth(2);
  g1syst->SetFillColor(kRed-8);
  g2syst->SetFillColor(kBlue-8);
  g2syst->SetFillStyle(1001);
  g1syst->SetFillStyle(1001);
  g2syst->DrawClone("2 same");
  g1syst->DrawClone("2 same");
  g1->DrawClone("psame");
  g2->DrawClone("psame");
  TLegend* l = new TLegend(0.6,0.80,0.98,0.98);
  l->SetFillColor(0);
  l->SetBorderSize(0);
  l->AddEntry(g1syst,"muon-tracks (LEE)");
  l->AddEntry(g2syst,"muon-tracks (LEW)");
  l->Draw("same");

  TLatex* latex = new TLatex();
  latex->SetTextSize(0.04);
  latex->SetTextFont(42);
  latex->SetTextAlign(21);
  latex->SetNDC();
  latex->DrawLatex(0.35,0.93,"ALICE");
  latex->DrawLatex(0.35,0.86,kSystemEnergy);
  latex->DrawLatex(0.35,0.79,"V0S mult: (0-20%) - (60-100%)");
  TFile* fout = new TFile("trk_bug_final.root","recreate");
  g1->Write("g_trk");
  g2->Write("g_bug");
  g1syst->Write("gsyst_trk");
  g2syst->Write("gsyst_bug");
  fout->Close();

  gPad->Print("eps/trk_vs_bug.eps");
  gPad->Print("png/trk_vs_bug.png");
}


void compare_effCorr_trk(TGraphErrors* g1, TGraphErrors* g2,TGraphErrors* g1syst, TGraphErrors* g2syst){
  TCanvas* c1 = new TCanvas("Efftrk","Efftrk",800,600);
  TH1F* fr = draw_frame(4,0.13);
  g1->SetMarkerStyle(kFullSquare);
  g2->SetMarkerStyle(kOpenCircle);
  g2->SetMarkerColor(kBlue);
  g1->SetMarkerColor(kRed);
  g2->SetLineColor(kBlue);
  g1->SetLineColor(kRed);
  g1->SetLineWidth(2);
  g2->SetLineWidth(2);
  g1syst->SetMarkerStyle(kFullSquare);
  g2syst->SetMarkerStyle(kOpenCircle);
  g2syst->SetMarkerColor(kBlue);
  g1syst->SetMarkerColor(kRed);
  g2syst->SetLineColor(kBlue);
  g1syst->SetLineColor(kRed);
  g1syst->SetLineWidth(2);
  g2syst->SetLineWidth(2);
  g1syst->SetFillColor(kRed-8);
  g2syst->SetFillColor(kBlue-8);
  g2syst->SetFillStyle(1001);
  g1syst->SetFillStyle(1001);
  //  g2syst->DrawClone("2 same");
  //  g1syst->DrawClone("2 same");
  g1->DrawClone("psame");
  g2->DrawClone("psame");
  TLegend* l = new TLegend(0.6,0.80,0.98,0.98);
  l->SetFillColor(0);
  l->SetBorderSize(0);
  l->AddEntry(g1syst,"muon-tracks with efficiency correction");
  l->AddEntry(g2syst,"muon-tracks w/o efficiency correction");
  l->Draw("same");

  TLatex* latex = new TLatex();
  latex->SetTextSize(0.04);
  latex->SetTextFont(42);
  latex->SetTextAlign(21);
  latex->SetNDC();
  latex->DrawLatex(0.35,0.93,"ALICE");
  latex->DrawLatex(0.35,0.86,kSystemEnergy);
  latex->DrawLatex(0.35,0.79,"V0S mult: (0-20%) - (60-100%)");


  gPad->Print("eps/Efftrk.eps");
  gPad->Print("png/Efftrk.png");
}
void compare_effCorr_tkl(TGraphErrors* g1, TGraphErrors* g2,TGraphErrors* g1syst, TGraphErrors* g2syst,TString period){
  TCanvas* c1 = new TCanvas(Form("Efftkl_%s",period.Data()),Form("Efftkl_%s",period.Data()),800,600);
  TH1F* fr = draw_frame(4,0.13);
  if(period.Contains("ratio")){
    fr->SetMinimum(0.9);
    fr->SetMaximum(2.5);
  }
  g1->Print();
  g1syst->Print();
  g1->SetMarkerStyle(kFullSquare);
  g2->SetMarkerStyle(kOpenCircle);
  g2->SetMarkerColor(kBlue);
  g1->SetMarkerColor(kRed);
  g2->SetLineColor(kBlue);
  g1->SetLineColor(kRed);
  g1->SetLineWidth(2);
  g2->SetLineWidth(2);
  g1syst->SetMarkerStyle(kFullSquare);
  g2syst->SetMarkerStyle(kOpenCircle);
  g2syst->SetMarkerColor(kBlue);
  g1syst->SetMarkerColor(kRed);
  g2syst->SetLineColor(kBlue);
  g1syst->SetLineColor(kRed);
  g1syst->SetLineWidth(2);
  g2syst->SetLineWidth(2);
  g1syst->SetFillColor(kRed-8);
  g2syst->SetFillColor(kBlue-8);
  g2syst->SetFillStyle(1001);
  g1syst->SetFillStyle(1001);
  g2syst->DrawClone("2 same");
  g1syst->DrawClone("2 same");
  g1->DrawClone("psame");
  g2->DrawClone("psame");
  TLegend* l = new TLegend(0.6,0.80,0.98,0.98);
  l->SetFillColor(0);
  l->SetBorderSize(0);
  l->AddEntry(g1syst,Form("muon-tracklet w/o efficiency correction - %s",period.Data()));
  l->AddEntry(g2syst,Form("muon-tracklet with efficiency correction - %s",period.Data()));
  l->Draw("same");

  TLatex* latex = new TLatex();
  latex->SetTextSize(0.04);
  latex->SetTextFont(42);
  latex->SetTextAlign(21);
  latex->SetNDC();
  latex->DrawLatex(0.35,0.93,"ALICE");
  latex->DrawLatex(0.35,0.86,kSystemEnergy);
  latex->DrawLatex(0.35,0.79,"V0S mult: (0-20%) - (60-100%)");


  gPad->Print(Form("eps/Efftkl_%s.eps",period.Data()));
  gPad->Print(Form("png/Efftkl_%s.png",period.Data()));
}

void draw_mu(TGraphErrors* g1,TGraphErrors* g1syst,TGraphErrors* g2,TGraphErrors* g2syst, TString method="V0S", TString opt=""){
  TCanvas* c1 = new TCanvas(Form("mu_tkl_%s",opt.Data()),Form("mu_tkl_%s",opt.Data()),800,600);
  TH1F* fr = draw_frame(4,0.125);
  g1->SetMarkerStyle(kFullSquare);
  g2->SetMarkerStyle(kOpenCircle);
  g1->SetMarkerColor(kRed);
  g2->SetMarkerColor(kBlue);
  g1->SetLineColor(kRed);
  g2->SetLineColor(kBlue);
  g1->SetLineWidth(2);
  g2->SetLineWidth(2);
  g1syst->SetMarkerStyle(kFullSquare);
  g2syst->SetMarkerStyle(kOpenCircle);
  g1syst->SetMarkerColor(kRed);
  g2syst->SetMarkerColor(kBlue);
  g1syst->SetLineColor(kRed);
  g2syst->SetLineColor(kBlue);
  g1syst->SetLineWidth(2);
  g2syst->SetLineWidth(2);
  g1syst->SetFillColor(kRed);
  g2syst->SetFillColor(kBlue);
  g1syst->SetFillStyle(3004);
  g2syst->SetFillStyle(3005);
  g1syst->DrawClone("2 same");
  g1->DrawClone("p same");
  g2syst->DrawClone("2 same");
  g2->DrawClone("p same");
  TLegend* l = new TLegend(0.7,0.84,0.98,0.98);
  l->SetFillColor(0);
  l->SetBorderSize(0);
  l->AddEntry(g2syst,"Pb-going");
  l->AddEntry(g1syst,"p-going");
  l->Draw("same");
  TLatex* latex = new TLatex();
  latex->SetTextSize(fontSize);
  latex->SetTextFont(42);
  latex->SetTextAlign(21);
  latex->SetNDC();
  latex->DrawLatex(0.40,0.93,"ALICE");
  latex->DrawLatex(0.40,0.86,kSystemEnergy);
  latex->DrawLatex(0.40,0.79,Form("%s mult: (0-20%) - (60-100%)",method.Data()));
  TFile* fout=new TFile(Form("v2final_mu_tkl%s.root",opt.Data()),"recreate");
  Printf("\nv2 bcde");
  g1->Print("all");
  Printf("\nv2 f");
  g2->Print("all");
  g1syst->Write("gsyst_b");
  g2syst->Write("gsyst_f");
  g1->Write("g_b");
  g2->Write("g_f");
  fout->Close();

  gPad->Print("eps/v2_mu_tkl_bcde_vs_f_with_syst.eps");
  gPad->Print("png/v2_mu_tkl_bcde_vs_f_with_syst.png");
}

void draw_tkl(TGraphErrors* g1,TGraphErrors* g1syst,TGraphErrors* g2,TGraphErrors* g2syst, TString method="V0S"){
  TCanvas* c1 = new TCanvas("tkl_tkl","tkl_tkl",800,600);
  TH1F* fr = draw_frame(6,0.15);
  fr->SetTitle(";#Delta#varphi (mrad);v_{2}{2PC,sub}");
  g1->SetMarkerStyle(kFullSquare);
  g2->SetMarkerStyle(kOpenCircle);
  g1->SetMarkerColor(kRed);
  g2->SetMarkerColor(kBlue);
  g1->SetLineColor(kRed);
  g2->SetLineColor(kBlue);
  g1->SetLineWidth(2);
  g2->SetLineWidth(2);
  if (g1syst) {
    g1syst->SetMarkerStyle(kFullSquare);
    g1syst->SetMarkerColor(kRed);
    g1syst->SetLineColor(kRed);
    g1syst->SetLineWidth(2);
    g1syst->SetFillColor(kRed);
    g1syst->SetFillStyle(3004);
    g1syst->DrawClone("2 same");
  }
  if (g2syst) {
    g2syst->SetMarkerStyle(kOpenCircle);
    g2syst->SetMarkerColor(kBlue);
    g2syst->SetLineColor(kBlue);
    g2syst->SetLineWidth(2);
    g2syst->SetFillColor(kBlue);
    g2syst->SetFillStyle(3005);
    g2syst->DrawClone("2 same");
  }
  g1->DrawClone("p same");
  g2->DrawClone("p same");
  TLegend* l = new TLegend(0.7,0.84,0.98,0.98);
  l->SetFillColor(0);
  l->SetBorderSize(0);
  l->AddEntry(g2syst ? g2syst : g2,"LHC13f");
  l->AddEntry(g1syst ? g1syst : g1,"LHC13bcde");
  l->Draw("same");
  TLatex* latex = new TLatex();
  latex->SetTextSize(fontSize);
  latex->SetTextFont(42);
  latex->SetTextAlign(21);
  latex->SetNDC();
  latex->DrawLatex(0.40,0.93,"ALICE");
  latex->DrawLatex(0.40,0.86,kSystemEnergy);
  latex->DrawLatex(0.40,0.79,Form("%s mult: (0-20%) - (60-100%)",method.Data()));

  gPad->Print("eps/v2_tkl_tkl_bcde_vs_f_with_syst.eps");
  gPad->Print("png/v2_tkl_tkl_bcde_vs_f_with_syst.png");
}


void draw_ratio(TGraphErrors* g1, TGraphErrors* g1syst, TString method="V0S", TString opt=""){
  TCanvas* c1 = new TCanvas(Form("ratio%s",opt.Data()),Form("ratio%s",opt.Data()),800,600);
  TH1F* fr = draw_frame(4,1.5);
  fr->SetMinimum(0.9);
  fr->SetMaximum(2.5);
  fr->GetYaxis()->SetTitle("v_{2}(Pb-going) / v_{2}(p-going)");
  TLine* l = new TLine(0,1,4,1);
  l->SetLineStyle(2);
  l->Draw();
  g1syst->SetMarkerStyle(kFullSquare);
  g1syst->SetMarkerColor(kBlue);
  g1syst->SetLineColor(kBlue);
  g1syst->SetLineWidth(2);
  g1syst->SetFillColor(kBlue);
  g1syst->SetFillStyle(3004);
  g1syst->Draw("2 same");

  g1->SetMarkerStyle(kFullSquare);
  g1->SetMarkerColor(kBlue);
  g1->SetLineColor(kBlue);
  g1->SetLineWidth(2);
  g1->Draw("p same");
  Printf("\nv2 ratio");
  g1->Print("all");

  TLatex* latex = new TLatex();
  latex->SetTextSize(fontSize);
  latex->SetTextFont(42);
  latex->SetTextAlign(21);
  latex->SetNDC();
  latex->DrawLatex(0.40,0.93,"ALICE");
  latex->DrawLatex(0.40,0.86,kSystemEnergy);
  latex->DrawLatex(0.40,0.79,Form("%s mult: (0-20%) - (60-100%)",method.Data()));
  TFile* fout=new TFile(Form("ratio%s.root",opt.Data()),"recreate");
  g1syst->Write("gsyst_ratio");
  g1->Write("g_ratio");
  fout->Close();
  gPad->Print("eps/v2_ratio.eps");
  gPad->Print("png/v2_ratio.png");

}


void draw_diff(TGraphErrors* g1, TGraphErrors* g1syst, TString method="V0S"){
  TCanvas* c1 = new TCanvas("diff","diff",800,600);
  TH1F* fr = draw_frame(4,0.04);
  //  fr->SetMinimum(-);
  //  fr->SetMaximum(2.5);
  fr->GetYaxis()->SetTitle("v_{2}(Pb-going) - v_{2}(p-going)");
  TLine* l = new TLine(0,1,4,1);
  l->SetLineStyle(2);
  l->Draw();
  g1syst->SetMarkerStyle(kFullSquare);
  g1syst->SetMarkerColor(kBlue);
  g1syst->SetLineColor(kBlue);
  g1syst->SetLineWidth(2);
  g1syst->SetFillColor(kBlue);
  g1syst->SetFillStyle(3004);
  g1syst->Draw("2 same");

  g1->SetMarkerStyle(kFullSquare);
  g1->SetMarkerColor(kBlue);
  g1->SetLineColor(kBlue);
  g1->SetLineWidth(2);
  g1->Draw("p same");
  Printf("\nv2 diff");
  g1->Print("all");

  TLatex* latex = new TLatex();
  latex->SetTextSize(fontSize);
  latex->SetTextFont(42);
  latex->SetTextAlign(21);
  latex->SetNDC();
  latex->DrawLatex(0.40,0.93,"ALICE");
  latex->DrawLatex(0.40,0.86,kSystemEnergy);
  latex->DrawLatex(0.40,0.79,Form("%s mult: (0-20%) - (60-100%)",method.Data()));
  TFile* fout=new TFile("diff.root","recreate");
  g1syst->Write("gsyst_diff");
  g1->Write("g_diff");
  fout->Close();
  gPad->Print("eps/v2_diff.eps");
  gPad->Print("png/v2_diff.png");

}


void compare_multiplicity_class(TGraphErrors* g1,TGraphErrors* g2,TGraphErrors* g3,TGraphErrors* g1syst,TGraphErrors* g2syst,TGraphErrors* g3syst,TString opt="",Bool_t isRatio=0){
  TCanvas* c1 = new TCanvas(Form("mu_tkl_mult_class_%s",opt.Data()),Form("mu_tkl_mult_class_%s",opt.Data()),800,600);
  TH1F* fr = draw_frame(4,isRatio?2.5:0.15);
  g1->SetMarkerStyle(kFullSquare);
  g2->SetMarkerStyle(kFullCircle);
  g3->SetMarkerStyle(kOpenCircle);
  g1->SetMarkerColor(kRed);
  g2->SetMarkerColor(kBlue);
  g3->SetMarkerColor(kGreen+2);
  g1->SetLineColor(kRed);
  g2->SetLineColor(kBlue);
  g3->SetLineColor(kGreen+2);
  g1->SetLineWidth(2);
  g2->SetLineWidth(2);
  g3->SetLineWidth(2);
  g1syst->SetMarkerStyle(kFullSquare);
  g2syst->SetMarkerStyle(kFullCircle);
  g3syst->SetMarkerStyle(kOpenCircle);
  g1syst->SetMarkerColor(kRed);
  g2syst->SetMarkerColor(kBlue);
  g3syst->SetMarkerColor(kGreen+2);
  g1syst->SetLineColor(kRed);
  g2syst->SetLineColor(kBlue);
  g3syst->SetLineColor(kGreen+2);
  g1syst->SetLineWidth(2);
  g2syst->SetLineWidth(2);
  g3syst->SetLineWidth(2);
  g1syst->SetFillColor(kRed-8);
  g2syst->SetFillColor(kBlue-8);
  g3syst->SetFillColor(kGreen-6);
  g1syst->SetFillStyle(1001);
  g2syst->SetFillStyle(1001);
  g3syst->SetFillStyle(1001);
  g3syst->DrawClone("2 same");
  g2syst->DrawClone("2 same");
  g1syst->DrawClone("2 same");
  g3->DrawClone("p same");
  g2->DrawClone("p same");
  g1->DrawClone("p same");
  TLegend* l = new TLegend(0.7,0.84,0.98,0.98);
  l->SetFillColor(0);
  l->SetBorderSize(0);
  l->AddEntry(g1syst,Form("0-20 %s",opt.Data()));
  l->AddEntry(g2syst,Form("20-40 %s",opt.Data()));
  l->AddEntry(g3syst,Form("40-60 %s",opt.Data()));
  l->Draw("same");
  TLatex* latex = new TLatex();
  latex->SetTextSize(fontSize);
  latex->SetTextFont(42);
  latex->SetTextAlign(21);
  latex->SetNDC();
  latex->DrawLatex(0.40,0.93,"ALICE");
  latex->DrawLatex(0.40,0.86,kSystemEnergy);
  gPad->Print(Form("eps/v2_mu_tkl_mult_class_%s.eps",opt.Data()));
  gPad->Print(Form("png/v2_mu_tkl_mult_class_%s.png",opt.Data()));
}

//void draw_vs_mult(){
//  Int_t fullMarkers[]={20,21,22,23,27};
//  Int_t openMarkers[]={24,25,26,32,33};
//
//  gStyle->SetOptTitle(0);
//  TFile* f1[5];
//  TFile* f2[5];
//  TFile* f3[5];
//  TFile* f4[5];
//  TGraphErrors* vg1[5][3]; // mu_tkl_b
//  TGraphErrors* vg2[5][3]; // mu_tkl_f
//  TGraphErrors* vg3[5][3]; // tkl_tkl_b
//  TGraphErrors* vg4[5][3]; // tkl_tkl_f
//  TH1D* hm1[5][3]; // mu_tkl_b
//  TH1D* hm2[5][3]; // mu_tkl_f
//  TH1D* hm3[5][3]; // tkl_tkl_b
//  TH1D* hm4[5][3]; // tkl_tkl_f
//
//  for (Int_t im=0;im<5;im++){
//    f1[im] = new TFile(Form("v2_mu_tkl_%s_bcde.root" ,centMethods[im].Data()));
//    f2[im] = new TFile(Form("v2_mu_tkl_%s_f.root"    ,centMethods[im].Data()));
//    f3[im] = new TFile(Form("v2_tkl_tkl_%s_bcde.root",centMethods[im].Data()));
//    f4[im] = new TFile(Form("v2_tkl_tkl_%s_f.root"   ,centMethods[im].Data()));
//
//    for(Int_t ic=0;ic<nc;ic++){
//      // TODO 4->0 for tracks
//      vg1[im][ic]=(TGraphErrors*)f1[im]->Get(Form("gv2_cent%d_ass4",ic));
//      vg2[im][ic]=(TGraphErrors*)f2[im]->Get(Form("gv2_cent%d_ass4",ic));
//      vg3[im][ic]=(TGraphErrors*)f3[im]->Get(Form("gv2_cent%d_ass4",ic));
//      vg4[im][ic]=(TGraphErrors*)f4[im]->Get(Form("gv2_cent%d_ass4",ic));
//      hm1[im][ic]=(TH1D*)f1[im]->Get(Form("TriggerStat_%d",ic));
//      hm2[im][ic]=(TH1D*)f2[im]->Get(Form("TriggerStat_%d",ic));
//      hm3[im][ic]=(TH1D*)f3[im]->Get(Form("TriggerStat_%d",ic));
//      hm4[im][ic]=(TH1D*)f4[im]->Get(Form("TriggerStat_%d",ic));
//    }
//  }
//
//  TGraphErrors* vgm1[nc][6]; // mu_tkl_b
//  TGraphErrors* vgm2[nc][6]; // mu_tkl_f
//  TGraphErrors* vgm3[nc][5]; // tkl_tkl_b
//  TGraphErrors* vgm4[nc][5]; // tkl_tkl_f
//
//  for(Int_t ic=0;ic<nc;ic++){
//    for(Int_t ip=0;ip<6;ip++){
//      vgm1[ic][ip]=new TGraphErrors(5);
//      vgm2[ic][ip]=new TGraphErrors(5);
//      vgm3[ic][ip]=new TGraphErrors(5);
//      vgm4[ic][ip]=new TGraphErrors(5);
//      vgm1[ic][ip]->SetTitle(Form("mu_tkl_b_Cent%d_pt%d",ic,ip));
//      vgm2[ic][ip]->SetTitle(Form("mu_tkl_f_Cent%d_pt%d",ic,ip));
//      vgm3[ic][ip]->SetTitle(Form("tkl_tkl_b_Cent%d_dphi%d",ic,ip+1));
//      vgm4[ic][ip]->SetTitle(Form("tkl_tkl_f_Cent%d_dphi%d",ic,ip+1));
//      for (Int_t im=0;im<5;im++){
//        if(0){
//          vgm1[ic][ip]->SetPoint(im,hm1[im][ic]->Integral(),vg1[im][ic]->GetY()[ip]);//total mult of mu
//          vgm2[ic][ip]->SetPoint(im,hm2[im][ic]->Integral(),vg2[im][ic]->GetY()[ip]);
//          vgm1[ic][ip]->GetXaxis()->SetTitle(Form("<#mu>, cent class #%.1d",ic));
//          vgm2[ic][ip]->GetXaxis()->SetTitle(Form("<#mu>, cent class #%.1d",ic));
//        }else{
//          vgm1[ic][ip]->SetPoint(im,hm3[im][ic]->Integral(),vg1[im][ic]->GetY()[ip]);//tkl mult with dphi=5
//          vgm2[ic][ip]->SetPoint(im,hm4[im][ic]->Integral(),vg2[im][ic]->GetY()[ip]);
//          vgm1[ic][ip]->GetXaxis()->SetTitle(Form("<tkl> (d#varphi_{cut}=5 mrad), cent class #%.1d",ic));
//          vgm2[ic][ip]->GetXaxis()->SetTitle(Form("<tkl> (d#varphi_{cut}=5 mrad), cent class #%.1d",ic));
//        }
//        vgm1[ic][ip]->SetPointError(im,0,vg1[im][ic]->GetEY()[ip]);
//        vgm1[ic][ip]->SetMarkerColor(ip+1);
//        vgm1[ic][ip]->SetLineColor(ip+1);
//        vgm1[ic][ip]->SetMarkerStyle(20+ip);
//        vgm2[ic][ip]->SetPointError(im,0,vg2[im][ic]->GetEY()[ip]);
//        vgm2[ic][ip]->SetMarkerColor(ip+1);
//        vgm2[ic][ip]->SetLineColor(ip+1);
//        vgm2[ic][ip]->SetMarkerStyle(20+ip);
//
//        vgm3[ic][ip]->SetPoint(im,hm3[im][ic]->Integral(),vg3[im][ic]->GetY()[ip]);
//        vgm3[ic][ip]->SetPointError(im,0,vg3[im][ic]->GetEY()[ip]);
//        vgm3[ic][ip]->SetMarkerColor(ip+1);
//        vgm3[ic][ip]->SetLineColor(ip+1);
//        vgm3[ic][ip]->SetMarkerStyle(20+ip);
//        vgm4[ic][ip]->SetPoint(im,hm4[im][ic]->Integral(),vg4[im][ic]->GetY()[ip]);
//        vgm4[ic][ip]->SetPointError(im,0,vg4[im][ic]->GetEY()[ip]);
//        vgm4[ic][ip]->SetMarkerColor(ip+1);
//        vgm4[ic][ip]->SetLineColor(ip+1);
//        vgm4[ic][ip]->SetMarkerStyle(20+ip);
//      }
//      vgm3[ic][ip]->GetXaxis()->SetTitle(Form("<tkl> (d#varphi_{cut}=5 mrad), cent class #%.1d",ic));
//      vgm4[ic][ip]->GetXaxis()->SetTitle(Form("<tkl> (d#varphi_{cut}=5 mrad), cent class #%.1d",ic));
//      vgm1[ic][ip]->GetYaxis()->SetTitle(Form("v2_{#mu}"));
//      vgm2[ic][ip]->GetYaxis()->SetTitle(Form("v2_{#mu}"));
//      vgm3[ic][ip]->GetYaxis()->SetTitle(Form("v2_{#tkl}"));
//      vgm4[ic][ip]->GetYaxis()->SetTitle(Form("v2_{#tkl}"));
//      vgm1[ic][ip]->SetFillColor(0);
//      vgm2[ic][ip]->SetFillColor(0);
//      vgm3[ic][ip]->SetFillColor(0);
//      vgm4[ic][ip]->SetFillColor(0);
//      vgm1[ic][ip]->GetYaxis()->SetTitleOffset(1.4);
//      vgm2[ic][ip]->GetYaxis()->SetTitleOffset(1.4);
//      vgm3[ic][ip]->GetYaxis()->SetTitleOffset(1.4);
//      vgm4[ic][ip]->GetYaxis()->SetTitleOffset(1.4);
//    }
//  }
//
//  for(Int_t ic=0;ic<1;ic++){
//    TCanvas*c=new TCanvas(Form("mu_tkl_Cent%d",ic),Form("mu_tkl_Cent%d",ic));
//    c->Divide(2,1);
//    for(Int_t ip=0;ip<6;ip++){
//      c->cd(1);
//      vgm1[ic][ip]->GetYaxis()->SetRangeUser(0,.15);
//      vgm1[ic][ip]->Draw((ip==0)?"ap":"psame");
//      c->cd(2);
//      vgm2[ic][ip]->GetYaxis()->SetRangeUser(0,.15);
//      vgm2[ic][ip]->Draw((ip==0)?"ap":"psame");
//    }
//    c->cd(1);
//    gPad->BuildLegend()->SetFillColor(0);
//    c->cd(2);
//    gPad->BuildLegend()->SetFillColor(0);
//  }
//  for(Int_t ic=0;ic<1;ic++){
//    TCanvas*c=new TCanvas(Form("tkl_tkl_Cent%d",ic),Form("tkl_tkl_Cent%d",ic));
//    c->Divide(2,1);
//    for(Int_t ip=0;ip<5;ip++){
//      c->cd(1);
//      vgm3[ic][ip]->GetYaxis()->SetRangeUser(0,.1);
//      vgm3[ic][ip]->Draw((ip==0)?"ap":"psame");
//      c->cd(2);
//      vgm4[ic][ip]->GetYaxis()->SetRangeUser(0,.1);
//      vgm4[ic][ip]->Draw((ip==0)?"ap":"psame");
//    }
//    c->cd(1);
//    gPad->BuildLegend()->SetFillColor(0);
//    c->cd(2);
//    gPad->BuildLegend()->SetFillColor(0);
//  }
//
//  TGraphErrors * hmInt1[5][3]; //integrated mult
//  TGraphErrors * hmInt2[5][3]; //integrated mult
//  TGraphErrors * hmInt3[5][3]; //integrated mult
//  TGraphErrors * hmInt4[5][3]; //integrated mult
//
//  for(Int_t ic=0;ic<1;ic++){
//    TCanvas *cmult=new TCanvas(Form("cmult_Cent%d",ic),Form("cmult_Cent%d",ic),1200,800);
//    cmult->Divide(2,1);
//    TCanvas *cmultInt=new TCanvas(Form("cmultInt_Cent%d",ic),Form("cmultInt_Cent%d",ic),1200,800);
//    cmultInt->Divide(2,1);
//    for (Int_t im=0;im<5;im++){
//      hmInt1[im][ic] = new TGraphErrors(1);
//      hmInt2[im][ic] = new TGraphErrors(1);
//      hmInt3[im][ic] = new TGraphErrors(1);
//      hmInt4[im][ic] = new TGraphErrors(1);
//      hmInt1[im][ic]->GetYaxis()->SetTitle(Form("<mu>, cent class #%.1d",ic));
//      hmInt2[im][ic]->GetYaxis()->SetTitle(Form("<mu>, cent class #%.1d",ic));
//      hmInt3[im][ic]->GetYaxis()->SetTitle(Form("<tkl> (d#varphi_{cut}=5 mrad) , cent class #%.1d",ic));
//      hmInt4[im][ic]->GetYaxis()->SetTitle(Form("<tkl> (d#varphi_{cut}=5 mrad) , cent class #%.1d",ic));
//      hmInt1[im][ic]->SetTitle(Form("mu_tkl_b_%s",centMethods[im].Data()));
//      hmInt2[im][ic]->SetTitle(Form("mu_tkl_f_%s",centMethods[im].Data()));
//      hmInt3[im][ic]->SetTitle(Form("tkl_tkl_b_%s",centMethods[im].Data()));
//      hmInt4[im][ic]->SetTitle(Form("tkl_tkl_f_%s",centMethods[im].Data()));
//      hmInt1[im][ic]->SetFillStyle(0);
//      hmInt2[im][ic]->SetFillStyle(0);
//      hmInt3[im][ic]->SetFillStyle(0);
//      hmInt4[im][ic]->SetFillStyle(0);
//      hmInt1[im][ic]->SetMarkerSize(2);
//      hmInt2[im][ic]->SetMarkerSize(2);
//      hmInt3[im][ic]->SetMarkerSize(2);
//      hmInt4[im][ic]->SetMarkerSize(2);
//      hmInt1[im][ic]->SetMarkerStyle(fullMarkers[im]);
//      hmInt2[im][ic]->SetMarkerStyle(openMarkers[im]);
//      hmInt3[im][ic]->SetMarkerStyle(fullMarkers[im]);
//      hmInt4[im][ic]->SetMarkerStyle(openMarkers[im]);
//      hmInt1[im][ic]->SetPoint(0,1.,hm1[im][ic]->Integral());
//      hmInt2[im][ic]->SetPoint(0,1.,hm2[im][ic]->Integral());
//      hmInt3[im][ic]->SetPoint(0,1.,hm3[im][ic]->Integral());
//      hmInt4[im][ic]->SetPoint(0,1.,hm4[im][ic]->Integral());
//
//      hm1[im][ic]->SetYTitle(Form("<mu>, cent class #%.1d",ic));
//      hm2[im][ic]->SetYTitle(Form("<mu>, cent class #%.1d",ic));
//      hm3[im][ic]->SetYTitle(Form("<tkl>, cent class #%.1d",ic));
//      hm4[im][ic]->SetYTitle(Form("<tkl>, cent class #%.1d",ic));
//      hm1[im][ic]->SetXTitle(Form("pT"));
//      hm2[im][ic]->SetXTitle(Form("pT"));
//      hm1[im][ic]->GetYaxis()->SetTitleOffset(1.6);
//      hm2[im][ic]->GetYaxis()->SetTitleOffset(1.6);
//      hm3[im][ic]->SetXTitle(Form("d#varphi value"));
//      hm4[im][ic]->SetXTitle(Form("d#varphi value"));
//      hm1[im][ic]->SetTitle(Form("mu_tkl_b_%s",centMethods[im].Data()));
//      hm2[im][ic]->SetTitle(Form("mu_tkl_f_%s",centMethods[im].Data()));
//      hm3[im][ic]->SetTitle(Form("tkl_tkl_b_%s",centMethods[im].Data()));
//      hm4[im][ic]->SetTitle(Form("tkl_tkl_f_%s",centMethods[im].Data()));
//      hm1[im][ic]->SetMarkerStyle(20+im);
//      hm2[im][ic]->SetMarkerStyle(24+im);
//      hm3[im][ic]->SetMarkerStyle(20+im);
//      hm4[im][ic]->SetMarkerStyle(24+im);
//      cmultInt->cd(1);
//      hmInt1[im][ic]->Draw((im==0)?"ap":"psame");
//      hmInt2[im][ic]->Draw("psame");
//      cmultInt->cd(2);
//      hmInt3[im][ic]->Draw((im==0)?"ap":"psame");
//      hmInt4[im][ic]->Draw("psame");
//      cmult->cd(1);
//      hm1[im][ic]->DrawCopy((im==0)?"p":"psame");
//      hm2[im][ic]->DrawCopy("psame");
//      cmult->cd(2);
//      hm3[im][ic]->DrawCopy((im==0)?"p":"psame");
//      hm4[im][ic]->DrawCopy("psame");
//    }
//    cmultInt->cd(1);
//    gPad->BuildLegend()->SetFillColor(0);
//    cmultInt->cd(2);
//    gPad->BuildLegend()->SetFillColor(0);
//    cmult->cd(1);
//    gPad->BuildLegend()->SetFillColor(0);
//    cmult->cd(2);
//    gPad->BuildLegend()->SetFillColor(0);
//  }
//}

//void draw_v2(Int_t imethod = 2, Int_t icentrality=0,TString opt="_LowMultScale"){
//void draw_v2(Int_t imethod = 2, Int_t icentrality=0,TString opt="_hist_integral"){
//void draw_v2(Int_t imethod = kV0S, Int_t icentrality=0,TString opt="_no_subtraction"){
//void draw_v2(Int_t imethod = kV0S, Int_t icentrality=0,TString opt="_constant_fit"){
void draw_v2(Int_t imethod = kV0S, Int_t icentrality=0,TString opt=""){
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetPadRightMargin(0.01);
  gStyle->SetPadTopMargin(0.01);

  // read resolution from Tim-trk
  TFile* fResTrk = new TFile("../Code_ResolutionMC/v2ratio_mutrk_p.root");
  TH1D* hAveTrk = (TH1D*) fResTrk->Get("hAve");
  TH1D* hMinTrk = (TH1D*) fResTrk->Get("hMin");
  TH1D* hMaxTrk = (TH1D*) fResTrk->Get("hMax");
  TFile* fResTkl_b = new TFile("../Code_ResolutionMC/v2ratio_mutrklet_p.root");
  TH1D* hAveTkl_b = (TH1D*) fResTkl_b->Get("hAve");
  TH1D* hMinTkl_b = (TH1D*) fResTkl_b->Get("hMin");
  TH1D* hMaxTkl_b = (TH1D*) fResTkl_b->Get("hMax");
  TFile* fResTkl_f = new TFile("../Code_ResolutionMC/v2ratio_mutrklet_Pb.root");
  TH1D* hAveTkl_f = (TH1D*) fResTkl_f->Get("hAve");
  TH1D* hMinTkl_f = (TH1D*) fResTkl_f->Get("hMin");
  TH1D* hMaxTkl_f = (TH1D*) fResTkl_f->Get("hMax");
  TFile* fResTkl_r = new TFile("../Code_ResolutionMC/RatioCorr.root");
  TH1D* hAveTkl_r = (TH1D*) fResTkl_r->Get("hAve");
  TH1D* hMinTkl_r = (TH1D*) fResTkl_r->Get("hMin");
  TH1D* hMaxTkl_r = (TH1D*) fResTkl_r->Get("hMax");
  new TCanvas();
  Int_t color=4;
  hMinTrk->SetLineColor(color);hMinTrk->SetMarkerColor(color);
  hAveTrk->SetLineColor(color);hAveTrk->SetMarkerColor(color);
  hMaxTrk->SetLineColor(color);  hMaxTrk->SetMarkerColor(color);
  color=1;
  hMinTkl_b->SetLineColor(color);hMinTkl_b->SetMarkerColor(color);
  hAveTkl_b->SetLineColor(color);hAveTkl_b->SetMarkerColor(color);
  hMaxTkl_b->SetLineColor(color);  hMaxTkl_b->SetMarkerColor(color);
  color=2;
  hMinTkl_f->SetLineColor(color);hMinTkl_f->SetMarkerColor(color);
  hAveTkl_f->SetLineColor(color);hAveTkl_f->SetMarkerColor(color);
  hMaxTkl_f->SetLineColor(color);  hMaxTkl_f->SetMarkerColor(color);
  color=7;
  hMinTkl_r->SetLineColor(color);hMinTkl_r->SetMarkerColor(color);
  hAveTkl_r->SetLineColor(color);hAveTkl_r->SetMarkerColor(color);
  hMaxTkl_r->SetLineColor(color);  hMaxTkl_r->SetMarkerColor(color);

  hMinTrk->Draw("l hist");
  hAveTrk->Draw("same");
  hMaxTrk->Draw("l hist same");
  hMinTkl_b->Draw("l hist same");
  hAveTkl_b->Draw("same");
  hMaxTkl_b->Draw("l hist same");
  hMinTkl_f->Draw("l hist same");
  hAveTkl_f->Draw("same");
  hMaxTkl_f->Draw("l hist same");
  hMinTkl_r->Draw("l hist same");
  hAveTkl_r->Draw("same");
  hMaxTkl_r->Draw("l hist same");

  TGraphErrors* gResTrk_b = new TGraphErrors(nbins_muon_trk);
  TGraphErrors* gResTrk_f    = new TGraphErrors(nbins_muon_trk);
  gResTrk_b->SetLineWidth(2);
  gResTrk_f->SetLineWidth(2);
  gResTrk_f->SetLineWidth(2);
  TGraphErrors* gResErrorTrk_b = (TGraphErrors*) gResTrk_b->Clone("gResErrorTrk_b");
  TGraphErrors* gResErrorTrk_f = (TGraphErrors*) gResTrk_f->Clone("gResErrorTrk_f");
  //tkl
  TGraphErrors* gResTkl_b = new TGraphErrors(nbins_muon_tkl);
  TGraphErrors* gResTkl_f    = new TGraphErrors(nbins_muon_tkl);
  TGraphErrors* gResTkl_r    = new TGraphErrors(nbins_muon_tkl);
  gResTkl_b->SetLineWidth(2);
  gResTkl_f->SetLineWidth(2);
  gResTkl_r->SetLineWidth(2);
  TGraphErrors* gResErrorTkl_b = (TGraphErrors*) gResTkl_b->Clone("gResErrorTkl_b");
  TGraphErrors* gResErrorTkl_f = (TGraphErrors*) gResTkl_f->Clone("gResErrorTkl_f");
  TGraphErrors* gResErrorTkl_r = (TGraphErrors*) gResTkl_r->Clone("gResErrorTkl_r");

  for (Int_t ibin=0;ibin<nbins_muon_trk;ibin++)binst_centersTrk[ibin]=(bins_muon_trk[ibin]+bins_muon_trk[ibin+1])/2.;
  for (Int_t ibin=0;ibin<nbins_muon_tkl;ibin++)binst_centersTkl[ibin]=(bins_muon_tkl[ibin]+bins_muon_tkl[ibin+1])/2.;

  // tracks
  for (Int_t ipt=0;ipt<nbins_muon_trk;ipt++){
    Float_t ave = hAveTrk->GetBinContent(hAveTrk->GetXaxis()->FindBin(binst_centersTrk[ipt]));
    Float_t min = hMinTrk->GetBinContent(hAveTrk->GetXaxis()->FindBin(binst_centersTrk[ipt]));
    Float_t max = hMaxTrk->GetBinContent(hAveTrk->GetXaxis()->FindBin(binst_centersTrk[ipt]));
    gResTrk_f->SetPoint(ipt,binst_centersTrk[ipt],ave);
    gResTrk_f->SetPointError(ipt,0,TMath::Max(max-ave,ave-min)/sqrt(2));
    gResTrk_b->SetPoint(ipt,binst_centersTrk[ipt],ave);
    gResTrk_b->SetPointError(ipt,0,TMath::Max(max-ave,ave-min)/sqrt(2));
    gResErrorTrk_b->SetPoint(ipt,binst_centersTrk[ipt],gResTrk_b->GetEY()[ipt]);
    gResErrorTrk_f->SetPoint(ipt,binst_centersTrk[ipt],gResTrk_f->GetEY()[ipt]);
  }
  //tracklets
  for (Int_t ipt=0;ipt<nbins_muon_tkl;ipt++){
    //b
    Float_t ave = hAveTkl_b->GetBinContent(hAveTkl_b->GetXaxis()->FindBin(binst_centersTkl[ipt]));
    Float_t min = hMinTkl_b->GetBinContent(hAveTkl_b->GetXaxis()->FindBin(binst_centersTkl[ipt]));
    Float_t max = hMaxTkl_b->GetBinContent(hAveTkl_b->GetXaxis()->FindBin(binst_centersTkl[ipt]));
    gResTkl_b->SetPoint(ipt,binst_centersTkl[ipt],ave);
    gResTkl_b->SetPointError(ipt,0,TMath::Max(max-ave,ave-min)/sqrt(2));
    gResErrorTkl_b->SetPoint(ipt,binst_centersTkl[ipt],gResTkl_b->GetEY()[ipt]);
    //f
    ave = hAveTkl_f->GetBinContent(hAveTkl_f->GetXaxis()->FindBin(binst_centersTkl[ipt]));
    min = hMinTkl_f->GetBinContent(hAveTkl_f->GetXaxis()->FindBin(binst_centersTkl[ipt]));
    max = hMaxTkl_f->GetBinContent(hAveTkl_f->GetXaxis()->FindBin(binst_centersTkl[ipt]));
    gResTkl_f->SetPoint(ipt,binst_centersTkl[ipt],ave);
    gResTkl_f->SetPointError(ipt,0,TMath::Max(max-ave,ave-min)/sqrt(2));
    gResErrorTkl_f->SetPoint(ipt,binst_centersTkl[ipt],gResTkl_f->GetEY()[ipt]);
    //r
    ave = hAveTkl_r->GetBinContent(hAveTkl_r->GetXaxis()->FindBin(binst_centersTkl[ipt]));
    min = hMinTkl_r->GetBinContent(hAveTkl_r->GetXaxis()->FindBin(binst_centersTkl[ipt]));
    max = hMaxTkl_r->GetBinContent(hAveTkl_r->GetXaxis()->FindBin(binst_centersTkl[ipt]));
    gResTkl_r->SetPoint(ipt,binst_centersTkl[ipt],ave);
    gResTkl_r->SetPointError(ipt,0,TMath::Max(max-ave,ave-min)/sqrt(2));
    gResErrorTkl_r->SetPoint(ipt,binst_centersTkl[ipt],gResTkl_r->GetEY()[ipt]);
  }

  const Int_t maxna = 5;
  TFile* f1[nm];
  TFile* f2[nm];
  TFile* f3[nm];
  TFile* f4[nm];
  TFile* fileSyst_f[nm];
  TFile* fileSyst_b[nm];
  TFile* fileSyst_r[nm];
  TFile* fileSyst_d[nm];
  TF1* fSyst_f[nm][3][nsyst];
  TF1* fSyst_b[nm][3][nsyst];
  TF1* fSyst_r[nm][3][nsyst];
  TF1* fSyst_d[nm][3][nsyst];
  TGraphErrors* gSyst_f[nm][3][nsyst];
  TGraphErrors* gSyst_b[nm][3][nsyst];
  TGraphErrors* gSyst_r[nm][3][nsyst];
  TGraphErrors* gSyst_d[nm][3][nsyst];

  TGraphErrors* vg1[nm][3][maxna]; // mu_tkl_b
  TGraphErrors* vg2[nm][3][maxna]; // mu_tkl_f
  TGraphErrors* vg3[nm][3][maxna]; // tkl_tkl_b
  TGraphErrors* vg4[nm][3][maxna]; // tkl_tkl_f
  TGraphErrors* vg5[nm][3][maxna]; // ratio
  TGraphErrors* vg6[nm][3][maxna]; // diff
  TGraphErrors* vg1s[nm][3];   // mu_tkl_b with syst error
  TGraphErrors* vg2s[nm][3];   // mu_tkl_f with syst error
  TGraphErrors* vg3s[nm][3];   // tkl_tkl_b with syst error
  TGraphErrors* vg4s[nm][3];   // tkl_tkl_f with syst error
  TGraphErrors* vg5s[nm][3];   // ratio with syst error
  TGraphErrors* vg6s[nm][3];   // diff with syst error

  TGraphErrors* vg1_nocorr[nm][3][maxna]; // mu_tkl_b without resolution correction for comparison with Tim
  TGraphErrors* vg1s_nocorr[nm][3];   // mu_tkl_b with syst error but without resolution
  TGraphErrors* vg2_nocorr[nm][3][maxna]; // mu_tkl_f without resolution correction for comparison with Tim
  TGraphErrors* vg2s_nocorr[nm][3];   // mu_tkl_f with syst error but without resolution

  TGraphErrors* fSyst_f_combined[nm][3];
  TGraphErrors* fSyst_b_combined[nm][3];
  TGraphErrors* fSyst_r_combined[nm][3];
  TGraphErrors* fSyst_d_combined[nm][3];
  TGraphErrors* fSyst_b_combined_nocorr[nm][3];
  TGraphErrors* fSyst_f_combined_nocorr[nm][3];

  TCanvas* csyst  = new TCanvas("csyst","syst detailed",1800,1000);
  csyst->Divide(4,nm,0.001,0.001);

  for (Int_t im=0;im<nm;im++){
    Printf("ALIVE %s",centMethods[im].Data());
        Bool_t isTkl=im>kNtrk && im<=kNtkl? 1 : 0;
    Bool_t isAMPTgen=im>kNtkl? 1 : 0;
    Int_t binAss = isTkl ? 4 : 0;
    const Int_t nbinsa = isTkl ? nbins_tkl      : nbins_trk;
    const Int_t nbinst = isTkl ? nbins_muon_tkl : isAMPTgen? nbins_muon_ampt_gen : nbins_muon_trk;
    Double_t* binsa    = isTkl ? bins_tkl       : bins_trk;
    Double_t* binst    = isTkl ? bins_muon_tkl  : isAMPTgen? bins_muon_ampt_gen : bins_muon_trk;
    for (Int_t ibin=0;ibin<nbinst;ibin++) {
      binst_centers[ibin]=(binst[ibin]+binst[ibin+1])/2.;
      binst_widths[ibin]=(binst[ibin+1]-binst[ibin])/2.;
    }
    f1[im] = new TFile(Form("v2/v2_mu_tkl_%s_bcde%s.root" ,centMethods[im].Data(),opt.Data()));
    f2[im] = new TFile(Form("v2/v2_mu_tkl_%s_f%s.root"    ,centMethods[im].Data(),opt.Data()));
    f3[im] = new TFile(Form("v2/v2_tkl_tkl_%s_bcde%s.root",centMethods[im].Data(),opt.Data()));
    f4[im] = new TFile(Form("v2/v2_tkl_tkl_%s_f%s.root"   ,centMethods[im].Data(),opt.Data()));
    fileSyst_b[im] = new TFile(Form("syst/syst_v2_mu_tkl_%s_bcde.root" ,centMethods[im].Data()));
    fileSyst_f[im] = new TFile(Form("syst/syst_v2_mu_tkl_%s_f.root"    ,centMethods[im].Data()));
    fileSyst_r[im] = new TFile(Form("syst/syst_v2_mu_tkl_%s_Ratio.root",centMethods[im].Data()));
    fileSyst_d[im] = new TFile(Form("syst/syst_v2_mu_tkl_%s_Diff.root" ,centMethods[im].Data()));
    for (Int_t ic=0;ic<nc;ic++){
      for (Int_t is=1;is<nsyst;is++){
        fSyst_f[im][ic][is] = (TF1*) fileSyst_f[im]->Get(Form("Unc_%s%s_Cent%i",centMethods[im].Data(),systCheck[is].Data(),ic));
        fSyst_b[im][ic][is] = (TF1*) fileSyst_b[im]->Get(Form("Unc_%s%s_Cent%i",centMethods[im].Data(),systCheck[is].Data(),ic));
        fSyst_r[im][ic][is] = (TF1*) fileSyst_r[im]->Get(Form("Unc_%s%s_Cent%i",centMethods[im].Data(),systCheck[is].Data(),ic));
        fSyst_d[im][ic][is] = (TF1*) fileSyst_d[im]->Get(Form("Unc_%s%s_Cent%i",centMethods[im].Data(),systCheck[is].Data(),ic));
        gSyst_f[im][ic][is] = (TGraphErrors*) fileSyst_f[im]->Get(Form("gUnc_%s%s_Cent%i",centMethods[im].Data(),systCheck[is].Data(),ic));
        gSyst_b[im][ic][is] = (TGraphErrors*) fileSyst_b[im]->Get(Form("gUnc_%s%s_Cent%i",centMethods[im].Data(),systCheck[is].Data(),ic));
        gSyst_r[im][ic][is] = (TGraphErrors*) fileSyst_r[im]->Get(Form("gUnc_%s%s_Cent%i",centMethods[im].Data(),systCheck[is].Data(),ic));
        gSyst_d[im][ic][is] = (TGraphErrors*) fileSyst_d[im]->Get(Form("gUnc_%s%s_Cent%i",centMethods[im].Data(),systCheck[is].Data(),ic));
      }
      for (Int_t ia=0;ia<nbinsa;ia++){
        vg1[im][ic][ia] = (TGraphErrors*) f1[im]->Get(Form("gv2_cent%i_ass%i",ic,ia));
        vg2[im][ic][ia] = (TGraphErrors*) f2[im]->Get(Form("gv2_cent%i_ass%i",ic,ia));
        vg3[im][ic][ia] = (TGraphErrors*) f3[im]->Get(Form("gv2_cent%i_ass%i",ic,ia));
        vg4[im][ic][ia] = (TGraphErrors*) f4[im]->Get(Form("gv2_cent%i_ass%i",ic,ia));
        vg1[im][ic][ia]->SetFillColor(0);
        vg2[im][ic][ia]->SetFillColor(0);
        vg3[im][ic][ia]->SetFillColor(0);
        vg4[im][ic][ia]->SetFillColor(0);

        vg1_nocorr[im][ic][ia] = (TGraphErrors*) vg1[im][ic][ia]->Clone(Form("gv2dphi_nocorr%i%i",ic,ia));
        vg2_nocorr[im][ic][ia] = (TGraphErrors*) vg2[im][ic][ia]->Clone(Form("gv2dphi_nocorr%i%i",ic,ia));

        // calc ratio and diff in binst
        vg5[im][ic][ia] =  (TGraphErrors*) vg1[im][ic][ia]->Clone();
        vg6[im][ic][ia] =  (TGraphErrors*) vg1[im][ic][ia]->Clone();
        for (Int_t ipt=0;ipt<nbinst;ipt++) {
          Float_t val1 = vg1[im][ic][ia]->GetY()[ipt];
          Float_t val2 = vg2[im][ic][ia]->GetY()[ipt];
          Float_t err1 = vg1[im][ic][ia]->GetEY()[ipt];
          Float_t err2 = vg2[im][ic][ia]->GetEY()[ipt];
          Float_t val  = val2/val1;
          Float_t err  = val*TMath::Sqrt(err1*err1/val1/val1+err2*err2/val2/val2);
          vg5[im][ic][ia]->GetY()[ipt]  = val;
          vg5[im][ic][ia]->GetEY()[ipt] = err;
          val  = val2-val1;
          err  = sqrt(err2*err2+err1*err1);
          vg6[im][ic][ia]->GetY()[ipt]  = val;
          vg6[im][ic][ia]->GetEY()[ipt] = err;
        } // trigger pt

        // resolution correction after calculation of the ratio
        // resolution correction
        if(im<kNtkl){//not applied for AMPT
          for (Int_t ipt=0;ipt<nbinst;ipt++) {
            if(im<=kNtrk){
              vg1[im][ic][ia]->GetY()[ipt] = vg1[im][ic][ia]->GetY()[ipt]/gResTrk_b->GetY()[ipt];
              vg1[im][ic][ia]->GetEY()[ipt] = vg1[im][ic][ia]->GetEY()[ipt]/gResTrk_b->GetY()[ipt];
              vg2[im][ic][ia]->GetY()[ipt] = vg2[im][ic][ia]->GetY()[ipt]/gResTrk_f->GetY()[ipt];
              vg2[im][ic][ia]->GetEY()[ipt] = vg2[im][ic][ia]->GetEY()[ipt]/gResTrk_f->GetY()[ipt];
            }else{
              vg1[im][ic][ia]->GetY()[ipt]  = vg1[im][ic][ia]->GetY()[ipt]/gResTkl_b->GetY()[ipt];
              vg1[im][ic][ia]->GetEY()[ipt] = vg1[im][ic][ia]->GetEY()[ipt]/gResTkl_b->GetY()[ipt];
              vg2[im][ic][ia]->GetY()[ipt]  = vg2[im][ic][ia]->GetY()[ipt]/gResTkl_f->GetY()[ipt];
              vg2[im][ic][ia]->GetEY()[ipt] = vg2[im][ic][ia]->GetEY()[ipt]/gResTkl_f->GetY()[ipt];
              vg5[im][ic][ia]->GetY()[ipt]  = vg5[im][ic][ia]->GetY()[ipt]/gResTkl_r->GetY()[ipt];
              vg5[im][ic][ia]->GetEY()[ipt] = vg5[im][ic][ia]->GetEY()[ipt]/gResTkl_r->GetY()[ipt];
            }
          }

        }//trigger pt
      } // associate
    } // centrality bin
  }

  for (Int_t im=0;im<nm;im++){
    Bool_t isTkl=im>kNtrk && im<=kNtkl? 1 : 0;
    Bool_t isAMPTgen=im>kNtkl? 1 : 0;
    Int_t binAss = isTkl ? 4 : 0;
    const Int_t nbinsa = isTkl ? nbins_tkl      : nbins_trk;
    const Int_t nbinst = isTkl ? nbins_muon_tkl : isAMPTgen? nbins_muon_ampt_gen : nbins_muon_trk;
    Double_t* binsa    = isTkl ? bins_tkl       : bins_trk;
    Double_t* binst    = isTkl ? bins_muon_tkl  : isAMPTgen? bins_muon_ampt_gen : bins_muon_trk;
    for (Int_t ibin=0;ibin<nbinst;ibin++) {
      binst_centers[ibin]=(binst[ibin]+binst[ibin+1])/2.;
      binst_widths[ibin]=(binst[ibin+1]-binst[ibin])/2.;
    }
    TLegend* leg;
    TCanvas* csyst2  = new TCanvas(Form("csyst2_%d",im),Form("syst detailed %d",im),1800,1000);
    csyst2->Divide(2,2,0.001,0.001);
    for(Int_t ic=0;ic<nc;ic++){
      fSyst_f_combined[im][ic] = (TGraphErrors*) vg2[im][ic][0]->Clone(Form("fSyst_f_combined%i%i",im,ic));
      fSyst_b_combined[im][ic] = (TGraphErrors*) vg1[im][ic][0]->Clone(Form("fSyst_b_combined%i%i",im,ic));
      fSyst_r_combined[im][ic] = (TGraphErrors*) vg2[im][ic][0]->Clone(Form("fSyst_r_combined%i%i",im,ic));
      fSyst_d_combined[im][ic] = (TGraphErrors*) vg2[im][ic][0]->Clone(Form("fSyst_d_combined%i%i",im,ic));
      Float_t combined_error_f[10]={0};
      Float_t combined_error_b[10]={0};
      Float_t combined_error_r[10]={0};
      Float_t combined_error_d[10]={0};
      Float_t combined_error_nocorr_b[10]={0};
      Float_t combined_error_nocorr_f[10]={0};

      vg1s[im][ic] = (TGraphErrors*) vg1[im][ic][binAss]->Clone(Form("vg1s%i%i",im,ic));
      vg2s[im][ic] = (TGraphErrors*) vg2[im][ic][binAss]->Clone(Form("vg2s%i%i",im,ic));
      vg5s[im][ic] = (TGraphErrors*) vg5[im][ic][binAss]->Clone(Form("vg5s%i%i",im,ic));
      vg6s[im][ic] = (TGraphErrors*) vg6[im][ic][binAss]->Clone(Form("vg6s%i%i",im,ic));
      vg1s_nocorr[im][ic] = (TGraphErrors*) vg1_nocorr[im][ic][binAss]->Clone(Form("vg1s_nocorr%i%i",im,ic));
      vg2s_nocorr[im][ic] = (TGraphErrors*) vg2_nocorr[im][ic][binAss]->Clone(Form("vg2s_nocorr%i%i",im,ic));
      for (Int_t ipt=0;ipt<fSyst_f_combined[im][ic]->GetN();ipt++){
        Float_t error_f[nsyst];
        Float_t error_b[nsyst];
        Float_t error_r[nsyst];
        Float_t error_d[nsyst];
        for (Int_t is=1;is<nsyst;is++){
          if(gSystfromfit){

            error_f[is] = TMath::Abs(fSyst_f[im][ic][is]->Eval(binst_centers[ipt])/TMath::Sqrt(2));
            error_b[is] = TMath::Abs(fSyst_b[im][ic][is]->Eval(binst_centers[ipt])/TMath::Sqrt(2));
            error_r[is] = TMath::Abs(fSyst_r[im][ic][is]->Eval(binst_centers[ipt])/TMath::Sqrt(2));
            error_d[is] = TMath::Abs(fSyst_d[im][ic][is]->Eval(binst_centers[ipt])/TMath::Sqrt(2));
            if (im==kLEE && is==k_LowMultScale) { 
              // we take systematics on remaining jet from tracklet analysis
              // because scaling factors are similar 
              error_b[is] = TMath::Abs(fSyst_b[kV0S][ic][is]->Eval(binst_centers[ipt])/TMath::Sqrt(2));
            }
          }else{
            error_f[is] = TMath::Abs(gSyst_f[im][ic][is]->GetY()[ipt]/TMath::Sqrt(2));
            error_b[is] = TMath::Abs(gSyst_b[im][ic][is]->GetY()[ipt]/TMath::Sqrt(2));
            error_r[is] = TMath::Abs(gSyst_r[im][ic][is]->GetY()[ipt]/TMath::Sqrt(2));
            error_d[is] = TMath::Abs(gSyst_d[im][ic][is]->GetY()[ipt]/TMath::Sqrt(2));
          }
        }
        // acceptance and efficiency correction !! only vertex cut
        combined_error_f[ipt]+=TMath::Power(error_f[k_vertex05],2);//barely significant!
        combined_error_b[ipt]+=TMath::Power(error_b[k_vertex05],2);//barely significant!
        combined_error_r[ipt]+=TMath::Power(error_r[k_vertex05],2);
        combined_error_d[ipt]+=TMath::Power(error_d[k_vertex05],2);
        //remaining jet
        combined_error_f[ipt]+=TMath::Power(TMath::Max(error_f[k_LowMultScale],TMath::Max(error_f[k_exclusion10],error_f[k_exclusion08])),2);
        combined_error_b[ipt]+=TMath::Power(TMath::Max(error_b[k_LowMultScale],TMath::Max(error_b[k_exclusion10],error_b[k_exclusion08])),2);
        combined_error_r[ipt]+=TMath::Power(TMath::Max(error_r[k_LowMultScale],TMath::Max(error_r[k_exclusion10],error_r[k_exclusion08])),2);
        combined_error_d[ipt]+=TMath::Power(TMath::Max(error_d[k_LowMultScale],TMath::Max(error_d[k_exclusion10],error_d[k_exclusion08])),2);
        // remaining ridge in low mult
        combined_error_f[ipt]+=TMath::Power(error_f[k_60_to_70],2);
        combined_error_b[ipt]+=TMath::Power(error_b[k_60_to_70],2);
        combined_error_r[ipt]+=TMath::Power(error_r[k_60_to_70],2);
        combined_error_d[ipt]+=TMath::Power(error_d[k_60_to_70],2);
        // v2 calculation !! NO "_hist_integral","_fit_opt"
        //11 June added constant_fit contribution
        combined_error_f[ipt]+=TMath::Power(TMath::MaxElement(5,&(error_f[k_constant_fit])),2);
        combined_error_b[ipt]+=TMath::Power(TMath::MaxElement(5,&(error_b[k_constant_fit])),2);
        combined_error_r[ipt]+=TMath::Power(TMath::MaxElement(5,&(error_r[k_constant_fit])),2);
        combined_error_d[ipt]+=TMath::Power(TMath::MaxElement(5,&(error_d[k_constant_fit])),2);
        //
        if(ic==0 && im==kV0S){//V0S
          Printf("\033[1;31m\n  %s pt=%f\033[m",centMethods[im].Data(),binst_centers[ipt]);
          Printf("f period");
          Printf("acceptance and efficiency : %.1f %%",100*error_f[k_vertex05]);
          Printf("            remaining jet : %.1f %%",100*TMath::Max(error_f[k_LowMultScale],TMath::Max(error_f[k_exclusion10],error_f[k_exclusion08])));
          Printf("          remaining ridge : %.1f %%",100*error_f[k_60_to_70]);
          Printf("           v2 calculation : %.1f %%",100*TMath::MaxElement(5,&(error_f[k_constant_fit])));
          Printf("b period");
          Printf("acceptance and efficiency : %.1f %%",100*error_b[k_vertex05]);
          Printf("            remaining jet : %.1f %%",100*TMath::Max(error_b[k_LowMultScale],TMath::Max(error_b[k_exclusion10],error_b[k_exclusion08])));
          Printf("          remaining ridge : %.1f %%",100*error_b[k_60_to_70]);
          Printf("           v2 calculation : %.1f %%",100*TMath::MaxElement(5,&(error_b[k_constant_fit])));
          Printf("ratio");
          Printf("acceptance and efficiency : %.1f %%",100*error_r[k_vertex05]);
          Printf("            remaining jet : %.1f %%",100*TMath::Max(error_r[k_LowMultScale],TMath::Max(error_r[k_exclusion10],error_r[k_exclusion08])));
          Printf("          remaining ridge : %.1f %%",100*error_r[k_60_to_70]);
          Printf("           v2 calculation : %.1f %%",100*TMath::MaxElement(5,&(error_r[k_constant_fit])));
          Printf("diff");
          Printf("acceptance and efficiency : %.1f %%",100*error_d[k_vertex05]);
          Printf("            remaining jet : %.1f %%",100*TMath::Max(error_d[k_LowMultScale],TMath::Max(error_d[k_exclusion10],error_d[k_exclusion08])));
          Printf("          remaining ridge : %.1f %%",100*error_d[k_60_to_70]);
          Printf("           v2 calculation : %.1f %%",100*TMath::MaxElement(5,&(error_d[k_constant_fit])));
        }
        if(ic==0 && im==kLEE){//LEE
          Printf("\033[1;31m\n  %s pt=%f\033[m",centMethods[im].Data(),binst_centers[ipt]);
          Printf("b period");
          Printf("acceptance and efficiency : %.1f %%",100*error_b[k_vertex05]);
          Printf("            remaining jet : %.1f %%",100*TMath::Max(error_b[k_LowMultScale],TMath::Max(error_b[k_exclusion10],error_b[k_exclusion08])));
          Printf("          remaining ridge : %.1f %%",100*error_b[k_60_to_70]);
          Printf("           v2 calculation : %.1f %%",100*TMath::MaxElement(4,&(error_b[k_baselineGausFit])));
        }
        combined_error_nocorr_b[ipt]=combined_error_b[ipt];
        combined_error_nocorr_f[ipt]=combined_error_f[ipt];
        // resolution
        if(im<=2){
          combined_error_f[ipt]+=TMath::Power(gResErrorTrk_f->GetY()[ipt],2);
          combined_error_b[ipt]+=TMath::Power(gResErrorTrk_b->GetY()[ipt],2);
          combined_error_r[ipt]+=0;
          combined_error_d[ipt]+=0;
        }else{
          combined_error_f[ipt]+=TMath::Power(gResErrorTkl_f->GetY()[ipt],2);
          combined_error_b[ipt]+=TMath::Power(gResErrorTkl_b->GetY()[ipt],2);
          combined_error_r[ipt]+=TMath::Power(gResErrorTkl_r->GetY()[ipt],2);
          combined_error_d[ipt]+=0;
        }
        // sqrt
        combined_error_f[ipt] =TMath::Sqrt(combined_error_f[ipt]);
        combined_error_b[ipt] =TMath::Sqrt(combined_error_b[ipt]);
        combined_error_r[ipt] =TMath::Sqrt(combined_error_r[ipt]);
        combined_error_d[ipt] =TMath::Sqrt(combined_error_d[ipt]);
        combined_error_nocorr_b[ipt] =TMath::Sqrt(combined_error_nocorr_b[ipt]);
        combined_error_nocorr_f[ipt] =TMath::Sqrt(combined_error_nocorr_f[ipt]);
        fSyst_f_combined[im][ic]->SetPoint(ipt,binst_centers[ipt],combined_error_f[ipt]);
        fSyst_b_combined[im][ic]->SetPoint(ipt,binst_centers[ipt],combined_error_b[ipt]);
        fSyst_r_combined[im][ic]->SetPoint(ipt,binst_centers[ipt],combined_error_r[ipt]);
        fSyst_d_combined[im][ic]->SetPoint(ipt,binst_centers[ipt],combined_error_d[ipt]);
        fSyst_f_combined[im][ic]->SetPointError(ipt,0,0);
        fSyst_b_combined[im][ic]->SetPointError(ipt,0,0);
        fSyst_r_combined[im][ic]->SetPointError(ipt,0,0);
        fSyst_d_combined[im][ic]->SetPointError(ipt,0,0);
        vg1s[im][ic]->SetPointError(ipt,0.1,vg1[im][ic][binAss]->GetY()[ipt]*combined_error_b[ipt]);
        vg2s[im][ic]->SetPointError(ipt,0.1,vg2[im][ic][binAss]->GetY()[ipt]*combined_error_f[ipt]);
        vg5s[im][ic]->SetPointError(ipt,0.1,vg5[im][ic][binAss]->GetY()[ipt]*combined_error_r[ipt]);

        // error on diff taking p-going and Pb-going uncertainties as uncorrelated
        Double_t e1 = vg1s[im][ic]->GetErrorY(ipt);
        Double_t e2 = vg2s[im][ic]->GetErrorY(ipt);
        vg6s[im][ic]->SetPointError(ipt,0.1,TMath::Sqrt(e1*e1+e2*e2));

        // error on diff taking into account correlations
        // vg6s[im][ic]->SetPointError(ipt,0.1,vg6[im][ic][binAss]->GetY()[ipt]*combined_error_d[ipt]);

        vg1s_nocorr[im][ic]->SetPointError(ipt,0.1,vg1_nocorr[im][ic][binAss]->GetY()[ipt]*combined_error_nocorr_b[ipt]);
        vg2s_nocorr[im][ic]->SetPointError(ipt,0.1,vg2_nocorr[im][ic][binAss]->GetY()[ipt]*combined_error_nocorr_f[ipt]);
      }

      leg = new TLegend(0.6,0.7,0.99,0.95);
      leg->AddEntry(fSyst_f_combined[0][icentrality],"Pb-going (f)","l");
      leg->AddEntry(fSyst_b_combined[0][icentrality],"p-going (bcde)","l");
      leg->AddEntry(fSyst_r_combined[0][icentrality],"ratio (f/bcde)","l");

      if (ic!=icentrality) continue;
      TLatex* l = new TLatex();
      l->SetTextSize(0.06);
      l->SetTextFont(42);
      l->SetTextAlign(21);
      l->SetNDC();
      csyst->cd(im*4+1);
      TH1F* fr_f = gPad->DrawFrame(0,-0.5,4,0.5);
      gResErrorTrk_f->SetLineStyle(2);
      gResErrorTkl_f->SetLineStyle(2);
      gResErrorTrk_b->SetLineStyle(2);
      gResErrorTkl_b->SetLineStyle(2);
      gResErrorTkl_r->SetLineStyle(2);
      if(im<=2)gResErrorTrk_f->Draw("l same");
      gResErrorTkl_f->Draw("l same");
      csyst->cd(im*4+2);
      TH1F* fr_b = gPad->DrawFrame(0,-0.5,4,0.5);
      if(im<=2)gResErrorTrk_b->Draw("l same");
      else gResErrorTkl_b->Draw("l same");
      csyst->cd(im*4+3);
      TH1F* fr_r = gPad->DrawFrame(0,-0.5,4,0.5);
      if(im>2)gResErrorTkl_r->Draw("l same");

      for (Int_t is=1;is<nsyst;is++){
        if(is==k_mixed_norm || is==k_hist_integral || is==k_parabolic_fit || is==k_fit_opt || is==k_v2_analytical || is==k_no_subtraction || is==k_vertex01)continue;
        csyst->cd(im*4+1);
        fSyst_f[im][ic][is]->Draw("l same");
        l->DrawLatex(0.25,0.9,Form("f %s",centMethods[im].Data()));
        csyst->cd(im*4+2);
        fSyst_b[im][ic][is]->Draw("l same");
        l->DrawLatex(0.25,0.9,Form("bcde %s",centMethods[im].Data()));
        csyst->cd(im*4+3);
        fSyst_r[im][ic][is]->Draw("l same");
        l->DrawLatex(0.25,0.9,Form("ratio %s",centMethods[im].Data()));
      }
      csyst2->cd(1);
      TH1F* fr_f = gPad->DrawFrame(0,-0.5,4,0.5);
      if(im<=2) gResErrorTrk_f->Draw("l same");
      else gResErrorTkl_f->Draw("l same");
      csyst2->cd(2);
      TH1F* fr_b = gPad->DrawFrame(0,-0.5,4,0.5);
      if(im<=2) gResErrorTrk_b->Draw("l same");
      else gResErrorTkl_b->Draw("l same");
      //      gResErrorTkl_b->Print("all");
      csyst2->cd(3);
      TH1F* fr_r = gPad->DrawFrame(0,-0.5,4,0.5);
      if(im>2)gResErrorTkl_r->Draw("l same");
      for (Int_t is=1;is<nsyst;is++){
        //if(is==2 || is==5 || is==6 || is==11 || is== 15 )continue; //not drawing syst not considered
        csyst2->cd(1);
        fSyst_f[im][ic][is]->Draw("l same");
        l->DrawLatex(0.25,0.9,Form("f %s",centMethods[im].Data()));
        csyst2->cd(2);
        fSyst_b[im][ic][is]->Draw("l same");
        l->DrawLatex(0.25,0.9,Form("bcde %s",centMethods[im].Data()));
        csyst2->cd(3);
        fSyst_r[im][ic][is]->Draw("l same");
        l->DrawLatex(0.25,0.9,Form("ratio %s",centMethods[im].Data()));
      }

    }
    csyst->cd(4*im+3+1);
    TH1F* fr = gPad->DrawFrame(0,0.,4,0.2);
    fSyst_f_combined[im][icentrality]->SetLineColor(kBlue);
    fSyst_b_combined[im][icentrality]->SetLineColor(kRed);
    fSyst_r_combined[im][icentrality]->SetLineWidth(2);
    fSyst_f_combined[im][icentrality]->SetLineWidth(2);
    fSyst_b_combined[im][icentrality]->SetLineWidth(2);
    fSyst_r_combined[im][icentrality]->Draw("l same");
    fSyst_f_combined[im][icentrality]->Draw("l same");
    fSyst_b_combined[im][icentrality]->Draw("l same");
    leg->Draw();
    csyst2->cd(3+1);
    TH1F* fr = gPad->DrawFrame(0,0.,4,0.2);
    fSyst_r_combined[im][icentrality]->Draw("l same");
    fSyst_f_combined[im][icentrality]->Draw("l same");
    fSyst_b_combined[im][icentrality]->Draw("l same");
    leg->Draw();
  }

  csyst->Print("eps/syst_summary.eps");
  csyst->Print("png/syst_summary.png");

  //tkl trk with resolution ocrrection
  //  compare_trk_tkl(vg1_nocorr[kLEE][icentrality][0],vg1_nocorr[kV0S][icentrality][4],vg1s_nocorr[kLEE][icentrality],vg1s_nocorr[kV0S][icentrality]);
  //  compare_trk_tkl(vg1[kLEE][icentrality][0],vg1[kV0S][icentrality][4],vg1s[kLEE][icentrality],vg1s[kV0S][icentrality]);
  //  compare_trk_bug(vg1_nocorr[kLEE][icentrality][0],vg1_nocorr[kLEW][icentrality][0],vg1s_nocorr[kLEE][icentrality],vg1s_nocorr[kLEW][icentrality]);
  //
  //  compare_effCorr_trk(vg1_nocorr[kTIM][icentrality][0],vg1_nocorr[kTIC][icentrality][0],vg1s_nocorr[kTIM][icentrality],vg1s_nocorr[kTIC][icentrality]);
  //
  //  //  draw_tkl (vg3[kV0S][icentrality][4],0x0,vg4[kV0S][icentrality][4],0x0,centMethods[kV0S].Data());
  //  draw_muon_factorization(vg1[kV0S][icentrality],"bcde",centMethods[kV0S].Data());
  //  draw_muon_factorization(vg2[kV0S][icentrality],"f"   ,centMethods[kV0S].Data());
  //  draw_tracklet_factorization(vg3[kV0S][icentrality],"bcde",centMethods[kV0S].Data());
  //  draw_tracklet_factorization(vg4[kV0S][icentrality],"f"   ,centMethods[kV0S].Data());
  //  draw_mu (vg1_nocorr[kV0S][icentrality][4],vg1s_nocorr[kV0S][icentrality],vg2_nocorr[kV0S][icentrality][4],vg2s_nocorr[kV0S][icentrality],centMethods[kV0S].Data(),"no_corr");
//  draw_mu (vg1[kDPM][icentrality][4],vg1s[kDPM][icentrality],vg2[kDPM][icentrality][4],vg2s[kDPM][icentrality],centMethods[kDPM].Data());
//  return;
  draw_mu (vg1[kV0S][icentrality][4],vg1s[kV0S][icentrality],vg2[kV0S][icentrality][4],vg2s[kV0S][icentrality],centMethods[kV0S].Data());
  draw_mu (vg1[kV0L][icentrality][4],vg1s[kV0L][icentrality],vg2[kV0L][icentrality][4],vg2s[kV0L][icentrality],centMethods[kV0L].Data(),"_abs");
  draw_ratio(vg5[kV0S][icentrality][4],vg5s[kV0S][icentrality],centMethods[kV0S].Data());
  draw_ratio(vg5[kV0L][icentrality][4],vg5s[kV0L][icentrality],centMethods[kV0L].Data(),"_abs");
  draw_diff(vg6[kV0S][icentrality][4],vg6s[kV0S][icentrality],centMethods[kV0S].Data());

  //  compare_eta(vg1[kV0S][icentrality][4],vg1s[kV0S][icentrality],vg1[kET1][icentrality][4],vg1s[kET1][icentrality],vg1[kET2][icentrality][4],vg1s[kET2][icentrality],vg1[kET3][icentrality][4],vg1s[kET3][icentrality],"bcde");
  //  compare_eta(vg2[kV0S][icentrality][4],vg2s[kV0S][icentrality],vg2[kET1][icentrality][4],vg2s[kET1][icentrality],vg2[kET2][icentrality][4],vg2s[kET2][icentrality],vg2[kET3][icentrality][4],vg2s[kET3][icentrality],"f");
  //  compare_eta(vg5[kV0S][icentrality][4],vg5s[kV0S][icentrality],vg5[kET1][icentrality][4],vg5s[kET1][icentrality],vg5[kET2][icentrality][4],vg5s[kET2][icentrality],vg5[kET3][icentrality][4],vg5s[kET3][icentrality],"ratio",1);
  //
  //  compare_ampt(vg1[kAMA][icentrality][0],vg1[kAMP][icentrality][0],vg1[kAMK][icentrality][0],vg1[kAMM][icentrality][0],vg1[kAMR][icentrality][0],"bcde");
  //  compare_ampt(vg2[kAMA][icentrality][0],vg2[kAMP][icentrality][0],vg2[kAMK][icentrality][0],vg2[kAMM][icentrality][0],vg2[kAMR][icentrality][0],"f");
  //  compare_ampt(vg5[kAMA][icentrality][0],vg5[kAMP][icentrality][0],vg5[kAMK][icentrality][0],vg5[kAMM][icentrality][0],vg5[kAMR][icentrality][0],"ratio",1);
  //
  //  compare_effCorr_tkl(vg1[kV0S][icentrality][4],vg1[kU0S][icentrality][4],vg1s[kV0S][icentrality],vg1s[kU0S][icentrality],"bcde");
  //  compare_effCorr_tkl(vg2[kV0S][icentrality][4],vg2[kU0S][icentrality][4],vg2s[kV0S][icentrality],vg2s[kU0S][icentrality],"f");
  //  compare_effCorr_tkl(vg5[kV0S][icentrality][4],vg5[kU0S][icentrality][4],vg5s[kV0S][icentrality],vg5s[kU0S][icentrality],"ratio");
  //
  //  compare_multiplicity_estimators(vg2[kV0S][icentrality][4],vg2[kV0D][icentrality][4],vg2[kV0M][icentrality][4],vg2[kZDC][icentrality][4],vg2[kCL1][icentrality][4],Form("f%s",opt.Data()),0);
  //  compare_multiplicity_estimators(vg1[kV0S][icentrality][4],vg1[kV0D][icentrality][4],vg1[kV0M][icentrality][4],vg1[kZDC][icentrality][4],vg1[kCL1][icentrality][4],Form("bcde%s",opt.Data()),0);
  //  compare_multiplicity_estimators(vg5[kV0S][icentrality][4],vg5[kV0D][icentrality][4],vg5[kV0M][icentrality][4],vg5[kZDC][icentrality][4],vg5[kCL1][icentrality][4],Form("ratio%s",opt.Data()),0,1);
  //
  //  compare_multiplicity_class(vg1[kV0S][0][4],vg1[kV0S][1][4],vg1[kV0S][2][4],vg1s[kV0S][0],vg1s[kV0S][1],vg1s[kV0S][2],"bcde");
  //  compare_multiplicity_class(vg2[kV0S][0][4],vg2[kV0S][1][4],vg2[kV0S][2][4],vg2s[kV0S][0],vg2s[kV0S][1],vg2s[kV0S][2],"f");
  //  compare_multiplicity_class(vg5[kV0S][0][4],vg5[kV0S][1][4],vg5[kV0S][2][4],vg5s[kV0S][0],vg5s[kV0S][1],vg5s[kV0S][2],"ratio",1);

}


