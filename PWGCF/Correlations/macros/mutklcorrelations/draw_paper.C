#include "env.h"

Float_t fontSize = 0.04;
const char* kSystemEnergy = "p-Pb #sqrt{s_{NN}} = 5.02 TeV";
Bool_t gShiftToMeanPt=1;

void draw_paper(){
  //  draw_projection(1);
  //  draw_projection(2);
  //  draw_projection(3);
  //  draw2D(1,0);
  //  draw2D(2,0);
  //  draw2D(3,0);
  //  draw2D(1,1);
  //  draw2D(2,1);
  //  draw2D(3,1);
  //  draw2D(1,2);
  //  draw2D(2,2);
  //  draw2D(3,2);

  draw_v2_final_ampt_muon();
//  draw_v2_final();
//  draw_v2_final_ampt();
//  draw_ratio_final();
//  draw_comparison();
  //
}

void draw_v2_final_ampt_muon(){
  TFile* f = new TFile("../Paper_ForwardMuCorr/macros/DPMJET_AMPT.root");
  f->ls();
  TH1D* hEffPi_b = (TH1D*) f->Get("AMPT_pPb_Pion");
  TH1D* hEffKa_b = (TH1D*) f->Get("AMPT_pPb_Kaon");
  TH1D* hEffMu_b = (TH1D*) f->Get("AMPT_pPb_Muon");
  TH1D* hEffPi_f = (TH1D*) f->Get("AMPT_Pbp_Pion");
  TH1D* hEffKa_f = (TH1D*) f->Get("AMPT_Pbp_Kaon");
  TH1D* hEffMu_f = (TH1D*) f->Get("AMPT_Pbp_Muon");
  
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadRightMargin(0.01);
  gStyle->SetPadTopMargin(0.01);
  fontSize = 0.048;
  
//  TFile *famp_b=new TFile("v2/v2_mu_tkl_AMP_bcde_LowMultScale.root");
//  TFile *famp_f=new TFile("v2/v2_mu_tkl_AMP_f_LowMultScale.root");
//  TFile *famk_b=new TFile("v2/v2_mu_tkl_AMK_bcde_LowMultScale.root");
//  TFile *famk_f=new TFile("v2/v2_mu_tkl_AMK_f_LowMultScale.root");
//  TFile *famm_b=new TFile("v2/v2_mu_tkl_AMM_bcde_LowMultScale.root");
//  TFile *famm_f=new TFile("v2/v2_mu_tkl_AMM_f_LowMultScale.root");
  TFile *famp_b=new TFile("v2/v2_mu_tkl_AMP_bcde.root");
  TFile *famp_f=new TFile("v2/v2_mu_tkl_AMP_f.root");
  TFile *famk_b=new TFile("v2/v2_mu_tkl_AMK_bcde.root");
  TFile *famk_f=new TFile("v2/v2_mu_tkl_AMK_f.root");
  TFile *famm_b=new TFile("v2/v2_mu_tkl_AMM_bcde.root");
  TFile *famm_f=new TFile("v2/v2_mu_tkl_AMM_f.root");
  TGraphErrors *gamp_b=famp_b->Get("gv2_cent0_ass0");
  TGraphErrors *gamp_f=famp_f->Get("gv2_cent0_ass0");
  TGraphErrors *gamk_b=famk_b->Get("gv2_cent0_ass0");
  TGraphErrors *gamk_f=famk_f->Get("gv2_cent0_ass0");
  TGraphErrors *gamm_b=famm_b->Get("gv2_cent0_ass0");
  TGraphErrors *gamm_f=famm_f->Get("gv2_cent0_ass0");
  
  gamp_b->SetMarkerStyle(kFullTriangleUp);
  gamp_f->SetMarkerStyle(kOpenTriangleUp);
  gamk_b->SetMarkerStyle(kFullTriangleDown);
  gamk_f->SetMarkerStyle(kOpenTriangleDown);
  gamm_b->SetMarkerStyle(kFullCircle);
  gamm_f->SetMarkerStyle(kOpenCircle);
  gamp_b->SetMarkerColor(kRed);
  gamp_f->SetMarkerColor(kRed);
  gamk_b->SetMarkerColor(kBlue);
  gamk_f->SetMarkerColor(kBlue);
  gamm_b->SetMarkerColor(kBlack);
  gamm_f->SetMarkerColor(kBlack);
  gamp_b->SetMarkerSize(1.4);
  gamp_f->SetMarkerSize(1.4);
  gamk_b->SetMarkerSize(1.4);
  gamk_f->SetMarkerSize(1.4);
  gamm_b->SetMarkerSize(1.4);
  gamm_f->SetMarkerSize(1.4);
  gamp_b->SetLineColor(kRed);
  gamp_b->SetLineWidth(3);
  gamp_b->SetLineStyle(0);
  gamp_b->SetFillColor(10);
  gamp_f->SetLineColor(kRed);
  gamp_f->SetLineWidth(3);
  gamp_f->SetLineStyle(2);
  gamp_f->SetFillColor(10);
  gamk_b->SetLineColor(kBlue);
  gamk_b->SetLineWidth(3);
  gamk_b->SetLineStyle(0);
  gamk_b->SetFillColor(10);
  gamk_f->SetLineColor(kBlue);
  gamk_f->SetLineWidth(3);
  gamk_f->SetLineStyle(2);
  gamk_f->SetFillColor(10);
  gamm_b->SetLineColor(kBlack);
  gamm_b->SetLineWidth(3);
  gamm_b->SetLineStyle(0);
  gamm_b->SetFillColor(10);
  gamm_f->SetLineColor(kBlack);
  gamm_f->SetLineWidth(3);
  gamm_f->SetLineStyle(2);
  gamm_f->SetFillColor(10);
  TGraphErrors *gama_b = (TGraphErrors*) gamm_b->Clone();
  TGraphErrors *gama_f = (TGraphErrors*) gamm_f->Clone();
  TGraphErrors *gratio = (TGraphErrors*) gamm_f->Clone();
  gama_b->SetLineColor(kGreen+2);
  gama_f->SetLineColor(kGreen+2);
  gama_b->SetMarkerColor(kGreen+2);
  gama_f->SetMarkerColor(kGreen+2);

  gama_b->Set(gamm_b->GetN());
  gama_f->Set(gamm_b->GetN());
  gratio->Set(gamm_b->GetN());

  for (Int_t i=0;i<gamm_b->GetN();i++){
    Float_t pt = gamm_b->GetX()[i];
    Float_t v2Pi_b = gamp_b->GetY()[i];
    Float_t v2Ka_b = gamk_b->GetY()[i];
    Float_t v2Mu_b = gamm_b->GetY()[i];
    Float_t v2Pi_f = gamp_f->GetY()[i];
    Float_t v2Ka_f = gamk_f->GetY()[i];
    Float_t v2Mu_f = gamm_f->GetY()[i];
    Float_t v2Pi_e_b = gamp_b->GetErrorY(i);
    Float_t v2Ka_e_b = gamk_b->GetErrorY(i);
    Float_t v2Mu_e_b = gamm_b->GetErrorY(i);
    Float_t v2Pi_e_f = gamp_f->GetErrorY(i);
    Float_t v2Ka_e_f = gamk_f->GetErrorY(i);
    Float_t v2Mu_e_f = gamm_f->GetErrorY(i);
    Float_t v2Mu_e_b = 0;
    Float_t v2Mu_e_f = 0;
    Float_t v2Mu_b = 0.0;
    Float_t v2Mu_f = 0.0;
    gamm_b->SetPoint(i,pt,0);
    gamm_f->SetPoint(i,pt,0);
    gamm_b->SetPointError(i,0,0);
    gamm_f->SetPointError(i,0,0);
    Float_t effPi_b = hEffPi_b->GetBinContent(hEffPi_b->FindBin(pt));
    Float_t effKa_b = hEffKa_b->GetBinContent(hEffKa_b->FindBin(pt));
    Float_t effMu_b = hEffMu_b->GetBinContent(hEffMu_b->FindBin(pt));
    Float_t effPi_f = hEffPi_f->GetBinContent(hEffPi_f->FindBin(pt));
    Float_t effKa_f = hEffKa_f->GetBinContent(hEffKa_f->FindBin(pt));
    Float_t effMu_f = hEffMu_f->GetBinContent(hEffMu_f->FindBin(pt));
    Float_t v2Al_b = v2Pi_b*effPi_b+v2Ka_b*effKa_b+v2Mu_b*effMu_b;
    Float_t v2Al_f = v2Pi_f*effPi_f+v2Ka_f*effKa_f+v2Mu_f*effMu_f;
    Float_t v2Al_e_b = v2Pi_e_b*effPi_b+v2Ka_e_b*effKa_b+v2Mu_e_b*effMu_b;
    Float_t v2Al_e_f = v2Pi_e_f*effPi_f+v2Ka_e_f*effKa_f+v2Mu_e_f*effMu_f;
    Float_t ratio = v2Al_f/v2Al_b;
    gama_b->SetPoint(i,pt,v2Al_b);
    gama_f->SetPoint(i,pt,v2Al_f);
    gama_b->SetPointError(i,0,v2Al_e_b);
    gama_f->SetPointError(i,0,v2Al_e_f);
    gratio->SetPoint(i,pt,ratio);
//    if (i!=gamm_b->GetN()-1) continue;
//    Int_t n=gamm_b->GetN();
//    gama_b->SetPoint(i+2,pt+0.3,gama_b->GetY()[n]+(gama_b->GetY()[n]-gama_b->GetY()[n-1])/(gama_b->GetX()[n]-gama_b->GetX()[n-1])*0.3);
//    gama_f->SetPoint(i+2,pt+0.3,gama_f->GetY()[n]+(gama_f->GetY()[n]-gama_f->GetY()[n-1])/(gama_f->GetX()[n]-gama_f->GetX()[n-1])*0.3);
//    gama_b->SetPointError(i+2,0,v2Al_e_b);
//    gama_f->SetPointError(i+2,0,v2Al_e_f);
//    gratio->SetPoint(i+2,pt+0.3,gratio->GetY()[n]+(gratio->GetY()[n]-gratio->GetY()[n-1])/(gratio->GetX()[n]-gratio->GetX()[n-1])*0.3);
  }
//  gama_b->SetPoint(0,0.0,gama_b->GetY()[1]+(gama_b->GetY()[2]-gama_b->GetY()[1])/(gama_b->GetX()[2]-gama_b->GetX()[1])*(0.-gama_b->GetX()[1]));
//  gama_f->SetPoint(0,0.0,gama_f->GetY()[1]+(gama_f->GetY()[2]-gama_f->GetY()[1])/(gama_f->GetX()[2]-gama_f->GetX()[1])*(0.-gama_f->GetX()[1]));
//  gratio->SetPoint(0,0.0,gratio->GetY()[1]+(gratio->GetY()[2]-gratio->GetY()[1])/(gratio->GetX()[2]-gratio->GetX()[1])*(0.-gratio->GetX()[1]));
//  gratio->SetPointError(0,0,0);
  
  TCanvas* c1 = new TCanvas("mu_tkl_ampt","mu_tkl_ampt",800,650);
  TH1F* fr = draw_frame(4,0.255);
  fr->SetMinimum(-0.02);
  fr->GetYaxis()->SetTitleSize(0.065);
  fr->GetYaxis()->SetTitle("v_{2}");
  fr->GetYaxis()->SetTitleOffset(1.0);
  fr->GetXaxis()->SetTitleSize(0.06);
  fr->GetXaxis()->SetTitleOffset(1.00);

  TLatex* latex = new TLatex();
  latex->SetTextSize(fontSize-0.003);
  latex->SetTextFont(42);
  latex->SetTextAlign(21);
  latex->SetNDC();
  latex->DrawLatex(0.35,0.93,"ALICE");
  latex->DrawLatex(0.35,0.87,kSystemEnergy);
  latex->DrawLatex(0.35,0.81,"V0S: (0-20%)-(60-100%)");

  gamp_b->Draw("lpsame");
  gamp_f->Draw("lpsame");
  gamk_b->Draw("lpsame");
  gamk_f->Draw("lpsame");
  gamm_b->Draw("lpsame");
  gamm_f->Draw("lpsame");
  gama_b->Draw("lpsame");
  gama_f->Draw("lpsame");
  TLegend* l = new TLegend(0.53,0.60,0.988,0.97);
  l->SetMargin(0.2);
  l->SetFillColor(0);
  l->SetBorderSize(0);
  l->AddEntry(gamp_b,"AMPT, #pi #rightarrow #mu, p-going");
  l->AddEntry(gamp_f,"AMPT, #pi #rightarrow #mu, Pb-going");
  l->AddEntry(gamk_b,"AMPT, K #rightarrow #mu, p-going");
  l->AddEntry(gamk_f,"AMPT, K #rightarrow #mu, Pb-going");
  l->AddEntry(gamm_b,"AMPT, HFM, p-going");
  l->AddEntry(gamm_f,"AMPT, HFM, Pb-going");
  l->AddEntry(gama_b,"AMPT, Mix, p-going");
  l->AddEntry(gama_f,"AMPT, Mix, Pb-going");
  l->Draw("same");

  gPad->Print("eps/v2_our_ampt.eps");
  gPad->Print("eps/v2_our_ampt.pdf");
  gPad->Print("png/v2_our_ampt.png");
  
//  TFile* f1 = new TFile("our_ampt_scaled.root","recreate");
  TFile* f1 = new TFile("our_ampt.root","recreate");
  gama_b->Write("gama_b");
  gama_f->Write("gama_f");
  gratio->Write("gratio");
  f1->Close();
}


void draw_v2_final_ampt(){
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadRightMargin(0.01);
  gStyle->SetPadTopMargin(0.01);
  fontSize = 0.048;
  TFile* f = new TFile("v2final_mu_tkl.root");
  f->ls();
  TGraphErrors* g_b     = (TGraphErrors*) f->Get("g_b");
  TGraphErrors* g_f     = (TGraphErrors*) f->Get("g_f");
  TGraphErrors* gsyst_b = (TGraphErrors*) f->Get("gsyst_b");
  TGraphErrors* gsyst_f = (TGraphErrors*) f->Get("gsyst_f");
  if(gShiftToMeanPt){
    ShiftToMeanPt(g_b,kV0S,"bcde");
    ShiftToMeanPt(g_f,kV0S,"f");
    ShiftToMeanPt(gsyst_b,kV0S,"bcde");
    ShiftToMeanPt(gsyst_f,kV0S,"f");
  }

  TCanvas* c1 = new TCanvas("mu_tkl_ampt","mu_tkl_ampt",800,650);
  TH1F* fr = draw_frame(4,0.255);
  fr->GetYaxis()->SetTitleSize(0.065);
  fr->GetYaxis()->SetTitle("v_{2}");
  fr->GetYaxis()->SetTitleOffset(1.0);
  fr->GetXaxis()->SetTitleSize(0.06);
  fr->GetXaxis()->SetTitleOffset(1.00);

  TLatex* latex = new TLatex();
  latex->SetTextSize(fontSize-0.003);
  latex->SetTextFont(42);
  latex->SetTextAlign(21);
  latex->SetNDC();
  latex->DrawLatex(0.35,0.93,"ALICE");
  latex->DrawLatex(0.35,0.87,kSystemEnergy);
  latex->DrawLatex(0.35,0.81,"V0S: (0-20%)-(60-100%)");
  //  latex->DrawLatex(0.40,0.79,Form("%s mult: (0-20%) - (60-100%)",method.Data()));
  TGraph* gHydro_f = Create_Hydro_f();
  TGraph* gHydro_b = Create_Hydro_b();
  TGraph* gAMPT_f = Create_AMPT_f();
  TGraph* gAMPT_b = Create_AMPT_b();
  gsyst_b->SetMarkerStyle(kOpenCircle);
  gsyst_f->SetMarkerStyle(kFullCircle);
  gsyst_b->SetMarkerSize(1.5);
  gsyst_f->SetMarkerSize(1.5);
  gsyst_b->SetMarkerColor(kRed);
  gsyst_f->SetMarkerColor(kRed);
  gsyst_b->SetLineColor(kRed-10);
  gsyst_f->SetLineColor(kRed);
  gsyst_b->SetFillColor(kRed-10);
  gsyst_f->SetFillColor(kRed);
  gsyst_b->SetFillStyle(1001);
  gsyst_f->SetFillStyle(1);
  gsyst_f->SetLineWidth(1);
  g_b->SetMarkerStyle(kOpenCircle);
  g_f->SetMarkerStyle(kFullCircle);
  g_b->SetMarkerSize(1.5);
  g_f->SetMarkerSize(1.5);
  g_b->SetMarkerColor(kRed);
  g_f->SetMarkerColor(kRed);
  g_b->SetLineColor(kRed);
  g_f->SetLineColor(kRed);
  g_b->SetLineWidth(3);
  g_f->SetLineWidth(3);
  gsyst_b->Draw("2 same");
  gsyst_f->DrawClone("2 same");
  gHydro_b->Draw("lsame");
  gHydro_f->Draw("lsame");
  gAMPT_b->Draw("lsame");
  gAMPT_f->Draw("lsame");
  g_b->Draw("pz same");
  g_f->DrawClone("pz same");

  TFile *fampt_b=new TFile("v2/v2_mu_tkl_AMA_bcde.root");
  TGraphErrors *gampt_b=fampt_b->Get("gv2_cent0_ass0");
  TFile *fampt_f=new TFile("v2/v2_mu_tkl_AMA_f.root");
  TGraphErrors *gampt_f=fampt_f->Get("gv2_cent0_ass0");
  gampt_b->SetLineColor(kGreen+2);
  gampt_b->SetLineWidth(3);
  gampt_b->SetLineStyle(0);
  gampt_b->SetFillColor(10);
  gampt_f->SetLineColor(kGreen+2);
  gampt_f->SetLineWidth(3);
  gampt_f->SetLineStyle(2);
  gampt_f->SetFillColor(10);
  gampt_b->Draw("lsame");
  gampt_f->Draw("lsame");

  //AMPT unsubtracted
  TFile *fampt_b_nosub=new TFile("v2/v2_mu_tkl_AMA_bcde_no_subtraction.root");
  TGraphErrors *gampt_b_nosub=fampt_b_nosub->Get("gv2_cent0_ass0");
  TFile *fampt_f_nosub=new TFile("v2/v2_mu_tkl_AMA_f_no_subtraction.root");
  TGraphErrors *gampt_f_nosub=fampt_f_nosub->Get("gv2_cent0_ass0");
  gampt_b_nosub->SetLineColor(kRed);
  gampt_b_nosub->SetLineWidth(3);
  gampt_b_nosub->SetLineStyle(0);
  gampt_b_nosub->SetFillColor(10);
  gampt_f_nosub->SetLineColor(kRed);
  gampt_f_nosub->SetLineWidth(3);
  gampt_f_nosub->SetLineStyle(2);
  gampt_f_nosub->SetFillColor(10);
  gampt_b_nosub->Draw("lsame");
  gampt_f_nosub->Draw("lsame");

  //AMPT unsubtracted
  TFile *fampt_b_lmscal=new TFile("v2/v2_mu_tkl_AMA_bcde_LowMultScale.root");
  TGraphErrors *gampt_b_lmscal=fampt_b_lmscal->Get("gv2_cent0_ass0");
  TFile *fampt_f_lmscal=new TFile("v2/v2_mu_tkl_AMA_f_LowMultScale.root");
  TGraphErrors *gampt_f_lmscal=fampt_f_lmscal->Get("gv2_cent0_ass0");
  gampt_b_lmscal->SetLineColor(kGray);
  gampt_b_lmscal->SetLineWidth(3);
  gampt_b_lmscal->SetLineStyle(0);
  gampt_b_lmscal->SetFillColor(10);
  gampt_f_lmscal->SetLineColor(kGray);
  gampt_f_lmscal->SetLineWidth(3);
  gampt_f_lmscal->SetLineStyle(2);
  gampt_f_lmscal->SetFillColor(10);
  gampt_b_lmscal->Draw("lsame");
  gampt_f_lmscal->Draw("lsame");


  TLegend* l = new TLegend(0.55,0.63,0.988,0.97);
  l->SetMargin(0.2);
  l->SetFillColor(0);
  l->SetBorderSize(0);
  l->AddEntry(gsyst_f,"v_{2}^{#mu}{2PC,sub}, Pb-going","PF");
  l->AddEntry(gsyst_b,"v_{2}^{#mu}{2PC,sub}, p-going"),"PF";
  l->AddEntry(gHydro_f,"Hydro, Pb-going");
  l->AddEntry(gHydro_b,"Hydro, p-going");
  l->AddEntry(gAMPT_f,"AMPT, Pb-going");
  l->AddEntry(gAMPT_b,"AMPT, p-going");
  l->AddEntry(gampt_f,"AMPT (our, subtracted), Pb-going","l");
  l->AddEntry(gampt_b,"AMPT (our, subtracted), p-going","l");
  l->AddEntry(gampt_f_nosub,"AMPT (our, unsubtracted), Pb-going","l");
  l->AddEntry(gampt_b_nosub,"AMPT (our, unsubtracted), p-going","l");
  l->AddEntry(gampt_f_lmscal,"AMPT (our, low-mult scaled), Pb-going","l");
  l->AddEntry(gampt_b_lmscal,"AMPT (our, low-mult scaled), p-going","l");
  l->Draw("same");
  gPad->Print("eps/v2_final_ampt.eps");
  gPad->Print("eps/v2_final_ampt.pdf");
  gPad->Print("png/v2_final_ampt.png");
}


void draw_projection(Int_t i=1){ 
  TString fileName;
  TString analysisType;
  TString assRange;
  TString syst;

  if (i==1){
    fileName     = "debug_mu_tkl_LEE_bcde.root";
    analysisType = "Assoc. tracks";
    assRange     = "0.5 < #it{p}_{T}^{a} (GeV/c) < 4";
    syst         = "p-Pb #sqrt{s_{NN}} = 5.02 TeV";
  } else if (i==2){
    fileName     = "debug_mu_tkl_V0S_bcde.root";
    analysisType = "Assoc. tracklets";
    assRange     = "0 < #Delta#varphi^{a} (mrad) < 5";
    syst         = "p-Pb #sqrt{s_{NN}} = 5.02 TeV";
  } else if (i==3){
    fileName     = "debug_mu_tkl_V0S_f.root";
    analysisType = "Assoc. tracklets";
    assRange     = "0 < #Delta#varphi^{a} (mrad) < 5";
    syst         = "Pb-p #sqrt{s_{NN}} = 5.02 TeV";
  }
  gStyle->SetPadRightMargin(0.01);
  gStyle->SetPadTopMargin(0.01);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.11);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);
  TFile* f1 = new TFile(fileName.Data());
  TH2D* hist1 = (TH2D*) f1->Get(i==1 ? "dphi_0_0_0_proj1x" : "dphi_0_4_0_proj1x");
  TF1* vn = (TF1*) f1->Get("vn");
  TF1* v1 = (TF1*) f1->Get("v1");
  TF1* v2 = (TF1*) f1->Get("v2");
  TF1* v3 = (TF1*) f1->Get("v3");
  Float_t fontSizeLocal = 0.05;
  hist1->GetYaxis()->SetTitle("Y_{sub} per #Delta#eta (rad^{-1})");
  hist1->GetXaxis()->SetTickLength(0.02);
  hist1->GetYaxis()->SetTickLength(0.02);
  hist1->GetYaxis()->SetDecimals();
  hist1->GetYaxis()->SetTitleSize(fontSizeLocal);
  hist1->GetXaxis()->SetTitleSize(fontSizeLocal);
  hist1->GetXaxis()->SetLabelSize(fontSizeLocal);
  hist1->GetYaxis()->SetLabelSize(fontSizeLocal);
  hist1->GetYaxis()->SetTitleOffset(1.55);
  TLatex* latex = new TLatex();
  latex->SetTextSize(0.044);
  latex->SetTextFont(42);
  latex->SetTextAlign(21);
  latex->SetNDC();
  TCanvas* c1 = new TCanvas(Form("c%i",i),Form("c%i",i),600*(i-1),20,600,600);
  hist1->SetLineWidth(3);
  hist1->SetLineColor(kBlack);
  hist1->SetMarkerSize(1.1);

  hist1->Draw("E0 X0");
  vn->Draw("same");
  v1->Draw("same");
  v2->Draw("same");
  v3->Draw("same");
  vn->SetLineWidth(3);
  v1->SetLineWidth(3);
  v2->SetLineWidth(3);
  v3->SetLineWidth(3);
  v1->SetLineStyle(7);
  v2->SetLineStyle(5);
  v3->SetLineStyle(3);
  v1->SetLineColor(kRed);
  v2->SetLineColor(kBlue);
  v3->SetLineStyle(3);
  latex->DrawLatex(0.40,0.94,"ALICE");
  latex->DrawLatex(0.40,0.88,syst);
  latex->DrawLatex(0.40,0.82,"V0S: (0-20%)-(60-100%)");
  //latex->DrawLatex(0.805,0.94,analysisType.Data());
  latex->DrawLatex(0.805,0.94,"0.5 < #it{p}_{T}^{t} (GeV/c) < 1");
  latex->DrawLatex(0.805,0.88,analysisType.Data());
  //latex->DrawLatex(0.805,0.82,assRange.Data());
  TLegend* legend1 = new TLegend(0.17, 0.707, 0.33, 0.78);
  TLegend* legend2 = new TLegend(0.35, 0.707, 0.95, 0.78);
  TLegend* legend3 = new TLegend(0.17, 0.60, 0.39, 0.68);
  TLegend* legend4 = new TLegend(0.47, 0.60, 0.69, 0.68);
  TLegend* legend5 = new TLegend(0.76, 0.60, 0.98, 0.68);
  legend1->AddEntry(hist1,"Data", "P");
  legend2->AddEntry(vn   , "#it{a}_{0} + #sum_{#it{n}=1}^{3} 2#it{a}_{#it{n}} cos(#it{n}#Delta#varphi) fit", "L");
  legend3->AddEntry(v1   ,"#it{n} = 1", "L");
  legend4->AddEntry(v2   ,"#it{n} = 2", "L");
  legend5->AddEntry(v3   ,"#it{n} = 3", "L");
  legend1->Draw();
  legend2->Draw();
  legend3->Draw();
  legend4->Draw();
  legend5->Draw();
  gPad->Print(Form("eps/proj%i.eps",i));
  gPad->Print(Form("eps/proj%i.pdf",i));
}


void draw2D(Int_t i=1, Int_t opt=0){
  fontSize = 0.05;
  gStyle->SetPadLeftMargin(0.01);
  gStyle->SetPadBottomMargin(0.00);
  gStyle->SetPadRightMargin(0.00);
  gStyle->SetPadTopMargin(0.00);
  gStyle->SetPalette(55);
  gStyle->SetNumberContours(50);
  TString fileName;
  TString analysisType;
  TString assRange;
  TString syst;
  Float_t zmin;
  Float_t zmax;

  if (i==1){
    fileName     = "debug_mu_tkl_LEE_bcde.root";
    analysisType = "Assoc. tracks";
    assRange     = "0.5 < #it{p}_{T}^{a} (GeV/c) < 4";
    syst         = "p-Pb #sqrt{s_{NN}} = 5.02 TeV";
    if (opt==0){
      zmin = 1.9301;
      zmax = 2.00;
    } else if (opt==1){
      zmin = 2.45;
      zmax = 2.55;
    } else if (opt==2){
      zmin = 0.495;
      zmax = 0.57;
    }
  } else if (i==2){
    fileName     = "debug_mu_tkl_V0S_bcde.root";
    analysisType = "Assoc. tracklets";
    assRange     = "0 < #Delta#varphi^{a} (mrad) < 5";
    syst         = "p-Pb #sqrt{s_{NN}} = 5.02 TeV";
    if (opt==0){
      zmin = 1.7601;
      zmax = 1.81;
    } else if (opt==1){
      zmin = 2.267;
      zmax = 2.36;
    } else if (opt==2){
      zmin = 0.495;
      zmax = 0.56;
    }
  } else if (i==3){
    fileName     = "debug_mu_tkl_V0S_f.root";
    analysisType = "Assoc. tracklets";
    assRange     = "0 < #Delta#varphi^{a} (mrad) < 5";
    syst         = "Pb-p #sqrt{s_{NN}} = 5.02 TeV";
    if (opt==0){
      zmin = 1.8401;
      zmax = 1.89;
    } else if (opt==1){
      zmin = 2.36;
      zmax = 2.46;
    } else if (opt==2){
      zmin = 0.5101;
      zmax = 0.56;
    }
  }

  TString cent;
  if (opt==0)      cent = "V0S: (0-20%)-(60-100%)";
  else if (opt==1) cent = "V0S: 0-20%";
  else if (opt==2) cent = "V0S: 60-100%";

  TFile* f1 = new TFile(fileName.Data());
  TH2D* hist1;
  if (opt==0) hist1 = (TH2D*) f1->Get("histSub");
  if (opt==1) hist1 = (TH2D*) f1->Get(i==1 ? "dphi_0_0_0" : "dphi_0_4_0");
  if (opt==2) hist1 = (TH2D*) f1->Get(i==1 ? "dphi_0_0_3" : "dphi_0_4_3");
  hist1->GetZaxis()->SetDecimals(kTRUE);
  hist1->GetYaxis()->SetNdivisions(6,kTRUE);
  hist1->GetZaxis()->SetNdivisions(6,kTRUE);
  hist1->GetXaxis()->SetTitleSize(fontSize+0.02);
  hist1->GetYaxis()->SetTitleSize(fontSize+0.02);
  hist1->GetZaxis()->SetTitleSize(fontSize);
  hist1->GetXaxis()->SetLabelSize(fontSize);
  hist1->GetYaxis()->SetLabelSize(fontSize);
  hist1->GetZaxis()->SetLabelSize(fontSize);
  hist1->GetXaxis()->CenterTitle();
  hist1->GetYaxis()->CenterTitle();
  hist1->GetZaxis()->CenterTitle();
  hist1->GetXaxis()->SetTitleOffset(1.2);
  hist1->GetYaxis()->SetTitleOffset(1.2);
  hist1->GetZaxis()->SetTitleOffset(1.7);
  hist1->GetZaxis()->SetTitle(opt==0 ? "Y_{sub} (rad^{-1})" : "Y (rad^{-1})");

  TLatex* latex = new TLatex();
  latex->SetTextSize(fontSize);
  latex->SetTextFont(42);
  latex->SetTextAlign(21);
  latex->SetNDC();
  TCanvas* c1 = new TCanvas(Form("c_subt%i%i",i,opt),Form("c_subt%i",i),600*(i-1),20,600,600);
  PadFor2DCorr();
  hist1->GetZaxis()->SetRangeUser(zmin,zmax);
  hist1->Draw("surf1 fb");
  latex->DrawLatex(0.25,0.95,"ALICE");
  latex->DrawLatex(0.25,0.89,syst.Data());
  latex->DrawLatex(0.25,0.83,cent.Data());
  //latex->DrawLatex(0.25,0.77,analysisType.Data());
  latex->DrawLatex(0.78,0.95,"0.5 < #it{p}_{T}^{t} (GeV/c) < 1");
  latex->DrawLatex(0.78,0.88,analysisType.Data());
  //latex->DrawLatex(0.78,0.88,assRange.Data());
  TString outName;
  if (opt==0) outName = "subt";
  if (opt==1) outName = "high";
  if (opt==2) outName = "low";
  gPad->Print(Form("eps/%s%i.eps",outName.Data(),i));
  gPad->Print(Form("eps/%s%i.pdf",outName.Data(),i));
}


void draw_comparison(){
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadRightMargin(0.01);
  gStyle->SetPadTopMargin(0.01);
  fontSize = 0.048;
  TFile* f = new TFile("trk_tkl_final.root");
  f->ls();
  TGraphErrors* g_trk     = (TGraphErrors*) f->Get("g_trk");
  TGraphErrors* g_tkl     = (TGraphErrors*) f->Get("g_tkl");
  TGraphErrors* gsyst_trk = (TGraphErrors*) f->Get("gsyst_trk");
  TGraphErrors* gsyst_tkl = (TGraphErrors*) f->Get("gsyst_tkl");
  if(gShiftToMeanPt){
    ShiftToMeanPt(g_trk,kLEE,"bcde");
    ShiftToMeanPt(g_tkl,kV0S,"bcde");
    ShiftToMeanPt(gsyst_trk,kLEE,"bcde");
    ShiftToMeanPt(gsyst_tkl,kV0S,"bcde");
  }

  TCanvas* c1 = new TCanvas("trk_tkl","trk_tkl",800,650);
  TH1F* fr = draw_frame(4,0.125);
  fr->GetYaxis()->SetTitleSize(0.065);
  fr->GetYaxis()->SetTitle("v_{2}^{#mu}{2PC,sub}");
  fr->GetXaxis()->SetTitleSize(0.06);
  fr->GetXaxis()->SetTitleOffset(1.00);
  fr->GetYaxis()->SetTitleOffset(1.);
  TLatex* latex = new TLatex();
  latex->SetTextSize(fontSize-0.003);
  latex->SetTextFont(42);
  latex->SetTextAlign(21);
  latex->SetNDC();
  latex->DrawLatex(0.35,0.93,"ALICE");
  latex->DrawLatex(0.35,0.87,kSystemEnergy);
  latex->DrawLatex(0.35,0.81,"V0S: (0-20%)-(60-100%)");
  gsyst_trk->SetMarkerStyle(kFullSquare);
  gsyst_tkl->SetMarkerStyle(kOpenCircle);
  gsyst_trk->SetMarkerSize(1.3);
  gsyst_tkl->SetMarkerSize(1.3);
  gsyst_trk->SetMarkerColor(kBlue);
  gsyst_tkl->SetMarkerColor(kRed);
  gsyst_trk->SetLineColor(kBlue);
  gsyst_tkl->SetLineColor(kRed-10);
  gsyst_trk->SetFillColor(kBlue);
  gsyst_tkl->SetFillColor(kRed-10);
  gsyst_trk->SetFillStyle(1);
  gsyst_tkl->SetFillStyle(1001);
  gsyst_tkl->SetLineWidth(1);
  gsyst_trk->SetLineWidth(1);
  g_trk->SetMarkerStyle(kFullSquare);
  g_tkl->SetMarkerStyle(kOpenCircle);
  g_trk->SetMarkerSize(1.3);
  g_tkl->SetMarkerSize(1.3);
  g_trk->SetMarkerColor(kBlue);
  g_tkl->SetMarkerColor(kRed);
  g_trk->SetLineColor(kBlue);
  g_tkl->SetLineColor(kRed);
  g_trk->SetLineWidth(3);
  g_tkl->SetLineWidth(3);
  gsyst_tkl->DrawClone("2 same");
  g_tkl->DrawClone("pz same");
  gsyst_trk->Draw("2 same");
  g_trk->Draw("pz same");
  TLegend* l = new TLegend(0.60,0.80,0.988,0.97);
  l->SetMargin(0.15);
  l->SetFillColor(0);
  l->SetBorderSize(0);
  l->AddEntry(gsyst_trk,"Assoc. tracks","PF");
  l->AddEntry(gsyst_tkl,"Assoc. tracklets","PF");
  l->Draw("same");
  gPad->Print("eps/trk_vs_tkl.eps");
  gPad->Print("eps/trk_vs_tkl.pdf");
  gPad->Print("png/trk_vs_tkl.png");
}


void draw_v2_final(){
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadRightMargin(0.01);
  gStyle->SetPadTopMargin(0.01);
  gStyle->SetHatchesLineWidth(1);
  gStyle->SetHatchesSpacing(0.6);
  fontSize = 0.048;
  TFile* f = new TFile("v2final_mu_tkl.root");
  f->ls();
  TGraphErrors* g_b     = (TGraphErrors*) f->Get("g_b");
  TGraphErrors* g_f     = (TGraphErrors*) f->Get("g_f");
  TGraphErrors* gsyst_b = (TGraphErrors*) f->Get("gsyst_b");
  TGraphErrors* gsyst_f = (TGraphErrors*) f->Get("gsyst_f");
  if(gShiftToMeanPt){
    ShiftToMeanPt(g_b,kV0S,"bcde");
    ShiftToMeanPt(g_f,kV0S,"f");
    ShiftToMeanPt(gsyst_b,kV0S,"bcde");
    ShiftToMeanPt(gsyst_f,kV0S,"f");
  }

  TCanvas* c1 = new TCanvas("mu_tkl","mu_tkl",800,650);
  TH1F* fr = draw_frame(4,0.13);
  // TH1F* fr = draw_frame(4,0.255);
  fr->GetYaxis()->SetTitleSize(0.065);
  // fr->GetYaxis()->SetTitle("v_{2}");
  fr->GetYaxis()->SetTitle("v_{2}^{#mu}{2PC,sub}");
  fr->GetYaxis()->SetTitleOffset(1.0);
  fr->GetXaxis()->SetTitleSize(0.06);
  fr->GetXaxis()->SetTitleOffset(1.00);

  TLatex* latex = new TLatex();
  latex->SetTextSize(fontSize-0.003);
  latex->SetTextFont(42);
  latex->SetTextAlign(21);
  latex->SetNDC();
  latex->DrawLatex(0.35,0.93,"ALICE");
  latex->DrawLatex(0.35,0.87,kSystemEnergy);
  latex->DrawLatex(0.35,0.81,"V0S: (0-20%)-(60-100%)");
  //  latex->DrawLatex(0.40,0.79,Form("%s mult: (0-20%) - (60-100%)",method.Data()));
  TGraph* gHydro_f = Create_Hydro_f();
  TGraph* gHydro_b = Create_Hydro_b();
  TGraph* gAMPT_f = Create_AMPT_f();
  TGraph* gAMPT_b = Create_AMPT_b();
  gsyst_b->SetMarkerStyle(kOpenCircle);
  gsyst_f->SetMarkerStyle(kFullCircle);
  gsyst_b->SetMarkerSize(1.5);
  gsyst_f->SetMarkerSize(1.5);
  gsyst_b->SetMarkerColor(kRed);
  gsyst_f->SetMarkerColor(kRed);
  gsyst_b->SetLineColor(kRed-10);
  gsyst_f->SetLineColor(kRed);
  gsyst_b->SetFillColor(kRed-10);
  gsyst_f->SetFillColor(kRed);
  gsyst_b->SetFillStyle(1001);
  gsyst_f->SetFillStyle(1);
  gsyst_f->SetLineWidth(1);
  g_b->SetMarkerStyle(kOpenCircle);
  g_f->SetMarkerStyle(kFullCircle);
  g_b->SetMarkerSize(1.5);
  g_f->SetMarkerSize(1.5);
  g_b->SetMarkerColor(kRed);
  g_f->SetMarkerColor(kRed);
  g_b->SetLineColor(kRed);
  g_f->SetLineColor(kRed);
  g_b->SetLineWidth(3);
  g_f->SetLineWidth(3);
//  gHydro_b->Draw("lsame");
//  gHydro_f->Draw("lsame");
//  gAMPT_b->Draw("lsame");
//  gAMPT_f->Draw("lsame");
  TFile* f1 = new TFile("our_ampt.root");
  TGraphErrors* gama_b = (TGraphErrors* ) f1->Get("gama_b");
  TGraphErrors* gama_f = (TGraphErrors* ) f1->Get("gama_f");
  gama_b->SetLineColor(kBlue);
  gama_f->SetLineColor(kMagenta);
  TFile* f1 = new TFile("our_ampt_scaled.root");
  TGraphErrors* gama_b_scaled = (TGraphErrors* ) f1->Get("gama_b");
  TGraphErrors* gama_f_scaled = (TGraphErrors* ) f1->Get("gama_f");
  gama_b_scaled->SetLineColor(kBlue);
  gama_f_scaled->SetLineColor(kMagenta);
  TGraph* gama_b_shaded = GetShade(gama_b,gama_b_scaled,kBlue); 
  TGraph* gama_f_shaded = GetShade(gama_f,gama_f_scaled,kGreen);
  gama_b_shaded->SetLineWidth(0);
  gama_f_shaded->SetLineWidth(0);
  gama_b_shaded->SetFillStyle(3345);
  gama_f_shaded->SetFillStyle(3354);
  gama_b_shaded->Draw("fsame");
  gama_f_shaded->Draw("fsame");
//  gama_b->Draw("lsame");
//  gama_f->Draw("lsame");
//  gama_b_scaled->Draw("lsame");
//  gama_f_scaled->Draw("lsame");
  
  gsyst_b->Draw("2 same");
  gsyst_f->DrawClone("2 same");
  g_b->Draw("pz same");
  g_f->DrawClone("pz same");

//  TLegend* l = new TLegend(0.53,0.60,0.988,0.97);
  TLegend* l = new TLegend(0.60,0.75,0.988,0.97);
  l->SetMargin(0.2);
  l->SetFillColor(0);
  l->SetBorderSize(0);
//  l->AddEntry(gsyst_f,"v_{2}^{#mu}{2PC,sub}, Pb-going","PF");
//  l->AddEntry(gsyst_b,"v_{2}^{#mu}{2PC,sub}, p-going"),"PF";
//  l->AddEntry(gHydro_f,"v_{2}^{h}, Hydro, Pb-going");
//  l->AddEntry(gHydro_b,"v_{2}^{h}, Hydro, p-going");
//  l->AddEntry(gAMPT_f,"v_{2}^{h}, AMPT, Pb-going");
//  l->AddEntry(gAMPT_b,"v_{2}^{h}, AMPT, p-going");
//  l->AddEntry(gama_f,"v_{2}^{#mu}{2PC,sub}, AMPT, Pb-going");
//  l->AddEntry(gama_b,"v_{2}^{#mu}{2PC,sub}, AMPT, p-going");
  l->AddEntry(gsyst_f,"ALICE, Pb-going","PF");
  l->AddEntry(gsyst_b,"ALICE, p-going"),"PF";
  l->AddEntry(gama_f_shaded,"AMPT, Pb-going","F");
  l->AddEntry(gama_b_shaded,"AMPT, p-going","F");
  l->Draw("same");
  gPad->Print("eps/v2_final.eps");
  gPad->Print("eps/v2_final.pdf");
  gPad->Print("png/v2_final.png");
}



void draw_ratio_final(){
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadRightMargin(0.01);
  gStyle->SetPadTopMargin(0.01);
  gStyle->SetHatchesLineWidth(1);
  gStyle->SetHatchesSpacing(0.6);
  fontSize = 0.048;
  TFile* f = new TFile("ratio.root");
  f->ls();
  TGraphErrors* g_ratio     = (TGraphErrors*) f->Get("g_ratio");
  TGraphErrors* gsyst_ratio = (TGraphErrors*) f->Get("gsyst_ratio");
  if(gShiftToMeanPt){
    //FIXME what distribution should we use for the ratio?
    ShiftToMeanPt(g_ratio,kV0S,"bcde");
    ShiftToMeanPt(gsyst_ratio,kV0S,"bcde");
  }

  gsyst_ratio->SetLineWidth(3);
  gsyst_ratio->SetMarkerColor(kRed);
  gsyst_ratio->SetLineColor(kRed-10);
  gsyst_ratio->SetMarkerSize(1.5);
  gsyst_ratio->SetFillStyle(1001);
  gsyst_ratio->SetLineWidth(1);
  g_ratio->SetLineWidth(3);
  g_ratio->SetMarkerColor(kRed);
  g_ratio->SetLineColor(kRed);
  g_ratio->SetMarkerSize(1.5);
  gsyst_ratio->SetFillColor(kRed-10);
  TCanvas* c1 = new TCanvas("ratio","ratio",800,650);
  TH1F* fr = draw_frame(4,1.5);
  fr->GetYaxis()->SetTitleSize(0.065);
  fr->GetYaxis()->SetTitleOffset(0.85);
  fr->GetXaxis()->SetTitleSize(0.06);
  fr->GetXaxis()->SetTitleOffset(1.00);

  fr->SetMinimum(0.9);
  fr->SetMaximum(2.3);
  fr->GetYaxis()->SetTitle("v_{2}(Pb-going) / v_{2}(p-going)");
  TLine* l = new TLine(0,1,4,1);
  l->SetLineStyle(2);
  l->Draw();

  TFile* f2 = new TFile("our_ampt_scaled.root");
  TGraphErrors* gratio_scaled = (TGraphErrors* ) f2->Get("gratio");
  //gratio_scaled->SetLineColor(kMagenta);
  //gratio_scaled->Draw("lsame");

  TFile* f1 = new TFile("our_ampt.root");
  TGraphErrors* gratio = (TGraphErrors* ) f1->Get("gratio");
  //gratio->SetLineColor(kMagenta);
  //gratio->Draw("lsame");

  TGraph* gratio_shaded = GetShade(gratio,gratio_scaled,kBlue); 
  gratio_shaded->SetLineWidth(0);
  gratio_shaded->SetFillStyle(3345);

  TLatex* latex = new TLatex();
  latex->SetTextSize(fontSize-0.003);
  latex->SetTextFont(42);
  latex->SetTextAlign(21);
  latex->SetNDC();
  latex->DrawLatex(0.35,0.93,"ALICE");
  latex->DrawLatex(0.35,0.87,kSystemEnergy);
  latex->DrawLatex(0.35,0.81,"V0S: (0-20%)-(60-100%)");
  //  latex->DrawLatex(0.40,0.79,Form("%s mult: (0-20%) - (60-100%)",method.Data()));
  TGraph* gHydro_ratio = Create_Hydro_ratio();
  TGraph* gAMPT_ratio = Create_AMPT_ratio();
  gratio_shaded->Draw("fsame");
  gsyst_ratio->Draw("2 same");
  //gHydro_ratio->Draw("lsame");
  //gAMPT_ratio->Draw("lsame");
  //gratio->Draw("lsame");
  g_ratio->Draw("pz same");
//  g_ratio->Fit("pol0");
//  TLegend* leg = new TLegend(0.55,0.75,0.85,0.97);
  TLegend* leg = new TLegend(0.60,0.82,0.85,0.97);
  leg->SetMargin(0.3);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(gsyst_ratio,"ALICE","FP");
//  leg->AddEntry(gHydro_ratio,"Hydro");
//  leg->AddEntry(gAMPT_ratio,"AMPT");
//  leg->AddEntry(gratio,"AMPT v_{2}^{#mu}{2PC,sub}"); 
  leg->AddEntry(gratio_shaded,"AMPT","F"); 
  leg->Draw("same");
  gPad->Print("eps/ratio_final.eps");
  gPad->Print("eps/ratio_final.pdf");
  gPad->Print("png/ratio_final.png");
}

TGraph* Create_Hydro_f(){
  const Int_t nPoints = 8;
  Double_t x[nPoints] = {0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5};
  Double_t y[nPoints] = {0.006888005,0.0249338,0.0449684,0.067962,0.09024495,0.1193485,0.149601,0.1836515};
  TGraph* g = new TGraph(nPoints,x,y);
  g->SetLineColor(kBlack);
  g->SetLineWidth(3);
  g->SetLineStyle(0);
  g->SetFillColor(10);
  return g;
}

TGraph* Create_Hydro_b(){
  const Int_t nPoints = 8;
  Double_t x[8] = {0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5};
  Double_t y[8] = {0.004445615,0.0171062,0.030959,0.04539025,0.05888625,0.0739857,0.088069,0.1039095};
  TGraph* g = new TGraph(nPoints,x,y);
  g->SetLineColor(kBlack);
  g->SetLineWidth(3);
  g->SetLineStyle(10);
  g->SetFillColor(10);
  return g;
}

TGraph* Create_Hydro_ratio(){
  const Int_t nPoints = 8;
  Double_t x[nPoints] = {0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5};
  Double_t y_f[nPoints] = {0.006888005,0.0249338,0.0449684,0.067962,0.09024495,0.1193485,0.149601,0.1836515};
  Double_t y_b[nPoints] = {0.004445615,0.0171062,0.030959,0.04539025,0.05888625,0.0739857,0.088069,0.1039095};
  Double_t y[nPoints]= {0};
  for (Int_t i=0;i<nPoints;i++) y[i]=y_f[i]/y_b[i];
  TGraph* g = new TGraph(nPoints,x,y);
  g->SetLineColor(kBlack);
  g->SetLineWidth(3);
  g->SetLineStyle(0);
  g->SetFillColor(10);
  return g;
}


TGraph* Create_AMPT_f(){
  const Int_t nPoints = 8;
  Double_t x[nPoints] = {0.3,0.75,1.25,1.75,2.25,2.75,3.5,4.5};
  Double_t y[nPoints] = {0.0178254,0.0598114,0.0963362,0.11249,0.119999,0.121476,0.113305,0.112206};
  TGraph* g = new TGraph(nPoints,x,y);
  g->SetLineColor(kBlue);
  g->SetLineWidth(3);
  g->SetLineStyle(9);
  g->SetFillColor(10);
  return g;
}

TGraph* Create_AMPT_b(){
  const Int_t nPoints = 8;
  Double_t x[nPoints] = {0.3,0.75,1.25,1.75,2.25,2.75,3.5,4.5};
  Double_t y[nPoints] = {0.0136546,0.0445013,0.0701439,0.0864519,0.0868313,0.0903688,0.0892526,0.0962324};
  TGraph* g = new TGraph(nPoints,x,y);
  g->SetLineColor(kBlue);
  g->SetLineWidth(3);
  g->SetLineStyle(7);
  g->SetFillColor(10);
  return g;
}


TGraph* Create_AMPT_ratio(){
  const Int_t nPoints = 8;
  Double_t x[nPoints] = {0.3,0.75,1.25,1.75,2.25,2.75,3.5,4.5};
  Double_t y[nPoints] = {1.30545,1.34404,1.37341,1.30119,1.38198,1.34422,1.26948,1.16599};
  TGraph* g = new TGraph(nPoints,x,y);
  g->SetLineColor(kBlue);
  g->SetLineWidth(3);
  g->SetLineStyle(9);
  g->SetFillColor(10);
  return g;
}


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

void PadFor2DCorr(){
  gPad->SetPad(0, 0, 1, 1);
  gPad->SetLeftMargin(0.17);
  gPad->SetTopMargin(0.15);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.05);
  gPad->SetTheta(55);
  gPad->SetPhi(45);
}

void ShiftToMeanPt(TGraphErrors* g, Int_t im,TString period){
  Printf("\033[1;31m Shifting graph %s to <pt> \033[m",g->GetName());
  Bool_t isTkl=im>kNtrk && im<=kNtkl? 1 : 0;
  Int_t nbinst = isTkl ? nbins_muon_tkl : nbins_muon_trk;
  Double_t* binst  = isTkl ? bins_muon_tkl  : bins_muon_trk;

  TFile *fin=new TFile("muon_pt_dist.root");
  //old different <pt> value for different periods
  //  TString hname;
  //  if(isTkl)hname=Form("h_muon_pt_dist_tkl_%s",period.Data());
  //  else hname=Form("h_muon_pt_dist_trk_%s",period.Data());
  //new sum of bcde+f
  TString hname_bcde;
  if(isTkl)hname_bcde=Form("h_muon_pt_dist_tkl_bcde");
  else hname_bcde=Form("h_muon_pt_dist_trk_bcde");
  TH1F * hpt_bcde=fin->Get(hname_bcde);
  TString hname_f;
  if(isTkl)hname_f=Form("h_muon_pt_dist_tkl_f");
  else hname_f=Form("h_muon_pt_dist_trk_f");
  TH1F * hpt_f=fin->Get(hname_f);

  TH1F* hpt=(TH1F*)hpt_bcde->Clone("hpt");
  hpt->Add(hpt_f);

  Double_t eps=0.0001;
  for (Int_t ibin=0;ibin<nbinst;ibin++){
    hpt->GetXaxis()->SetRangeUser(binst[ibin]+eps,binst[ibin+1]-eps);
    Printf("%f [%f,%f] ->  %f",g->GetX()[ibin],binst[ibin],binst[ibin+1],hpt->GetMean());
    g->GetX()[ibin]=hpt->GetMean();
  }
}

TGraph* GetShade(TGraph* gmin, TGraph* gmax, Color_t color){
  Int_t n = gmin->GetN();
  Double_t xmin = 0.0001;
  Double_t xmax = 0.1;
  Double_t dx = (xmax-xmin)/n;
  TGraph *gshade = new TGraph(2*n);
  for (int i=0;i<n;i++){
    double x    = gmin->GetX()[i];
    double ymin = gmin->GetY()[i];
    double ymax = gmax->GetY()[i];
    gshade->SetPoint(i      ,x,ymin);
    gshade->SetPoint(2*n-i-1,x,ymax);
  }
  gshade->SetFillColor(color);
  gshade->SetLineColor(color);
  return gshade;
}
