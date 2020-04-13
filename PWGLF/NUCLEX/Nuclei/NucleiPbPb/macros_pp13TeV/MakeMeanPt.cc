#include <iostream>
#include "TSystem.h" //includi sempre
#include "TRandom.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THnSparse.h"
#include "TF1.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TPaveStats.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TFile.h"
#include "TList.h"
#include "TMath.h"
#include "TROOT.h"
#include "TStyle.h"
#include "src/Common.h"

bool useAntideuterons = true;

void MakeMeanPt() {

  //pp
  TFile *fileDeut_pp=new TFile(Form("%soMeanPt_ddbar.root",kBaseOutputDir.c_str()));
  TFile *fileProt_pp=new TFile(Form("%soGetYieldMeanProtonNicolo.root",kBaseOutputDir.c_str()));

  TGraphErrors **grDeut_pp=new TGraphErrors*[2];
  TGraphErrors **grProt_pp=new TGraphErrors*[3];
  grDeut_pp[0]=(TGraphErrors *)fileDeut_pp->Get(Form("ptmean_stat_%s%s","d","#bar{d}"));
  grDeut_pp[1]=(TGraphErrors *)fileDeut_pp->Get(Form("ptmean_sys_%s%s","d","#bar{d}"));

  grProt_pp[0]=(TGraphErrors *)fileProt_pp->Get(Form("ptmean_stat_p_nicolo"));
  grProt_pp[1]=(TGraphErrors *)fileProt_pp->Get(Form("ptmean_uncor_p_nicolo"));//Total systematic (cor+uncor, despite the name)
  grProt_pp[2]=(TGraphErrors *)fileProt_pp->Get(Form("ptmean_uncorReal_p_nicolo"));//uncor systematic

  TGraphErrors *grDeutCoal_pp=(TGraphErrors *)fileDeut_pp->Get(Form("meanPt_%s%s_coalescence","d","#bar{d}"));

  TFile *fileDeutBW_pp=new TFile(Form("%sptmean_blastwave_deuteron_vsColls_rms.root",kBaseOutputDir.c_str()));
  TGraphAsymmErrors *grDeutBW_pp=(TGraphAsymmErrors *)fileDeutBW_pp->Get(Form("ptMeanBW"))->Clone();

  TFile *fileProtBW_pp=new TFile(Form("%sptmean_blastwave_proton_vsColls_rms.root",kBaseOutputDir.c_str()));
  TGraphAsymmErrors *grProtBW_pp=(TGraphAsymmErrors *)fileProtBW_pp->Get(Form("ptMeanBW"))->Clone();

  //Rimuovo lo 0-100%:
  for(Int_t j=0;j<3;j++) {
    grProt_pp[j]->RemovePoint(3);
  }

  //pPb
  TFile *file_pPb=new TFile(Form("%smeanPtDeuteron.root",kBaseOutputDir.c_str()));
  TGraphErrors **grDeut_pPb=new TGraphErrors*[3];
  TGraphErrors **grProt_pPb=new TGraphErrors*[3];
  grDeut_pPb[0]=(TGraphErrors *)((TCanvas *)file_pPb->Get("cCanv"))->GetListOfPrimitives()->At(18);//stat
  grDeut_pPb[1]=(TGraphErrors *)((TCanvas *)file_pPb->Get("cCanv"))->GetListOfPrimitives()->At(16);//sys tot
  grDeut_pPb[2]=(TGraphErrors *)((TCanvas *)file_pPb->Get("cCanv"))->GetListOfPrimitives()->At(17);//sys uncor

  grProt_pPb[0]=(TGraphErrors *)((TCanvas *)file_pPb->Get("cCanv"))->GetListOfPrimitives()->At(13);//stat
  grProt_pPb[1]=(TGraphErrors *)((TCanvas *)file_pPb->Get("cCanv"))->GetListOfPrimitives()->At(5);//sys tot
  grProt_pPb[2]=(TGraphErrors *)((TCanvas *)file_pPb->Get("cCanv"))->GetListOfPrimitives()->At(9);//sys uncor

  for(Int_t i=0;i<grDeut_pPb[0]->GetN();i++) {
    grDeut_pPb[0]->SetPointError(i,0,grDeut_pPb[0]->GetErrorY(i));
  }
  TGraphErrors *grDeutCoal_pPb=(TGraphErrors *)((TCanvas *)file_pPb->Get("cCanv"))->FindObject("grDeuteronCoal");

  TGraphErrors *grDeutBW_pPb=(TGraphErrors *)((TCanvas *)file_pPb->Get("cCanv"))->FindObject("grDeuteronBW");
  TGraphErrors *grProtBW_pPb=(TGraphErrors *)((TCanvas *)file_pPb->Get("cCanv"))->FindObject("grProtonBW");

  //pp 13 TeV
  const char*  particlename = (useAntideuterons) ? "Antideuterons" : "Deuterons";
  TFile *fileDeut_pp13=new TFile(Form("%s%s_mean_pt.root",kBaseOutputDir.c_str(),particlename));
  TFile *fileProt_pp13=new TFile(Form("%sWellDoneProtonPt.root",kBaseOutputDir.c_str()));

  TGraphErrors **grDeut_pp13=new TGraphErrors*[2];
  TGraphErrors **grProt_pp13=new TGraphErrors*[3];
  grDeut_pp13[0]=(TGraphErrors *)fileDeut_pp13->Get(Form("%s_meanpt_stat",particlename));
  grDeut_pp13[1]=(TGraphErrors *)fileDeut_pp13->Get(Form("%s_meanpt_syst",particlename));
  grProt_pp13[0]=(TGraphErrors *)fileProt_pp13->Get("obj4");
  grProt_pp13[0]->SetName("proton_meanpt_pp13tev_stat");
  grProt_pp13[1]=(TGraphErrors *)fileProt_pp13->Get("obj1");
  grProt_pp13[1]->SetName("proton_meanpt_pp13tev_syst");

  Int_t lineWidth=1;

  //pp
  for(Int_t j=0;j<2;j++) {
    grDeut_pp[j]->SetMarkerStyle(kFullCircle);
    grDeut_pp[j]->SetMarkerColor(kGreen+3);
    grDeut_pp[j]->SetLineColor(grDeut_pp[j]->GetMarkerColor());
    grDeut_pp[j]->SetMarkerSize(1.7);
    grDeut_pp[j]->SetLineWidth(lineWidth);
  }

  for(Int_t j=0;j<3;j++) {
    grProt_pp[j]->SetMarkerStyle(kFullCross);
    grProt_pp[j]->SetMarkerColor(kGreen+3);
    grProt_pp[j]->SetLineColor(grProt_pp[j]->GetMarkerColor());
    grProt_pp[j]->SetMarkerSize(1.7);
    grProt_pp[j]->SetLineWidth(lineWidth);
  }
  grProt_pp[2]->SetFillStyle(3002);

  grDeutCoal_pp->SetFillColor(kGreen-6);

  grDeutCoal_pp->SetMarkerColor(grDeutCoal_pp->GetFillColor());
  grDeutCoal_pp->SetLineColor(grDeutCoal_pp->GetFillColor());

  grDeutCoal_pp->SetFillStyle(1001);

  grDeutBW_pp->SetFillColor(kGreen+2);
  grDeutBW_pp->SetMarkerColor(grDeutBW_pp->GetFillColor());
  grDeutBW_pp->SetLineColor(grDeutBW_pp->GetFillColor());
  grDeutBW_pp->SetFillStyle(3006);

  grProtBW_pp->SetFillColor(kGreen+2);
  grProtBW_pp->SetMarkerColor(grProtBW_pp->GetFillColor());
  grProtBW_pp->SetLineColor(grProtBW_pp->GetFillColor());
  grProtBW_pp->SetFillStyle(grDeutBW_pp->GetFillStyle());//3012);

  //pPb
  for(Int_t j=0;j<3;j++) {
    grDeut_pPb[j]->SetMarkerStyle(kOpenCircle);
    grDeut_pPb[j]->SetMarkerColor(kBlue+1);
    grDeut_pPb[j]->SetLineColor(grDeut_pPb[j]->GetMarkerColor());
    grDeut_pPb[j]->SetMarkerSize(1.6);
    grDeut_pPb[j]->SetLineWidth(lineWidth);
  }
  grDeut_pPb[2]->SetFillStyle(3002);

  for(Int_t j=0;j<3;j++) {
    grProt_pPb[j]->SetMarkerStyle(kOpenCross);
    grProt_pPb[j]->SetMarkerColor(kBlue+1);
    grProt_pPb[j]->SetLineColor(grProt_pPb[j]->GetMarkerColor());
    grProt_pPb[j]->SetFillColor(grProt_pPb[j]->GetMarkerColor());
    grProt_pPb[j]->SetMarkerSize(1.7);
    grProt_pPb[j]->SetLineWidth(lineWidth);
  }
  grProt_pPb[2]->SetFillStyle(3002);

  grDeutCoal_pPb->SetFillStyle(3020);

  grProtBW_pPb->SetFillStyle(grDeutBW_pPb->GetFillStyle());//3005);
  grDeutBW_pPb->SetFillColor(grDeut_pPb[1]->GetMarkerColor());
  grDeutBW_pPb->SetMarkerColor(grDeutBW_pPb->GetFillColor());
  grDeutBW_pPb->SetLineColor(grDeutBW_pPb->GetFillColor());

  grProtBW_pPb->SetFillStyle(grDeutBW_pPb->GetFillStyle());//3544);
  grProtBW_pPb->SetFillColor(grProt_pPb[1]->GetMarkerColor());
  grProtBW_pPb->SetMarkerColor(grProtBW_pPb->GetFillColor());
  grProtBW_pPb->SetLineColor(grProtBW_pPb->GetFillColor());

  //pp
  for(Int_t j=0;j<2;j++) {
    grDeut_pp13[j]->SetMarkerStyle(kFullCircle);
    grDeut_pp13[j]->SetMarkerColor(kOrange-3);
    grDeut_pp13[j]->SetLineColor(grDeut_pp13[j]->GetMarkerColor());
    grDeut_pp13[j]->SetMarkerSize(1.7);
    grDeut_pp13[j]->SetLineWidth(lineWidth);
  }

  for(Int_t j=0;j<2;j++) {
    grProt_pp13[j]->SetMarkerStyle(kFullCross);
    grProt_pp13[j]->SetMarkerColor(kOrange-3);
    grProt_pp13[j]->SetLineColor(grProt_pp13[j]->GetMarkerColor());
    grProt_pp13[j]->SetMarkerSize(1.7);
    grProt_pp13[j]->SetLineWidth(lineWidth);
  }

  TH1F *h=new TH1F("h","h",10000,0,100);
  h->GetXaxis()->SetTitle("#LTd#it{N}_{ch}/d#it{#eta}_{lab}#GT_{|#it{#eta}_{lab}| < 0.5}");
  h->GetXaxis()->SetRangeUser(1,49);
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetXaxis()->SetTitleOffset(1.2);

  h->GetYaxis()->SetTitle("#LT#it{p}_{T}#GT (GeV/#it{c})");
  h->GetYaxis()->SetRangeUser(0.15,2.15);//0.25,2.45
  h->GetYaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleOffset(0.9);

  const Int_t Nlg=3;
  TLegend *lg[Nlg];
  for(Int_t i=0;i<Nlg;i++) {
    if(i==0) lg[i]=new TLegend(0.75,0.62,0.95,0.79,"","brNDC");
    else if(i==1) lg[i]=new TLegend(0.75,0.32,0.93,0.49,"","brNDC");
    else if(i==2) lg[i]=new TLegend(0.75,0.10,0.95,0.19,"","brNDC");
    lg[i]->SetNColumns(1);
    lg[i]->SetBorderSize(0);
    lg[i]->SetFillStyle(0);//0(empty) //1001(full)
    lg[i]->SetFillColor(kWhite);
    lg[i]->SetTextFont(42);
    lg[i]->SetTextSize(0.028);
  }
  
  //pp
  lg[0]->AddEntry(grProt_pp[1],"p+#bar{p}","P");
  lg[0]->AddEntry(grDeut_pp[1],"d+#bar{d}","P");
  //lg[0]->AddEntry(grProtBW_pp,"Blast-Wave p+#bar{p}","F");
  lg[0]->AddEntry(grDeutBW_pp,"Blast-Wave d+#bar{d}","F");
  lg[0]->AddEntry(grDeutCoal_pp,"Coalescence d+#bar{d}","F");

  //pPb
  lg[1]->AddEntry(grProt_pPb[1],"p+#bar{p}","P");
  lg[1]->AddEntry(grDeut_pPb[1],"d+#bar{d}","P");
  //lg[1]->AddEntry(grProtBW_pPb,"Blast-Wave p+#bar{p}","F");
  lg[1]->AddEntry(grDeutBW_pPb,"Blast-Wave d+#bar{d}","F");
  lg[1]->AddEntry(grDeutCoal_pPb,"Coalescence d+#bar{d}","F");

  //pp13TeV
  if(useAntideuterons){
    lg[2]->AddEntry(grDeut_pp13[1],"#bar{d}","P");
  } else {
    lg[2]->AddEntry(grDeut_pp13[1],"d","P");
  }
  lg[2]->AddEntry(grProt_pp13[1],"p+#bar{p}","P");

  printf("ciao\n");
  

  const Int_t Npv=7;
  TPaveText *pv[Npv];
  for(Int_t i=0;i<Npv;i++) {

    if(i==0) pv[i]=new TPaveText(0.755,0.89,0.85,0.96,"brNDC");//0.17,0.78,0.45,0.85

    else if(i==1) pv[i]=new TPaveText(0.75,0.83,0.95,0.84,"brNDC");//0.23
    else if(i==4) pv[i]=new TPaveText(0.75,0.79,0.95,0.81,"brNDC");//0.23

    else if(i==2) pv[i]=new TPaveText(0.75,0.54,0.93,0.55,"brNDC");//0.23*/
    else if(i==5) pv[i]=new TPaveText(0.75,0.50,0.93,0.52,"brNDC");//0.23*/

    else if(i==3) pv[i]=new TPaveText(0.75,0.25,0.95,0.26,"brNDC");//0.23
    else if(i==6) pv[i]=new TPaveText(0.75,0.21,0.93,0.22,"brNDC");//0.23*/
    
    pv[i]->SetFillStyle(0);
    pv[i]->SetBorderSize(0);
    pv[i]->SetTextFont(42);
    pv[i]->SetTextAlign(12);
    if(i!=0)pv[i]->SetTextSize(0.030);
    if(i==3 || i==4) pv[i]->SetTextSize(0.025);
  }
  pv[0]->AddText("ALICE");
  //pv[0]->AddText("#bf{ALICE Preliminary}");
  pv[1]->AddText("pp, #sqrt{#it{s}} = 7 TeV");//"ALICE");
  pv[2]->AddText("p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV");//"ALICE");
  pv[3]->AddText("pp, #sqrt{#it{s}} = 13 TeV");//"ALICE");

  pv[4]->AddText("V0M Multiplicity Classes");
  pv[5]->AddText("V0A Multiplicity Classes (Pb-side)");
  pv[6]->AddText("V0M Multiplicity Classes");
  

  const Int_t Ncv=1;
  TCanvas *cv[Ncv];
  for(Int_t i=0;i<Ncv;i++) {
    cv[i]=new TCanvas();
    cv[i]->SetWindowSize(1050,700);
    cv[i]->ToggleEditor();
    cv[i]->ToggleToolBar();
    //cv[i]->Divide(1,1);
    cv[i]->SetTicks();
    cv[i]->SetRightMargin(0.25);
    cv[i]->SetLeftMargin(0.1);
    cv[i]->SetBottomMargin(0.14);
    cv[i]->SetTopMargin(0.05);
  }
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  h->Draw();

  grDeutBW_pp->Draw("3,same");
  grDeutBW_pPb->Draw("3,same");
  //grProtBW_pp->Draw("3,same");
  //grProtBW_pPb->Draw("3,same");

  grDeutCoal_pp->Draw("3,same");
  grDeutCoal_pPb->Draw("3,same");

  grDeut_pp[1]->Draw("E2,same");
  grDeut_pp[0]->Draw("PZ,same");
  grProt_pp[2]->Draw("E2,same");
  grProt_pp[1]->Draw("E2,same");
  grProt_pp[0]->Draw("PZ,same");

  grDeut_pPb[2]->Draw("E2,same");
  grDeut_pPb[1]->Draw("E2,same");
  grDeut_pPb[0]->Draw("PZ,same");
  grProt_pPb[2]->Draw("E2,same");
  grProt_pPb[1]->Draw("E2,same");
  grProt_pPb[0]->Draw("PZ,same");

  grDeut_pp13[1]->Draw("E2,same");
  grDeut_pp13[0]->Draw("PZ,same");
  grProt_pp13[1]->Draw("E2,same");
  grProt_pp13[0]->Draw("PZ,same");

  if(1) {
    lg[0]->Draw("same");
    lg[1]->Draw("same");
    lg[2]->Draw("same");
    pv[0]->Draw("same");
    pv[1]->Draw("same");
    pv[2]->Draw("same");
    pv[3]->Draw("same");
    pv[4]->Draw("same");
    pv[5]->Draw("same");
    pv[6]->Draw("same");
  }

  // TPad *pd=new TPad("pd","pd",0.34,0.16,0.72,0.52);
  // pd->SetFillStyle(0);
  // pd->SetRightMargin(0.02);
  // pd->SetLeftMargin(0.14);
  // pd->SetBottomMargin(0.24);
  // pd->SetTopMargin(0.05);
  // pd->SetTicks();
  // pd->Draw("same");
  // pd->cd();

  // TH1F *h2=(TH1F *)h->Clone();
  // //h2->GetXaxis()->SetRangeUser(1,49);
  // h2->GetXaxis()->SetTitleSize(0.09);
  // h2->GetXaxis()->SetLabelSize(0.09);
  // h2->GetXaxis()->SetTitleOffset(1.2);
  // h2->GetXaxis()->SetNdivisions(505);

  // h2->GetYaxis()->SetRangeUser(0.55,1.35);
  // h2->GetYaxis()->SetTitleSize(0.09);
  // h2->GetYaxis()->SetLabelSize(0.09);
  // h2->GetYaxis()->SetTitleOffset(0.7);
  // h2->GetYaxis()->SetNdivisions(505);

  // TGraphErrors **grProt_pp_2=new TGraphErrors*[3];
  // TGraphErrors **grProt_pPb_2=new TGraphErrors*[3];
  // for(Int_t i=0;i<3;i++) {
  //   grProt_pp_2[i]=(TGraphErrors *)grProt_pp[i]->Clone();
  //   grProt_pPb_2[i]=(TGraphErrors *)grProt_pPb[i]->Clone();
  // }
  // TGraphAsymmErrors *grProtBW_pp_2=(TGraphAsymmErrors *)grProtBW_pp->Clone();
  // TGraphErrors *grProtBW_pPb_2=(TGraphErrors *)grProtBW_pPb->Clone();

  // const Int_t Nlg2=1;
  // TLegend *lg2[Nlg2];
  // for(Int_t i=0;i<Nlg2;i++) {
  //   lg2[i]=new TLegend(0.4,0.3,0.7,0.55,"","brNDC");
  //   lg2[i]->SetNColumns(1);
  //   lg2[i]->SetBorderSize(0);
  //   lg2[i]->SetFillStyle(0);//0(empty) //1001(full)
  //   lg2[i]->SetFillColor(kWhite);
  //   lg2[i]->SetTextFont(42);
  //   lg2[i]->SetTextSize(0.080);
  // }

  // lg2[0]->AddEntry(grProtBW_pp,"pp, Blast-Wave p+#bar{p}","F");
  // lg2[0]->AddEntry(grProtBW_pPb,"p-Pb, Blast-Wave p+#bar{p}","F");

  // const Int_t Npv2=1;
  // TPaveText *pv2[Npv2];
  // for(Int_t i=0;i<Npv2;i++) {
  //   pv2[i]=new TPaveText(0.17,0.78,0.42,0.9,"brNDC");//0.17,0.78,0.45,0.85
  //   pv2[i]->SetFillStyle(0);
  //   pv2[i]->SetBorderSize(0);
  //   pv2[i]->SetTextFont(42);
  //   pv2[i]->SetTextAlign(12);
  //   pv2[i]->SetTextSize(0.080);
  // }
  // pv2[0]->AddText("(anti-)protons");

  // h2->Draw("same");

  // grProtBW_pPb_2->Draw("3,same");
  // grProtBW_pp_2->Draw("3,same");

  // grProt_pp_2[2]->Draw("E2,same");
  // grProt_pp_2[1]->Draw("E2,same");
  // grProt_pp_2[0]->Draw("PZ,same");

  // grProt_pPb_2[2]->Draw("E2,same");
  // grProt_pPb_2[1]->Draw("E2,same");
  // grProt_pPb_2[0]->Draw("PZ,same");
  // lg2[0]->Draw("same");
  // pv2[0]->Draw("same");

  cv[0]->cd(1)->Update();
  cv[0]->Print(Form("%sptmean_%s_all.png",kBaseOutputDir.data(),"d"));
  TFile file_out(Form("%scustom_mean_pt.root",kBaseOutputDir.data()),"recreate");
  cv[0]->Write();
  //cv[0]->Print(Form("ptmean_%s_all.root","d"));

  return;

}
