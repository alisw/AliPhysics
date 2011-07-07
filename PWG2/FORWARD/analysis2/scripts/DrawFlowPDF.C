#include <TFile.h>
#include <TProfile.h>
#include <TList.h>
#include <TError.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TMath.h>
#include <TLatex.h>

// Data members
const char* pdfName = "Flow.pdf";

TCanvas* SetupCanvas(TString name)
{
  TCanvas* c = new TCanvas("c","c",640,960);

  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetBorderSize(0);
  c->cd();

  if (name.Contains("Monitoring")) return c;

  TPad* p1 = new TPad("p1", "p1", 0, 0.5, 1.0, 1.0, 0, 0);
  p1->SetTopMargin(0.1);
  p1->SetBottomMargin(0.1);
  p1->SetRightMargin(0.05);
  p1->SetGridx();
  p1->SetTicks(1, 1);
  p1->Draw();
  p1->cd();

  TProfile* p = new TProfile(name, "p", 48, -6, 6);
  p->GetYaxis()->SetRangeUser(0, 0.083);
  p->GetXaxis()->SetTitleFont(22);
  p->GetXaxis()->SetLabelFont(22);
  p->GetYaxis()->SetTitleFont(22);
  p->GetYaxis()->SetLabelFont(22);
  p->SetXTitle("#eta");
  p->SetYTitle("v_{2}");
  p->Draw();

  TLatex* pt = new TLatex(.12, .87, "FMD Preliminary");
  pt->SetNDC();
  pt->SetTextFont(22);
  pt->SetTextSize(0.07);
  pt->SetTextColor(kRed+1);
  pt->SetTextAlign(13);
  pt->Draw();

  c->cd();

  TPad* p2 = new TPad("p2", "p2", 0, 0.25, 1.0, 0.5, 0, 0);
  p2->SetTopMargin(0.);
  p2->SetRightMargin(0.05);
  p2->SetBottomMargin(0.2);
  p2->SetGridx();
  p2->SetGridy();
  p2->SetTicks(1,1);
  p2->Draw();
  p2->cd();

  TProfile* pp = new TProfile(Form("%s_diff", name.Data()), "p", 48, -6, 6);
  pp->GetXaxis()->SetTitleFont(22);
  pp->GetXaxis()->SetLabelFont(22);
  pp->GetYaxis()->SetTitleFont(22);
  pp->GetYaxis()->SetLabelFont(22);
  pp->SetXTitle("#eta");
  if ( name.Contains("MC")) {
    pp->SetYTitle("measured / truth");
    pp->GetYaxis()->SetRangeUser(0.51, 1.15);
  }
  if (!name.Contains("MC")) {
    pp->SetYTitle("v_{2}{2} / v_{2}{4}");
    pp->GetYaxis()->SetRangeUser(0.98, 1.17);
  }
  pp->GetXaxis()->SetTitleSize(2.*p->GetXaxis()->GetTitleSize());
  pp->GetXaxis()->SetLabelSize(2.*p->GetXaxis()->GetLabelSize());
  pp->GetYaxis()->SetTitleSize(2.*p->GetYaxis()->GetTitleSize());
  pp->GetYaxis()->SetLabelSize(2.*p->GetYaxis()->GetLabelSize());
  pp->GetYaxis()->SetTitleOffset(0.5);
  pp->Draw();

  TF1* oneLine = new TF1("oneLine","1",-6,6);
  oneLine->SetLineWidth(0.5);
  oneLine->Draw("same");

  if (!name.Contains("MC")) return c;

  c->cd();

  TPad* p3 = new TPad("p3", "p3", 0, 0, 1.0, 0.25, 0, 0);
  p3->SetTopMargin(0.);
  p3->SetRightMargin(0.05);
  p3->SetBottomMargin(0.2);
  p3->SetGridx();
  p3->SetGridy();
  p3->SetTicks(1,1);
  p3->Draw();
  p3->cd();

  TProfile* ppp = new TProfile(Form("%s_trdiff", name.Data()), "p", 48, -6, 6);
  ppp->GetXaxis()->SetTitleFont(22);
  ppp->GetXaxis()->SetLabelFont(22);
  ppp->GetYaxis()->SetTitleFont(22);
  ppp->GetYaxis()->SetLabelFont(22);
  ppp->GetYaxis()->SetRangeUser(0.66, 1.19);
  ppp->SetXTitle("#eta");
  ppp->SetYTitle("measured / track ref");
  ppp->GetXaxis()->SetTitleSize(2.*p->GetXaxis()->GetTitleSize());
  ppp->GetXaxis()->SetLabelSize(2.*p->GetXaxis()->GetLabelSize());
  ppp->GetYaxis()->SetTitleSize(2.*p->GetYaxis()->GetTitleSize());
  ppp->GetYaxis()->SetLabelSize(2.*p->GetYaxis()->GetLabelSize());
  ppp->GetYaxis()->SetTitleOffset(0.5);
  ppp->Draw();

  oneLine->Draw("same");

  return c;
}

void MakeFmdAndSpdPlots(TFile* f) 
{
  TList* qList = static_cast<TList*>(f->Get("FlowResults/QCumulants"));
  TList* aList = static_cast<TList*>(qList->FindObject("FMD"));
  TList* fmdList = static_cast<TList*>(aList->FindObject("v2"));
  aList = static_cast<TList*>(qList->FindObject("SPD"));
  TList* spdList = static_cast<TList*>(aList->FindObject("v2"));

  TPad*     pad = 0;
   
  Int_t y[11] = { 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 100 };
  for (Int_t c = -2; c < 10; c++) {
    TProfile* fmd2DiffHist = 0;
    TProfile* spd2DiffHist = 0;
    TProfile* fmd4DiffHist = 0;
    TProfile* spd4DiffHist = 0;
    TH1D*       fmdDiff = 0;
    TH1D*       spdDiff = 0;
    TH1D*       fmd = 0;
    TH1D*       spd = 0;
 
     if (c == -2) {
      fmd2DiffHist = (TProfile*)fmdList->FindObject("hQ2Cumulant2DiffFlowFMD_mb");
      spd2DiffHist = (TProfile*)spdList->FindObject("hQ2Cumulant2DiffFlowSPD_mb");
      fmd4DiffHist = (TProfile*)fmdList->FindObject("hQ2Cumulant4DiffFlowFMD_mb");
      spd4DiffHist = (TProfile*)spdList->FindObject("hQ2Cumulant4DiffFlowSPD_mb");
    }
    else if (c == -1) {
      fmd2DiffHist = (TProfile*)fmdList->FindObject("hQ2Cumulant2DiffFlowFMD_0_40");
      spd2DiffHist = (TProfile*)spdList->FindObject("hQ2Cumulant2DiffFlowSPD_0_40");
      fmd4DiffHist = (TProfile*)fmdList->FindObject("hQ2Cumulant4DiffFlowFMD_0_40");
      spd4DiffHist = (TProfile*)spdList->FindObject("hQ2Cumulant4DiffFlowSPD_0_40");
    }
    else {
      fmd2DiffHist = (TProfile*)fmdList->FindObject(Form("hQ2Cumulant2DiffFlowFMD_%d_%d", y[c], y[c+1]));
      spd2DiffHist = (TProfile*)spdList->FindObject(Form("hQ2Cumulant2DiffFlowSPD_%d_%d", y[c], y[c+1]));
      fmd4DiffHist = (TProfile*)fmdList->FindObject(Form("hQ2Cumulant4DiffFlowFMD_%d_%d", y[c], y[c+1]));
      spd4DiffHist = (TProfile*)spdList->FindObject(Form("hQ2Cumulant4DiffFlowSPD_%d_%d", y[c], y[c+1]));
    }
    
    TString name = "v_{2}{2} and v_{2}{4} ";
    if (c == -2)      name += "Minimum Bias";
    else if (c == -1) name += "0-40% Centrality"; 
    else {
                      name += TString::Format("%d-%d", y[c], y[c+1]);
                      name += "% Centrality";
    }

    TCanvas* c1 = SetupCanvas(name);
    pad = (TPad*)c1->GetPrimitive("p1");
    pad->cd();
    
    TLatex* tit = new TLatex(0.10, 0.95, name);
    tit->SetNDC();
    tit->SetTextFont(132);
    tit->SetTextSize(0.07);
    tit->Draw();

    Int_t nEvs = (Int_t)fmd2DiffHist->GetBinContent(fmd2DiffHist->GetNbinsX()+1);
    TLatex* et = new TLatex(.93, .87, Form("%d events", nEvs));
    et->SetNDC();
    et->SetTextFont(132);
    et->SetTextAlign(33);
    et->Draw();

    fmd2DiffHist->SetMarkerStyle(22);
    spd2DiffHist->SetMarkerStyle(23);
    fmd4DiffHist->SetMarkerStyle(22);
    spd4DiffHist->SetMarkerStyle(23);

    fmd2DiffHist->SetMarkerColor(kRed);
    spd2DiffHist->SetMarkerColor(kRed);
    fmd4DiffHist->SetMarkerColor(kBlue);
    spd4DiffHist->SetMarkerColor(kBlue);

    fmd2DiffHist->SetLineColor(kRed);
    spd2DiffHist->SetLineColor(kRed);
    fmd4DiffHist->SetLineColor(kBlue);
    spd4DiffHist->SetLineColor(kBlue);

    spd2DiffHist->Draw("same e1");
    fmd2DiffHist->Draw("same e1");
    spd4DiffHist->Draw("same e1");
    fmd4DiffHist->Draw("same e1");

    TLegend* l = new TLegend(.12,.40,.35,.80);
//    l->SetNColumns(2);
    l->SetFillColor(0);
    l->SetFillStyle(0);
    l->SetBorderSize(0);
    l->SetTextFont(132);
    l->AddEntry(fmd2DiffHist, "FMD v_{2}{2}", "lpe");
    l->AddEntry(spd2DiffHist, "SPD v_{2}{2}", "lpe");
    l->AddEntry(fmd4DiffHist, "FMD v_{2}{4}", "lpe");
    l->AddEntry(spd4DiffHist, "SPD v_{2}{4}", "lpe");
    l->Draw();

    pad = (TPad*)c1->GetPrimitive("p2");
    pad->cd();

    fmdDiff = fmd2DiffHist->ProjectionX("fmdDiff", "e");
    fmd = fmd4DiffHist->ProjectionX("fmd4Diff", "e");
    fmdDiff->Divide(fmd);
    fmdDiff->SetLineColor(kBlue);
    fmdDiff->Draw("same e1");

    spdDiff = spd2DiffHist->ProjectionX("spdDiff", "e");
    spd = spd4DiffHist->ProjectionX("spd4Diff", "e");
    spdDiff->Divide(spd);
    spdDiff->SetLineColor(kRed);
    spdDiff->Draw("same e1");

    name.Prepend("Title:");
    cout<<name<<endl;
    c1->Print(pdfName, name);

    delete fmd;
    delete spd;
    delete fmdDiff;
    delete spdDiff;
    delete c1;
   } // end of c
}

void Make2ParticlePlots(TFile* f) 
{
  TList* qList = static_cast<TList*>(f->Get("FlowResults/QCumulants"));
  TList* aList = static_cast<TList*>(qList->FindObject("MC"));
  TList* mcList = static_cast<TList*>(aList->FindObject("v2"));
  aList = static_cast<TList*>(qList->FindObject("FMDTR"));
  TList* fmdtrList = static_cast<TList*>(aList->FindObject("v2"));
  aList = static_cast<TList*>(qList->FindObject("SPDTR"));
  TList* spdtrList = static_cast<TList*>(aList->FindObject("v2"));
  aList = static_cast<TList*>(qList->FindObject("FMD"));
  TList* fmdList = static_cast<TList*>(aList->FindObject("v2"));
  aList = static_cast<TList*>(qList->FindObject("SPD"));
  TList* spdList = static_cast<TList*>(aList->FindObject("v2"));

  TProfile2D* mcProfile = (TProfile2D*)qList->FindObject("pMCTruth");
  TPad*       pad = 0;

  Int_t y[11] = { 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 100 };
  for (Int_t c = -2; c < 10; c++) {
    TProfile*   mcDiffHist = 0;
    TProfile*   fmdtrDiffHist = 0;
    TProfile*   spdtrDiffHist = 0;
    TProfile*   fmdDiffHist = 0;
    TProfile*   spdDiffHist = 0;
    TH1D*       mcTruth = 0;
    TH1D*       fmdDiff = 0;
    TH1D*       spdDiff = 0;
    TH1D*       fmdTrDiff = 0;
    TH1D*       spdTrDiff = 0;

    if (c == -2) {
      mcDiffHist = (TProfile*)mcList->FindObject("hQ2Cumulant2DiffFlowMC_mb");
      fmdtrDiffHist = (TProfile*)fmdtrList->FindObject("hQ2Cumulant2DiffFlowFMDTR_mb");
      spdtrDiffHist = (TProfile*)spdtrList->FindObject("hQ2Cumulant2DiffFlowSPDTR_mb");
      fmdDiffHist = (TProfile*)fmdList->FindObject("hQ2Cumulant2DiffFlowFMD_mb");
      spdDiffHist = (TProfile*)spdList->FindObject("hQ2Cumulant2DiffFlowSPD_mb");
      mcTruth = (TH1D*)mcProfile->ProjectionX("mcTruth", 1, 9, "E");
      mcTruth->Scale(1./9.);
    }
    else if (c == -1) {
      mcDiffHist = (TProfile*)mcList->FindObject("hQ2Cumulant2DiffFlowMC_0_40");
      fmdtrDiffHist = (TProfile*)fmdtrList->FindObject("hQ2Cumulant2DiffFlowFMDTR_0_40");
      spdtrDiffHist = (TProfile*)spdtrList->FindObject("hQ2Cumulant2DiffFlowSPDTR_0_40");
      fmdDiffHist = (TProfile*)fmdList->FindObject("hQ2Cumulant2DiffFlowFMD_0_40");
      spdDiffHist = (TProfile*)spdList->FindObject("hQ2Cumulant2DiffFlowSPD_0_40");
      mcTruth = (TH1D*)mcProfile->ProjectionX("mcTruth", 1, 5, "E");
      mcTruth->Scale(1./5.);
    }
    else {
      mcDiffHist = (TProfile*)mcList->FindObject(Form("hQ2Cumulant2DiffFlowMC_%d_%d", y[c], y[c+1]));
      fmdtrDiffHist = (TProfile*)fmdtrList->FindObject(Form("hQ2Cumulant2DiffFlowFMDTR_%d_%d", y[c], y[c+1]));
      spdtrDiffHist = (TProfile*)spdtrList->FindObject(Form("hQ2Cumulant2DiffFlowSPDTR_%d_%d", y[c], y[c+1]));
      fmdDiffHist = (TProfile*)fmdList->FindObject(Form("hQ2Cumulant2DiffFlowFMD_%d_%d", y[c], y[c+1]));
      spdDiffHist = (TProfile*)spdList->FindObject(Form("hQ2Cumulant2DiffFlowSPD_%d_%d", y[c], y[c+1]));
      mcTruth = (TH1D*)mcProfile->ProjectionX("mcTruth", c+1, c+1, "E");
    }

    TString name = "v_{2}{2} MC ";
    if (c == -2)      name += "Minimum Bias";
    else if (c == -1) name += "0-40% Centrality"; 
    else {
                      name += TString::Format("%d-%d", y[c], y[c+1]);
                      name += "% Centrality";
    }
    
    TCanvas* c1 = SetupCanvas(name);
    pad = (TPad*)c1->GetPrimitive("p1");
    pad->cd();

    TLatex* tit = new TLatex(0.10, 0.95, name);
    tit->SetNDC();
    tit->SetTextFont(132);
    tit->SetTextSize(0.07);
    tit->Draw();

    Int_t nEvs = (Int_t)fmdDiffHist->GetBinContent(fmdDiffHist->GetNbinsX()+1);
    TLatex* et = new TLatex(.93, .87, Form("%d events", nEvs));
    et->SetNDC();
    et->SetTextFont(132);
    et->SetTextAlign(33);
    et->Draw();

    mcDiffHist->SetMarkerStyle(21);
    fmdDiffHist->SetLineColor(kBlue);
    spdDiffHist->SetLineColor(kRed);
    fmdtrDiffHist->SetLineColor(kGreen+1);
    spdtrDiffHist->SetLineColor(kCyan+1);
    mcTruth->SetFillColor(kGray);

    mcTruth->Draw("same e3");
    mcDiffHist->Draw("same e1");
    fmdDiffHist->Draw("same e1");
    spdDiffHist->Draw("same e1");
    fmdtrDiffHist->Draw("same e1");
    spdtrDiffHist->Draw("same e1");
 
    TLegend* l = new TLegend(.12,.40,.35,.80);
    l->SetFillColor(0);
    l->SetFillStyle(0);
    l->SetBorderSize(0);
    l->SetTextFont(132);
    l->AddEntry(fmdDiffHist, "FMD", "lpe");
    l->AddEntry(spdDiffHist, "SPD", "lpe");
    l->AddEntry(fmdtrDiffHist, "FMD TrackRefs", "lpe");
    l->AddEntry(spdtrDiffHist, "SPD TrackRefs", "lpe");
    l->AddEntry(mcDiffHist, "MC Truth", "lpe");
    l->AddEntry(mcTruth, "cos(2(#varphi-#Psi_{rp}))", "f");
    l->Draw();

    pad = (TPad*)c1->GetPrimitive("p2");
    pad->cd();

    fmdDiff = fmdDiffHist->ProjectionX("fmdDiff", "e");
    mcDiff1 = mcDiffHist->ProjectionX("mcDiff1", "e");
    fmdDiff->Divide(mcDiff1);
    fmdDiff->SetLineColor(kBlue);
    fmdDiff->Draw("same e1");

    spdDiff = spdDiffHist->ProjectionX("spdDiff", "e");
    mcDiff2 = mcDiffHist->ProjectionX("mcDiff2", "e");
    spdDiff->Divide(mcDiff2);
    spdDiff->SetLineColor(kRed);
    spdDiff->Draw("same e1");

    pad = (TPad*)c1->GetPrimitive("p3");
    pad->cd();

    fmdTrDiff = fmdDiffHist->ProjectionX("fmdTrDiff", "e");
    mcDiff3 = fmdtrDiffHist->ProjectionX("mcDiff3", "e");
    fmdTrDiff->Divide(mcDiff3);
    fmdTrDiff->SetLineColor(kBlue);
    fmdTrDiff->Draw("same e1");

    spdTrDiff = spdDiffHist->ProjectionX("spdTrDiff", "e");
    mcDiff4 = spdtrDiffHist->ProjectionX("mcDiff4", "e");
    spdTrDiff->Divide(mcDiff4);
    spdTrDiff->SetLineColor(kRed);
    spdTrDiff->Draw("same e1");

    name.Prepend("Title:");
    cout<<name<<endl;
    c1->Print(pdfName, name);
  
    delete mcDiff1;
    delete mcDiff2;
    delete mcDiff3;
    delete mcDiff4; 
    delete fmdDiff;
    delete spdDiff;
    delete fmdTrDiff;
    delete spdTrDiff;
    delete c1;
   } // end of c
}

void Make4ParticlePlots(TFile* f) 
{
  TList* qList = static_cast<TList*>(f->Get("FlowResults/QCumulants"));
  TList* aList = static_cast<TList*>(qList->FindObject("MC"));
  TList* mcList = static_cast<TList*>(aList->FindObject("v2"));
  aList = static_cast<TList*>(qList->FindObject("FMDTR"));
  TList* fmdtrList = static_cast<TList*>(aList->FindObject("v2"));
  aList = static_cast<TList*>(qList->FindObject("SPDTR"));
  TList* spdtrList = static_cast<TList*>(aList->FindObject("v2"));
  aList = static_cast<TList*>(qList->FindObject("FMD"));
  TList* fmdList = static_cast<TList*>(aList->FindObject("v2"));
  aList = static_cast<TList*>(qList->FindObject("SPD"));
  TList* spdList = static_cast<TList*>(aList->FindObject("v2"));

  TProfile2D* mcProfile = (TProfile2D*)qList->FindObject("pMCTruth");
  TPad*       pad = 0;
 
  Int_t y[11] = { 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 100 };
  for (Int_t c = -2; c < 10; c++) {
    TProfile* mcDiffHist = 0;
    TProfile* fmdtrDiffHist = 0;
    TProfile* spdtrDiffHist = 0;
    TProfile* fmdDiffHist = 0;
    TProfile* spdDiffHist = 0;
    TH1D*       mcTruth = 0;
    TH1D*       fmdDiff = 0;
    TH1D*       spdDiff = 0;
    TH1D*       fmdTrDiff = 0;
    TH1D*       spdTrDiff = 0;
    TH1D*       mcDiff1 = 0;
    TH1D*       mcDiff2 = 0;
    TH1D*       mcDiff3 = 0;
    TH1D*       mcDiff4 = 0;

    if (c == -2) {
      mcDiffHist = (TProfile*)mcList->FindObject("hQ2Cumulant4DiffFlowMC_mb");
      fmdtrDiffHist = (TProfile*)fmdtrList->FindObject("hQ2Cumulant4DiffFlowFMDTR_mb");
      spdtrDiffHist = (TProfile*)spdtrList->FindObject("hQ2Cumulant4DiffFlowSPDTR_mb");
      fmdDiffHist = (TProfile*)fmdList->FindObject("hQ2Cumulant4DiffFlowFMD_mb");
      spdDiffHist = (TProfile*)spdList->FindObject("hQ2Cumulant4DiffFlowSPD_mb");
      mcTruth = (TH1D*)mcProfile->ProjectionX("mcTruth", 1, 9, "E");
      mcTruth->Scale(1./9.);
    }
    else if (c == -1) {
      mcDiffHist = (TProfile*)mcList->FindObject("hQ2Cumulant4DiffFlowMC_0_40");
      fmdtrDiffHist = (TProfile*)fmdtrList->FindObject("hQ2Cumulant4DiffFlowFMDTR_0_40");
      spdtrDiffHist = (TProfile*)spdtrList->FindObject("hQ2Cumulant4DiffFlowSPDTR_0_40");
      fmdDiffHist = (TProfile*)fmdList->FindObject("hQ2Cumulant4DiffFlowFMD_0_40");
      spdDiffHist = (TProfile*)spdList->FindObject("hQ2Cumulant4DiffFlowSPD_0_40");
      mcTruth = (TH1D*)mcProfile->ProjectionX("mcTruth", 1, 5, "E");
      mcTruth->Scale(1./5.);
    }
    else {
      mcDiffHist = (TProfile*)mcList->FindObject(Form("hQ2Cumulant4DiffFlowMC_%d_%d", y[c], y[c+1]));
      fmdtrDiffHist = (TProfile*)fmdtrList->FindObject(Form("hQ2Cumulant4DiffFlowFMDTR_%d_%d", y[c], y[c+1]));
      spdtrDiffHist = (TProfile*)spdtrList->FindObject(Form("hQ2Cumulant4DiffFlowSPDTR_%d_%d", y[c], y[c+1]));
      fmdDiffHist = (TProfile*)fmdList->FindObject(Form("hQ2Cumulant4DiffFlowFMD_%d_%d", y[c], y[c+1]));
      spdDiffHist = (TProfile*)spdList->FindObject(Form("hQ2Cumulant4DiffFlowSPD_%d_%d", y[c], y[c+1]));
      mcTruth = (TH1D*)mcProfile->ProjectionX("mcTruth", c+1, c+1, "E");
    }
 
    TString name = "v_{2}{4} MC ";
    if (c == -2)      name += "Minimum Bias";
    else if (c == -1) name += "0-40% Centrality"; 
    else {
                      name += TString::Format("%d-%d", y[c], y[c+1]);
                      name += "% Centrality";
    }    
     
    TCanvas* c1 = SetupCanvas(name);
    pad = (TPad*)c1->GetPrimitive("p1");
    pad->cd();

    TLatex* tit = new TLatex(0.10, 0.95, name);
    tit->SetNDC();
    tit->SetTextFont(132);
    tit->SetTextSize(0.07);
    tit->Draw();

    Int_t nEvs = (Int_t)fmdDiffHist->GetBinContent(fmdDiffHist->GetNbinsX()+1);
    TLatex* et = new TLatex(.93, .87, Form("%d events", nEvs));
    et->SetNDC();
    et->SetTextFont(132);
    et->SetTextAlign(33);
    et->Draw();

    mcDiffHist->SetMarkerStyle(21);
    fmdDiffHist->SetLineColor(kBlue);
    spdDiffHist->SetLineColor(kRed);
    fmdtrDiffHist->SetLineColor(kGreen+1);
    spdtrDiffHist->SetLineColor(kCyan+1);
    mcTruth->SetFillColor(kGray);

    mcTruth->Draw("same e3");
    mcDiffHist->Draw("same e1");
    fmdDiffHist->Draw("same e1");
    spdDiffHist->Draw("same e1");
    fmdtrDiffHist->Draw("same e1");
    spdtrDiffHist->Draw("same e1");
 
    TLegend* l = new TLegend(.12,.40,.35,.80);
//    l->SetNColumns(2);
    l->SetFillColor(0);
    l->SetFillStyle(0);
    l->SetBorderSize(0);
    l->SetTextFont(132);
    l->AddEntry(fmdDiffHist, "FMD", "lpe");
    l->AddEntry(spdDiffHist, "SPD", "lpe");
    l->AddEntry(fmdtrDiffHist, "FMD TrackRefs", "lpe");
    l->AddEntry(spdtrDiffHist, "SPD TrackRefs", "lpe");
    l->AddEntry(mcDiffHist, "MC Truth", "lpe");
    l->AddEntry(mcTruth, "cos(2(#varphi-#Psi_{rp}))", "f");
    l->Draw();

    pad = (TPad*)c1->GetPrimitive("p2");
    pad->cd();

    fmdDiff = fmdDiffHist->ProjectionX("fmdDiff", "e");
    mcDiff1 = mcDiffHist->ProjectionX("mcDiff1", "e");
    fmdDiff->Divide(mcDiff1);
    fmdDiff->SetLineColor(kBlue);
    fmdDiff->Draw("same e1");

    spdDiff = spdDiffHist->ProjectionX("spdDiff", "e");
    mcDiff2 = mcDiffHist->ProjectionX("mcDiff2", "e");
    spdDiff->Divide(mcDiff2);
    spdDiff->SetLineColor(kRed);
    spdDiff->Draw("same e1");

    pad = (TPad*)c1->GetPrimitive("p3");
    pad->cd();

    fmdTrDiff = fmdDiffHist->ProjectionX("fmdTrDiff", "e");
    mcDiff3 = fmdtrDiffHist->ProjectionX("mcDiff3", "e");
    fmdTrDiff->Divide(mcDiff3);
    fmdTrDiff->SetLineColor(kBlue);
    fmdTrDiff->Draw("same e1");

    spdTrDiff = spdDiffHist->ProjectionX("spdTrDiff", "e");
    mcDiff4 = spdtrDiffHist->ProjectionX("mcDiff4", "e");
    spdTrDiff->Divide(mcDiff4);
    spdTrDiff->SetLineColor(kRed);
    spdTrDiff->Draw("same e1");

    name.Prepend("Title:");
    cout<<name<<endl;
    c1->Print(pdfName, name);
  
    delete mcDiff1;
    delete mcDiff2;
    delete mcDiff3;
    delete mcDiff4;
    delete fmdDiff;
    delete spdDiff;
    delete fmdTrDiff;
    delete spdTrDiff;
    delete c1;
  } // end of c
}

void MakeMonitoringPlots(TFile* f) 
{
  TList* qList = static_cast<TList*>(f->Get("FlowResults/QCumulants"));

  TH1D* cent = (TH1D*)qList->FindObject("Centralities");
  TH2D* vert = (TH2D*)qList->FindObject("CoverageVsVertex");

  TCanvas* c1 = SetupCanvas("Monitoring");
  c1->Divide(1,2);
  c1->cd(1);
  cent->GetXaxis()->SetTitleFont(22);
  cent->GetXaxis()->SetLabelFont(22);
  cent->GetYaxis()->SetTitleFont(22);
  cent->GetYaxis()->SetLabelFont(22);
  cent->SetXTitle("Centrality %");
  cent->SetYTitle("# of events");
  cent->Draw();
  c1->cd(2);
  vert->GetXaxis()->SetTitleFont(22);
  vert->GetXaxis()->SetLabelFont(22);
  vert->GetYaxis()->SetTitleFont(22);
  vert->GetYaxis()->SetLabelFont(22);
  vert->SetXTitle("#eta");
  vert->SetYTitle("Vertex z coordinate [cm]");
  vert->Draw("colz");

  c1->Print(pdfName, "Title:Monitoring Plots");
}


void DrawFlowPDF(char* file = "AnalysisResults.root") 
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TFile* f = TFile::Open(file, "READ");
  TCanvas* c0 = SetupCanvas("empty");
  c0->SetName("empty");
  c0->Print(Form("%s[", pdfName));
  MakeFmdAndSpdPlots(f);
  Make2ParticlePlots(f);
  Make4ParticlePlots(f);
  MakeMonitoringPlots(f);
  c0->Print(Form("%s]", pdfName));
}
