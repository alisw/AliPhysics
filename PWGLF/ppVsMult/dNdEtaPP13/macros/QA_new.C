#if !defined (__CINT__) || (defined(__MAKECINT__))
#include <iostream>
#include "TFile.h"
#include "TList.h"
#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "TProfile.h"
#include "TStyle.h"

#endif

#include "style.C"
#include "TObjString.h"
#include "TObjArray.h"
// TODO: add run/MCperiod to the output files
void QAtracklets(const Char_t *fdata, const Char_t *fmc);
void QAtrackletsRatio(const Char_t *fdata, const Char_t *fmc);
void QAvertex(const Char_t *fdata, const Char_t *fmc);
void QAcentrality(const Char_t *fdata);
void QAoccupancy(const Char_t *fdata, const Char_t *fmc);
void QAClusters0(const Char_t *fdata, const Char_t *fmc);
void QA_V0(const Char_t *fdata, const Char_t *fmc);
void QAChi2(const Char_t *fdata, const Char_t *fmc);
void QA_tracklets(const Char_t *fdata, const Char_t *fmc);


TString canvasPrefix = "QA_"; //225705  ../QA_RBR/MB/QA_
void QA_new(TString fdata = "data.root",
TString fmc   = "mc.root")

{

  //  qacentrality(fdata);
  QAcentrality(fdata);
  QAoccupancy(fdata, fmc);
  QAtracklets(fdata, fmc);
  QAvertex(fdata, fmc);
  QAtrackletsRatio(fdata,fmc);
  QAClusters0(fdata,fmc);
  QA_V0(fdata, fmc);
  QAChi2(fdata,fmc);
  QA_tracklets(fdata,fmc);


}

void QAtracklets(const Char_t *fdata, const Char_t *fmc)
{

  style();

  TFile *fdtin = TFile::Open(fdata);
  TList *ldtin = (TList *)fdtin->Get("clist");
  TH2 *hdtin = (TH2 *)ldtin->FindObject("etaphiTracklets");
  TH1 *pdtin = (TH1 *)hdtin->ProjectionY("pdtin_tracklets");
  pdtin->SetMarkerStyle(20);
  pdtin->SetMarkerSize(1.);
  pdtin->SetMarkerColor(kAzure-3);
  hdtin->Scale(1. / hdtin->GetEntries());
  pdtin->Scale(1. / hdtin->GetEntries());

  TFile *fmcin = TFile::Open(fmc);
  TList *lmcin = (TList *)fmcin->Get("clist");

  if(!lmcin) {
      std::cout << "NOLIST" << std::endl;

  }
  lmcin->ls();
  TH2 *hmcin = (TH2 *)lmcin->FindObject("etaphiTracklets");
  if(!hmcin) {
    std::cout << "NO H!! etaphiTracklets" << std::endl;

  }
  TH1 *pmcin = (TH1 *)hmcin->ProjectionY("pmcin_tracklets");
  pmcin->SetLineColor(kRed+1);
  pmcin->SetFillStyle(1001);
  pmcin->SetFillColorAlpha(kRed+1, 0.1);
  hmcin->Scale(1. / hmcin->GetEntries());
  pmcin->Scale(1. / hmcin->GetEntries());

  TCanvas *cData = new TCanvas("cTrackletData", "cTrackletData", 800, 800);
  //  cData->SetLogz();
  TH1 * hfr = cData->DrawFrame(-2.5, 0., 2.5, 2. * TMath::Pi());
  hfr->SetTitle(";#eta;#varphi");
  hdtin->Draw("same,col");
  cData->SaveAs(canvasPrefix+"trackletData.pdf");

  TCanvas *cMC = new TCanvas("cTrackletMC", "cTrackletMC", 800, 800);
  //  cMC->SetLogz();
  hfr = cMC->DrawFrame(-2.5, 0., 2.5, 2. * TMath::Pi());
  hfr->SetTitle(";#eta;#varphi");
  hmcin->Draw("same,col");
  cMC->SaveAs(canvasPrefix+"trackletMC.pdf");

  TCanvas *cPhi = new TCanvas("cTrackletPhi", "cTrackletPhi", 800, 800);
  //  cPhi->SetLogy();
  hfr = cPhi->DrawFrame(0., 0., 2. * TMath::Pi(), 0.013);
  hfr->SetTitle(";#varphi;");
  pdtin->DrawCopy("same");
  pmcin->DrawCopy("same,histo");
  TLegend *legend = new TLegend(0.20, 0.21+0.63, 0.40, 0.31+0.63);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.04);
  legend->AddEntry(pdtin, "data", "pl");
  legend->AddEntry(pmcin, "MC (LHC15g3c3)", "l");
  legend->Draw("same");
  cPhi->SaveAs(canvasPrefix+"trackletPhi.pdf");

  TCanvas *cPhir = new TCanvas("cTrackletPhir", "cTrackletPhir", 800, 800);
  //  cPhi->SetLogy();
  hfr = cPhir->DrawFrame(0., 0.5, 2. * TMath::Pi(), 1.5);
  hfr->SetTitle(";#varphi;data / Monte Carlo");
  pdtin->Divide(pmcin);
  pdtin->Draw("same");
  cPhir->SaveAs(canvasPrefix+"trackletPhir.pdf");
}

void QAClusters0(const Char_t *fdata, const Char_t *fmc)
{

  style();

  TFile *fdtin = TFile::Open(fdata);
  TList *ldtin = (TList *)fdtin->Get("clist");
  TH2 *hdtin = (TH2 *)ldtin->FindObject("cl0InfoUsed");
  TH1 *pdtin = (TH1 *)hdtin->ProjectionY("pdtin_cls0");
  pdtin->SetMarkerStyle(20);
  pdtin->SetMarkerSize(1.);
  pdtin->SetMarkerColor(kAzure-3);
  hdtin->Scale(1. / hdtin->GetEntries());
  pdtin->Scale(1. / hdtin->GetEntries());

  TFile *fmcin = TFile::Open(fmc);
  TList *lmcin = (TList *)fmcin->Get("clist");

  if(!lmcin) {
      std::cout << "NOLIST" << std::endl;

  }
  lmcin->ls();
  TH2 *hmcin = (TH2 *)lmcin->FindObject("cl0InfoUsed");
  if(!hmcin) {
    std::cout << "NO H!! clusters0" << std::endl;

  }
  TH1 *pmcin = (TH1 *)hmcin->ProjectionY("pmcin_cls0");
  pmcin->SetLineColor(kRed+1);
  pmcin->SetFillStyle(1001);
  pmcin->SetFillColorAlpha(kRed+1, 0.1);
  hmcin->Scale(1. / hmcin->GetEntries());
  pmcin->Scale(1. / hmcin->GetEntries());


  TCanvas *cData = new TCanvas("cClusters0Data", "cClusters0Data", 800, 800);
  //  cData->SetLogz();
  TH1 * hfr = cData->DrawFrame(-20, 0., 20, 2. * TMath::Pi());
  hfr->SetTitle(";Z_{vertex}; #varphi");
  hdtin->Draw("same,col");
  cData->SaveAs(canvasPrefix+"Clusters0Data_Used.pdf");

  TCanvas *cMC = new TCanvas("cClusters0MC", "cClusters0MC", 800, 800);
  //  cMC->SetLogz();
  hfr = cMC->DrawFrame(-20, 0., 20, 2. * TMath::Pi());
  hfr->SetTitle(";Z_{vertex}; #varphi");
  hmcin->Draw("same,col");
  cMC->SaveAs(canvasPrefix+"Clusters0MC_Used.pdf");

  TCanvas *cPhi = new TCanvas("cClusters0Phi", "cClusters0Phi", 800, 800);
  //  cPhi->SetLogy();
  hfr = cPhi->DrawFrame(0., 0., 2. * TMath::Pi(), 0.028);
  hfr->SetTitle(";#varphi;");
  pdtin->DrawCopy("same");
  pmcin->DrawCopy("same,histo");
  TLegend *legend = new TLegend(0.20, 0.21+0.63, 0.40, 0.31+0.63);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.04);
  legend->AddEntry(pdtin, "data", "pl");
  legend->AddEntry(pmcin, "MC (LHC15g3c3)", "l");
  legend->Draw("same");
  cPhi->SaveAs(canvasPrefix+"Clusters0Phi_Used.pdf");

  TCanvas *cPhir = new TCanvas("cTrackletPhir", "cClusters0Phir", 800, 800);
  //  cPhi->SetLogy();
  hfr = cPhir->DrawFrame(0., 0.5, 2. * TMath::Pi(), 1.5);
  hfr->SetTitle(";#varphi;data / Monte Carlo");
  pdtin->Divide(pmcin);
  pdtin->Draw("same");
  cPhir->SaveAs(canvasPrefix+"Clusters0Phir_Used.pdf");
}


void QAtrackletsRatio(const Char_t *fdata, const Char_t *fmc)
{

  style();

  TFile *fdtin = TFile::Open(fdata);
  TList *ldtin = (TList *)fdtin->Get("clist");
  TH2 *hdtin = (TH2 *)ldtin->FindObject("etaphiTracklets");

  hdtin->Scale(1. / hdtin->GetEntries());

  TFile *fmcin = TFile::Open(fmc);
  TList *lmcin = (TList *)fmcin->Get("clist");

  if(!lmcin) {
      std::cout << "NOLIST" << std::endl;

  }
  lmcin->ls();
  TH2 *hmcin = (TH2 *)lmcin->FindObject("etaphiTracklets");
  if(!hmcin) {
    std::cout << "NO H!! etaphiTracklets" << std::endl;

  }

  hmcin->Scale(1. / hmcin->GetEntries());

  TCanvas *cPhir = new TCanvas("cTrackletPhir", "cTrackletPhir", 800, 800);

  hdtin->Divide(hmcin);

  hdtin->SetMinimum(0.8);
  hdtin->SetMaximum(1.2);
  gPad->SetRightMargin(0.13);
  hdtin->GetXaxis()->SetRangeUser(-1.9,1.9);
  hdtin->Draw("colz");
  cPhir->SaveAs(canvasPrefix+"trackletPhiEtaratio.pdf");
}


void QAoccupancy(const Char_t *fdata, const Char_t *fmc)
{
  style();

  TFile *fdtin = TFile::Open(fdata);
  TList *ldtin = (TList *)fdtin->Get("clist");
  TH2 *hdtin = (TH2 *)ldtin->FindObject("NClustersSPD2");
  TProfile *pdtin = hdtin->ProfileY("pdtin_clusters");
  pdtin->SetMarkerStyle(20);
  pdtin->SetMarkerSize(2);
  pdtin->SetMarkerColor(kAzure-3);

  TFile *fmcin = TFile::Open(fmc);
  TList *lmcin = (TList *)fmcin->Get("clist");
  TH2 *hmcin = (TH2 *)lmcin->FindObject("NClustersSPD2");
  TProfile *pmcin = hmcin->ProfileY("pmcin_clusters");
  pmcin->SetMarkerStyle(25);
  pmcin->SetMarkerSize(2);
  pmcin->SetMarkerColor(kRed+1);

  TCanvas *c = new TCanvas("cOccupancy", "cOccupancy", 800, 800);
  c->SetLogy();
  TH1 * hfr = c->DrawFrame(-1.5, 2., 10.5, 500.);
  DrawBinLabelsX(hfr, kTRUE);
  hfr->SetTitle(";;#LT#it{N}_{clusters,SPD-1}#GT");
  pdtin->DrawCopy("same");
  pmcin->DrawCopy("same");
  TLegend *legend = new TLegend(0.20, 0.18, 0.50, 0.30);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.04);
  legend->AddEntry(pdtin, "data", "pl");
  legend->AddEntry(pmcin, "MC (LHC15g3c3)", "pl");
  legend->Draw("same");
  c->SaveAs(canvasPrefix+"occupancy.pdf");
  return;
  TCanvas *cr = new TCanvas("cOccupancyr", "cOccupancyr", 800, 800);
  // hfr = cr->DrawFrame(-0.5, 0.75, 10.5, 1.25);
  // DrawBinLabelsX(hfr, kTRUE);
  // hfr->SetTitle(";;#LT#it{N}_{clusters,SPD-1}#GT ratio");
  pdtin->SetLineColor(kAzure-3);
  pdtin->SetLineWidth(3);
  pdtin->Divide(pmcin);
  pdtin->Draw("same,histo");
  legend = new TLegend(0.505025, 0.760673, 0.805276, 0.930142);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.04);
  legend->AddEntry(pdtin, "data / Monte Carlo", "l");
  legend->Draw("same");
  cr->SaveAs(canvasPrefix+"occupancyr.pdf");

}

void QAcentrality(const Char_t *fdata)
{
  style();

  TFile *fin = TFile::Open(fdata);
  TList *lin = (TList *)fin->Get("clist");
  lin->ls();
  TH1 *hin = (TH1 *)lin->FindObject("EvCentrDist");
  Float_t sum = 1.2 * hin->Integral(hin->FindBin(0.1), hin->FindBin(79.9));
  //Float_t sum =  hin->Integral(hin->FindBin(0.1), hin->FindBin(99.9));

  hin->Scale(1. / sum);
  SetHistoStyle(hin, 20, kRed+1);
  TCanvas *c = new TCanvas("cQAcentrality", "cQAcentrality", 800, 800);
  TH1 * hfr = c->DrawFrame(0., 0.005, 100., 0.015);
  hfr->SetTitle(";centrality percentile;events");
  hin->Draw("same");
  c->SaveAs(canvasPrefix+"centrality.pdf");

  TH2 *hinv0 = (TH2 *)lin->FindObject("V0");
  TCanvas *cv0 = new TCanvas("cQAcentralityV0", "cQAcentralityV0", 800, 800);
  cv0->SetLogx();
  cv0->SetLogz();
  //  TH1 * hfrv0 = cv0->DrawFrame(100., -0.5, 50000., 10.5);
  // DrawBinLabelsY(hfrv0, kTRUE);
  // hfrv0->SetTitle(";V0 signal;");
  //hinv0->Draw("same,col");
  hinv0->Draw("col");
  cv0->SaveAs(canvasPrefix+"centralityV0.pdf");
}

void QAvertex(const Char_t *fdata, const Char_t *fmc)
{
  style();

  TFile *fdtin = TFile::Open(fdata);
  TList *ldtin = (TList *)fdtin->Get("clist");
  TH2 *hdtin = (TH2 *)ldtin->FindObject("zv");
  TH1 *hdt = (TH1 *)ldtin->FindObject("zvNoSel");
  SetHistoStyle(hdt, 20, kRed+1);
  hdt->Scale(1. / hdt->Integral());

  TH1 *hdt0010 = hdtin->ProjectionX("hdt0010", 1, 4);
  SetHistoStyle(hdt0010, 20, kRed+1);
  hdt0010->Scale(1. / hdt0010->Integral());

  TH1 *hdt70100 = hdtin->ProjectionX("hdt70100", 10, 11);
  SetHistoStyle(hdt70100, 25, kAzure-3);
  hdt70100->Scale(1. / hdt70100->Integral());

  TFile *fmcin = TFile::Open(fmc);
  TList *lmcin = (TList *)fmcin->Get("clist");
  TH1 *hmc = (TH1 *)lmcin->FindObject("zvNoSel");
  SetHistoStyle(hmc, 25, kAzure-3);
  hmc->Scale(1. / hmc->Integral());

  TCanvas *c = new TCanvas("cVertex", "cVertex", 800, 800);
  TH1 * hfr = c->DrawFrame(-20., 0., 20., 0.115);
  hfr->SetTitle(";#it{z}_{vtx};");
  hdt0010->Draw("same");
  hdt70100->Draw("same");
  TLegend *legend = new TLegend(0.20, 0.18+0.60, 0.40, 0.30+0.60);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.04);
  legend->AddEntry(hdt0010, "0-10%", "p");
  legend->AddEntry(hdt70100, "70-100%", "p");
  legend->Draw("same");
  c->SaveAs(canvasPrefix+"Data_vertex.pdf");

  TCanvas *c1 = new TCanvas("cVertexDataMC", "cVertexDataMC", 800, 800);
  hfr = c1->DrawFrame(-20., 0., 20., 0.115);
  hfr->SetTitle(";#it{z}_{vtx};");
  hdt->Draw("same");
  hmc->Draw("same");
  legend = new TLegend(0.20, 0.18+0.60, 0.40, 0.30+0.64);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.04);
  legend->AddEntry(hdt, "data", "p");
  legend->AddEntry(hmc, "MC (LHC15g3c3)", "p");
  legend->Draw("same");
  c1->SaveAs(canvasPrefix+"vertexDataMC.pdf");

  //return 0;
}

void QA_V0(const Char_t *fdata, const Char_t *fmc)
{

  style();

  TFile *fdtin = TFile::Open(fdata);
  TList *ldtin = (TList *)fdtin->Get("clist");
  TH2 *hdtin = (TH2 *)ldtin->FindObject("V0");
  TH1 *hldtin1 = hdtin->ProjectionX("V0withPS1", 1, 2);
  hldtin1->Scale(1. / hldtin1->GetEntries(),"width");
  hldtin1->SetMarkerColor(kBlue);
  hldtin1->SetMarkerStyle(kFullCircle);
  hldtin1->SetMarkerSize(0.5);
  hldtin1->GetYaxis()->SetTitle("Events");
 // hldtin1->GetYaxis()->SetLabelSize(0.1);

  TFile *fmcin = TFile::Open(fmc);
  TList *lmcin = (TList *)fmcin->Get("clist");
  TH2 *hmcin = (TH2 *)lmcin->FindObject("V0");
  TH1 *hldtin1mc = hmcin->ProjectionX("V0withPS1mc", 1, 2);
  hldtin1mc->Scale(1. / hldtin1mc->GetEntries(),"width");
  hldtin1mc->SetMarkerColor(kRed);
  hldtin1mc->SetMarkerStyle(kFullCircle);
  hldtin1mc->SetMarkerSize(0.5);

  TCanvas *c1_v0 = new TCanvas("c1_v0", "c1_v0", 1200, 800);
  c1_v0->Divide(1,1);
  c1_v0->cd(1);
  // Upper plot will be in pad1a
  TPad *pad1a = new TPad("pad1a", "pad1a", 0, 0.3, 1, 1.0);
  pad1a->SetBottomMargin(0); // Upper and lower plot are joined
  pad1a->SetGridx();
  pad1a->SetGridy();
  pad1a->SetLogy();         // Vertical grid
  // Vertical grid
  pad1a->Draw();             // Draw the upper pad: pad1a
  pad1a->cd();               // pad1a becomes the current pad
  //    h1->SetStats(0);          // No statistics on upper plot
  hldtin1->Draw();               // Draw h1
  hldtin1mc->Draw("same");         // Draw h2 on top of h1

  TLegend * l1 = new TLegend(0.6,0.8,0.85,0.9);
  l1->AddEntry(hldtin1, " (Data)", "p");
  l1->AddEntry(hldtin1mc, "(MC)", "p");
  //l->AddEntry(hldtin3, "GRL_3 (INT7) (Data)", "p");
  l1->Draw("same");

  c1_v0->cd(1);          // Go back to the main canvas before defining pad2a
  TPad *pad2a = new TPad("pad2a", "pad2a", 0, 0.05, 1, 0.3);
  pad2a->SetTopMargin(0);
  pad2a->SetBottomMargin(0.2);
  pad2a->SetGridx(); // vertical grid
  pad2a->SetLogy();         // Vertical grid

  pad2a->Draw();
  pad2a->cd();       // pad2a becomes the current pad

  // Define the ratio plot
  TH1F *h1 = (TH1F*)hldtin1->Clone("h1");
  h1->SetMarkerColor(kBlack);
  h1->GetYaxis()->SetTitle("Data/MC ");
  h1->SetStats(0);      // No statistics on lower plot
  h1->Divide(hldtin1mc);
  h1->SetMarkerStyle(21);
  h1->Draw("ep");       // Draw the ratio plot

  c1_v0->SaveAs(canvasPrefix+"Sel_V0Signal.pdf");
}

void QAChi2(const Char_t *fdata, const Char_t *fmc)
{

  style();

  TFile *fdtin = TFile::Open(fdata);
  TList *ldtin = (TList *)fdtin->Get("clist");
  TH2 *hfdtin = (TH2 *)ldtin->FindObject("b0_TrData_WDvsEta");
  TH1 *hdtin = (TH1 *)hfdtin->ProjectionY("b0_WDist");
  hfdtin->Scale(1. / hfdtin->GetEntries());
  hdtin->Scale(1. / hfdtin->GetEntries());




  hdtin->SetMarkerColor(kBlue);
  hdtin->SetMarkerStyle(kFullCircle);
  hdtin->SetMarkerSize(0.5);
  hdtin->GetYaxis()->SetTitle("Normalized Events");
  hdtin->GetXaxis()->SetTitle("#chi^{2}");
  hdtin->GetXaxis()->SetRangeUser(0.,2.0);



  TFile *fmcin = TFile::Open(fmc);
  TList *hlmcin = (TList *)fmcin->Get("clist");
    TH2 *hfmcin = (TH2 *)hlmcin->FindObject("b0_TrData_WDvsEta");
    TH1 *hmcin = (TH1 *)hfmcin->ProjectionY("b0_WDist");

    hfmcin->Scale(1. / hfmcin->GetEntries());
   hmcin->Scale(1. / hfmcin->GetEntries());


  if(!hmcin) {
      std::cout << "NOLIST" << std::endl;

  }
  hmcin->ls();

    hmcin->SetMarkerColor(kRed);
    hmcin->SetMarkerStyle(kFullCircle);
    hmcin->SetMarkerSize(0.5);


  if(!hmcin) {
    std::cout << "NO H!! etaphiTracklets" << std::endl;

  }



  TCanvas *Chi2 = new TCanvas("Chi2", "Chi2", 800, 800);
  Chi2->SetLogy();
  //TH1 * hfr = Chi2->DrawFrame(0., 0., 2.0, 2);
  hdtin->Draw("same");
  hmcin->Draw("same");
  TLegend *legend = new TLegend(0.20, 0.18+0.60, 0.40, 0.30+0.60);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.04);
  legend->AddEntry(hdtin, "data", "p");
  legend->AddEntry(hmcin, "MC (LHC15g3c3)", "p");
  legend->Draw("same");
  Chi2->SaveAs(canvasPrefix+"Chi2.pdf");


  TCanvas *Chi2r = new TCanvas("Chi2r", "Chi2r", 800, 800);

  hdtin->Divide(hmcin);

 // hdtin->SetMinimum(0.8);
 // hdtin->SetMaximum(1.2);
  gPad->SetRightMargin(0.13);
  hdtin->GetXaxis()->SetRangeUser(0.,2.0);
  hdtin->GetYaxis()->SetRangeUser(0.,2.0);
  hdtin->GetYaxis()->SetTitle("DA/MC");
  hdtin->GetXaxis()->SetTitle("#chi^{2}");
  hdtin->Draw();
  Chi2r->SaveAs(canvasPrefix+"Chi2ratio.pdf");
}

void QA_tracklets(const Char_t *fdata, const Char_t *fmc)
{

  style();

  TFile *fdtin = TFile::Open(fdata);
  TList *ldtin = (TList *)fdtin->Get("clist");
  TH1 *hldtin1 = (TH1 *)ldtin->FindObject("mltTracklets");
  //TH1 *hldtin1 = hdtin->ProjectionX("V0withPS1", 1, 2);
  hldtin1->Scale(1. / hldtin1->GetEntries(),"width");
  hldtin1->SetMarkerColor(kBlue);
  hldtin1->SetMarkerStyle(kFullCircle);
  hldtin1->SetMarkerSize(0.5);
  hldtin1->GetYaxis()->SetTitle("P(N_{ch}) (SPD Tracklets)");
  hldtin1->GetXaxis()->SetRangeUser(0.,150.0);

 // hldtin1->GetYaxis()->SetLabelSize(0.1);

  TFile *fmcin = TFile::Open(fmc);
  TList *lmcin = (TList *)fmcin->Get("clist");
 // TH2 *hmcin = (TH2 *)lmcin->FindObject("V0");
  TH1 *hmcin1mc = (TH1 *)lmcin->FindObject("mltTracklets");
 // TH1 *hldtin1mc = hmcin->ProjectionX("V0withPS1mc", 1, 2);
  hmcin1mc->Scale(1. / hmcin1mc->GetEntries(),"width");
  hmcin1mc->SetMarkerColor(kRed);
  hmcin1mc->SetMarkerStyle(kFullCircle);
  hmcin1mc->SetMarkerSize(0.5);

  TCanvas *c1_tracklets = new TCanvas("c1_tracklets", "c1_tracklets", 1200, 800);
  c1_tracklets->Divide(1,1);
  c1_tracklets->cd(1);
  // Upper plot will be in pad1a
  TPad *pad1a = new TPad("pad1a", "pad1a", 0, 0.3, 1, 1.0);
  pad1a->SetBottomMargin(0); // Upper and lower plot are joined
  pad1a->SetGridx();
  pad1a->SetGridy();
  pad1a->SetLogy();         // Vertical grid
  // Vertical grid
  pad1a->Draw();             // Draw the upper pad: pad1a
  pad1a->cd();               // pad1a becomes the current pad
  //    h1->SetStats(0);          // No statistics on upper plot
  hldtin1->Draw();               // Draw h1
  hmcin1mc->Draw("same");         // Draw h2 on top of h1

  TLegend * l1 = new TLegend(0.6,0.8,0.85,0.9);
  l1->AddEntry(hldtin1, " (Data)", "p");
  l1->AddEntry(hmcin1mc, "(MC)", "p");
  //l->AddEntry(hldtin3, "GRL_3 (INT7) (Data)", "p");
  l1->Draw("same");

  c1_tracklets->cd(1);          // Go back to the main canvas before defining pad2a
  TPad *pad2a = new TPad("pad2a", "pad2a", 0, 0.05, 1, 0.3);
  pad2a->SetTopMargin(0);
  pad2a->SetBottomMargin(0.2);
  pad2a->SetGridx(); // vertical grid
  pad2a->SetLogy();         // Vertical grid

  pad2a->Draw();
  pad2a->cd();       // pad2a becomes the current pad

  // Define the ratio plot
  TH1F *h1 = (TH1F*)hldtin1->Clone("h1");
  h1->SetMarkerColor(kBlack);
  h1->GetYaxis()->SetTitle("Data/MC ");
  h1->GetXaxis()->SetTitle(" SPD Tracklets Multiplicity ");

  h1->SetStats(0);      // No statistics on lower plot
  h1->Divide(hmcin1mc);
  h1->SetMarkerStyle(21);
  h1->Draw("ep");       // Draw the ratio plot

  c1_tracklets->SaveAs(canvasPrefix+"Mult_Tracklets.pdf");
}
