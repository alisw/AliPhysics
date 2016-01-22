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
void QAvertex(const Char_t *fdata, const Char_t *fmc);
void QAcentrality(const Char_t *fdata);
void QAoccupancy(const Char_t *fdata, const Char_t *fmc);

TString canvasPrefix = "QA_";

void QA(TString fdata = "trmult.data.root",
   TString fmc   = "trmult.mc.root")
{


  // The following code assumes that we have a filenme with the format
  // "../output/LHC15f_pass2/226062/AnalysisResults.226062.root
  // and uses this assumption to extract MC and Data taks + runname
  TObjArray * tokensData = fdata.Tokenize("/");
  TObjArray * tokensMC   = fmc.Tokenize("/");
  
  // Check assuptions on data and MC formats, the second token should be "output"
  if (((TObjString*)tokensData->At(1))->String() == "output" && 
      ((TObjString*)tokensMC  ->At(1))->String() == "output"
      ) {
    
    canvasPrefix =
      canvasPrefix +
      ((TObjString*)tokensData->At(2))->String() + "_" +// period data
      ((TObjString*)tokensMC  ->At(2))->String() + "_" +// period mc
      ((TObjString*)tokensData->At(3))->String() + "_" ;// run number
      
      
  }
  std::cout << "Prefix: " << canvasPrefix.Data() << std::endl;
  
  //  qacentrality(fdata);
  QAoccupancy(fdata, fmc);
  QAtracklets(fdata, fmc);
  QAvertex(fdata, fmc);
  
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
  

  /*  
  pdtin->Scale(pmcin->Integral(pmcin->FindBin(0. + 0.001),
			       pmcin->FindBin(1. - 0.001))
	       / pdtin->Integral(pdtin->FindBin(0. + 0.001),
				 pdtin->FindBin(1. - 0.001)));
  */  
  
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
  hfr = cPhi->DrawFrame(0., 0., 2. * TMath::Pi(), 0.01);
  hfr->SetTitle(";#varphi;");
  pdtin->DrawCopy("same");
  pmcin->DrawCopy("same,histo");
  TLegend *legend = new TLegend(0.20, 0.18+0.63, 0.50, 0.30+0.63);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.04);
  legend->AddEntry(pdtin, "data", "pl");
  legend->AddEntry(pmcin, "Monte Carlo", "l");
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
  TH1 * hfr = c->DrawFrame(-0.5, 2., 10.5, 500.);
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
  legend->AddEntry(pmcin, "Monte Carlo", "pl");
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

  TH1 *hdt7080 = hdtin->ProjectionX("hdt7080", 11, 11);
  SetHistoStyle(hdt7080, 25, kAzure-3);
  hdt7080->Scale(1. / hdt7080->Integral());
  
  TFile *fmcin = TFile::Open(fmc);
  TList *lmcin = (TList *)fmcin->Get("clist");
  TH1 *hmc = (TH1 *)lmcin->FindObject("zvNoSel");
  SetHistoStyle(hmc, 25, kAzure-3);
  hmc->Scale(1. / hmc->Integral());

  TCanvas *c = new TCanvas("cVertex", "cVertex", 800, 800);
  TH1 * hfr = c->DrawFrame(-20., 0., 20., 0.1);
  hfr->SetTitle(";#it{z}_{vtx};");
  hdt0010->Draw("same");
  hdt7080->Draw("same");
  TLegend *legend = new TLegend(0.20, 0.18+0.60, 0.50, 0.30+0.60);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.04);
  legend->AddEntry(hdt0010, "0-10%", "p");
  legend->AddEntry(hdt7080, "70-80%", "p");
  legend->Draw("same");
  c->SaveAs(canvasPrefix+"vertex.pdf");
  
  TCanvas *c1 = new TCanvas("cVertexDataMC", "cVertexDataMC", 800, 800);
  hfr = c1->DrawFrame(-20., 0., 20., 0.1);
  hfr->SetTitle(";#it{z}_{vtx};");
  hdt->Draw("same");
  hmc->Draw("same");
  legend = new TLegend(0.20, 0.18+0.60, 0.50, 0.30+0.60);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.04);
  legend->AddEntry(hdt, "data", "p");
  legend->AddEntry(hmc, "Monte Carlo", "p");
  legend->Draw("same");
  c1->SaveAs(canvasPrefix+"vertexDataMC.pdf");
  
  //return 0;  
}

