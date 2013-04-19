#include "TCanvas.h"
#include <TH1.h>
#include <TF1.h>
#include <TFile.h>
#include <Rtypes.h>
#include <TLegend.h>
#include <TString.h>
#include <TAttMarker.h>
#include <RtypesCint.h>
#include <TNamed.h>
#include "TDirectoryFile.h"
#include <TPRegexp.h>
#include <X3DDefs.h>

namespace RawProduction {
  class Output;
}

enum EffVersion { LHC10h_1234_Apr_10 };
EffVersion effVersion = LHC10h_1234_Apr_10;
bool canvHalfWidth=false;

TF1* GetEfficency(const TString& trigger, int fromCent, int toCent, const TString& pid, const TString& methode)
{
  // Dmitri LHC10h Apr. 10 Efficiencies
  // avalible pids: Allcore Disp2core CPVcore Both2core Disp2 All CPV Both2
  if( effVersion == LHC10h_1234_Apr_10 ) {
    if( ! methode.Contains("yr1") ) {
      Printf("ERROR:GetEfficency: only pol1 efficiancies avalible");
      return 0x0;
    }
    
    // Cent bin defintion as extracted from PHOS_embedding/
    int cent = -1;
    if( 0 == fromCent && 5 == toCent ) cent = 0;
    if( 5 == fromCent && 10 == toCent ) cent = 1;
    if( 10 == fromCent && 20 == toCent ) cent = 2;
    if( 20 == fromCent && 40 == toCent ) cent = 3;
    if( 40 == fromCent && 60 == toCent ) cent = 4;
    if( 60 == fromCent && 80 == toCent ) cent = 5;
    if( cent < 0 ) {
      Printf("ERROR:GetEfficency not avalible for centrality [%i,%i)", fromCent, toCent);
      return 0x0;
    }
    
    // determine name and return efficiancy function
    TDirectory* pastDir = gDirectory;
    TFile* file = TFile::Open("PHOS_eff_Full_PbPb_1234.root", "READ");
    char funcName[256]; 
    if ( methode.Contains("int") ) 
      sprintf(funcName, "eff_int_Pi0_Gaus_PbPb_%s_cen%i", pid.Data(), cent);
    else 
      sprintf(funcName, "eff_int_Pi0_Gaus_PbPb_%s_cen%i", pid.Data(), cent);
    TF1* func = dynamic_cast<TF1*> ( file->Get(funcName) );
    pastDir->cd();
    
    if( ! func )
      Printf("ERROR:GetEfficency: efficiancy function %s does not exist", funcName);
    return func;
  }
  return 0x0; // should not be reached
}

TH1* GetRawProduction(const RawProduction::Output& rawOutput, const TString& trigger="kMB", int fromCent=0, int toCent=10,
		      const TString& pid="All", const TString& graphName="yr1", Color_t color=kBlack, Style_t style=kFullDotSmall)
{
  TString newName = Form("raw_%s_%02i-%02i_%s_%s", trigger.Data(), fromCent, toCent, pid.Data(), graphName.Data());
  TH1* hist = dynamic_cast<TH1*> ( gDirectory->Get(newName.Data()) );

  if( ! hist ) {
    TString oldName(Form("%s/c%02i-%02i/%s/%s", trigger.Data(), fromCent, toCent, pid.Data(), graphName.Data()));
    TH1* hist = rawOutput.GetHistogram(oldName.Data());
    hist->SetName(newName.Data());
    hist->SetTitle(Form("%s, %02i-%02i%%, %s", trigger.Data(), fromCent, toCent, graphName.Data()));
    hist->GetXaxis()->SetTitle("p_{T}");
  }
  hist->SetLineColor(color);
  hist->SetMarkerColor(color);
  hist->SetMarkerStyle(style);
  return hist;
}

TH1* MakeProduction(const RawProduction::Output& rawOutput, const TString& trigger="kMB", int fromCent=0, int toCent=10,
		      const TString& pid="All", const TString& graphName="yr1", Color_t color=kBlack, Style_t style=kFullDotSmall)
{
  // First, check if production histogram allready exist in cd.
  TString name = Form("prod_%s_%02i-%02i_%s_%s", trigger.Data(), fromCent, toCent, pid.Data(), graphName.Data());
  TH1* hist = dynamic_cast<TH1*> ( gDirectory->Get(name.Data()) );

  if( ! hist ) { // else clone raw and correct for efficiancy
    const TH1* rawHist = GetRawProduction(rawOutput, trigger, fromCent, toCent, pid, graphName, color, style);
    hist = (TH1*) rawHist->Clone(name.Data());
    hist->Divide(GetEfficency(trigger, fromCent, toCent, pid, graphName));
    hist->GetYaxis()->SetTitle("#frac{d^{2}N_{#pi^{0}}}{p_{T}dp_{T}dy N_{ev}}");
  }
  hist->SetLineColor(color);
  hist->SetMarkerColor(color);
  hist->SetMarkerStyle(style);
  return hist;
}

TH1* MakeRatio(TH1* h1, TH1* h2, const TString& title ="")
{
  TString name = Form("%s_%s", h1->GetName(), h2->GetName());
  TH1* hist = dynamic_cast<TH1*> ( gDirectory->Get(name.Data()) );
  if( hist )
    return hist;

  hist = (TH1*)h1->Clone(name.Data());
  hist->Divide(h2);
  hist->SetTitle(title.Data());
  hist->GetYaxis()->SetTitle("Ratio");
  return hist;
}


TCanvas* DrawPIDProductionWithRatios(const RawProduction::Output& rawOutput, const TString& trigger,
		       int fromCent, int toCent, const TString& methode="yr1",
		       const TString& pids = TString("Allcore CPVcore Disp2core Both2core"), bool raw = false)
{
  const int capacity = 8;
  const int nHists = TMath::Min( pids.CountChar(' ')+1, capacity);
  TStringToken pidst = TStringToken(pids, " ");
  char pidsa[capacity][64] ={""};
  
  TH1* hists[capacity] = {0x0};
  const Style_t markers[capacity] = {22, 22, 23, 33, 24, 26, 32, 27};
  const Color_t colors[capacity] = {kBlack, kRed, kBlue, kGreen, kGray, kMagenta, kCyan, kOrange};
  
  int index = 0;
  while(pidst.NextToken() && index < capacity) {
    sprintf(pidsa[index], "%s", pidst.Data());
    if(raw)
      hists[index] = GetRawProduction(rawOutput, trigger, fromCent, toCent, pidst, methode, colors[index], markers[index]);
    else
      hists[index] = MakeProduction(rawOutput, trigger, fromCent, toCent, pidst, methode, colors[index], markers[index]);
    index++;
  }
  
  //TString pids_underscore = TString(pids); pids_underscore.ReplaceAll(" ", "_");
  TString key = Form("PIDRatios_%s_c%02i-%02i_%s_%s_raw%i_h%i", trigger.Data(), fromCent, toCent, methode.Data(), TString(pids).ReplaceAll(" ", "_").Data(), raw, canvHalfWidth);
  //TString key = Form("PIDRatios_%s_c%02i-%02i_%s_raw%i", trigger.Data(), fromCent, toCent, methode.Data(), raw);
  
  TCanvas* canv = 0x0;
  if( canvHalfWidth) 
    canv = new TCanvas(key.Data(), key.Data(), 1024/2, 768);
  else 
    canv = new TCanvas(key.Data(), key.Data(), 1024, 768);

  // Direct
  canv->cd();
  TPad *pad1 = new TPad("pad1","pad1",0,0.4,1,1);
  pad1->SetBottomMargin(0);
  pad1->SetLogy();
  pad1->Draw();
  pad1->cd();
  if( canvHalfWidth ) {
    hists[0]->GetYaxis()->SetLabelFont(63); //font in pixels
    hists[0]->GetYaxis()->SetLabelSize(10); //in pixels
    hists[0]->GetYaxis()->SetTitleSize(0.03);
    hists[0]->GetYaxis()->SetTitleOffset(1.2);
  } else {
    hists[0]->GetYaxis()->SetLabelFont(63); //font in pixels
    hists[0]->GetYaxis()->SetLabelSize(20); //in pixels
    hists[0]->GetYaxis()->SetTitleSize(0.055);
    hists[0]->GetYaxis()->SetTitleOffset(0.7);
  }
  if( raw )
    hists[0]->GetYaxis()->SetRangeUser(1.e-7, 1.e1);
  else
    hists[0]->GetYaxis()->SetRangeUser(1.e-5, 1.e3);
  // Draw
  hists[0]->DrawCopy();
  for(int i=1;i<nHists;++i)
    hists[i]->DrawCopy("same");
  // Legend
  TLegend* leg1 = new TLegend(0.7,0.6,0.85,0.88);
  for(int i=0;i<nHists;++i)
    leg1->AddEntry(hists[i], pidsa[i] , "lep");
  leg1->Draw();
  
  // Ratios
  canv->cd();
  TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.4);
  pad2->SetTopMargin(0);
  pad2->Draw();
  pad2->cd();
  pad2->SetGridy();
  TH1* firstRatio = MakeRatio(hists[1], hists[0]);
  if ( raw )
    firstRatio->GetYaxis()->SetRangeUser(0., 2.);
  else
    firstRatio->GetYaxis()->SetRangeUser(0.6, 1.4);
  firstRatio->GetYaxis()->SetTitle("");
  firstRatio->GetYaxis()->SetLabelFont(63); 
  firstRatio->GetYaxis()->SetLabelSize(25);
  firstRatio->GetXaxis()->SetLabelFont(63); 
  firstRatio->GetXaxis()->SetLabelSize(20);
  firstRatio->SetTitle("Ratio");
  // Draw
  firstRatio->DrawCopy("AXIS");
  firstRatio->DrawCopy("AXIGsame");
  firstRatio->DrawCopy("same");
  for(int i=2;i<nHists;++i)
    MakeRatio(hists[i],hists[0])->DrawCopy("same");
  // Ratios
  TLegend* leg2 = new TLegend(0.7,0.63,0.85,0.98);
  for(int i=1;i<nHists;++i)
    leg2->AddEntry(MakeRatio(hists[i],hists[0]), Form("%s/%s", pidsa[i], pidsa[0]), "lep");
  leg2->Draw();
  
  canv->SaveAs(Form("imgs/%s.pdf", key.Data()));
  canv->SaveAs(Form("imgs/%s.png", key.Data()));
  
  return canv;
  //delete canv;
}


void DrawPIDRatios(const RawProduction::Output& rawOutput)
{
  const int nCent = 2;
  int centBins[nCent][2] = {{0,5}, {5,10}/*, {10,20}, {20,40}, {40,60}, {60,80}*/};
  TStringToken methodes("yr1 yr1int", " ");
  while( methodes.NextToken() ) {
    for(int ic=0; ic<nCent; ++ic) {
      DrawPIDProductionWithRatios(rawOutput, "kCentral", centBins[ic][0], centBins[ic][1], methodes.Data(), "Allcore CPVcore Disp2core Both2core All CPV Disp2 Both2", true );
      DrawPIDProductionWithRatios(rawOutput, "kCentral", centBins[ic][0], centBins[ic][1], methodes.Data(), "Allcore CPVcore Disp2core Both2core All CPV Disp2 Both2", false );
      DrawPIDProductionWithRatios(rawOutput, "kCentral", centBins[ic][0], centBins[ic][1], methodes.Data(), "Allcore CPVcore Disp2core Both2core", true );
      DrawPIDProductionWithRatios(rawOutput, "kCentral", centBins[ic][0], centBins[ic][1], methodes.Data(), "Allcore CPVcore Disp2core Both2core", false );
      DrawPIDProductionWithRatios(rawOutput, "kCentral", centBins[ic][0], centBins[ic][1], methodes.Data(), "All Disp2 CPV Both2 Allcore", true );
      DrawPIDProductionWithRatios(rawOutput, "kCentral", centBins[ic][0], centBins[ic][1], methodes.Data(), "All Disp2 CPV Both2 Allcore", false );
    }
  }
}

void DrawProduction()
{
  gROOT->LoadMacro("MakeRawProduction.C+g");
  RawProduction::Output rawOutput;
  gStyle->SetOptStat(0);
  
  DrawPIDRatios(rawOutput);
}
