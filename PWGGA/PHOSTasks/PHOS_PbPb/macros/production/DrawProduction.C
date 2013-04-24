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
#include <cfloat>

namespace RawProduction {
  class Output;
}

enum EffVersion { LHC10h_1234_Apr_10 };
EffVersion effVersion = LHC10h_1234_Apr_10;
bool canvHalfWidth=false;
int maxFailedCombined = 0;

TH1* MakeCombinedMethodeProduction(const RawProduction::Output& rawOutput, const TString& trigger="kMB", int fromCent=0, int toCent=10, const TString& pid="All", const TString& graphName="yr1", Color_t color=kBlack, Style_t style=kFullDotSmall);
TH1* MakeCombinedProduction(const RawProduction::Output& rawOutput, const TString& trigger="kMB", int fromCent=0, int toCent=10, const TString& pid="All", const TString& graphName="yr1", Color_t color=kBlack, Style_t style=kFullDotSmall);

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
		      const TString& pid="All", const TString& methode="yr1", Color_t color=kBlack, Style_t style=kFullDotSmall)
{  
  // First, check if production histogram allready exist in cd.
  TString name = Form("prod_%s_%02i-%02i_%s_%s", trigger.Data(), fromCent, toCent, pid.Data(), methode.Data());
  TH1* hist = dynamic_cast<TH1*> ( gDirectory->Get(name.Data()) );
  
  // Check if combined.
  if( ! hist ) hist = MakeCombinedProduction(rawOutput, trigger, fromCent, toCent, pid, methode, color, style);
  if( ! hist ) hist = MakeCombinedMethodeProduction(rawOutput, trigger, fromCent, toCent, pid, methode, color, style);
  

  // else clone raw and correct for efficiancy
  if( ! hist ) { 
    const TH1* rawHist = GetRawProduction(rawOutput, trigger, fromCent, toCent, pid, methode, color, style);
    hist = (TH1*) rawHist->Clone(name.Data());
    hist->Divide(GetEfficency(trigger, fromCent, toCent, pid, methode));
    hist->GetYaxis()->SetTitle("#frac{d^{2}N_{#pi^{0}}}{p_{T}dp_{T}dy N_{ev}}");
  }
  hist->SetLineColor(color);
  hist->SetMarkerColor(color);
  hist->SetMarkerStyle(style);
  return hist;
}

TH1* CombinePID(const RawProduction::Output& rawOutput, const TString& trigger, int fromCent, int toCent, const TString& combinedPidName, const TString& pids, const TString& methode, Color_t color, Style_t style)
{
  // First, check if production histogram allready exist in cd.
  TString name = Form("prod_%s_%02i-%02i_%s_%s", trigger.Data(), fromCent, toCent, combinedPidName.Data(), methode.Data());
  TH1* hist = dynamic_cast<TH1*> ( gDirectory->Get(name.Data()) );
  if( hist ) return hist;
  
  const int capacity = 16;
  const int nHists = TMath::Min( pids.CountChar(' ')+1, capacity);
  TH1* hists[capacity];
  TStringToken pidst(pids, " ");
  pidst.NextToken();
  const TString firstPID = pidst;
  hists[0] = MakeProduction(rawOutput, trigger, fromCent, toCent, firstPID, methode, color, style);
  int index = 1;
  while( pidst.NextToken() && index < capacity ) {
    hists[index] = MakeProduction(rawOutput, trigger, fromCent, toCent, pidst, methode, color, style);
    ++index;
  }
  TH1* combHist = hists[0]->Clone(name.Data());
  combHist->SetTitle(TString(combHist->GetTitle()).ReplaceAll(firstPID.Data(), combinedPidName.Data()));
  
  for(int ptBin=0; ptBin<combHist->GetNbinsX(); ++ptBin){
    double wmeansum = 0.;
    double sumw = 0.;
    double mins = DBL_MAX;
    int nFailed = 0;
    for(int i=0; i<nHists; ++i) {
      const double x = hists[i]->GetBinContent(ptBin);
      const double s = hists[i]->GetBinError(ptBin);
      if( 0.==x || 0.==s) {
	nFailed++;
	continue;
      }
      const double w = 1./(s*s);
      wmeansum += w*x;
      sumw +=w;
      if( mins > s )
	mins = s;
    }
    if(nFailed > maxFailedCombined) {
      combHist->SetBinContent(ptBin, 0.);
      combHist->SetBinError(ptBin, 0.);
    }
    else {
      const double wmean = wmeansum/sumw;
      combHist->SetBinContent(ptBin, wmean);
      combHist->SetBinError(ptBin, mins);
    }
  }

  return combHist;
}


TH1* MakeCombinedProduction(const RawProduction::Output& rawOutput, const TString& trigger, int fromCent, int toCent, const TString& combinedPidName, const TString& methode, Color_t color, Style_t style)
{
  TString pids;
  if( combinedPidName.Contains("Combcore") )
    pids = "Allcore CPVcore Disp2core Both2core";
  else
    return 0x0;
  
  // First, check if production histogram allready exist in cd.
  TString name = Form("prod_%s_%02i-%02i_%s_%s", trigger.Data(), fromCent, toCent, combinedPidName.Data(), methode.Data());
  TH1* hist = dynamic_cast<TH1*> ( gDirectory->Get(name.Data()) );
  if( hist ) 
    return hist;

  return CombinePID(rawOutput, trigger, fromCent, toCent, combinedPidName, pids, methode, color, style);
}

TH1* MergeMethodes(const RawProduction::Output& rawOutput, const TString& trigger, int fromCent, int toCent, const TString& pid, const TString& combMethodeName, const TString& methodes, Color_t color, Style_t style)
{
  // First, check if production histogram allready exist in cd.
  TString name = Form("prod_%s_%02i-%02i_%s_%s", trigger.Data(), fromCent, toCent, pid.Data(), combMethodeName.Data());
  TH1* hist = dynamic_cast<TH1*> ( gDirectory->Get(name.Data()) );
  if( hist ) return hist;
  
  const int capacity = 16;
  const int nHists = TMath::Min( methodes.CountChar(' ')+1, capacity);
  TH1* hists[capacity];
  TStringToken mst(methodes, " ");
  mst.NextToken();
  const TString firstMethode = mst;
  hists[0] = MakeProduction(rawOutput, trigger, fromCent, toCent, pid, firstMethode, color, style);
  int index = 1;
  while( mst.NextToken() && index < capacity ) {
    hists[index] = MakeProduction(rawOutput, trigger, fromCent, toCent, pid, mst, color, style);
    ++index;
  }
  TH1* combHist = hists[0]->Clone(name.Data());
  combHist->SetTitle(TString(combHist->GetTitle()).ReplaceAll(firstMethode.Data(),combMethodeName.Data()));
  
  for(int ptBin=0; ptBin<combHist->GetNbinsX(); ++ptBin){
    double wmeansum = 0.;
    double sumw = 0.;
    double mins = DBL_MAX;
    int nFailed = 0;
    for(int i=0; i<nHists; ++i) {
      const double x = hists[i]->GetBinContent(ptBin);
      const double s = hists[i]->GetBinError(ptBin);
      if( 0.==x || 0.==s) {
	nFailed++;
	continue;
      }
      const double w = 1./(s*s);
      wmeansum += w*x;
      sumw +=w;
      if( mins > s )
	mins = s;
    }
    if(nFailed > maxFailedCombined) {
      combHist->SetBinContent(ptBin, 0.);
      combHist->SetBinError(ptBin, 0.);
    }
    else {
      const double wmean = wmeansum/sumw;
      combHist->SetBinContent(ptBin, wmean);
      combHist->SetBinError(ptBin, mins);
    }
  }

  return combHist;
}

TH1* MakeCombinedMethodeProduction(const RawProduction::Output& rawOutput, const TString& trigger, int fromCent, int toCent, const TString& pid, const TString& combMethodeName, Color_t color, Style_t style)
{  
  TString methodes;
  if( combMethodeName.EqualTo("yr1comb") )
    methodes = "yr1 yr1int";
  else 
    return 0x0;
 
  // First, check if production histogram allready exist in cd.
  TString name = Form("prod_%s_%02i-%02i_%s_%s", trigger.Data(), fromCent, toCent, pid.Data(), combMethodeName.Data());
  TH1* hist = dynamic_cast<TH1*> ( gDirectory->Get(name.Data()) );
  if( hist ) 
    return hist;

  return MergeMethodes(rawOutput, trigger, fromCent, toCent, pid, combMethodeName, methodes, color, style);
}

TH1* MakeRatio(TH1* h1, TH1* h2, const TString& title ="")
{
  TString name = Form("%s_%s", h1->GetName(), h2->GetName());
  TH1* hist = dynamic_cast<TH1*> ( gDirectory->Get(name.Data()) );
  if( hist )
    return hist;

  hist = (TH1*)h1->Clone(name.Data());
  hist->Divide(h1, h2, 1, 1, "B");
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


TCanvas* DrawMethodeProductionWithRatios(const RawProduction::Output& rawOutput, const TString& trigger,
		       int fromCent, int toCent, const TString& methodes="yr1comb yr1 yr1int",
		       const TString& pid = "Allcore", bool raw = false)
{
  const int capacity = 8;
  const int nHists = TMath::Min( methodes.CountChar(' ')+1, capacity);
  TStringToken methodest = TStringToken(methodes, " ");
  char methodesa[capacity][64] ={""};
  
  TH1* hists[capacity] = {0x0};
  const Style_t markers[capacity] = {22, 22, 23, 33, 24, 26, 32, 27};
  const Color_t colors[capacity] = {kBlack, kRed, kBlue, kGreen, kGray, kMagenta, kCyan, kOrange};
  
  int index = 0;
  while(methodest.NextToken() && index < capacity) {
    sprintf(methodesa[index], "%s", methodest.Data());
    if(raw)
      hists[index] = GetRawProduction(rawOutput, trigger, fromCent, toCent, pid, methodest, colors[index], markers[index]);
    else
      hists[index] = MakeProduction(rawOutput, trigger, fromCent, toCent, pid, methodest, colors[index], markers[index]);
    index++;
  }
  
  TString key = Form("MethodeRatios_%s_c%02i-%02i_%s_%s_raw%i_h%i", trigger.Data(), fromCent, toCent, TString(methodes).ReplaceAll(" ", "_").Data(), pid.Data(), raw, canvHalfWidth);
  
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
    leg1->AddEntry(hists[i], methodesa[i] , "lep");
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
    leg2->AddEntry(MakeRatio(hists[i],hists[0]), Form("%s/%s", methodesa[i], methodesa[0]), "lep");
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

void DrawPIDRatiosCombined(const RawProduction::Output& rawOutput)
{
  const int nCent = 2;
  int centBins[nCent][2] = {{0,5}, {5,10}/*, {10,20}, {20,40}, {40,60}, {60,80}*/};
  TStringToken methodes("yr1comb", " ");
  while( methodes.NextToken() ) {
    for(int ic=0; ic<nCent; ++ic) {
      DrawPIDProductionWithRatios(rawOutput, "kCentral", centBins[ic][0], centBins[ic][1], methodes.Data(), "Allcore CPVcore Disp2core Both2core All CPV Disp2 Both2", false );
      DrawPIDProductionWithRatios(rawOutput, "kCentral", centBins[ic][0], centBins[ic][1], methodes.Data(), "Allcore CPVcore Disp2core Both2core", false );
      DrawPIDProductionWithRatios(rawOutput, "kCentral", centBins[ic][0], centBins[ic][1], methodes.Data(), "Combcore Allcore CPVcore Disp2core Both2core", false );
      DrawPIDProductionWithRatios(rawOutput, "kCentral", centBins[ic][0], centBins[ic][1], methodes.Data(), "All Disp2 CPV Both2 Allcore", false );
    }
  }
}

void DrawMethodeRatios(const RawProduction::Output& rawOutput)
{
  const int nCent = 2;
  int centBins[nCent][2] = {{0,5}, {5,10}/*, {10,20}, {20,40}, {40,60}, {60,80}*/};
  TStringToken pids("Allcore CPVcore Disp2core Both2core", " ");
  while( pids.NextToken() ) {
    for(int ic=0; ic<nCent; ++ic) {
      DrawMethodeProductionWithRatios(rawOutput, "kCentral", centBins[ic][0], centBins[ic][1], "yr1comb yr1 yr1int", pids, false );
    }
  }
}

void DrawProduction()
{
  gROOT->LoadMacro("MakeRawProduction.C+g");
  RawProduction::Output rawOutput;
  gStyle->SetOptStat(0);
  
  DrawMethodeRatios(rawOutput);
  DrawPIDRatiosCombined(rawOutput);
  
  canvHalfWidth = true;
  DrawPIDRatios(rawOutput);
}
