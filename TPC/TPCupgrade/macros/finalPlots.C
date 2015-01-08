/*
 *

 .L $ALICE_ROOT/TPC/Upgrade/macros/finalPlots.C+g
 finalPlots("eps10/medStat/*.debug.root","eps20/medStat/*.debug.root")

 finalPlots("eps10/medStat/*.debug.root","eps20/medStat/*.debug.root","/data/Work/software/svncern/papers/TDR/08-Monitoring_calib/figs")

 
*/

#include "TLegend.h"
#include "TCanvas.h"
#include "TString.h"
#include "TStyle.h"
#include "TColor.h"
#include "TPaveStats.h"
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TSystem.h"
#include "AliToyMCReconstruction.h"

void SetStyle();
TCanvas *GetCanvas(TString name, TString title, Float_t nx=1, Float_t ny=1);
void SaveCanvas(TCanvas *c);
void DrawOnTop(TPad *c, TObjArray &arrHists, Bool_t stats);

TString fSaveDir;

//void finalPlots(const char* filesEps10, const char* filesEps20, TString saveDir="")
void finalPlots(const char* filesEps20, TString saveDir="")
{
  fSaveDir=saveDir;
  
  TString idealUndistorted("t1_0_0_130_10.");
  TString idealDistorted("t1_1_3_130_10.");
  TString distorted("t0_1_0_130_10.");
  TString realTracking("t0_1_2_130_10.");
//  TString realTrackingPreT0("t0_1_4_130_10.");
  
  SetStyle();
  
//   TTree *tEps10=AliToyMCReconstruction::ConnectTrees(filesEps10);
  TTree *tEps20=AliToyMCReconstruction::ConnectTrees(filesEps20);

  TString test(tEps20->GetCurrentFile()->GetName());
  if (!test.Contains("0_0_0_")) {
    printf("ERROR: default file is not '0_0_0'\n");
    return;
  }

  if (/*tEps10->GetListOfFriends()->GetEntries()!=5 ||*/ tEps20->GetListOfFriends()->GetEntries()!=5) {
    printf("ERROR: wrong number of entries in the friends, not default\n");
    return;
  }

  TString drawStr;
  
  //
  // T0seed resolution
  //
  TCanvas *cT0res=GetCanvas("T0seedResolution","T0 seed resolution");
  //ideal undistorted
  TH1F *hT0resI = new TH1F("hT0resI","T0 resolution;(#it{t}_{0}^{seed}-#it{t}_{0}) #upoint #it{v}_{drift};#tracks",100,-50.1,50.1);
  drawStr=Form("(%sfTime0-t0)*vDrift",idealUndistorted.Data());
  tEps20->Draw(drawStr+">>hT0resI","","goff");
  //fully distorted
  TH1F *hT0resD = new TH1F("hT0resD","T0 resolution;(#it{t}_{0}^{seed}-#it{t}_{0}) #upoint #it{v}_{drift};#tracks",100,-50.1,50.1);
  hT0resD->SetLineColor(kRed);
  drawStr=Form("(%sfTime0-t0)*vDrift",distorted.Data());
  tEps20->Draw(drawStr+">>hT0resD","","goff");
  //distorted and average correction
  TH1F *hT0resDC = new TH1F("hT0resDC","T0 resolution;(#it{t}_{0}^{seed}-#it{t}_{0}) #upoint #it{v}_{d} (cm);#tracks",100,-50.1,50.1);
//   hT0resDC->SetLineColor(kGreen+2);
  drawStr=Form("(%sfTime0-t0)*vDrift",realTracking.Data());
  tEps20->Draw(drawStr+">>hT0resDC","","goff");

//  TH1F *hT0resDCPreT0 = new TH1F("hT0resDCPreT0","T0 resolution;(#it{t}_{0}^{seed}-#it{t}_{0}) #upoint #it{v}_{drift};#tracks",100,-50.1,50.1);
//  drawStr=Form("(%sfTime0-t0)*vDrift",realTrackingPreT0.Data());
//  tEps20->Draw(drawStr+">>hT0resDCPreT0","","goff");
  
  //   hT0resI->Draw();
  //   hT0resD->Draw("same");
  hT0resDC->Draw(/*"same"*/);
  //hT0resDCPreT0->Draw(/*"same"*/);
  
  SaveCanvas(cT0res);

  //
  // Track parameter resolution (y) with ideal clusters at the ITS and inner wall of the TPC
  //
  TCanvas *cYresComparison=GetCanvas("YresComparison","Comparison of Yres for ideal clusters");
  //ideal clusters at the ITS outermost point
  TH1F *hYresITS = new TH1F("hYresITS",";#it{y}_{TPC}-#it{y}_{ITS} (cm);#tracks",100,-0.21,0.21);
  drawStr=Form("%strackITS.fP[0]-%stRealITS.fP[0]",idealUndistorted.Data(),idealUndistorted.Data());
  tEps20->Draw(drawStr+">>hYresITS","","goff");
  hYresITS->SetLineColor(kRed);

  TH1F *hYresTPC = new TH1F("hYresTPC",";#it{y}_{TPC}-#it{y}_{ITS} (cm);#tracks",100,-0.21,0.21);
  drawStr=Form("%strackITS2.fP[0]-%stRealITS2.fP[0]",idealUndistorted.Data(),idealUndistorted.Data());
  tEps20->Draw(drawStr+">>hYresTPC","","goff");

  hYresTPC->Draw();
  hYresITS->Draw("same");

  SaveCanvas(cYresComparison);

  //
  //  Track parameter resolution (y) with fully distorted clusters at the inner wall of the TPC
  //

  TCanvas *cYresDistorted=GetCanvas("YresDistorted","Yres for fully distorted clusters");
  //ideal clusters at the ITS outermost point
  TH1F *hYresDist = new TH1F("hYresDist",";#it{y}_{TPC}-#it{y}_{ITS} (cm);#tracks",100,-15.5,5.5);
  drawStr=Form("%strackITS2.fP[0]-%stRealITS2.fP[0]",distorted.Data(),distorted.Data());
  tEps20->Draw(drawStr+">>hYresDist","","goff");
  
  hYresDist->Draw();
  
  SaveCanvas(cYresDistorted);
  
  //
  //  Track parameter resolution (y) with fully distorted and corrected clusters (Tzero seed)
  //  at the inner wall of the TPC
  //
  
  TCanvas *cYresDistCorrTzeroSeed=GetCanvas("YresDistCorrTzeroSeed","Yres for fully distorted/corrected clusters (Tzero seed)");
  //ideal clusters at the ITS outermost point
  TH1F *hYresDistCorrTzeroSeed = new TH1F("hYresDistCorrTzeroSeed",";#it{y}_{TPC}-#it{y}_{ITS} (cm);#tracks",100,-.85,1.15);
  drawStr=Form("%strackITS2.fP[0]-%stRealITS2.fP[0]",realTracking.Data(),realTracking.Data());
  tEps20->Draw(drawStr+">>hYresDistCorrTzeroSeed","","goff");
  
  hYresDistCorrTzeroSeed->Draw();
  
  SaveCanvas(cYresDistCorrTzeroSeed);

  //
  //  Track parameter resolution (y) with fully distorted and corrected clusters (Tzero seed)
  //  at the inner wall of the TPC
  //
  
  TCanvas *cYresDistCorrTzero=GetCanvas("YresDistCorrTzero","Yres for fully distorted/corrected clusters (Tzero)");
  //ideal clusters at the ITS outermost point
  TH1F *hYresDistCorrTzero = new TH1F("hYresDistCorrTzero",";#it{y}_{TPC}-#it{y}_{ITS} (cm);#tracks",100,-.5,.85);
  drawStr=Form("%strackITS2.fP[0]-%stRealITS2.fP[0]",idealDistorted.Data(),idealDistorted.Data());
  tEps20->Draw(drawStr+">>hYresDistCorrTzero","","goff");
  
  hYresDistCorrTzero->Draw();
  
  SaveCanvas(cYresDistCorrTzero);


  //
  // plot all params
  //

  TString titles[5]={"#it{y}_{TPC}-#it{y}_{ITS} (cm)","#it{z}_{TPC}-#it{z}_{ITS} (cm)","sin(#it{#alpha})_{TPC}-sin(#it{#alpha})_{ITS}","tan(#it{#lambda})_{TPC}-tan(#it{#lambda})_{ITS}","1/#it{p}_{T TPC}-1/#it{p}_{T ITS} ((GeV/#it{c})^{-1})"};
  //Double_t min[5]={-.85,-15,-.009,-.005,-.05};
  //Double_t max[5]={ .85, 15, .009, .005, .05};
  Double_t min[5]={-.85,-15,-.015,-.005,-.075};
  Double_t max[5]={1.15, 15, .015, .005, .075};
  TString type[3]={idealUndistorted,idealDistorted,realTracking};
  Int_t colors[3]={kBlack,kGreen-2,kRed};

  TLegend *leg=new TLegend(.1,.55,.95,.95);
  leg->SetTextSize(0.075);
  leg->SetBorderSize(1);
  leg->SetFillColor(10);
  TCanvas *cResParams=GetCanvas("ResParams","Resolution of parameters",1.2,1.8);
  cResParams->Divide(2,3);
  for (Int_t i=0;i<5;++i){
    TPad *pad=(TPad*)cResParams->cd(i+1);
    TObjArray arr;
    for (Int_t it=0; it<3; ++it) {
      TH1F *hResParams=new TH1F(Form("hResParams_%d_%d",i,it),
                                Form(";%s;#tracks",titles[i].Data()),
                                100,min[i],max[i]);
      drawStr=Form("%strackITS2.fP[%d]-tRealITS2.fP[%d]",type[it].Data(),i,i);
      tEps20->Draw(drawStr+Form(">>hResParams_%d_%d",i,it),"","goff");
      hResParams->SetLineColor(colors[it]);
      arr.Add(hResParams);
    }
    if (i==0) {
      leg->AddEntry(arr.At(0),"no distortions (ideal)","l");
      leg->AddEntry(arr.At(1),"distorted/corrected (t_{0})","l");
      leg->AddEntry(arr.At(2),"distorted/corrected (t_{0}^{seed})","l");
    }
    DrawOnTop(pad,arr,kTRUE);
  }

  cResParams->cd(6);
  leg->Draw();
  SaveCanvas(cResParams);

  TCanvas *cResRPhi=GetCanvas("ResRPhi","Resolution of rPhi");
  TLegend *leg2=new TLegend(.12,.7,.48,.95);
  TObjArray arr;
  Int_t i=0;
  for (Int_t it=0; it<3; ++it) {
    TH1F *hResParams=new TH1F(Form("hResParams_%d_%d",i,it),
                              Form(";%s;#tracks",titles[i].Data()),
                              100,min[i],max[i]);
    drawStr=Form("%strackITS2.fP[%d]-tRealITS2.fP[%d]",type[it].Data(),i,i);
    tEps20->Draw(drawStr+Form(">>hResParams_%d_%d",i,it),"","goff");
    hResParams->SetLineColor(colors[it]);
    arr.Add(hResParams);
  }
  leg2->AddEntry(arr.At(0),"no distortions (ideal)","l");
  leg2->AddEntry(arr.At(1),"distorted/corrected (t_{0})","l");
  leg2->AddEntry(arr.At(2),"distorted/corrected (t_{0}^{seed})","l");
  DrawOnTop(cResRPhi,arr,kTRUE);
  leg2->Draw("same");
  SaveCanvas(cResRPhi);
  
}


TCanvas *GetCanvas(TString name, TString title, Float_t nx, Float_t ny)
{
  TCanvas *c=(TCanvas*)gROOT->GetListOfCanvases()->FindObject(name.Data());
  if (!c) c=new TCanvas(name,title,nx*700,ny*500);
  c->Clear();
  c->cd();
  return c;
}

void SaveCanvas(TCanvas *c)
{
  //
  //
  //
  
  if (fSaveDir.IsNull()) return;
  
  c->SaveAs(Form("/tmp/%s.eps",c->GetName()));
  c->SaveAs(Form("%s/%s.png",fSaveDir.Data(),c->GetName()));
  gSystem->Exec(Form("ps2pdf -dEPSCrop /tmp/%s.eps %s/%s.pdf",c->GetName(),fSaveDir.Data(),c->GetName()));
}

void DrawOnTop(TPad *c, TObjArray &arrHists, Bool_t /*stats*/)
{
  Double_t min=0,max=0;
  Double_t ystatMax=gStyle->GetStatY();
  Double_t ystatH  =gStyle->GetStatH()*2./3.;
  const Int_t nHists=arrHists.GetEntriesFast();
  for (Int_t iHist=0; iHist<nHists; ++iHist) {
    TH1 *h=(TH1*)arrHists.UncheckedAt(iHist);
    TPad *pad = new TPad(Form("%s_%d",c->GetName(),iHist+1),"",0,0,1,1);
    pad->SetFillStyle(4000);
    pad->SetFrameFillStyle(0);
    pad->Draw();
    pad->cd();
    if (iHist>0) {
      h->SetMinimum(min);
      h->SetMaximum(max);
      pad->SetTicky(0);
      pad->SetTickx(0);
    }
    h->Draw((iHist==0)?"":"AH");
    pad->Update();
    if (iHist==0){
      min=pad->GetUymin();
      max=pad->GetUymax();
      printf("min: %.2f %.2f\n",min,max);
    }
    TPaveStats *ps = (TPaveStats*)h->GetListOfFunctions()->FindObject("stats");
    if (!ps) printf("shitttt %s\n",h->GetName());
    else {
    ps->SetTextColor(h->GetLineColor());
    ps->SetY2NDC(ystatMax-iHist*ystatH);
    ps->SetY1NDC(ystatMax-(iHist+1)*ystatH);}
  }
}

void SetStyle()
{
  const Int_t NCont=255;
  //const Int_t NCont=50;
  TH1::AddDirectory();
  TStyle *st = new TStyle("mystyle","mystyle");
  gROOT->GetStyle("Plain")->Copy((*st));
  st->SetTitleX(0.1);
  st->SetTitleW(0.8);
  st->SetTitleH(0.08);
  st->SetStatX(.95);
  st->SetStatY(.95);
  st->SetStatW(.25);
  st->SetStatH(.25);
  st->SetNumberContours(NCont);
  st->SetPalette(1,0);
  st->SetOptStat("rm");
  st->SetOptTitle(0);
  st->SetOptFit(0);
  st->SetGridColor(kGray+1);
//   st->SetPadGridX(kTRUE);
//   st->SetPadGridY(kTRUE);
  st->SetPadTickX(kTRUE);
  st->SetPadTickY(kTRUE);
  st->SetMarkerStyle(20);
  st->SetMarkerSize(.5);

  st->SetPadLeftMargin(0.12);
  st->SetPadBottomMargin(0.12);
  st->SetPadRightMargin(0.05);
  st->SetPadTopMargin(0.05);
  st->cd();
  
  Int_t nimTPCFont=42; //or 62 for sans serif font
  //default definitions
  st->SetTextFont(nimTPCFont);
  st->SetTitleFont(nimTPCFont, "T");
  st->SetTitleFont(nimTPCFont, "XYZ");
  st->SetLabelFont(nimTPCFont,"XYZ");
  st->SetLabelSize(0.045,"XYZ");
  st->SetTitleSize(0.05,"XYZ");
  st->SetTitleOffset(1.1,"XZ");
  st->SetTitleOffset(1.3,"Y");
  st->SetStatFont(nimTPCFont);
  st->SetOptTitle(0);
  st->SetPalette(1,0);
  st->SetStatBorderSize(1);
  new TColor(2001,1,1,1);
  st->SetFillColor(2001);
  st->SetTickLength(gStyle->GetTickLength()/696.*472.,"y");
  
  const Int_t NRGBs = 5;
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  //grey
  //  Double_t       stops[5] = {0.00, 0.34, 0.61, 0.84, 1.00};
  //  Double_t         red[5]   = {1.00, 0.84, 0.61, 0.34, 0.00};
  //  Double_t         green[5]   = {1.00, 0.84, 0.61, 0.34, 0.00};
  //  Double_t         blue[5]   = {1.00, 0.84, 0.61, 0.34, 0.00};
  
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  
  st->cd();
}

