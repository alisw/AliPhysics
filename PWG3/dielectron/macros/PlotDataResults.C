#include <TROOT.h>
#include <TStyle.h>
#include <TH1D.h>
#include <TObjArray.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TAxis.h>

#include "AliDielectronSignalExt.h"
#include "AliDielectronCFdraw.h"
#include "AliDielectron.h"

/*
gSystem->AddIncludePath("-I/data/Work/software/svngsi/dielectron/trunk/dielectron");
.L PlotDataResults.C+g
PlotDataResuts("/data/Work/train/V005.data/2010-10-29_2144.3552/mergedPeriods/data/7TeV/LHC10c.pass2/jpsi.root");
PlotDataResuts("/data/Work/train/V005.data/2010-10-21_2342.3445/mergedPeriods/data/7TeV/LHC10pass2/jpsi.root")
PlotDataResuts("/data/Work/train/V005.data/2010-10-21_2342.3445/mergedPeriods/data/7TeV/LHC10d.pass1/wiechula_jpsi.root");
*/

AliDielectronSignalBase* GetSignalLS(AliDielectronCFdraw &d, Int_t step, const char* nameAdd);
AliDielectronSignalBase* GetSignalRot(AliDielectronCFdraw &d, Int_t step, const char* nameAdd);
void DrawSpectra(AliDielectronSignalBase *sig, const char* cname);

//_______________________________________
void PlotDataResuts(const char* filename)
{
  AliDielectronCFdraw d(filename);
  
  Int_t stepFirst=0, stepAny=1;
  
  gStyle->SetOptStat(0);
  //Set common Ranges
//   d.SetRangeUser("Leg1_Pt",1.,1000.);
//   d.SetRangeUser("Leg2_Pt",1.,1000.);
  
  //============================
  //SPD first
  //

  //--- Like sign subtraction
  AliDielectronSignalBase *sigFirst=GetSignalLS(d,stepFirst,"First");
  DrawSpectra(sigFirst,"cFirst");
  //--- Rotation subtraction
  AliDielectronSignalBase *sigFirstRot=GetSignalRot(d,stepFirst,"FirstRot");
  DrawSpectra(sigFirstRot,"cFirstRot");
  
  //============================
  //SPD any
  //
  AliDielectronSignalBase *sigAny=GetSignalLS(d,stepAny,"Any");
  DrawSpectra(sigAny,"cAny");
  //--- Rotation subtraction
  AliDielectronSignalBase *sigAnyRot=GetSignalRot(d,stepAny,"AnyRot");
  DrawSpectra(sigAnyRot,"cAnyRot");
  
}


//_______________________________________
AliDielectronSignalBase *GetSignalLS(AliDielectronCFdraw &d, Int_t step, const char* nameAdd)
{
  //
  // Get Extracted signal from likesign method
  //
  
  TObjArray *arr=new TObjArray;
  arr->SetOwner();

  for (Int_t iType=0;iType<3;++iType){
    d.SetRangeUser("PairType",iType,iType);
    arr->AddAt(d.Project("M",step),iType);
  }

  AliDielectronSignalExt *sig=new AliDielectronSignalExt;
  sig->SetIntegralRange(2.9,3.15);
  sig->SetMethod(AliDielectronSignalBase::kLikeSign);
  sig->Process(arr);
  
  delete arr;
  return sig;
}

//_______________________________________
AliDielectronSignalBase *GetSignalRot(AliDielectronCFdraw &d, Int_t step, const char* nameAdd)
{
  //
  // Get Extracted signal from likesign method
  //
  
  TObjArray *arr=new TObjArray;
  arr->SetOwner();

  Int_t iType=AliDielectron::kEv1PM;
  d.SetRangeUser("PairType",iType,iType);
  arr->AddAt(d.Project("M",step),iType);
  
  iType=AliDielectron::kEv1PMRot;
  d.SetRangeUser("PairType",iType,iType);
  arr->AddAt(d.Project("M",step),iType);
  
  AliDielectronSignalExt *sig=new AliDielectronSignalExt;
  sig->SetIntegralRange(2.9,3.15);
  sig->SetMethod(AliDielectronSignalBase::kRotation);
  sig->Process(arr);
  
  delete arr;
  return sig;
}

//_______________________________________
void DrawSpectra(AliDielectronSignalBase *sig, const char* cname)
{
  //
  //
  //
  gStyle->SetOptTitle(0);
  TCanvas *c=(TCanvas*)gROOT->FindObject(cname);
  if (!c) c=new TCanvas(cname,cname,400,600);
  c->Clear();
  c->Divide(1,2,0,0);
  c->cd(1);
  gPad->SetTopMargin(0.01);
  gPad->SetRightMargin(0.01);
  
  gPad->SetBottomMargin(0);
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetTickx();
  gPad->SetTicky();

  TH1 *hUS=sig->GetUnlikeSignHistogram();
  hUS->SetMarkerStyle(20);
  hUS->SetMarkerSize(0.7);
  hUS->SetMarkerColor(kRed);
  hUS->SetLineColor(kRed);
  hUS->SetStats(0);
  
  TH1* hBackground=sig->GetBackgroundHistogram();
  hBackground->SetMarkerStyle(24);
  hBackground->SetMarkerSize(0.7);
  hBackground->SetStats(0);
  hBackground->SetMarkerColor(kBlue);
  hBackground->SetLineColor(kBlue);
  
  hUS->Draw();
  hBackground->Draw("same");
  
  c->cd(2);
  gPad->SetRightMargin(0.01);
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetTopMargin(0);

  TH1* hSignal=sig->GetSignalHistogram();
  hSignal->SetMarkerStyle(20);
  hSignal->SetMarkerSize(0.7);
  hSignal->SetMarkerColor(kRed);
  hSignal->SetLineColor(kRed);
  hSignal->Draw();
}

