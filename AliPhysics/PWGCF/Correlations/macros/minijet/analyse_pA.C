/* $Id:  $ */
//--------------------------------------------------
//
// post-processing macro for the minijet analysis 
// uses input of analysis class AliAnalysisTaskPhiCorrelation
//
//  Author : Emilia Leogrande (Utrecht University)
//           emilia.leogrande@cern.ch
//
//-------------------------------------------------

#include <TChain.h>
#include <TList.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TCanvas.h>
#include "TRandom.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TF1.h"
#include "TMath.h"
#include "TDirectory.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TColor.h"

#include <iostream>
using namespace std;

TFile *file = TFile::Open("dphi_corr_V0A_Train445_0707_wingsCorrected.root");
TFile *output = TFile::Open("minijets_results.root","recreate");

const int centralityClasses = 20;
const Double_t Deta_max = 1.8;
const Double_t Deta_gap = 1.2;
const Double_t Dphi_binwidth = TMath::TwoPi()/72;
const Double_t Dphi_zyam_notsub = 1.1;
const Double_t Dphi_ns_range = 1.45;

TH2D *correlations[centralityClasses];
TH1D *triggers;
TH1D *events;
TH2D *correlationsME[centralityClasses];
TH1D *fDeltaPhiForZyam[centralityClasses];
TH1D *fLongRangeForZyamNotSub[centralityClasses];
TH1D *fDeltaPhiLRSub[centralityClasses];
TH1D *fDeltaPhiLRNotSub[centralityClasses];
TH1D *fFinalLRSub[centralityClasses];
TH1D *fFinalLRNotSub[centralityClasses];

Double_t zyam = 0, err_zyam = 0;
Double_t zyamNoSub = 0, err_zyamNoSub = 0;
Double_t ns_yield = 0, err_ns = 0;
Double_t nsNoSub_yield = 0, err_nsNoSub = 0;
Double_t as_yield = 0, err_as = 0;
Double_t asNoSub_yield = 0, err_asNoSub = 0;

TGraphErrors *NearSide = new TGraphErrors();
TGraphErrors *NearSideNotSub = new TGraphErrors();
TGraphErrors *AwaySide = new TGraphErrors();
TGraphErrors *AwaySideNotSub = new TGraphErrors();
TGraphErrors *AverageTriggers = new TGraphErrors();
TGraphErrors *UncorrelatedSeeds = new TGraphErrors();

TH2D *MERemoval(TH2D*, Int_t, TF1*);
TF1 *fTriangle(Double_t *, Double_t *);
TH1D *ZyamLRsubtraction(TH2D *, Int_t, Double_t, Bool_t);
TH1D *LRsubtraction(TH2D *, Int_t, Double_t, Double_t, Bool_t);
Double_t GetZyam(TH1D *, Bool_t, const Double_t);
Double_t GetZyamError(TH1D *, Bool_t, const Double_t);
TH1D *BaselineSubtraction(TH1D *, Int_t, Bool_t, Double_t, Double_t, Double_t);
Double_t YieldNS(TH1D *, const Double_t, Bool_t);
Double_t YieldAS(TH1D *, const Double_t, Bool_t, Bool_t, Double_t, Double_t);
void GraphProperties(TGraphErrors *, const char *, Int_t, Int_t);
void WriteIntoFile(TFile *, const char *, TGraphErrors *);

void analyse_pA(){
  
  for(Int_t i=0; i<centralityClasses; i++){
    correlations[i] = (TH2D*)file->Get(Form("dphi_0_0_%i",i));
  }
  triggers = (TH1D*)file->Get("triggers_0");
  events = (TH1D*)file->Get("events");

  TH1D *p0 = (TH1D*)triggers->Clone();
  TH1D *p1 = (TH1D*)events->Clone();
  p0->Sumw2();
  p1->Sumw2();
  p0->Divide(p0,p1,1,1,"B");

  TF1 *mixedEvent = new TF1("mixedEvent", fTriangle, -Deta_max, Deta_max, 4); 
  mixedEvent->FixParameter(0,1);
  mixedEvent->FixParameter(1,1./Deta_max);
  mixedEvent->FixParameter(2,1);
  mixedEvent->FixParameter(3,-1./Deta_max);

  Double_t scale_ns = mixedEvent->Integral(-Deta_gap,Deta_gap)/(mixedEvent->Integral(-Deta_max,-Deta_gap)+mixedEvent->Integral(Deta_gap,Deta_max));
  Double_t scale_as = mixedEvent->Integral(-Deta_max,Deta_max)/(mixedEvent->Integral(-Deta_max,-Deta_gap)+mixedEvent->Integral(Deta_gap,Deta_max));
  Int_t points = 0;

  for(Int_t i = 0; i < centralityClasses; i++){
    correlationsME[i] = MERemoval(correlations[i], i, mixedEvent);
    correlationsME[i]->Scale(1./triggers->GetBinContent(i+1));

    // Long Range subtraction
    fDeltaPhiForZyam[i] = ZyamLRsubtraction(correlationsME[i], i, scale_as, kTRUE);
    zyam = GetZyam(fDeltaPhiForZyam[i], kTRUE, Dphi_zyam_notsub);
    err_zyam = GetZyamError(fDeltaPhiForZyam[i], kTRUE, Dphi_zyam_notsub);
    fDeltaPhiLRSub[i] = LRsubtraction(correlationsME[i], i, scale_ns, scale_as, kTRUE);
    fFinalLRSub[i] = BaselineSubtraction(fDeltaPhiLRSub[i], i, kTRUE, zyam, err_zyam, scale_ns);
    ns_yield = YieldNS(fFinalLRSub[i], Dphi_ns_range, kFALSE);
    err_ns = YieldNS(fFinalLRSub[i], Dphi_ns_range, kTRUE);
    as_yield = YieldAS(fFinalLRSub[i], Dphi_ns_range, kFALSE, kTRUE, scale_ns, scale_as);
    err_as = YieldAS(fFinalLRSub[i], Dphi_ns_range, kTRUE, kTRUE, scale_ns, scale_as);

    NearSide->SetPoint(points, i*5+2.5, ns_yield);
    NearSide->SetPointError(points, 0.5, err_ns);
    GraphProperties(NearSide, "Near-side yield", kGreen+2, 21);

    AwaySide->SetPoint(points, i*5+2.5, as_yield);
    AwaySide->SetPointError(points, 0.5, err_as);
    GraphProperties(AwaySide, "Away-side yield", kBlue, 21);

    // no Long Range subtraction
    fLongRangeForZyamNotSub[i] = ZyamLRsubtraction(correlationsME[i], i, scale_as, kFALSE);
    zyamNoSub = GetZyam(fLongRangeForZyamNotSub[i], kFALSE, Dphi_zyam_notsub);
    err_zyamNoSub = GetZyamError(fLongRangeForZyamNotSub[i], kFALSE, Dphi_zyam_notsub);
    fDeltaPhiLRNotSub[i] = LRsubtraction(correlationsME[i], i, scale_ns, scale_as, kFALSE);
    fFinalLRNotSub[i] = BaselineSubtraction(fDeltaPhiLRNotSub[i], i, kFALSE, zyamNoSub, err_zyamNoSub, scale_ns);
    nsNoSub_yield = YieldNS(fFinalLRNotSub[i], Dphi_ns_range, kFALSE);
    err_nsNoSub = YieldNS(fFinalLRNotSub[i], Dphi_ns_range, kTRUE);
    asNoSub_yield = YieldAS(fFinalLRNotSub[i], Dphi_ns_range, kFALSE, kFALSE, scale_ns, scale_as);
    err_asNoSub = YieldAS(fFinalLRNotSub[i], Dphi_ns_range, kTRUE, kFALSE, scale_ns, scale_as);

    NearSideNotSub->SetPoint(points, i*5+2.5, nsNoSub_yield);
    NearSideNotSub->SetPointError(points, 0.5, err_nsNoSub);
    GraphProperties(NearSideNotSub, "Near-side yield without LR subtraction", kGreen, 21);

    AwaySideNotSub->SetPoint(points, i*5+2.5, asNoSub_yield);
    AwaySideNotSub->SetPointError(points, 0.5, err_asNoSub);
    GraphProperties(AwaySideNotSub, "Away-side yield without LR subtraction", kCyan+1, 21);
    
    //<triggers> and <uncorrelated_seeds>
    AverageTriggers->SetPoint(points, i*5+2.5, p0->GetBinContent(i+1));
    AverageTriggers->SetPointError(points, 0.5, p0->GetBinError(i+1));
    GraphProperties(AverageTriggers, "Average triggers", kBlack, 21);

    UncorrelatedSeeds->SetPoint(points, i*5+2.5, p0->GetBinContent(i+1)/(1+ ns_yield + as_yield));
    UncorrelatedSeeds->SetPointError(points, 0.5, TMath::Sqrt( pow(p0->GetBinError(i+1),2)/pow(1+ns_yield+as_yield,2) + pow(p0->GetBinContent(i+1),2)/pow((1+ns_yield+as_yield),4) * (pow(err_ns,2)+pow(err_as,2)))); //sqrt{(DA)^2/B^2 + A^2/B^4 * [(DB1)^2 + (DB2)^2]}
    GraphProperties(UncorrelatedSeeds, "Uncorrelated seeds", kRed, 21);

    points++;
  }


  WriteIntoFile(output, NearSide);
  WriteIntoFile(output, AwaySide);		
  WriteIntoFile(output, NearSideNotSub);
  WriteIntoFile(output, AwaySideNotSub);
  WriteIntoFile(output, AverageTriggers);
  WriteIntoFile(output, UncorrelatedSeeds);
  output->Close();
}


TH2D *MERemoval(TH2D* correlations, Int_t loopindex, TF1* mixedEvent){

  TH2D *correlationsME = (TH2D*)correlations->Clone(Form("%i",loopindex));
  correlationsME->Reset();
  correlationsME->Sumw2();

  TH1D *fEtaProjection = (TH1D*)correlations->ProjectionY(Form("_py_%i",loopindex));

  for(Int_t ix = 0; ix < correlations->GetNbinsX(); ix++){
    for(Int_t iy = 0; iy < correlations->GetNbinsY(); iy++){
      correlationsME->SetBinContent(ix+1, iy+1, correlations->GetBinContent(ix+1,iy+1) * mixedEvent->Eval(fEtaProjection->GetBinCenter(iy+1)));
      correlationsME->SetBinError(ix+1, iy+1, correlations->GetBinError(ix+1,iy+1) * mixedEvent->Eval(fEtaProjection->GetBinCenter(iy+1)));
    }
  }

  return correlationsME;

}


TF1 *fTriangle(Double_t *x, Double_t *par){
    
  if(x[0]>-Deta_max && x[0]<=0){
    return par[0]+par[1]*x[0];
  }
  else if(x[0]>0 && x[0]<Deta_max){
    return par[2]+par[3]*x[0];
  }
  else
    return 0;
}

TH1D *ZyamLRsubtraction(TH2D *correlationsME, Int_t loopindex, Double_t scale_as, Bool_t sub){

  //no eta gap
  TH1D *fDeltaPhi = (TH1D*)correlationsME->ProjectionX(Form("_px_%i",loopindex),1,correlationsME->GetNbinsY(),"e");
  TH1D *fDeltaEta = (TH1D*)correlationsME->ProjectionY(Form("_py_%i",loopindex),1,correlationsME->GetNbinsX(),"e");

  //eta gap for long range
  TH1D *fLongRangeNeg = (TH1D*)correlationsME->ProjectionX(Form("longrangeneg_px_%i",loopindex),fDeltaEta->FindBin(-Deta_max+0.0001),fDeltaEta->FindBin(-Deta_gap-0.0001),"e");

  TH1D *fLongRangePos = (TH1D*)correlationsME->ProjectionX(Form("longrangepos_px_%i",loopindex),fDeltaEta->FindBin(Deta_gap+0.0001),fDeltaEta->FindBin(Deta_max-0.0001),"e");

  TH1D *fLongRange = (TH1D*)fLongRangeNeg->Clone(Form("longrange_px_%i",loopindex));
  fLongRange->Reset();
  fLongRange->Sumw2();
  fLongRange->Add(fLongRangeNeg,fLongRangePos);

  //scale by _as factor (in the subtraction A - B where A has no eta gap at all --fDeltaPhi)

  Double_t k = 1;
  if(sub) k = scale_as;

  for(Int_t i = 0; i < fLongRange->GetNbinsX(); i++){
    fLongRange->SetBinContent(i+1, fLongRange->GetBinContent(i+1) * k); 
    fLongRange->SetBinError(i+1, fLongRange->GetBinError(i+1) * k);
  }

  TH1D *fTotMinusLong = (TH1D*)fDeltaPhi->Clone(Form("dphi-long_%i",loopindex));
  fTotMinusLong->Reset();
  fTotMinusLong->Sumw2();
  fTotMinusLong->Add(fDeltaPhi,fLongRange,1,-1);

  if(sub) return fTotMinusLong;
  else    return fLongRange;
}

TH1D *LRsubtraction(TH2D *correlationsME, Int_t loopindex, Double_t scale_ns, Double_t scale_as, Bool_t sub){

  //no eta gap
  TH1D *fDeltaPhi = (TH1D*)correlationsME->ProjectionX(Form("_px_phi_%i",loopindex),1,correlationsME->GetNbinsY(),"e");
  TH1D *fDeltaEta = (TH1D*)correlationsME->ProjectionY(Form("_py_eta_%i",loopindex),1,correlationsME->GetNbinsX(),"e");

  TH1D *fSignalNS = (TH1D*)correlationsME->ProjectionX(Form("_px_ns_%i",loopindex),fDeltaEta->FindBin(-Deta_gap+0.0001),fDeltaEta->FindBin(Deta_gap-0.0001),"e");
  TH1D *fSignalAS = (TH1D*)correlationsME->ProjectionX(Form("_px_as_%i",loopindex),fDeltaEta->FindBin(-Deta_max+0.0001),fDeltaEta->FindBin(Deta_max-0.0001),"e");

  TH1D *fSignal = (TH1D*)fDeltaPhi->Clone(Form("_px_ns+as_%i",loopindex));
  fSignal->Reset();
  fSignal->Sumw2();
  
  Double_t k = 1;
  if(!sub) k = scale_ns/scale_as;

  for(Int_t i = 0; i < fSignal->GetNbinsX()/2; i++){
    fSignal->SetBinContent(i+1, fSignalNS->GetBinContent(i+1));
    fSignal->SetBinError(i+1, fSignalNS->GetBinError(i+1));
  }
  for(Int_t i = fSignal->GetNbinsX()/2; i < fSignal->GetNbinsX(); i++){
    fSignal->SetBinContent(i+1, fSignalAS->GetBinContent(i+1)*k);
    fSignal->SetBinError(i+1, fSignalAS->GetBinError(i+1)*k);
  }
  
  //eta gap for long range
  TH1D *fLongRangeNeg = (TH1D*)correlationsME->ProjectionX(Form("longrangeneg_px_%i",loopindex),fDeltaEta->FindBin(-Deta_max+0.0001),fDeltaEta->FindBin(-Deta_gap-0.0001),"e");

  TH1D *fLongRangePos = (TH1D*)correlationsME->ProjectionX(Form("longrangepos_px_%i",loopindex),fDeltaEta->FindBin(Deta_gap+0.0001),fDeltaEta->FindBin(Deta_max-0.0001),"e");

  TH1D *fLongRange = (TH1D*)fLongRangeNeg->Clone(Form("longrange_px_%i",loopindex));
  fLongRange->Reset();
  fLongRange->Sumw2();
  fLongRange->Add(fLongRangeNeg,fLongRangePos);

  //mirror ns into as
  
  for(Int_t i = fLongRange->GetNbinsX()/2+1;i <= fLongRange->GetNbinsX(); i++){
    fLongRange->SetBinContent(i, fLongRange->GetBinContent(fLongRange->GetNbinsX()+1-i));
    fLongRange->SetBinError(i, fLongRange->GetBinError(fLongRange->GetNbinsX()+1-i));
  }
  
  //scale by _ns and _as factor 

  for(Int_t i = 0; i < fLongRange->GetNbinsX()/2; i++){
    fLongRange->SetBinContent(i+1, fLongRange->GetBinContent(i+1) * scale_ns); 
    fLongRange->SetBinError(i+1, fLongRange->GetBinError(i+1) * scale_ns);
  }
  for(Int_t i = fLongRange->GetNbinsX()/2; i < fLongRange->GetNbinsX(); i++){
    fLongRange->SetBinContent(i+1, fLongRange->GetBinContent(i+1) * scale_as); 
    fLongRange->SetBinError(i+1, fLongRange->GetBinError(i+1) * scale_as);
  }

  TH1D *fFinal = (TH1D*)fDeltaPhi->Clone(Form("_px_phi_sub_%i",loopindex));
  fFinal->Reset();
  fFinal->Sumw2();

  fFinal->Add(fSignal,fLongRange,1,-1);
  
  if(sub) return fFinal;
  else    return fSignal;
}

Double_t GetZyam(TH1D *fHisto, Bool_t sub, const Double_t Dphi_zyam_notsub){

   if(sub){
     TF1 *fit = new TF1("fit","pol0(0)", fHisto->GetBinCenter(fHisto->GetNbinsX()/2), fHisto->GetBinCenter(fHisto->GetNbinsX()));
     fHisto->Fit("fit","R0");
     Double_t zyam = fit->GetParameter(0);
   }
   else{
     Double_t zyam = ( fHisto->GetBinContent(fHisto->FindBin(-Dphi_zyam_notsub) -1) + fHisto->GetBinContent(fHisto->FindBin(-Dphi_zyam_notsub)) + fHisto->GetBinContent(fHisto->FindBin(Dphi_zyam_notsub)) + fHisto->GetBinContent(fHisto->FindBin(Dphi_zyam_notsub) +1) ) / 4;
   }
   
   return zyam;
}

Double_t GetZyamError(TH1D *fHisto, Bool_t sub, const Double_t Dphi_zyam_notsub){

   if(sub){
     TF1 *fit = new TF1("fit","pol0(0)", fHisto->GetBinCenter(fHisto->GetNbinsX()/2), fHisto->GetBinCenter(fHisto->GetNbinsX()));
     fHisto->Fit("fit","R0");
     Double_t zyam_err = fit->GetParError(0);
   }
   else{
     Double_t zyam_err = TMath::Sqrt( pow(fHisto->GetBinError(fHisto->FindBin(-Dphi_zyam_notsub) -1),2) + pow(fHisto->GetBinError(fHisto->FindBin(-Dphi_zyam_notsub)),2) + pow(fHisto->GetBinError(fHisto->FindBin(Dphi_zyam_notsub)),2) + pow( fHisto->GetBinError(fHisto->FindBin(Dphi_zyam_notsub) +1),2) ) / 4; 
   }
   
   return zyam_err;		      
}

TH1D *BaselineSubtraction(TH1D *fHisto, Int_t loopindex, Bool_t sub, Double_t zyam, Double_t zyam_err, Double_t scale_ns){

  TH1D *fHistoSub = (TH1D*)fHisto->Clone(Form("_baseSub_%i",loopindex));
  fHistoSub->Reset();
  fHistoSub->Sumw2();

  for(Int_t i = 0; i < fHistoSub->GetNbinsX(); i++){
    if(sub){
      fHistoSub->SetBinContent(i+1, fHisto->GetBinContent(i+1) - zyam);
      fHistoSub->SetBinError(i+1, TMath::Sqrt(pow(fHisto->GetBinError(i+1),2) + pow(zyam_err,2)));
    } 
    else{
      fHistoSub->SetBinContent(i+1, fHisto->GetBinContent(i+1) - zyam*scale_ns);
      fHistoSub->SetBinError(i+1, TMath::Sqrt(pow(fHisto->GetBinError(i+1),2) + pow(zyam_err,2)*pow(scale_ns,2)));
    }
  }

  return fHistoSub;
}

Double_t YieldNS(TH1D *fHisto, const Double_t Dphi_ns_range, Bool_t err){

  Double_t error = 0;
  Double_t ns_result = fHisto->IntegralAndError(fHisto->FindBin(-Dphi_ns_range), fHisto->FindBin(Dphi_ns_range), error, "width");

  if(err) return error;
  else return ns_result;

}

Double_t YieldAS(TH1D *fHisto, const Double_t Dphi_ns_range, Bool_t err, Bool_t sub, Double_t scale_ns, Double_t scale_as){

  Double_t error1 = 0, error2 = 0, error = 0;
  Double_t as_result = fHisto->IntegralAndError(1, fHisto->FindBin(-Dphi_ns_range)-1, error1, "width") + fHisto->IntegralAndError(fHisto->FindBin(Dphi_ns_range)+1, fHisto->GetNbinsX(), error2, "width");
  error = TMath::Sqrt(pow(error1,2) + pow(error2,2));

  if(!sub){
    as_result *= scale_as;
    as_result /= scale_ns;
    error *= scale_as;
    error /= scale_ns;
  }

  if(err) return error;
  else return as_result;

}

void GraphProperties(TGraphErrors *graph, const char *name, Int_t color, Int_t style){

  graph->SetName(name);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(style);
  graph->SetLineWidth(1);

}

void WriteIntoFile(TFile *file, TGraphErrors *graph){

  graph->Write();
}
