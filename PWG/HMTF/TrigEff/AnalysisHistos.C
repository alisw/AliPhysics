#if !defined (__CINT__) || (defined(__MAKECINT__))
#include <iostream>
#include "TFile.h"
#include "AliAnalysisTaskTrigEff.h"
#include "TCanvas.h"
#include "TString.h"
#include "TROOT.h"
#include "TLatex.h"
#endif

TH1 *  PlotEff(TH1 * hNum, TH1* hDenum, TCanvas * c = 0, Int_t npad = 0);
Double_t IntegrateTracklets(TH1 * h) ;


// This is used to rebin the V0M histograms so that they match the multiplicity classes used in the analysis.
Int_t nbinsNew = 11;
const Double_t bins[]= {0, 1, 5, 10, 15, 20, 30, 40, 50, 70, 100, 110};


void AnalysisHistos(const char * filename = "out_170413.root", Bool_t rebin = kTRUE) {

  TFile *fileIn = new TFile(filename);
  TH1F * hV0M[kNHist];
  TH1F *hTrk70100[kNHist];
  TList * l = (TList*) fileIn->Get("taskout");
  const char * histoNames[] = {"kHistoINT11", "kHistoINT7", "kHistoINT7Offline", "kHistoINT11EvSel", "kHistoINT7EvSel", "kHistoINT7OfflineEvSel"};
  for(Int_t ihist = 0; ihist < kNHist; ihist++){
    hV0M[ihist] = (TH1F*) l->FindObject(histoNames[ihist]);
    hTrk70100[ihist] = (TH1F*) l->FindObject(TString(histoNames[ihist])+"_ntrk70100");
  }
  TH1I * hEv=(TH1I*) l->FindObject("fHistEvCount");
  Int_t nev = hEv->GetBinContent(hEv->FindBin(kEvInelGT0)); // Normalize to histos passing the INEL > 0 Selection

  new TCanvas("cEvSel", "cEvSel");
  hEv->Draw();

  TCanvas * c = new TCanvas("cV0M","cV0M", 900,600);
  c->Divide(3,2);

  for(Int_t ihist = 0; ihist < kNHist; ihist++){
    c->cd(ihist+1);
    if(rebin)
      hV0M[ihist] = (TH1F*) hV0M[ihist]->Rebin(nbinsNew, TString(hV0M[ihist]->GetName())+"_rebin", bins);
    hV0M[ihist]->Scale(1./nev, "width");
    hV0M[ihist]->Draw("E");    
  }
  
  TCanvas * cEff = new TCanvas("cEff","cEff");
  cEff->Divide(2,2);

  PlotEff(hV0M[kHistoINT7]             , hV0M[kHistoINT11]      , cEff , 1);
  PlotEff(hV0M[kHistoINT7Offline]      , hV0M[kHistoINT11]      , cEff , 2);
  PlotEff(hV0M[kHistoINT7EvSel]        , hV0M[kHistoINT11EvSel] , cEff , 3);
  PlotEff(hV0M[kHistoINT7OfflineEvSel] , hV0M[kHistoINT11EvSel] , cEff , 4);    

  TCanvas * cTrk70100 = new TCanvas("cTrk70100","cTrk70100", 900,450);
  cTrk70100->Divide(2,1);
  cTrk70100->cd(1);
  hTrk70100[kHistoINT11EvSel]       ->Draw();
  hTrk70100[kHistoINT7OfflineEvSel] ->Draw("same");
  gROOT->ProcessLine("NewLegend(\"\", \"l\", 0,1)" );
  PlotEff(hTrk70100[kHistoINT7OfflineEvSel], hTrk70100[kHistoINT11EvSel], cTrk70100,2);
  Double_t ntrkAcc = IntegrateTracklets(hTrk70100[kHistoINT7OfflineEvSel]);
  Double_t ntrkAll = IntegrateTracklets(hTrk70100[kHistoINT11EvSel]);
  std::cout << "Tracklet Trigger Efficiency "<< ntrkAcc/ntrkAll << std::endl;
  TLatex * ltex = new TLatex(0.35,0.2,Form("Tracklet Trig. Eff. = %2.2f", ntrkAcc/ntrkAll));
  ltex->SetNDC();
  ltex->Draw();

  
}

TH1 *  PlotEff(TH1 * hNum, TH1* hDenum, TCanvas * c, Int_t npad) {
  if(c) c->cd(npad);
  TH1 * hEff = (TH1*) hNum->Clone(TString("hNun_")+"Eff");
  hEff->Divide(hNum, hDenum, 1, 1, "B");
  hEff->Draw();
  return hEff;
}

Double_t IntegrateTracklets(TH1 * h) {

  Int_t nbin = h->GetNbinsX();
  Double_t integral = 0;
  for(Int_t ibin = 1; ibin <= nbin; ibin++){
    integral += h->GetBinContent(ibin) * h->GetBinCenter(ibin);
    
  }
  return integral;
  
}
