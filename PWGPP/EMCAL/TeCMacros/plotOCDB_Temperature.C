#include <TCanvas.h>
#include <TDatime.h>
#include <TFile.h>
#include <TFile.h>
#include <TGrid.h>
#include <TH2F.h>
#include <TMap.h>
#include <TNtuple.h>
#include <TProfile.h>
#include <TLegend.h>
#include <stdio.h>
#include "TInfo.h"

TH2 *getT2D(TObjArray *arr, Int_t type=1)
{
  if (!arr) 
    return 0;
  const Int_t rns=arr->GetEntries();
  TInfo f; 
  TString title=f.Type(type);
  TH2F *ret = new TH2F(Form("h2f_%d",type),Form("%s;run idx;sensor idx",title.Data()),rns,0,rns,159,0,159);
  ret->SetDirectory(0);
  for (Int_t i=0;i<rns;++i) {
    TInfo *tinfo = dynamic_cast<TInfo*>(arr->At(i));
    if (!tinfo)
      continue;
    for (Int_t j=0;j<160;++j) {
      if (!tinfo->IsValid(j)) 
	continue;
      Double_t val = tinfo->T(j,type);
      ret->SetBinContent(ret->FindBin(i,j),val);
    }
  }
  return ret;
}

TH1 *getT(Int_t ns1, Int_t ns2, TObjArray *arr, Int_t type=1) 
{
  if (!arr) 
    return 0;
  const Int_t rns=arr->GetEntries();
  TProfile *ret = new TProfile(Form("h%d%d%d",ns1,ns2,type),Form(";run idx"),rns,0,rns);
  ret->SetDirectory(0);
  for (Int_t i=0;i<rns;++i) {
    TInfo *tinfo = dynamic_cast<TInfo*>(arr->At(i));
    if (!tinfo)
      continue;
    for (Int_t j=ns1;j<=ns2;++j) {
      if (!tinfo->IsValid(j)) 
	continue;
      Double_t val = 0;
      if (type==1) 
	val = tinfo->MinT(j);
      else if (type==2) 
	val = tinfo->MaxT(j);
      else 
	val = tinfo->Diff(j);
      ret->Fill(i,val);
    }
  }
  return ret;
}

void plot_LHC18d() 
{
  TFile *inf = TFile::Open("temperatures.root");
  TObjArray *arr = dynamic_cast<TObjArray*>(inf->Get("temperatures_lhc18d"));
  if (!arr) 
    return;
  //arr->Print();

  if (0) {
    TCanvas *c1 = new TCanvas("cMinT18d","cMinT18d");
    TH2 *hminT = getT2D(arr,1);
    hminT->SetStats(0);
    hminT->Draw("colz");
    c1->Print(Form("%s.pdf",c1->GetName()));
    TCanvas *c2 = new TCanvas("cMaxT18d","cMaxT18d");
    TH2 *hmaxT = getT2D(arr,2);
    hmaxT->SetStats(0);
    hmaxT->Draw("colz");
    c2->Print(Form("%s.pdf",c2->GetName()));
    TCanvas *c3 = new TCanvas("cDiffT18d","cDiffT18d");
    TH2 *hdiffT = getT2D(arr,0);
    hdiffT->SetStats(0);
    hdiffT->Draw("colz");
    c3->Print(Form("%s.pdf",c3->GetName()));
  }
  if (0) {
    TCanvas *c1 = new TCanvas("cTemp","cTemp");
    Int_t type=1;
    TH2 *h2f = new TH2F("h2f",";run idx;T",1,0,60,1,18,26);
    h2f->SetStats(0);
    h2f->Draw();
    TLegend *leg = new TLegend(0.5,0.5,0.9,0.9);
    for (Int_t i=0,n1=0,n2=7;i<20;++i) {
      TH1 *h=getT(n1,n2,arr,type);
      n1+=8;
      n2+=8;
      h->SetLineColor(i+1);
      h->SetLineWidth(3);
      h->SetName(Form("SM%d",i));
      h->Draw("same, hist");
      leg->AddEntry(h,h->GetName(),"l");
    }
    leg->Draw();
  }
  if (1) {
    TCanvas *ci = new TCanvas("cLHC18d","cLHC18d",1200,800);
    ci->Print(Form("%s.pdf[",ci->GetName()));
    const Int_t rns=arr->GetEntries();
    for (Int_t i=0;i<rns;++i) {
      TInfo *tinfo = dynamic_cast<TInfo*>(arr->At(i));
      if (!tinfo)
	continue;
      TH2 *h=tinfo->GetHist(3);
      h->Draw("colz text");
      ci->Print(Form("%s.pdf",ci->GetName()));
    }
    ci->Print(Form("%s.pdf]",ci->GetName()));
  }
}

void plotOCDB_Temperature()
{
  plot_LHC18d();
}

