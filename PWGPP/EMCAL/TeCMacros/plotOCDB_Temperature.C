#include <TCanvas.h>
#include <TDatime.h>
#include <TFile.h>
#include <TFile.h>
#include <TGrid.h>
#include <TH2F.h>
#include <TMap.h>
#include <TNtuple.h>
#include <stdio.h>
#include "TInfo.C"

TH2 *getT2D(TObjArray *arr, Int_t type=1)
{
  if (!arr) 
    return 0;
  const Int_t rns=arr->GetEntries();
  TH2F *ret = new TH2F("",Form(";run idx;sensor idx"),rns,0,rns,159,0,159);
  ret->SetDirectory(0);
  for (Int_t i=0;i<rns;++i) {
    TInfo *tinfo = dynamic_cast<TInfo*>(arr->At(i));
    if (!tinfo)
      continue;
    for (Int_t j=0;j<159;++j) {
      if (!tinfo->IsValid(j)) 
	continue;
      Double_t val = 0;
      if (type==1) 
	val = tinfo->MinT(j);
      else if (type==2) 
	val = tinfo->MaxT(j);
      else 
	val = tinfo->Diff(j);
      ret->SetBinContent(ret->FindBin(i,j),val);
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

void plotOCDB_Temperature()
{
  plot_LHC18d();
}

