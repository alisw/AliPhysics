#if !defined(__CINT__) || defined(__MAKECINT__)
#include <AliEMCALGeometry.h>
#include <TCanvas.h>
#include <TDatime.h>
#include <TFile.h>
#include <TGrid.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TMap.h>
#include <TNtuple.h>
#include <TProfile.h>
#include "LInfo.h"

class LDraw : public TNamed {
 public:
  LDraw(const char *name, const char *fname="tempinfo.root"); 
  virtual ~LDraw() {;}
  void       Compute();
  Int_t      GetNRuns()                    const { return fArr->GetEntries(); }
  void       Print(Option_t *opt="")       const {};
  void       SetPrint(Bool_t b=1)                { fDoPrint=b; }
  void       DrawAll();
  TCanvas   *DrawFrac(Int_t type)          const;
  TH1       *GetFracRun(Int_t sm, Int_t t) const;
  Int_t      GetRunNo(Int_t run)           const { return (static_cast<LInfo*>(fArr->At(run)))->GetRunNo(); }
  TObjArray *fArr;     // array with info
  Bool_t     fDoPrint; // if true then print canvases  
  ClassDef(LDraw, 1); // Led draw class
};
#endif

LDraw::LDraw(const char *name, const char *fname) : TNamed(name,fname), fArr(0), fDoPrint(0)
{
  TFile *inf = TFile::Open("ledinfo.root");
  fArr = dynamic_cast<TObjArray*>(inf->Get(Form("led_%s",name)));
  delete inf;
}

void LDraw::Compute()
{
  const Int_t rns=fArr->GetEntries();
  for (Int_t i=0;i<rns;++i) {
    LInfo *linfo = dynamic_cast<LInfo*>(fArr->At(i));
    if (!linfo)
      continue;
    linfo->Compute();
    if (fDoPrint) {
      linfo->Print();
      cout << "fraction good strips ";
      for (Int_t i=0;i<20;++i) 
	cout << linfo->FracStrips(i) << " ";
      cout << endl;
      cout << "fraction good towers ";
      for (Int_t i=0;i<20;++i) 
	cout << linfo->FracLeds(i) << " ";
      cout << endl;
    }
  }
}

void LDraw::DrawAll()
{
  TString cname(Form("cLedAll_%s",GetName()));
  TCanvas *c1 = DrawFrac(0);
  c1->Print(Form("%s.pdf[",cname.Data()));
  c1->Print(Form("%s.pdf",cname.Data()));
  TCanvas *c2=DrawFrac(1);
  c2->Print(Form("%s.pdf",cname.Data()));
  TCanvas *c = 0;
  const Int_t rns=fArr->GetEntries();
  for (Int_t i=0;i<rns;++i) {
    LInfo *linfo = dynamic_cast<LInfo*>(fArr->At(i));
    if (!linfo)
      continue;
    c = linfo->DrawHist(5,1);
    c->Print(Form("%s.pdf",cname.Data()));
    delete c;
    c = linfo->DrawHist(6,1);
    c->Print(Form("%s.pdf",cname.Data()));
    delete c;
    c = linfo->DrawHist(0,1);
    c->Print(Form("%s.pdf",cname.Data()));
    delete c;
  }
  c2->Print(Form("%s.pdf]",cname.Data()));
}

TCanvas *LDraw::DrawFrac(Int_t type) const
{
  const char *tname="Led";
  if (type==0) 
    tname = "LedMon";
  TString lab(Form("cLedFrac_%s_%s",tname,GetName()));
  TCanvas *c = new TCanvas(lab,lab,1200,800);
  TLegend *leg = new TLegend(0.92,0.1,0.99,0.99);
  TObjArray arr;
  Double_t min=1e10,max=-1e10;
  for (Int_t i=0;i<20;++i) {
    TH1 *h=GetFracRun(i,type);
    Int_t col=i+1;
    switch (i) {
    case  9: col=kOrange+2; break;
    case 10: col=kOrange+10; break;
    case 11: col=kMagenta+2; break;
    case 12: col=kCyan+2; break;
    case 13: col=kYellow+2; break;
    case 14: col=kGray+2; break;
    case 15: col=kOrange-2; break;
    case 16: col=kViolet+2; break;
    case 17: col=kRed-2; break;
    case 18: col=kGreen+2; break;
    case 19: col=kGray+2; break;
    }
    h->SetMarkerStyle(20+(i%2));
    h->SetMarkerSize(1.2);
    h->SetMarkerColor(col);
    h->SetLineColor(col);
    h->SetLineWidth(3);
    h->SetName(Form("SM%d",i));
    arr.Add(h);
    leg->AddEntry(h,h->GetName(),"p");
    min=TMath::Min(min,h->GetMinimum());
    max=TMath::Max(max,h->GetMaximum());
  }
  TH2 *h2f = new TH2F("h2f",";run idx;T",1,0,GetNRuns(),1,-0.05,1.05);
  h2f->SetTitle(Form("Fraction %s per SM: run %d to %d, min=%.1f, max=%.1f",tname, GetRunNo(0), GetRunNo(GetNRuns()-1), min, max));
  h2f->SetStats(0);
  h2f->Draw();
  leg->Draw();
  for (Int_t i=0;i<arr.GetEntries();++i)
    arr.At(i)->Draw("same, p, hist");
  c->SetGridx(1);
  c->SetGridy(1);
  if (fDoPrint)
    c->Print(Form("%s.pdf",c->GetName()));
  return c;
}

TH1 *LDraw::GetFracRun(Int_t sm, Int_t t) const
{
  const Int_t rns=fArr->GetEntries();
  TH1F *ret = new TH1F(Form("hFracRun_%d_%d",sm,t),Form(";run idx;fraction"),rns,0,rns);
  if (t==0)
    ret->SetTitle(Form("Fraction of strips per run: run %d to %d",GetRunNo(0), GetRunNo(GetNRuns()-1)));
  else
    ret->SetTitle(Form("Fraction of leds per run: run %d to %d",GetRunNo(0), GetRunNo(GetNRuns()-1)));
  ret->SetDirectory(0);
  ret->SetStats(0);
  ret->SetMarkerStyle(20);
  ret->SetMarkerSize(1);
  Int_t r1,r2=0;
  for (Int_t i=0;i<rns;++i) {
    LInfo *linfo = dynamic_cast<LInfo*>(fArr->At(i));
    if (!linfo)
      continue;
    Double_t val = 0;
    if (t==0) 
      val = linfo->FracStrips(sm);
    else
      val = linfo->FracLeds(sm);
    ret->SetBinContent(i+1,val);
  }
  return ret;
}

void plotL_period(const char *period, Bool_t doprint=0) 
{
  LDraw d(period);
  d.SetPrint(doprint);
  d.Compute();
  d.Print();
  if (0) {
    d.DrawFrac(0);
    d.DrawFrac(1);
  } else 
    d.DrawAll();
}

void plotOCDB_LED(const char *period="lhc18d")
{
  plotL_period(period);
}

#if 0
void test_geo()
{
  AliEMCALGeometry *g=AliEMCALGeometry::GetInstance("EMCAL_COMPLETE12SMV1_DCAL_8SM");
  Int_t kSM=g->GetNumberOfSuperModules();
  cout << "Number of SM: " << kSM << endl;
  for (Int_t i=0;i<kSM;++i) {
    Int_t np = g->GetNumberOfCellsInPhiDirection(i);
    Int_t ne = g->GetNumberOfCellsInEtaDirection(i);
    cout << i << ": " << np << " " << ne << endl;
  }
}
#endif


