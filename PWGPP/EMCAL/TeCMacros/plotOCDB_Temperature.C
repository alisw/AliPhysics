#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TCanvas.h>
#include <TDatime.h>
#include <TFile.h>
#include <TGrid.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMap.h>
#include <TNtuple.h>
#include <TObjArray.h>
#include <TProfile.h>
#include "TInfo.h"

class TDraw : public TNamed {
 public:
  TDraw(const char *name, const char *fname="tempinfo.root");
  virtual ~TDraw() {;}
  TObjArray *GetArray()                       { return fArr; }
  void       Compute();
  void       DrawAll();
  TCanvas   *DrawTSensorsPerSM(Int_t type, Int_t sm) const;
  TCanvas   *DrawTPerSM(Int_t type=3)        const;
  TCanvas   *DrawT2D(Int_t type=3)      const;
  TCanvas   *DrawOccRun()               const;
  TCanvas   *DrawOccSensor2D()          const;
  TCanvas   *DrawOcc2D()                const;
  Double_t   GetMaxT(Int_t type=3)      const;
  Double_t   GetMinT(Int_t type=3)      const;
  Double_t   GetMaxTperS(Int_t ns, Int_t type=3) const;
  Double_t   GetMinTperS(Int_t ns, Int_t type=3) const;
  Double_t   GetFraction(Int_t ns)      const { return (Double_t)GetNRunsValid(ns)/GetNRuns(); }
  Int_t      GetBad(Int_t ns)           const { return fBad.At(ns); }
  Int_t      GetNBad()                  const;
  Int_t      GetNRuns()                 const { return fArr->GetEntries(); }
  Int_t      GetNRunsValid(Int_t ns)    const;
  Int_t      GetRunNo(Int_t run)        const { return (static_cast<TInfo*>(fArr->At(run)))->GetRunNo(); }
  TH1       *GetOccRun()                const;
  TH1       *GetOccSensor()             const;
  TH2       *GetOccSensor2D()           const;
  TH2       *GetOcc2D()                 const;
  TH1       *GetT(Int_t nsensor, Int_t type=3) const;
  TH1       *GetTSM(Int_t sm, Int_t type=3) const;
  TH2       *GetT2D(Int_t type=3)       const;
  Bool_t     IsGood(Int_t ns)           const { return (fBad.At(ns)==0); }
  void       Print(Option_t *opt="")    const;
  void       SetPrint(Bool_t b=1)        { fDoPrint=b; }
  void       SetBad(Int_t ns, Bool_t b)  { fBad.SetAt(ns,b); }
 protected:
  TH2       *GetMask()                  const;
  TObjArray *fArr;     // array with info
  Bool_t     fDoPrint; // if true then print canvases
  TArrayI    fBad;     // list of bad channels
  Double_t   fMinFrac; // minimum fraction
  ClassDef(TDraw, 1); // Temperature draw class
};
#endif

TDraw::TDraw(const char *name, const char *fname) : TNamed(name,fname), fArr(0), fDoPrint(0), fBad(TInfo::NSensors()), fMinFrac(0.8)
{
  TFile *inf = TFile::Open("tempinfo.root");
  fArr = dynamic_cast<TObjArray*>(inf->Get(Form("temperatures_%s",name)));
  delete inf;
  fBad.SetAt(1,26);
  fBad.SetAt(1,40);
  fBad.SetAt(1,43);
  fBad.SetAt(1,44);
  fBad.SetAt(1,45);
  fBad.SetAt(1,84);
  fBad.SetAt(1,85);
  fBad.SetAt(1,86);
  fBad.SetAt(1,87);
  fBad.SetAt(1,92);
  fBad.SetAt(1,93);
  fBad.SetAt(1,94);
  fBad.SetAt(1,95);
  fBad.SetAt(1,117);
  fBad.SetAt(1,135);
  fBad.SetAt(1,140);
  fBad.SetAt(1,147);
  fBad.SetAt(1,148);
  fBad.SetAt(1,149);
  fBad.SetAt(1,150);
  fBad.SetAt(1,151);
  fBad.SetAt(1,156);
  fBad.SetAt(1,157);
  fBad.SetAt(1,158);
  fBad.SetAt(1,159);
}

void TDraw::Compute()
{
  const Int_t nbad = GetNBad();
  cout << "Known bad sensors: " << nbad << " (";
  for (Int_t ns=0,first=0;ns<TInfo::NSensors();++ns) {
    if (IsGood(ns))
      continue;
    if (first==0) {
      cout << ns;
      first = 1;
    } else
      cout << ", " << ns;
  }
  cout << ")" << endl;

  const Int_t rns=fArr->GetEntries();
  for (Int_t i=0;i<rns;++i) {
    TInfo *tinfo = dynamic_cast<TInfo*>(fArr->At(i));
    if (!tinfo)
      continue;
    if (tinfo->GetNFaulty()!=nbad) {
      cout << "Setting nfaulty to " << nbad << " for run " << i << endl;
      tinfo->SetNFaulty(nbad);
    }
    Double_t frac = tinfo->Fraction();
    if (frac<fMinFrac)
      cout << "Run " << i << " with run number " << tinfo->GetRunNo() << " below threshold: " << frac << endl;
  }

  for (Int_t ns=0;ns<TInfo::NSensors();++ns) {
    Double_t frac=GetFraction(ns);
    if (!IsGood(ns)) {
      if (frac>0) {
        cout << "Sensor " << ns << " marked bad, but has fraction: " << frac << endl;
      }
      continue;
    }
    if (frac<fMinFrac)
      cout << "Sensor " << ns << " below threshold: " << frac << endl;
  }
}

void TDraw::DrawAll()
{
  TString cname(Form("cTempAll_%s",GetName()));
  TCanvas *c1 = DrawOccRun();
  c1->Print(Form("%s.pdf[",cname.Data()));
  c1->Print(Form("%s.pdf",cname.Data()));
  TCanvas *c2=DrawOccSensor2D();
  c2->Print(Form("%s.pdf",cname.Data()));
  TCanvas *c3=DrawTPerSM(3);
  c3->Print(Form("%s.pdf",cname.Data()));
  for (Int_t i=0; i < 20; i++){
    TCanvas* c5 = DrawTSensorsPerSM(3,i);
    c5->Print(Form("%s.pdf",cname.Data()));
    delete c5;
  }

  TString lab(Form("cTemp2D_%s_%s",TInfo::Type(3),GetName()));
  TCanvas *c4 = new TCanvas(lab,lab,1200,800);
  const Int_t rns=fArr->GetEntries();
  Double_t min=-1;
  Double_t max=30;
  for (Int_t i=0;i<rns;++i) {
    TInfo *tinfo = dynamic_cast<TInfo*>(fArr->At(i));
    if (!tinfo)
      continue;
    TH2 *h=tinfo->GetHist(3);
    h->SetMinimum(min);
    h->SetMaximum(max);
    h->Draw("colz,text");
    TH2 *hm = GetMask();
    hm->DrawCopy("box,same");
    c4->Print(Form("%s.pdf",cname.Data()));
  }
  c4->Print(Form("%s.pdf]",cname.Data()));
}

TCanvas *TDraw::DrawOccRun() const
{
  TString lab(Form("cOccRun_%s",GetName()));
  TCanvas *c = new TCanvas(lab,lab,1200,800);
  TH1 *h = GetOccRun();
  h->SetMinimum(0);
  h->SetMaximum(1);
  h->Draw("P");
  TLine *l = new TLine(0,fMinFrac,GetNRuns(),fMinFrac);
  l->SetLineColor(2);
  l->SetLineWidth(2);
  l->SetLineStyle(9);
  l->Draw();
  if (fDoPrint)
    c->Print(Form("%s.pdf",c->GetName()));
  return c;
}

TCanvas *TDraw::DrawOccSensor2D() const
{
  TString lab(Form("cOccSensor2D_%s",GetName()));
  TCanvas *c = new TCanvas(lab,lab,1200,800);
  TH2 *h2 = GetOccSensor2D();
  h2->Draw("colz,text");
  if (fDoPrint)
    c->Print(Form("%s.pdf",c->GetName()));
  return c;
}

TCanvas *TDraw::DrawOcc2D() const
{
  TString lab(Form("cOcc2D_%s",GetName()));
  TCanvas *c = new TCanvas(lab,lab);
  TH2 *h = GetOcc2D();
  h->Draw("colz");
  if (fDoPrint)
    c->Print(Form("%s.pdf",c->GetName()));
  return c;
}

TCanvas *TDraw::DrawTSensorsPerSM(Int_t type, Int_t sm) const
{
  TString lab(Form("cTemp_%s_%s_%d",TInfo::Type(type),GetName(),sm));
  TCanvas *c = new TCanvas(lab,lab,1200,800);
  Double_t max=GetMaxT(type),maxp=-100;
  Double_t min=GetMinT(type),minp=+100;
  TLegend *leg = new TLegend(0.92,0.1,0.99,0.99);
  TObjArray arr;

  TH1 *hAv=GetTSM(sm,type);
  hAv->SetMarkerStyle(20);
  hAv->SetMarkerSize(1.4);
  hAv->SetMarkerColor(kBlack);
  hAv->SetLineColor(kBlack);
  hAv->SetLineWidth(3);
  hAv->SetName(Form("sm%d",sm));
  arr.Add(hAv);
  leg->AddEntry(hAv,hAv->GetName(),"p");

  for (Int_t i=0;i<8;++i) {
    TH1 *h=GetT(sm*8+i,type);
    Int_t col=i+2;
    h->SetMarkerStyle(24+(i%2));
    h->SetMarkerSize(1.2);
    h->SetMarkerColor(col);
    h->SetLineColor(col);
    h->SetLineWidth(3);
    h->SetName(Form("sensor%d",i));
    //h->Draw("same, p, hist");

    if (h->GetMaximum()>-1){
      arr.Add(h);
      Double_t minh=h->GetMinimum(-1);
      if (minh<minp) minp=minh;
      Double_t maxh=h->GetMaximum();
      if (maxh>maxp) maxp=maxh;
      leg->AddEntry(h,h->GetName(),"p");
    }
  }
  minp=0.9*minp; maxp=1.1*maxp;
  TH2 *h2f = new TH2F("h2f",";run idx;T",1,0,GetNRuns(),1,minp,maxp);
  h2f->SetTitle(Form("T per sensor in SM%d: run %d to %d, min=%.1f, max=%.1f",sm,GetRunNo(0), GetRunNo(GetNRuns()-1), min, max));
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

TCanvas *TDraw::DrawTPerSM(Int_t type) const
{
  TString lab(Form("cTemp_%s_%s",TInfo::Type(type),GetName()));
  TCanvas *c = new TCanvas(lab,lab,1200,800);
  Double_t max=GetMaxT(type),maxp=-100;
  Double_t min=GetMinT(type),minp=+100;
  TLegend *leg = new TLegend(0.92,0.1,0.99,0.99);
  TObjArray arr;
  for (Int_t i=0;i<20;++i) {
    TH1 *h=GetTSM(i,type);
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
    //h->Draw("same, p, hist");
    arr.Add(h);
    Double_t minh=h->GetMinimum();
    if (minh<minp) minp=minh;
    Double_t maxh=h->GetMaximum();
    if (maxh>maxp) maxp=maxh;
    leg->AddEntry(h,h->GetName(),"p");
  }
  minp=0.9*minp; maxp=1.1*maxp;
  TH2 *h2f = new TH2F("h2f",";run idx;T",1,0,GetNRuns(),1,minp,maxp);
  h2f->SetTitle(Form("Average T per SM: run %d to %d, min=%.1f, max=%.1f",GetRunNo(0), GetRunNo(GetNRuns()-1), min, max));
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

TCanvas *TDraw::DrawT2D(Int_t type) const
{
  TString lab(Form("cTemp2D%s_%s",TInfo::Type(type),GetName()));
  TCanvas *c = new TCanvas(lab,lab);
  TH2 *h = GetT2D(type);
  h->Draw("colz");
  if (fDoPrint)
    c->Print(Form("%s.pdf",c->GetName()));
  return c;
}

TH2 *TDraw::GetMask() const
{
  static TH2F *h = 0;
  if (h==0) {
    h = new TH2F("hBadMask",";col;rows",8,-0.5,7.5,20,-0.5,19.5);
    h->SetDirectory(0);
    h->SetStats(0);
    h->SetLineColor(0);
    h->SetMaximum(0);
    h->SetMinimum(-1);
    for (Int_t i=0;i<TInfo::NSensors();++i) {
      Double_t bin = TInfo::GetBin(i);
      if (!IsGood(i))
        h->SetBinContent(bin,-1);
      else
        h->SetBinContent(bin,-100);
    }
  }
  return h;
}

TH1 *TDraw::GetOccRun() const
{
  const Int_t rns=fArr->GetEntries();
  TH1F *ret = new TH1F("hOccRun",Form(";run idx;fraction"),rns,0,rns);
  ret->SetTitle(Form("Fraction of all sensors per run: run %d to %d, min.fraction=%.1f",GetRunNo(0), GetRunNo(GetNRuns()-1), fMinFrac));
  ret->SetDirectory(0);
  ret->SetStats(0);
  ret->SetMarkerStyle(20);
  ret->SetMarkerSize(1);
  Int_t r1,r2=0;
  for (Int_t i=0;i<rns;++i) {
    TInfo *tinfo = dynamic_cast<TInfo*>(fArr->At(i));
    if (!tinfo)
      continue;
    Double_t val=tinfo->Fraction();
    ret->SetBinContent(i+1,val);
  }
  return ret;
}

TH1 *TDraw::GetOccSensor() const
{
  const Int_t rns=fArr->GetEntries();
  TProfile *ret = new TProfile("hOccSensor",Form(";sensor id;fraction"),TInfo::NSensors(),0,TInfo::NSensors());
  ret->SetDirectory(0);
  ret->SetStats(0);
  ret->SetMarkerStyle(20);
  ret->SetMarkerSize(1);
  for (Int_t i=0;i<rns;++i) {
    TInfo *tinfo = dynamic_cast<TInfo*>(fArr->At(i));
    if (!tinfo)
      continue;
    for (Int_t j=0;j<TInfo::NSensors();++j) {
      Int_t val=tinfo->IsValid(j);
      ret->Fill(j,val);
    }
  }
  return ret;
}

TH2 *TDraw::GetOccSensor2D() const
{
  TH2F *hOcc2D = new TH2F("hOccSensor2D",";col;rows",8,-0.5,7.5,20,-0.5,19.5);
  hOcc2D->SetTitle(Form("Fraction of all runs per sensor: run %d to %d, min.fraction=%.1f",GetRunNo(0), GetRunNo(GetNRuns()-1), fMinFrac));
  hOcc2D->SetDirectory(0);
  hOcc2D->SetStats(0);
  hOcc2D->SetMinimum(-1);
  hOcc2D->SetMaximum(+1);
  for (Int_t ns=0;ns<TInfo::NSensors();++ns) {
    Double_t frac = GetFraction(ns);
    Int_t bin = TInfo::GetBin(ns);
    if (!IsGood(ns))
      frac = -1;
    hOcc2D->SetBinContent(bin,frac);
  }
  return hOcc2D;
}

TH2 *TDraw::GetOcc2D() const
{
  const Int_t rns=fArr->GetEntries();
  TH2F *ret = new TH2F("h2focc",Form("%s;run idx;sensor idx",GetName()),rns,0,rns,TInfo::NSensors(),0,TInfo::NSensors());
  ret->SetDirectory(0);
  ret->SetStats(0);
  for (Int_t i=0;i<rns;++i) {
    TInfo *tinfo = dynamic_cast<TInfo*>(fArr->At(i));
    if (!tinfo)
      continue;
    for (Int_t j=0;j<TInfo::NSensors();++j) {
      Int_t val=tinfo->IsValid(j);
      ret->SetBinContent(ret->FindBin(i,j),val);
    }
  }
  return ret;
}

Double_t TDraw::GetMaxT(Int_t type) const
{
  Double_t ret=-1e9;
  const Int_t rns=fArr->GetEntries();
  for (Int_t i=0;i<rns;++i) {
    TInfo *tinfo = dynamic_cast<TInfo*>(fArr->At(i));
    if (!tinfo)
      continue;
    Double_t v=tinfo->AbsMaxT(type);
    if (v>ret)
      ret=v;
  }
  return ret;
}

Double_t TDraw::GetMaxTperS(Int_t ns, Int_t type) const
{
  Double_t ret=-1e9;
  const Int_t rns=fArr->GetEntries();
  for (Int_t i=0;i<rns;++i) {
    TInfo *tinfo = dynamic_cast<TInfo*>(fArr->At(i));
    if (!tinfo)
      continue;
    if (!tinfo->IsValid(ns))
      continue;
    Double_t v=tinfo->T(ns,type);
    if (v>ret)
      ret=v;
  }
  return ret;
}

Double_t TDraw::GetMinT(Int_t type) const
{
  Double_t ret=1e9;
  const Int_t rns=fArr->GetEntries();
  for (Int_t i=0;i<rns;++i) {
    TInfo *tinfo = dynamic_cast<TInfo*>(fArr->At(i));
    if (!tinfo)
      continue;
    Double_t v=tinfo->AbsMinT(type);
    if (v<ret)
      ret=v;
  }
  return ret;
}

Double_t TDraw::GetMinTperS(Int_t ns, Int_t type) const
{
  Double_t ret=1e9;
  const Int_t rns=fArr->GetEntries();
  for (Int_t i=0;i<rns;++i) {
    TInfo *tinfo = dynamic_cast<TInfo*>(fArr->At(i));
    if (!tinfo)
      continue;
    if (!tinfo->IsValid(ns))
      continue;
    Double_t v=tinfo->T(ns,type);
    if (v<ret)
      ret=v;
  }
  return ret;
}

Int_t TDraw::GetNBad() const
{
  Int_t ret=0;
  const Int_t nbad = fBad.GetSize();
  for (Int_t ns=0,first=0;ns<nbad;++ns) {
    if (IsGood(ns))
      continue;
    ++ret;
  }
  return ret;
}

TH1 *TDraw::GetTSM(Int_t sm, Int_t type) const
{
  const Int_t rns=fArr->GetEntries();
  TProfile *ret = new TProfile(Form("h%d%d",sm,type),Form(";run idx"),rns,0,rns);
  ret->SetDirectory(0);
  ret->SetStats(0);
  for (Int_t i=0;i<rns;++i) {
    TInfo *tinfo = dynamic_cast<TInfo*>(fArr->At(i));
    if (!tinfo)
      continue;
    Double_t val = tinfo->AvgTempSM(sm);
    cout <<  i << "\t" << sm << "\t" << val << endl;
    ret->Fill(i,val);
  }
  return ret;
}

TH1 *TDraw::GetT(Int_t nsensor, Int_t type) const
{
  const Int_t rns=fArr->GetEntries();
  TProfile *ret = new TProfile(Form("h%d%d",nsensor,type),Form(";run idx"),rns,0,rns);
  ret->SetDirectory(0);
  ret->SetStats(0);
  for (Int_t i=0;i<rns;++i) {
    TInfo *tinfo = dynamic_cast<TInfo*>(fArr->At(i));
    if (!tinfo)
      continue;
    Double_t val = -1;
    if (tinfo->IsValid(nsensor))
      val = tinfo->T(nsensor,type);
    ret->Fill(i,val);
  }
  return ret;
}

TH2 *TDraw::GetT2D(Int_t type) const
{
  const Int_t rns=fArr->GetEntries();
  TString title=TInfo::Type(type);
  TH2F *ret = new TH2F(Form("h2f_%d",type),Form("%s: %s;run idx;sensor idx",GetName(),title.Data()),rns,0,rns,TInfo::NSensors(),0,TInfo::NSensors());
  ret->SetDirectory(0);
  ret->SetStats(0);
  for (Int_t i=0;i<rns;++i) {
    TInfo *tinfo = dynamic_cast<TInfo*>(fArr->At(i));
    if (!tinfo)
      continue;
    for (Int_t j=0;j<TInfo::NSensors();++j) {
      if (!tinfo->IsValid(j))
	continue;
      Double_t val = tinfo->T(j,type);
      ret->SetBinContent(ret->FindBin(i,j),val);
    }
  }
  return ret;
}

Int_t TDraw::GetNRunsValid(Int_t ns) const
{
  Int_t ret = 0;
  const Int_t rns=fArr->GetEntries();
  for (Int_t i=0;i<rns;++i) {
    TInfo *tinfo = dynamic_cast<TInfo*>(fArr->At(i));
    if (!tinfo)
      continue;
    ret += tinfo->IsValid(ns);
  }
  return ret;
}

void TDraw::Print(Option_t *opt) const
{
  const Int_t rns=fArr->GetEntries();
  cout << "Period: " << GetName() << " with " << rns << " runs" << endl;
  if (strlen(opt)>0)
    for (Int_t ns=0;ns<TInfo::NSensors();++ns)
      cout << "Channel " << ns << " " << GetFraction(ns) << endl;
  for (Int_t i=0;i<rns;++i) {
    TInfo *tinfo = dynamic_cast<TInfo*>(fArr->At(i));
    if (!tinfo)
      continue;
    cout << "-> " << i << ": run=" << tinfo->GetRunNo() << " frac=" << tinfo->Fraction() << " minT=" << tinfo->AbsMinT()<< " maxT=" << tinfo->AbsMaxT() << endl;
    if (strlen(opt)>0)
      tinfo->Print();
  }
}

void plotT_period(const char *period, Bool_t doprint=0)
{
  TDraw d(period);
  d.SetPrint(doprint);
  d.Compute();
  d.Print();
  if (0) {
    d.DrawOccRun();
    d.DrawOccSensor2D();
    d.DrawT2D(3);
    d.DrawTPerSM(3);
    for (Int_t i = 0; i < 20; i++)
      d.DrawTSensorsPerSM(3,i);
  } else
    d.DrawAll();
}

void plotOCDB_Temperature(const char *period="lhc18d")
{
  plotT_period(period);
}

