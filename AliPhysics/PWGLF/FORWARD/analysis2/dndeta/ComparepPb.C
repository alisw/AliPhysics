#include <TMultiGraph.h>
#include <TString.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TError.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TH1.h>
#include <TAxis.h>
#include <THStack.h>

//____________________________________________________________________
/** 
 * Get the marker style 
 * 
 * @param sys   Collision system 
 * @param style Base style 
 * 
 * @return Marker style 
 */
Int_t Marker(Int_t sys, Int_t style)
{
  return (sys == 4 ? (style == 23 ? 9 : 4) : 0) + style;
}

//____________________________________________________________________
/** 
 * Get the data for forward 
 * 
 * @param sys   System
 * @param meth  Centrality method 
 * @param style Base style 
 * 
 * @return TMultiGraph or null
 */
TMultiGraph*
GetOneFwd(UShort_t sys, const char* meth, Int_t style)
{

  
  TString what(meth);
  what.ReplaceAll("X", (sys == 3 ? "A" : "C"));
  TString sysN(sys == 3 ? "pPb" : "Pbp");

  TString cmd;
  cmd = Form("Drawer::GetStack(0,\"%s\",5023,\"CENT%s\")",
	     sysN.Data(), meth);

  Long_t ret = gROOT->ProcessLine(cmd);
  if (!ret) {
    Warning("GetOne", "didn't get anything back for %s", what.Data());
    return 0;
  }
  THStack* stack = reinterpret_cast<THStack*>(ret); 
  if (!stack) {
    Warning("GetOne", "didn't get anything back for %s", what.Data());
    return 0;
  }
  Int_t sty = Marker(sys, style);
  TMultiGraph* mg = new TMultiGraph;
  TIter nextH(stack->GetHists());
  TH1*  hist = 0;
  while ((hist = static_cast<TH1*>(nextH()))) {
    ret = gROOT->ProcessLine(Form("Drawer::H2G((TH1*)%p,0)",hist));
    if (!ret) {
      Warning("GetOne", "Couldn't convert %s", hist->GetName());
      continue;
    }
    TGraphErrors* g = reinterpret_cast<TGraphErrors*>(ret);
    g->SetMarkerStyle(sty);
    mg->Add(g);
  }

  Info("", "Got %d graphs", mg->GetListOfGraphs()->GetEntries());
  mg->GetListOfGraphs()->Sort();
  mg->SetName(Form("fwd_%s_%s", sysN.Data(), meth));
  return mg;
}
//____________________________________________________________________
/** 
 * Get the data for central 
 * 
 * @param sys   System
 * @param meth  Centrality method 
 * @param style Base style 
 * 
 * @return TMultiGraph or null
 */
TMultiGraph*
GetOneCen(UShort_t sys, const char* meth, Int_t style)
{
  TString what(meth);
  what.ReplaceAll("X", (sys == 3 ? "A" : "C"));
  UShort_t cnt = 0;
  if      (what.Contains("V0M"))   cnt = (0x01 << 4);
  else if (what.Contains("V0A"))   cnt = (0x02 << 4);
  else if (what.Contains("ZNA"))   cnt = (0x04 << 4);
  else if (what.Contains("ZNC"))   cnt = (0x08 << 4);
  else if (what.Contains("V0C"))   cnt = (0x10 << 4);
  else if (what.Contains("CL1"))   cnt = (0x20 << 4);
  Info("GetOne", "%d/%s -> 0x%04x", sys, meth, cnt);

  Long_t ret =
    gROOT->ProcessLine(Form("RefData::GetData(%d,5023,0x%x,0,100,0xf);",
			    sys, cnt));
  if (!ret) {
    Warning("GetOne", "didn't get anything back for %s", what.Data());
    return 0;
  }
  TMultiGraph* mg = reinterpret_cast<TMultiGraph*>(ret); 
  if (!mg) {
    Warning("GetOne", "didn't get anything back for %s", what.Data());
    return 0;
  }

  Int_t sty = Marker(sys, style);
  TIter next(mg->GetListOfGraphs());
  TGraph* g = 0;
  while ((g = static_cast<TGraph*>(next())))
    g->SetMarkerStyle(sty);
  Info("", "Got %d graphs", mg->GetListOfGraphs()->GetEntries());
  mg->GetListOfGraphs()->Sort();
  mg->SetName(Form("fwd_%s_%s", (sys == 3 ? "pPb" : "Pbp"), meth));
  return mg;
}

//____________________________________________________________________
/** 
 * Flip a graph round @f$ \eta=0@f$
 * 
 * @param g Graph 
 */
void
FlipG(TGraphErrors* g)
{
  Int_t n = g->GetN();
  for (Int_t i = 0; i < n/2; i++) {
    Int_t    i1   = i;
    Int_t    i2   = n-i-1;
    if (i1 == i2) continue;
    Double_t x1   = g->GetX()[i1];
    Double_t x2   = g->GetX()[i2];
    Double_t y1   = g->GetY()[i1];
    Double_t y2   = g->GetY()[i2];
    Double_t ex1 = g->GetEX()[i1];
    Double_t ex2 = g->GetEX()[i2];
    Double_t ey1 = g->GetEY()[i1];
    Double_t ey2 = g->GetEY()[i2];
    g->SetPoint(i1, -x2, y2);
    g->SetPointError(i1, ex2, ey2);
    g->SetPoint(i2, -x1, y1);
    g->SetPointError(i2, ex1, ey1);
  }
}
//____________________________________________________________________
/** 
 * Flip a graph round @f$ \eta=0@f$
 * 
 * @param g Graph 
 */
void
FlipG(TGraphAsymmErrors* g)
{
  Int_t n = g->GetN();
  for (Int_t i = 0; i < n/2; i++) {
    Int_t    i1   = i;
    Int_t    i2   = n-i-1;
    if (i1 == i2) continue;
    Double_t x1   = g->GetX()[i1];
    Double_t x2   = g->GetX()[i2];
    Double_t y1   = g->GetY()[i1];
    Double_t y2   = g->GetY()[i2];
    Double_t exl1 = g->GetEXlow()[i1];
    Double_t exl2 = g->GetEXlow()[i2];
    Double_t exh1 = g->GetEXhigh()[i1];
    Double_t exh2 = g->GetEXhigh()[i2];    
    Double_t eyl1 = g->GetEYlow()[i1];
    Double_t eyl2 = g->GetEYlow()[i2];
    Double_t eyh1 = g->GetEYhigh()[i1];
    Double_t eyh2 = g->GetEYhigh()[i2];
    g->SetPoint(i1, -x2, y2);
    g->SetPointError(i1, exl2, exh2, eyl2, eyh2);
    g->SetPoint(i2, -x1, y1);
    g->SetPointError(i2, exl1, exh1, eyl1, eyh1);
  }
}
//____________________________________________________________________
/** 
 * Flip a graph round @f$ \eta=0@f$
 * 
 * @param g Graph 
 */
void
FlipG(TGraph* g)
{
  if (g->IsA()->InheritsFrom(TGraphAsymmErrors::Class()))
    FlipG(static_cast<TGraphAsymmErrors*>(g));
  else if (g->IsA()->InheritsFrom(TGraphErrors::Class()))
    FlipG(static_cast<TGraphErrors*>(g));
  else
    Warning("FlipG", "Don't know how to flip a %s", g->ClassName());
}
  
//____________________________________________________________________
/** 
 * Flip all graphs in a multi-graph
 * 
 * @param mg Mulitgraph
 * 
 * @return New multi-graph 
 */
TMultiGraph*
FlipMG(TMultiGraph* mg)
{
  TMultiGraph* ret = static_cast<TMultiGraph*>(mg->Clone());
  TIter next(ret->GetListOfGraphs());
  TGraph* g = 0;
  while ((g = static_cast<TGraph*>(next())))
    FlipG(g);
  return ret;
}

//____________________________________________________________________
/** 
 * Take the ratio of two graphs 
 * 
 * @param num    Numerator - on return the ratio 
 * @param denom  Denominator 
 */
void
RatioGG(TGraphErrors* num, TGraph* denom)
{
  Info("RatioGG(GE)", "%s (%s:%d)/ %s (%s:%d)",
       num->GetName(), num->ClassName(), num->GetN(),
       denom->GetName(), denom->ClassName(), denom->GetN());
  Int_t n = num->GetN();
  Double_t dh = denom->GetX()[n-1]+denom->GetEX()[n-1];
  Double_t dl = denom->GetX()[0]  +denom->GetEX()[0];
  TBits clean(n);
  for (Int_t i = 0; i < n; i++) {
    Double_t x   = num->GetX()[i];
    if (x > dh || x < dl) {
      clean.SetBitNumber(i);
      continue;
    }
    Double_t d   = denom->Eval(x);
    Double_t y   = num->GetY()[i];
    Double_t ex  = num->GetEX()[i];
    Double_t ey  = num->GetEY()[i];
    num->SetPoint(i, x, y/d);
    num->SetPointError(i, ex, ey/d);
  }
  for (Int_t i = clean.GetNbits()-1; i >= 0; i--) {
    if (clean.TestBitNumber(i)) num->RemovePoint(i);
  }
}
//____________________________________________________________________
/** 
 * Take the ratio of two graphs 
 * 
 * @param num    Numerator - on return the ratio 
 * @param denom  Denominator 
 */
void
RatioGG(TGraphAsymmErrors* num, TGraph* denom)
{
  Info("RatioGG(GA)", "%s (%s:%d)/ %s (%s:%d)",
       num->GetName(), num->ClassName(), num->GetN(),
       denom->GetName(), denom->ClassName(), denom->GetN());
  Int_t n = num->GetN();
  for (Int_t i = 0; i < n; i++) {
    Double_t x   = num->GetX()[i];
    Double_t d   = denom->Eval(x);
    Double_t y   = num->GetY()[i];
    Double_t exl = num->GetEXlow()[i];
    Double_t exh = num->GetEXhigh()[i];
    Double_t eyl = num->GetEYlow()[i];
    Double_t eyh = num->GetEYhigh()[i];
    num->SetPoint(i, x, y/d);
    num->SetPointError(i, exl, exh, eyl/d, eyh/d);
  }  
}
//____________________________________________________________________
/** 
 * Take the ratio of two graphs 
 * 
 * @param num    Numerator - on return the ratio 
 * @param denom  Denominator 
 */
void
RatioGG(TGraph* num, TGraph* denom)
{
  if (num->IsA()->InheritsFrom(TGraphAsymmErrors::Class()))
    RatioGG(static_cast<TGraphAsymmErrors*>(num),denom);
  else if (num->IsA()->InheritsFrom(TGraphErrors::Class()))
    RatioGG(static_cast<TGraphErrors*>(num),denom);
  else
    Warning("RatioGG", "Don't know how to ratio a %s to a %s",
	    num->ClassName(), denom->ClassName());
}


//____________________________________________________________________
/** 
 * Take the ratio of all graphs in two multi-graphs 
 * 
 * @param num    Multi-graph of numerators 
 * @param denom  Multi-graph of denominators 
 * 
 * @return new multi-graph of ratios 
 */
TMultiGraph*
RatioMG(TMultiGraph* num, TMultiGraph* denom)
{
  TMultiGraph* ret = static_cast<TMultiGraph*>(num->Clone());
  TIter nextN(ret->GetListOfGraphs());
  TIter nextD(denom->GetListOfGraphs());
  TGraph* n = 0;
  TGraph* d = 0;
  while ((n = static_cast<TGraph*>(nextN())) &&
	 (d = static_cast<TGraph*>(nextD())))
    RatioGG(n, d);
  return ret;
}

//____________________________________________________________________
/** 
 * Fix the frame of a multi-graph
 * 
 * @param o      The frame 
 * @param hasFwd Whether we have forward data or not 
 * @param fac    Scale factor for y axis
 * @param ytitle Title on y axis 
 */
void
FixFrame(TObject* o,
	 Bool_t   hasFwd,
	 Double_t fac=1.2,
	 const char* ytitle="1/#it{N} d#it{N}_{ch}/d#it{#eta}")
{
  static Int_t cnt = 0;
  TMultiGraph* mg = static_cast<TMultiGraph*>(o);
  TH1* h = mg->GetHistogram();
  mg->SetTitle("");
  if (!h) return;

  
  Double_t etaMin = (hasFwd ? -6 : -2.5);
  Double_t etaMax = (hasFwd ? +6 : +2.5);
  h->GetXaxis()->Set((etaMax-etaMin)/0.5, etaMin, etaMax);
  h->Rebuild();
  h->SetName(Form("%s_%d", o->GetName(), cnt++));
  h->SetTitle("");
  h->SetXTitle("#it{#eta}");
  h->SetYTitle(ytitle);
  h->GetXaxis()->SetNdivisions(210);
  h->GetYaxis()->SetNdivisions(210);
  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetXaxis()->SetTitleFont(42);
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetTitleOffset(0.7);
  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleOffset(0.7);
  h->SetMaximum(fac*h->GetMaximum());
  h->SetMinimum(0.01);
}

//____________________________________________________________________
/** 
 * Extract the centrality string from a name 
 * 
 * @param name The name 
 * 
 * @return Centrality string 
 */
const char*
ExtractCent(const char* name)
{
  TString n(name);
  TObjArray*  a   = n.Tokenize("_");
  TObjString* oc1 = static_cast<TObjString*>(a->At(1));
  TObjString* oc2 = static_cast<TObjString*>(a->At(2));
  TString     sc1 = oc1->String().Strip(TString::kLeading, '0');
  TString     sc2 = oc2->String().Strip(TString::kLeading, '0');
  delete a;
  return Form("%3d-%3d%%", sc1.Atoi(), sc2.Atoi());
}

//____________________________________________________________________
/** 
 * Compare one set of centrality graphs 
 * 
 * @param meth   The method to use 
 * @param style  The base style 
 */
void
CompareOne(const char* meth, Int_t style)
{
  TMultiGraph* fwdpPb = GetOneFwd(3, meth, style);
  TMultiGraph* fwdPbp = GetOneFwd(4, meth, style);
  Bool_t       hasFwd = true;
  if (!fwdpPb || !fwdPbp) {
    Warning("CompareOne", "FWD pPb=%p Pbp=%p", fwdpPb, fwdPbp);
    hasFwd = false;
  }
  TMultiGraph* cenpPb = GetOneCen(3, meth, style);
  TMultiGraph* cenPbp = GetOneCen(4, meth, style);
  Bool_t       hasCen = true;
  if (!cenpPb || !cenPbp) {
    Warning("CompareOne", "CEN pPb=%p Pbp=%p", cenpPb, cenPbp);
    hasCen = false;
  }
  if (!hasCen && !hasFwd) return;
  
  TCanvas* c = new TCanvas(meth, meth, 900, 1000);
  c->SetTopMargin(0.01);
  c->SetRightMargin(0.01);
  c->cd();
  c->Divide(1,3,0,0);
  Double_t lx = 0.8;
  c->GetPad(1)->SetRightMargin(1-lx+.01);
  c->GetPad(2)->SetRightMargin(1-lx+.01);
  c->GetPad(3)->SetRightMargin(1-lx+.01);
  c->GetPad(1)->SetGridx();
  c->GetPad(2)->SetGridx();
  c->GetPad(3)->SetGridx();
  // c->GetPad(2)->SetBottomMargin(0.1);

  c->cd(1);
  TMultiGraph* frame = 0;
  if (hasFwd) {
    fwdpPb->Draw("pa");
    fwdPbp->Draw("p");
    frame = fwdpPb;
  }
  if (hasCen) {
    cenpPb->Draw(frame ? "p" : "pa");
    cenPbp->Draw("p");
    if (!frame) frame = cenpPb;
  }    
  FixFrame(frame, hasFwd, 1.2);

  TLegend* l = new TLegend(lx, 0.01, 0.99, 0.99,
			   Form("%s estimator", meth));
  l->SetFillStyle(0);
  l->SetBorderSize(0);
  l->SetTextFont(42);
  l->SetNColumns(1);
  TLegendEntry* e = l->AddEntry("", "p-Pb", "p");
  e->SetMarkerStyle(Marker(3, style));
  e = l->AddEntry("", "Pb-p", "p");
  e->SetMarkerStyle(Marker(4, style));
  l->Draw();

  c->cd(2);
  TMultiGraph* ffwdPbp = 0;
  if (hasFwd) {
    ffwdPbp = FlipMG(fwdPbp);
    ffwdPbp->SetName(Form("fwd_%s_fPbp", meth));
    fwdpPb->Draw("ap");
    ffwdPbp->Draw("p");
  }
  TMultiGraph* fcenPbp = 0;
  if (hasCen) {
    fcenPbp = FlipMG(cenPbp);
    fcenPbp->SetName(Form("cen_%s_fPbp", meth));
    cenpPb->Draw(hasFwd ? "p" : "ap");
    fcenPbp->Draw("p");
  }

  TLegend* ll = new TLegend(lx, 0.01, 0.99, 0.9,
			    "Pb-p flipped around #it{#eta}=0");
  ll->SetFillStyle(0);
  ll->SetBorderSize(0);
  ll->SetTextFont(42);
  ll->SetNColumns(1);
  TIter next(frame->GetListOfGraphs());
  TGraph* g = 0;
  while ((g = static_cast<TGraph*>(next()))) {
    e = ll->AddEntry("", ExtractCent(g->GetName()), "f");
    e->SetFillColor(g->GetMarkerColor());
    e->SetFillStyle(1001);
    e->SetLineColor(g->GetMarkerColor());
  }
  ll->Draw();
  
  c->cd(3);
  frame = 0;

  Double_t minEta = (hasFwd ? -3.5 : -2.);
  Double_t maxEta = (hasFwd ? +5.0 : +2.);
  TGraphErrors* band = new TGraphErrors(2);
  band->SetPoint(0, minEta, 1);
  band->SetPoint(1, maxEta, 1);
  band->SetPointError(0, 0, .05);
  band->SetPointError(1, 0, .05);
  band->SetFillColor(kYellow-8);
  band->SetFillStyle(3001);
  band->SetLineColor(kYellow-8);
  
  if (hasFwd) {
    TMultiGraph* fwdRat = RatioMG(ffwdPbp, fwdpPb);
    fwdRat->SetName(Form("fwd_%s_rat", meth));
    fwdRat->GetListOfGraphs()->AddFirst(band, "3");
    frame = fwdRat;
    fwdRat->Draw("ap");
  }
  if (hasCen) {
    TMultiGraph* cenRat = RatioMG(fcenPbp, cenpPb);
    cenRat->SetName(Form("cen_%s_rat", meth));
    if (!frame) {
      cenRat->GetListOfGraphs()->AddFirst(band, "3");
      frame = cenRat;
    }
    cenRat->Draw(hasFwd ? "p" : "ap");
  }  
  FixFrame(frame, hasFwd, 1, "Flipped Pb-p over p-Pb");
  if (meth[0] != 'C') {
    frame->GetHistogram()->SetMinimum(0.74);
    frame->GetHistogram()->SetMaximum(1.26);
  }
  TLegend* lll = new TLegend(lx, 0.4, 0.99, 0.5);
  lll->SetFillStyle(0);
  lll->SetBorderSize(0);
  lll->SetTextFont(42);
  lll->SetNColumns(1);
  lll->AddEntry(band, "5% error band", "F");
  lll->Draw();
  
  c->Print(Form("%s.pdf", meth));
}

//____________________________________________________________________
void
ComparepPb()
{
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/FORWARD/analysis2/dndeta/Drawer.C+g");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/FORWARD/analysis2/OtherData.C+g");

  const char* meths[]  = { "V0M", "CL1", "V0X", "ZNX", 0 };
  Int_t       styles[] = { 20,    21,    22,    23,    0 };

  const char** pmeth   = meths;
  Int_t*       pstyle  = styles;
  
  while (*pmeth && *pstyle) {
    CompareOne(*pmeth, *pstyle);
    pmeth++;
    pstyle++;
  }
}
//____________________________________________________________________
//
// EOF
// 
