#ifndef DRAWER_C
#define DRAWER_C 
#ifndef __CINT__
# include <TROOT.h>
# include <TSystem.h>
# include <TInterpreter.h>
# include <TString.h>
# include <THStack.h>
# include <TMultiGraph.h>
# include <TGraph.h>
# include <TGraphErrors.h>
# include <TGraphAsymmErrors.h>
# include <TList.h>
# include <TError.h>
# include <TMath.h>
# include <TLegend.h>
# include <TLegendEntry.h>
# include <TCanvas.h>
# include <TClass.h>
# include <TFile.h>
# include <TParameter.h>
# include <TColor.h>
# include <TStyle.h>
# include <TLatex.h>
# include <TMap.h>
# include <map>
# include <vector>
# include "Combiner.C"
#else 
class TLatex;
class TCanvas;
class THStack;
class TMultiGraph;
class TGraph;
class TGraphErrors;
class TGraphAsymmErrors;
class TLegend;
class TH1;
class TPair;
class TString;
#endif

/**
 * Utility class 
 * 
 * @deprecated Use new GSE based drawing 
 */
struct Drawer {
  static const char* PlotPrefix() { return "plots"; }
  /** 
   * Make sure we have loaded the RefData class 
   * 
   */
  static void LoadOther() {
#if 0
    if (gROOT->GetClass("RefData")) return;
    const char* fwd = 0;
    if (gSystem->Getenv("FWD"))
      fwd = gSystem->Getenv("FWD");
    else 
      fwd = gSystem->ExpandPathName("$ALICE_PHYSICS/PWGLF/FORWARD/analysis2");
    gROOT->SetMacroPath(Form("%s:%s", gROOT->GetMacroPath(),fwd));
    TString path(gROOT->GetMacroPath());
    // if (!path.Contains(fwd)) {
    path.Append(Form(":%s", fwd));
    gROOT->SetMacroPath(path);
    // }
    gROOT->LoadMacro("OtherData.C+");
#endif 
  }
  /** 
   * Get other data 
   * 
   * @param system   System
   * @param sNN      Center of mass energy 
   * @param trigger  Trigger type 
   * 
   * @return Pointer to TMultiGraph or null
   */
  static TMultiGraph* GetOther(UShort_t       sys,
			       UShort_t       sNN,
			       UShort_t       trg,
			       UShort_t       exps=0xf,
			       Int_t          verbose=3)
  {
    return 0;
#if 0
    UShort_t c1  = 0;
    UShort_t c2  = 0;
    switch (trg) {
    case 0x004: if (sys == 3) exps = 0x4; break;
    case 0x010:
    case 0x020:
    case 0x040:
    case 0x080:
    case 0x100:
    case 0x200: c2 = 100; break;
    }

    LoadOther();
    gROOT->ProcessLine(Form("RefData::Verbose(%d);", verbose));
    ::Info("", "sys=%d sNN=%d trg=0x%x c1=%d c2=%d exp=0x%x",
	   sys, sNN, trg, c1, c2, exps);
    Long_t   ret   = 
      gROOT->ProcessLine(Form("RefData::GetData(%d,%d,%d,%d,%d,%d);",
                              sys,sNN,trg,c1,c2,exps));
    if (!ret) {
      if (verbose)
	 Warning("", "RefData::GetData(%d,%d,0x%x,%d,%d,0x%x); failed",
		 sys,sNN,trg,c1,c2,exps);
      return 0;
    }

    TMultiGraph* mg = reinterpret_cast<TMultiGraph*>(ret);  

    return mg;
#endif
  }
  /** 
   * Get Work-in-Progress data from other PWGs
   * 
   * @param system 
   * @param sNN 
   * @param trigger 
   * 
   * @return 
   */
  static TMultiGraph* GetOther(const TString& system, 
			       UShort_t       sNN, 
			       const TString& trigger,
			       Int_t          verbose=0)
  {
    return GetOther(system, sNN, trigger, "ALICE WIP", verbose);
  }
  /** 
   * Get other data 
   * 
   * @param system   System
   * @param sNN      Center of mass energy 
   * @param trigger  Trigger type 
   * 
   * @return Pointer to TMultiGraph or null
   */
  static TMultiGraph* GetOther(const TString& system, 
			       UShort_t       sNN, 
			       const TString& trigger,
			       const TString& exps,
			       Int_t          verbose=0)
  {
    TString s(system); 
    s.ToUpper();

    UShort_t sys = 0;
    if      (s.EqualTo("PP"))   sys = 1;
    else if (s.EqualTo("PPB"))  sys = 3;
    else if (s.EqualTo("PBPB")) sys = 2;
    else if (s.EqualTo("PBP"))  sys = 4;

    TString e(exps);
    e.ToUpper();
    UShort_t exp = 0x0;
    if (e.Contains("UA5"))   exp |= 0x01;
    if (e.Contains("CMS"))   exp |= 0x02;
    if (e.Contains("ALICE")) exp |= 0x04;
    if (e.Contains("WIP"))   exp |= 0x08;
    if (e.Contains("EG"))    exp |= 0x10;
    
    TString t(trigger); 
    t.ToUpper();
    UShort_t trg = 0;
    if      (t.EqualTo("INEL"))                            trg = 0x1;
    else if (t.EqualTo("INELGT0") || t.EqualTo("INEL>0"))  trg = 0x2;
    else if (t.EqualTo("NSD") || t.EqualTo("V0AND"))       trg = 0x4;
    else if (t.BeginsWith("CENT")) { 
      if      (t.EndsWith("V0M")) trg |= 0x010;
      else if (t.EndsWith("V0A")) trg |= 0x020;
      else if (t.EndsWith("V0X")) trg |= (sys == 3 ? 0x020 : 0x100);
      else if (t.EndsWith("ZNA")) trg |= 0x040;
      else if (t.EndsWith("ZNC")) trg |= 0x080;
      else if (t.EndsWith("V0C")) trg |= 0x100;
      else if (t.EndsWith("ZNX")) trg |= (sys == 3 ? 0x040 : 0x80);
      else if (t.EndsWith("MB"))  { trg =  0x4; exp = 0x4; }
    }

    ::Info("", "sys=%d sNN=%d trg=0x%x exp=0x%x", sys, sNN, trg, exp);
    TMultiGraph* mg = GetOther(sys, sNN, trg, exp, verbose);
    if (!mg) 
      Warning("GetOthers", 
	      "No other data for %s %s %d (%s)", 
	      system.Data(), trigger.Data(), sNN, exps.Data());
    return mg;
  }
  /** 
   * Get ALICE Blue
   * 
   * @return color number 
   */
  static Int_t AliceBlue()
  {
    return TColor::GetColor(40,   58, 68);
  }
  /** 
   * Get ALICE Red
   * 
   * @return color number 
   */
  static Int_t AliceRed()
  {
    return TColor::GetColor(226,   0, 26);
  }
  /** 
   * Get ALICE Purple 
   * 
   * @return color number 
   */
  static Int_t AlicePurple()
  {
    return TColor::GetColor(202,  71, 67);
  }
  /** 
   * Get ALICE Yellow 
   * 
   * @return color number 
   */
  static Int_t AliceYellow()
  {
    return TColor::GetColor(238, 125, 17);
  }
  static Color_t Brighten(Color_t origNum, Int_t nTimes=2)
  {
    Int_t origR, origG, origB;
    TColor::Pixel2RGB(TColor::Number2Pixel(origNum), origR, origG, origB);
    Int_t off    = nTimes*0x33;
    Int_t newR   = TMath::Min((origR+off),0xff);
    Int_t newG   = TMath::Min((origG+off),0xff);
    Int_t newB   = TMath::Min((origB+off),0xff);
    Int_t newNum = TColor::GetColor(newR, newG, newB);
    return newNum;
  }
  static TString SNNString(const TString& system, UShort_t sNN)
  {
    TString e(system.EqualTo("pp", TString::kIgnoreCase)
	      ? "#sqrt{s}=" : "#sqrt{s_{NN}}=");    
    if      (sNN < 1000)       e.Append(Form("%dGeV", sNN));
    else if (sNN % 1000 == 0)  e.Append(Form("%dTeV", sNN/1000));
    else                       e.Append(Form("%.2fTeV", float(sNN)/1000));
    return e;
  }
  
  /** 
   * Convert a histogram to a graph 
   * 
   * @param h    Histogram
   * @param xOff X offset 
   * @param sign if negative, mirror around this point
   * 
   * @return Newly created graph 
   */
  static TGraphErrors* H2G(TH1* h, Double_t xOff, Int_t sign=1)
  {
    TGraphErrors* g = new TGraphErrors();
    g->SetMarkerStyle(h->GetMarkerStyle());
    g->SetMarkerColor(h->GetMarkerColor());
    g->SetLineColor(h->GetLineColor());
    g->SetName(h->GetName());
    g->SetTitle(h->GetTitle());
  
    Int_t iP = 0;
    Int_t nX = h->GetNbinsX();
    for (Int_t iX = 1; iX <= nX; iX++) { 
      Double_t y = h->GetBinContent(iX);
      if (y < 1e-6) continue;
      Double_t x  = h->GetXaxis()->GetBinCenter(iX) - xOff;
      Double_t ex = h->GetXaxis()->GetBinWidth(iX)/2;
      Double_t ey = h->GetBinError(iX);

      g->SetPoint(iP, sign * x, y);
      g->SetPointError(iP++, ex, ey);
    }
    return g;
  }

  /** 
   * Convert a histogram to a graph 
   * 
   * @param h    Histogram
   * @param xOff X offset 
   * @param sign if negative, mirror around this point
   * 
   * @return Newly created graph 
   */
  static TGraphAsymmErrors* H2GA(TH1* h, TH1* sA, TH1* sC)
  {
    TGraphAsymmErrors* g = new TGraphAsymmErrors();
    g->SetMarkerStyle(h->GetMarkerStyle());
    g->SetMarkerColor(h->GetMarkerColor());
    g->SetLineColor(h->GetLineColor());
    g->SetName(h->GetName());
    g->SetTitle(h->GetTitle());
  
    Int_t iP = 0;
    Int_t nX = h->GetNbinsX();
    for (Int_t iX = 1; iX <= nX; iX++) { 
      Double_t y = h->GetBinContent(iX);
      if (y < 1e-6) continue;
      Double_t x  = h->GetXaxis()->GetBinCenter(iX);
      Double_t ex = h->GetXaxis()->GetBinWidth(iX)/2;
      Double_t ec = 0;
      Double_t ee = 0;
      if (x < 0) {
	ec = sC->GetBinContent(iX);
	ee = sC->GetBinError(iX);
      }
      else {
	ec = sA->GetBinContent(iX);
	ee = sA->GetBinError(iX);
      }
      Double_t eyl = y-(ec-ee);
      Double_t eyh = (ec+ee)-y;
      
      g->SetPoint(iP, x, y);
      g->SetPointError(iP++, ex, ex, eyl, eyh);
    }
    return g;
  }
  /** 
   * Get error bands from vanilla graph
   *
   * @param g 
   * @param low
   * @param up
   */
  static void ErrorGraphs(const TGraph* g, TGraph*& low, TGraph*& up)
  {
    Warning("ErrorGraphs", "Called with vanila TGraph (%s)",
	    g->IsA()->GetName());
    Int_t n = g->GetN();
    low = new TGraph(n);
    up  = new TGraph(n);
    for (Int_t i = 0; i < n; i++) {
      low->SetPoint(i, 0, 0);
      up->SetPoint(i, 0, 0);
    }
  }
  /** 
   * Get error bands 
   */
  static void ErrorGraphs(const TGraphErrors* g, TGraph*& low, TGraph*& up)
  {
    // Info("ErrorGraphs", "Called with TGraphErrors");
    Int_t n = g->GetN();
    low = new TGraph(n);
    up  = new TGraph(n);
    for (Int_t i = 0; i < n; i++) {
      low->SetPoint(i, g->GetX()[i], g->GetY()[i] - g->GetEY()[i]);
      up->SetPoint(i, g->GetX()[i], g->GetY()[i] + g->GetEY()[i]);
    }
  }
  /** 
   * Get error bands 
   */
  static void ErrorGraphs(const TGraphAsymmErrors* g, TGraph*& low, TGraph*& up)
  {
    // Info("ErrorGraphs", "Called with TGraphAsymmErrors");
    Int_t n = g->GetN();
    low = new TGraph(n);
    up  = new TGraph(n);
    for (Int_t i = 0; i < n; i++) {
      low->SetPoint(i, g->GetX()[i], g->GetY()[i] - g->GetEYlow()[i]);
      up->SetPoint(i, g->GetX()[i], g->GetY()[i] + g->GetEYhigh()[i]);
    }
  }
  /** 
   * Make a histogram of systematic errors with a fixed value 
   * 
   * @param h 
   * @param aSide 
   * @param factor 
   * 
   * @return 
   */
  static TH1* ErrorHist(const TH1* h, Bool_t aSide, Double_t factor=0.076)
  {
    TString nm(h->GetName());
    nm.ReplaceAll("Forward", "SysError");
    nm.Append(aSide ? "_a" : "_c");
    TH1* ret = static_cast<TH1*>(h->Clone(nm.Data()));
    ret->SetDirectory(0);
    for (Int_t i = 1; i <= ret->GetNbinsX(); i++) {
      Double_t x = ret->GetXaxis()->GetBinCenter(i);
      if (aSide && x < 0) {
	ret->SetBinContent(i,0);
	ret->SetBinError(i,0);
	continue;
      }
      if (!aSide && x >= 0) {
	ret->SetBinContent(i,0);
	ret->SetBinError(i,0);
      }
      ret->SetBinError(i, ret->GetBinContent(i)*factor);
    }
    return ret;
  }
      
				 
  /** 
   * Calculate and return ratio of a histogram to a graph
   * 
   * @param h Histogram
   * @param g Graph
   * @param l Lower error band on g 
   * @param u Upper error band on g
   * 
   * @return Ratio
   */
  static TH1* HOverG(const TH1* h, 
		     const TGraph* g, 
		     const TGraph* l, 
		     const TGraph* u,
		     Bool_t err=false) 
  {
    TH1* ret = static_cast<TH1*>(h->Clone());
    
    for (Int_t i = 1; i <= ret->GetNbinsX(); i++) { 
      Double_t c = h->GetBinContent(i);
      if (c < 1e-6) continue;
      
      Double_t e  = ret->GetBinError(i);
      Double_t x  = ret->GetXaxis()->GetBinCenter(i);;
      Double_t y  = g->Eval(x);
      Double_t yl = l->Eval(x);
      Double_t yu = u->Eval(x);
      Double_t ye = TMath::Max(y-yl,yu-yl);
      if (y < 1e-6) { 
	c = 0;
	e = 0;
      }
      else if (err) {
	c /= y;
	e /= y;
      }
      else { 
	c /= y;
	e = TMath::Sqrt(e*e*y*y+ye*ye*c*c);
      }
      ret->SetBinContent(i,c);
      ret->SetBinError(i,e);
    }
    return ret;
  }
  /** 
   * Calculate the ratio of two graphs 
   * 
   * @param num Numerator 
   * @param den Denominator 
   * 
   * @return newly allocated graph 
   */
  static TGraph* GOverG(const TGraph* num, 
			const TGraph* den,
			const TGraph* dlow, 
			const TGraph* dup,
			Double_t      dx=0)
  {
    TGraph* g = static_cast<TGraph*>(num->Clone());
    g->SetMarkerStyle(den->GetMarkerStyle());
    g->SetMarkerColor(num->GetMarkerColor());
    g->SetLineColor(num->GetLineColor());
    g->SetName(num->GetName());
    g->SetTitle(num->GetTitle());
    g->Set(0);
  
    const TGraphAsymmErrors* nga = 0;
    const TGraphErrors* nge = 0;
    TGraphAsymmErrors* ga = 0;
    TGraphErrors* ge = 0;
    if (g->IsA()->InheritsFrom(TGraphAsymmErrors::Class())) {
      nga = static_cast<const TGraphAsymmErrors*>(num);
      ga  = static_cast<TGraphAsymmErrors*>(g);
    }
    else if (g->IsA()->InheritsFrom(TGraphErrors::Class())) {
      nge = static_cast<const TGraphErrors*>(num);
      ge  = static_cast<TGraphErrors*>(g);
    }
    
    TGraph* nlow = 0;
    TGraph* nup = 0;
    // Info("", "Get other numerator error bands");
    if      (nga) ErrorGraphs(nga, nlow, nup);
    else if (nge) ErrorGraphs(nge, nlow, nup);

    Double_t dMin = den->GetX()[0];
    Double_t dMax = den->GetX()[den->GetN()-1];
    Int_t iP = 0;
    Int_t n = num->GetN();
    for (Int_t i = 0; i < n; i++) { 
      Double_t x   = num->GetX()[i];
      if (x < dMin || x > dMax) continue;

      Double_t yn  = num->GetY()[i];
      Double_t ynl = (nlow ? yn-nlow->Eval(x) : 0);
      Double_t ynu = (nup  ? nup->Eval(x)-yn  : 0);
      Double_t yd  = den->Eval(x-dx);
      Double_t ydl = (dlow ? yd-dlow->Eval(x-dx) : 0);
      Double_t ydu = (dup  ? dup->Eval(x-dx)-yd  : 0);
      if (yd < 1e-6 || yn < 1e-6) continue;

      g->SetPoint(iP, x, yn/yd);
      if (ga && nga) {
	// Printf("iP=%d yn=%f ynl=%f ynu=%f yd=%f ydl=%f ydu=%f",
	//        iP, yn, ynl, ynu, yd, ydl, ydu);
	Double_t el = TMath::Sqrt(ydl*ydl*yn*yn+ynl*ynl*yd*yd) / yd;
	Double_t eu = TMath::Sqrt(ydu*ydu*yn*yn+ynu*ynu*yd*yd) / yd;
	el = ynl / yd;
	eu = ynu / yd;
	ga->SetPointError(iP, nga->GetEXlow()[i], nga->GetEXhigh()[i], el, eu);
      }
      else if (ge && nge) {
	Double_t yne = TMath::Max(ynl,ynu);
	Double_t yde = TMath::Max(ydl,ydu);
	// Printf("iP=%d yn=%f yne=%f yd=%f yde=%f", iP, yn, yne, yd, yde);
	Double_t e   = TMath::Sqrt(yde*yde*yn*yn+yne*yne*yd*yd) / yd; 
	e = yne / yd;
	ge->SetPointError(iP, nge->GetEX()[i], e);
      }
      iP++;
    }
    delete nlow;
    delete nup;

    return g;
  }
  /** 
   * Get a result stack 
   * 
   * @param legend          Optional legend to fill
   * @param system          Collision system
   * @param sNN             Collision energy [GeV]
   * @param trigger         Trigger 
   * @param rebinned        If to get rebinned result
   * @param empirical       If to get empirical result 
   * 
   * @return Stack or null
   */
  static THStack* GetStack(TLegend*       legend, 
			   const TString& system, 
			   UShort_t       sNN, 
			   const TString& trigger, 
			   Bool_t         rebinned=true, 
			   Bool_t         empirical=true,
			   Int_t          marker=20)
  {
    TString sys(system); sys.ToUpper();
    TString trg(trigger); trg.ToUpper();
    if (trg.EqualTo("NSD")) trg = "V0AND";

    if (sys.EqualTo("PP")) { 
      sys = "pp";
      if      (trg.EqualTo("INEL"))                             trg = "INEL";
      else if (trg.EqualTo("INELGT0") || trg.EqualTo("INEL>0")) trg = "INELGt0";
      else if (trg.EqualTo("V0AND"))                            trg = "V0AND";
      else                                                      trg = ""; 
    }
    else if (sys.EqualTo("PPB")) { 
      sys = "pPb";
      if      (trg.EqualTo("CENTMB") || trg.EqualTo("V0AND")) trg = "CENTMB";
      else if (trg.EqualTo("CENTV0A"))                        trg = "CENTV0A";
      else if (trg.EqualTo("CENTV0X"))                        trg = "CENTV0A";
      else if (trg.EqualTo("CENTV0M"))                        trg = "CENTV0M";
      else if (trg.EqualTo("CENTZNA"))                        trg = "CENTZNA";
      else if (trg.EqualTo("CENTZNX"))                        trg = "CENTZNA";
      else                                                    trg = "";
    }
    else if (sys.EqualTo("PBP")) {
      sys = "Pbp";
      if      (trg.EqualTo("CENTMB") || trg.EqualTo("V0AND")) trg = "CENTMB";
      else if (trg.EqualTo("CENTV0C"))                        trg = "CENTV0C";
      else if (trg.EqualTo("CENTV0X"))                        trg = "CENTV0C";
      else if (trg.EqualTo("CENTV0M"))                        trg = "CENTV0M";
      else if (trg.EqualTo("CENTZNC"))                        trg = "CENTZNC";
      else if (trg.EqualTo("CENTZNX"))                        trg = "CENTZNC";
      else                                                    trg = "";
    }
    else if (sys.EqualTo("PBPB")) {
      sys = "PbPb";
      if (trg.EqualTo("CENT"))                                trg = "CENT";
      else if (trg.EqualTo("CENTV0M"))                        trg = "CENT";
      else                                                    trg = "";
    }
    else { 
      Error("", "Unknown system: %s", system.Data());
      return 0;
    }
    if (trg.IsNull()) { 
      Error("", "Unknown trigger %s for system %s", 
	    trigger.Data(), system.Data());
      return 0;
    }

    TString path(Form("%s/%s/%04d/%s/%s.C",
		      (empirical ? "nosec" : "normal"),
		      sys.Data(), sNN, trg.Data(), 
		      (rebinned ? "rebin" : "full")));
    if (gSystem->AccessPathName(path.Data())) { 
      Error("", "Script %s does not exist", path.Data());
      Error("", "sys=%s, sNN=%d, trg=%s, %d", 
	    sys.Data(), sNN, trg.Data(), rebinned);
      return 0;
    }

    THStack* stack = new THStack("stack", "");
    if (legend) 
      gROOT->Macro(Form("%s((THStack*)%p,(TLegend*)%p,%d)",
			path.Data(), stack, legend, marker));
    else 
      gROOT->Macro(Form("%s((THStack*)%p,0,%d)",
			path.Data(), stack, marker));

    gInterpreter->UnloadFile(path.Data());
    if (!stack->GetHists() || stack->GetHists()->GetEntries() <= 0) {
      Error("", "Got no histograms in stack from %s", 
	    path.Data());
      delete stack;
      return 0;
    }

    return stack;
  }
  /** 
   * Make a canvas 
   * 
   * @param system    System 
   * @param sNN       Energy
   * @param trigger   Trigger 
   * @param rebinned  Rebin 
   * @param empirical Empirical 
   * 
   * @return 
   */
  static TCanvas* MakeCanvas(const TString&  system, 
			     UShort_t        sNN, 
			     const TString&  trigger,
			     Bool_t          rebinned=true, 
			     Bool_t          empirical=true)
  {
    TString e = SNNString(system, sNN);
    TString t(Form("%s %s %s", system.Data(), e.Data(), trigger.Data()));
    TCanvas* canvas = new TCanvas(Form("%s_%04d_%s_%s_%s",
				       system.Data(), sNN, trigger.Data(), 
				       (rebinned ? "coarse" : "full"), 
				       (empirical ? "empirical" : "mc")), 
				  t.Data(),
				  1200, 1200);
    canvas->SetTopMargin(0.01);
    canvas->SetRightMargin(0.01);
    canvas->SetFillColor(0);
    canvas->SetFillStyle(0);
    canvas->SetBorderMode(0);
    canvas->SetBorderSize(0);
    
    return canvas;
  }
  /** 
   * Print canvas
   * 
   * @param canvas canvas
   * @param types  file types 
   */
  static void PrintCanvas(TCanvas* canvas, const TString& types)
  {
    if (types.IsNull()) return;

    canvas->Modified();
    canvas->Update();
    canvas->cd();

    TObjArray*  tokens = types.Tokenize(" ,");
    TObject*    token  = 0;
    TIter       next(tokens);
    while ((token = next())) {
      canvas->Print(Form("%s/%s.%s", PlotPrefix(), canvas->GetName(), 
			 token->GetName()));
    }
    tokens->Delete();
    delete tokens;
  }
  static Double_t* FixTriggerEff(const TString&  sys, 
				 UShort_t        sNN, 
				 TString&        trigger)
  {
    if (sNN != 8000 || !sys.EqualTo("pp", TString::kIgnoreCase))
      return 0;
    
    Double_t* ret = new Double_t[2];
    ret[1] = 0;
    if (trigger.EqualTo("INEL",TString::kIgnoreCase)) { 
      ret[0] = 0.85;
      return ret;
    }
    if (trigger.EqualTo("NSD", TString::kIgnoreCase) ||
	trigger.EqualTo("V0AND", TString::kIgnoreCase)) { 
      trigger = "NSD";
      ret[0]  = 0.93;
      return ret;
    }
    return 0;
  }
			    
  /** 
   * Draw a result
   * 
   * @param system          Collision system
   * @param sNN             Collision energy [GeV]
   * @param trigger         Trigger 
   * @param rebinned        If to get rebinned result
   * @param empirical       If to get empirical result 
   */
  static void Draw(const TString&  system, 
		   UShort_t        sNN, 
		   const TString&  trigger,
		   const Option_t* option="e3",
		   Bool_t          rebinned=true, 
		   Bool_t          empirical=true)
  {
    TCanvas* canvas = MakeCanvas(system,sNN,trigger,rebinned, empirical);

    TString     trg     = trigger; 
    Double_t*   effs    = FixTriggerEff(system,sNN,trg);
    const char* trigs[] = { trigger, 0 };
    const char* exps[]  = { "ALICE", "WIP", 0 };

    TLegend* l = 0;
    TObjArray u;
    TPair* dataOther = GetDataOther(l, u, system, sNN, trigs, exps, option,
				    rebinned, empirical, effs);
    if (effs) delete [] effs;
    if (!dataOther || !dataOther->Key()) {
      Error("", "No data found %s", canvas->GetTitle());
      return;
    }
    THStack*     data  = static_cast<THStack*>(dataOther->Key());
    TMultiGraph* other = static_cast<TMultiGraph*>(dataOther->Value());
    
    data->Draw("nostack");
    if (other) other->Draw("p");
    
    canvas->cd();
    data->SetTitle(Form("%s %s %s", system.Data(), data->GetTitle(), trigs[0]));
    data->SetMaximum(data->GetMaximum("nostack")*1.2);
    data->GetXaxis()->SetTitle("#eta");
    data->GetYaxis()->SetTitle("1/#it{N} d#it{N}_{ch}/d#it{#eta}");
    
    PrintCanvas(canvas, "pdf png");
  }
  /** 
   * Combine some graphs together 
   * 
   * @param mg 
   * 
   * @return 
   */
  static TGraphAsymmErrors* Combine(TMultiGraph* mg)
  {
    typedef std::map<double,std::vector<double> > ValueMap;
    ValueMap m;
    mg->GetListOfGraphs()->ls();
    TIter next(mg->GetListOfGraphs());
    TGraph* g = 0;
    while ((g = static_cast<TGraph*>(next()))) {
      TGraphErrors*      ge = 0;
      TGraphAsymmErrors* ga = 0;
      if (g->IsA()->InheritsFrom(TGraphErrors::Class()))
	ge = static_cast<TGraphErrors*>(g);
      else if (g->IsA()->InheritsFrom(TGraphAsymmErrors::Class()))
	ga = static_cast<TGraphAsymmErrors*>(g);

      Int_t n = g->GetN();
      for (Int_t i = 0; i < n; i++) {
	Double_t x   = g->GetX()[i];
	Double_t y   = g->GetY()[i];
	Double_t eyl = (ge ? ge->GetEY()[i] :
			ga ? ga->GetEYlow()[i] : 0);
	Double_t eyh = (ge ? ge->GetEY()[i] :
			ga ? ga->GetEYhigh()[i] : 0);
	// Double_t data[] = { y, eyl, eyh };

	ValueMap::iterator j = m.begin();
	for (; j != m.end(); ++j) {
	  ValueMap::iterator k = j;
	  k++;
	  if (TMath::Abs(x-j->first)<1e-1 && (k != m.end() && x < k->first)) {
	    break;
	  }
	}
	if (j != m.end()) {
	  j->second.push_back(y);
	  j->second.push_back(eyl);
	  j->second.push_back(eyh);
	}
	else {
	  m[x].push_back(y);
	  m[x].push_back(eyl);
	  m[x].push_back(eyh);
	}
      }
    }
    LinearSigmaCombiner sc;
    
    TGraphAsymmErrors* ret = new TGraphAsymmErrors(m.size());
    
    ValueMap::iterator j = m.begin();
    Int_t i = 0;
    for (; j != m.end(); ++j) {
      // printf("%7f (%d):", j->first, int(j->second.size()));
      if (j->second.size() == 3) {
	// Signular point
	ret->SetPoint(i, j->first, j->second[0]);
	ret->SetPointError(i, 0, 0, j->second[1], j->second[2]);
	i++;
	continue;
      }
      Combiner::List l;
      for (size_t k = 0; k < j->second.size(); k += 3) {
	l.Add(j->second[k],j->second[k+1],j->second[k+2]);
	// std::cout << (k/3)+1 << ": " << l.fData.back() << std::endl;
      }
      Combiner::const_iterator b = l.begin();
      Combiner::const_iterator e = l.end();
       
      Combiner::Final f = sc.Calculate(b,e);
      // std::cout << "Final: " << f << std::endl;
      ret->SetPoint(i, j->first, f.fX);
      ret->SetPointError(i, 0, 0, f.fEl, f.fEh);
      i++;
      // printf("\n");
    }
    return ret;
  }
  static Int_t FindCentBin(Color_t c)
  {
    // colors in decreasing cenrtality 
    Color_t cols[] = { 96, 86, 76, 66, 59, 55, 52, 0 };
    Color_t *pc = cols;
    Int_t idx = 0;
    while (*pc) {
      if (c == *pc) return idx;
      pc++;
      idx++;
    }
    return  -1;
  }
    
  /** 
   * Export a result
   * 
   * @param system          Collision system
   * @param sNN             Collision energy [GeV]
   * @param trigger         Trigger 
   * @param rebinned        If to get rebinned result
   * @param empirical       If to get empirical result 
   */
  static void Export(const TString& system, 
		     UShort_t       sNN, 
		     const TString& trigger,
		     Bool_t         rebinned=false, 
		     Bool_t         empirical=true)
  {
    TString e = SNNString(system, sNN);
    TString n(Form("%s_%04d_%s_%s_%s.root",
		   system.Data(), sNN, trigger.Data(), 
		   (rebinned ? "coarse" : "full"), 
		   (empirical ? "empirical" : "mc")));
    TFile* f = TFile::Open(n, "RECREATE");
    THStack* s = GetStack(0, system, sNN, trigger, rebinned, empirical);

    TString trg = trigger;
    Double_t* effs = FixTriggerEff(system, sNN, trg);
    Double_t  eff  = (effs? effs[0] : 1);
    TDirectory* fmd = f->mkdir("fmd");

    TMultiGraph* gf = new TMultiGraph("all","All Graphs");
    TList* listH = s->GetHists();
    TIter nextH(listH);
    TH1* oH = 0;
    TObjArray* byCent = (trigger.Contains("CENT") ? new TObjArray : 0);
    while ((oH = static_cast<TH1*>(nextH()))) {
      TString nm(oH->GetName());
      if (nm.Contains("syserror", TString::kIgnoreCase)) continue;
      // if (!nm.EndsWith("_all")) continue;

      oH->Scale(eff);
      TString sys(nm);
      sys.ReplaceAll("Forward", "SysError");
      TH1* sA = static_cast<TH1*>(listH->FindObject(Form("%s_a", sys.Data())));
      TH1* sC = static_cast<TH1*>(listH->FindObject(Form("%s_c", sys.Data())));
      if (!sA) {
	//Printf("W: Export - no A-side systematic error for %s",nm.Data());
	sA = ErrorHist(oH, true);
      }
      else
	sA->Scale(eff);
      if (!sC) {
	//Printf("W: Export - no C-side systematic error for %s",nm.Data());
	sC = ErrorHist(oH, false);
      }
      else
	sC->Scale(eff);
      
      TGraph* g = H2GA(oH, sA, sC);
      fmd->cd();
      g->Write();
      gf->Add(g, "p1");

      if (!byCent) continue;
      Int_t cBin = FindCentBin(g->GetMarkerColor());
      if (cBin < 0) {
	Printf("W: Export - didn't find centrality bin for %s (%d)",
	       oH->GetName(), g->GetMarkerColor());
	continue;
      }
      TMultiGraph* cGraph = static_cast<TMultiGraph*>(byCent->At(cBin));
      if (!cGraph) {
	nm.ReplaceAll("dndetForward_","");
	cGraph = new TMultiGraph;
	cGraph->SetName(nm); // Form("c%03d", cBin));
	byCent->AddAtAndExpand(cGraph, cBin);
      }
      cGraph->Add(g, "p1");
    }
    // s->Write();
    
    
    TMultiGraph* other = GetOther(system, sNN, trigger);
    if (other) { 
      TDirectory* spd = f->mkdir("spd");
      spd->cd();
      TIter nextG(other->GetListOfGraphs());
      TGraph* oG = 0;
      while ((oG = static_cast<TGraph*>(nextG()))) {
	spd->cd();
	oG->Write();
	gf->Add(oG, "p1");

	if (!byCent) continue;
	Int_t cBin = FindCentBin(oG->GetMarkerColor());
	if (cBin < 0) {
	  Printf("W: Export - didn't find centrality bin for %s (%d)",
		 oG->GetName(), oG->GetMarkerColor());
	  continue;
	}
	TMultiGraph* cGraph = static_cast<TMultiGraph*>(byCent->At(cBin));
	if (!cGraph) {
	  cGraph = new TMultiGraph;
	  byCent->AddAtAndExpand(cGraph, cBin);
	}
	cGraph->Add(oG, "p1");
      }
      // other->Write();
    }
    else 
      Warning("", "No other data for %s,%d,%s", 
	      system.Data(),sNN,trigger.Data());
    f->cd();
    gf->Write();
    gf->Draw("alp");

    if (!byCent) {
      // TMultiGraph* tmp = new TMultiGraph();
      // tmp->Add(gf);
      // if (other) tmp->Add(other);
      TGraphAsymmErrors* cf = Combine(gf);
      cf->SetName("combined");
      cf->SetTitle("Combined");
      cf->SetMarkerStyle(24);
      cf->SetMarkerColor(kGreen+2);
      gf->Add(cf, "p");
      cf->Write();
    } else {
      TMultiGraph* com = new TMultiGraph;
      com->SetName("combined");
      TIter next(byCent);
      TMultiGraph* mg = 0;
      while ((mg = static_cast<TMultiGraph*>(next()))) {
	mg->ls();
	TGraphAsymmErrors* cf = Combine(mg);
	cf->SetName(Form("combined_%s",mg->GetName()));
	cf->SetTitle("Combined");
	cf->SetMarkerStyle(24);
	cf->SetMarkerColor(static_cast<TGraph*>(mg->GetListOfGraphs()->At(0))
			   ->GetMarkerColor());
	gf->Add(cf, "p");
	com->Add(cf, "p");
	cf->Write();
      }
      byCent->Write("bycent", TObject::kSingleKey);
      com->Write();
    }

    
    f->Write();
    Info("Export", "Exported data to %s", f->GetName());
  }
  /** 
   * Add the systematics 
   * 
   * @param stack Stack
   * @param sys   Systematic error 
   */
  static void AddSystematics(THStack* stack, Double_t sys)
  {
    TIter next(stack->GetHists());
    TH1*  hist = 0;
    while ((hist = static_cast<TH1*>(next()))) {
      Int_t nBins = hist->GetNbinsX();
      for (Int_t bin = 1; bin <= nBins; bin++) {
	Double_t c = hist->GetBinContent(bin);
	if (c < 1e-6) continue;
	hist->SetBinError(bin, sys*c);
      }
    }
  }
  /**
   * Symmetrice between two stacks 
   * 
   * @param s1 Stack 
   * @param s2 Stack 
   * 
   * @return Return new stack with symmetriced histograms
   */
  static THStack* Symmetrice(THStack* s1, THStack* s2) 
  {

    THStack* res = new THStack("Result", "Result");

    Int_t nHist = s1->GetHists()->GetEntries();
    for (Int_t i = 0; i < nHist; i++) { 
      TH1*     h1 = static_cast<TH1*>(s1->GetHists()->At(i));
      TH1*     h2 = static_cast<TH1*>(s2->GetHists()->At(i));
      Int_t    nb = h1->GetNbinsX();
      Double_t x1 = h1->GetXaxis()->GetXmin();
      Double_t x2 = h1->GetXaxis()->GetXmax();
      Double_t dx = (x2-x1)/nb;
      
      Int_t    hb = (x2+x2) / dx;
      TH1*     h  = new TH1D(h1->GetName(), h1->GetTitle(), hb, -x2, x2);
      h->SetMarkerStyle(h1->GetMarkerStyle());
      h->SetMarkerColor(h1->GetMarkerColor());
      h->SetMarkerSize(h1->GetMarkerSize());
      h->SetFillStyle(kGray); // h1->GetFillStyle());
      h->SetFillColor(1001); // h1->GetFillColor());
      h->SetLineColor(h1->GetLineColor());

      for (Int_t j = 1; j <= hb; j++) { 
	Double_t x   =  h->GetXaxis()->GetBinCenter(j);
	Int_t    b1  =  h1->GetXaxis()->FindBin(x);
	Int_t    b2  =  h2->GetXaxis()->FindBin(-x);
	Double_t xx1 =  h1->GetXaxis()->GetBinCenter(b1);
	Double_t xx2 = -h2->GetXaxis()->GetBinCenter(b2);
	Double_t c1  =  h1->GetBinContent(b1);
	Double_t c2  =  h2->GetBinContent(b2);
	Double_t e1  =  0.08*c1; // h1->GetBinError(b1);
	Double_t e2  =  0.08*c2; // h2->GetBinError(b2);
	
	if (TMath::Abs(xx1+xx2) < 1e-6) { 
	  Double_t cc = .5 * (c1+c2);
	  Double_t ee = TMath::Sqrt(e1*e1+e2*e2);
	  // Info("", "@ %f, fill with mean %f+/-%f", x, cc, ee);
	  h->SetBinContent(j, cc);
	  h->SetBinError(j, ee);
	}
	else if (c1 > 0) {
	  // Info("", "@ %f fill with  %f+/-%f (1)", x, c1, e1);
	  h->SetBinContent(j, c1);
	  h->SetBinError(j, e1);
	}
	else if (c2 > 0) {
	  // Info("", "@ %f fill with  %f+/-%f (2)", x, c2, e2);
	  h->SetBinContent(j, c2);
	  h->SetBinError(j, e2);
	}
      } // for j
      res->Add(h);
    } // for i 

    return res;
  }
  /** 
   * Symmetrice histograms in a stack 
   * 
   * @param s Stack 
   * 
   * @return Return new stack with symmetriced histograms
   */
  static THStack* Symmetrice(THStack* s) 
  {
    THStack* res = new THStack("res", s->GetTitle());

    Int_t nHist = s->GetHists()->GetEntries();
    for (Int_t i = 0; i < nHist; i++) { 
      TH1*     h1 = static_cast<TH1*>(s->GetHists()->At(i));
      Int_t    nb = h1->GetNbinsX();
      Double_t x1 = h1->GetXaxis()->GetXmin();
      Double_t x2 = h1->GetXaxis()->GetXmax();
      Double_t dx = (x2-x1)/nb;
      
      Int_t    hb = (x2+x2) / dx;
      TH1*     h  = new TH1D(h1->GetName(), h1->GetTitle(), hb, -x2, x2);
      h->SetMarkerStyle(h1->GetMarkerStyle());
      h->SetMarkerColor(h1->GetMarkerColor());
      h->SetMarkerSize(h1->GetMarkerSize());
      h->SetFillStyle(kGray); // h1->GetFillStyle());
      h->SetFillColor(1001); // h1->GetFillColor());
      h->SetLineColor(h1->GetLineColor());

      for (Int_t j = 1; j <= hb; j++) { 
	Double_t x   =  h->GetXaxis()->GetBinCenter(j);
	Int_t    b1  =  h1->GetXaxis()->FindBin(x);
	Int_t    b2  =  h1->GetXaxis()->FindBin(-x);
	Double_t c1  =  h1->GetBinContent(b1);
	Double_t c2  =  h1->GetBinContent(b2);
	Double_t e1  =  0.08*c1; // h1->GetBinError(b1);
	Double_t e2  =  0.08*c2; // h2->GetBinError(b2);
	
	if (c1 > 0 && c2 > 0) { 
	  Double_t cc = .5 * (c1+c2);
	  Double_t ee = TMath::Sqrt(e1*e1+e2*e2);
	  // Info("", "@ %f, fill with mean %f+/-%f", x, cc, ee);
	  h->SetBinContent(j, cc);
	  h->SetBinError(j, ee);
	}
	else if (c1 > 0) {
	  // Info("", "@ %f fill with  %f+/-%f (1)", x, c1, e1);
	  h->SetBinContent(j, c1);
	  h->SetBinError(j, e1);
	}
	else if (c2 > 0) {
	  // Info("", "@ %f fill with  %f+/-%f (2)", x, c2, e2);
	  h->SetBinContent(j, c2);
	  h->SetBinError(j, e2);
	}
      } // for j
      res->Add(h);
    } // for i  
    return res;
  }
  /** 
   * Draw pPb symmetriced with Pbp 
   * 
   * @param trigger 
   * @param option 
   * @param rebinned 
   * @param empirical 
   */
  static void pPbSym(const TString& trigger,
		     const Option_t* option="e3",
		     Bool_t         rebinned=true, 
		     Bool_t         empirical=true) 
  {
    TString  system("pPb");
    TString  system2("Pbp");
    UShort_t sNN  = 5023;
    TString  e    = SNNString(system, sNN);
    TCanvas* canvas = MakeCanvas(system, sNN, trigger);

    TString t1(trigger); t1.ReplaceAll("X", "A");
    TString t2(trigger); t2.ReplaceAll("X", "C");
    THStack* s1 = GetStack(0, system,  sNN, t1, rebinned, empirical);
    THStack* s2 = GetStack(0, system2, sNN, t2, rebinned, empirical);
  
    if (!s1 || !s2) {
      delete canvas;
      return; 
    }
    s1->SetTitle(canvas->GetTitle());
    s2->SetTitle(canvas->GetTitle());

    THStack* res = Symmetrice(s1, s2); 
    res->Draw(Form("nostack %s", option));

    TMultiGraph* other = GetOther(system, sNN, t1);
    if (other) { 
      TIter next(other->GetListOfGraphs());
      TGraph* g = 0;
      while ((g = static_cast<TGraph*>(next()))) {
	g->SetFillColor(kGray);
	g->SetFillStyle(3001);
	g->Draw(Form("p same %s", option));
      }
    }
    else 
      Warning("", "No other data for %s,%d,%s", 
	      system.Data(),sNN,trigger.Data());


    PrintCanvas(canvas, "pdf png");
  }
  /** 
   * Draw a result scaled by pp result 
   * 
   * @param system          Collision system
   * @param sNN             Collision energy [GeV]
   * @param ppsNN           pp Collision energy [GeV]
   * @param trigger         Trigger 
   * @param ppTrigger       pp Trigger 
   * @param etaShift        Possible eta shift 
   * @param rebinned        If to get rebinned result
   * @param empirical       If to get empirical result 
   */
  static void ScaleBypp(const TString& system, 
			UShort_t       sNN, 
			UShort_t       ppsNN, 
			const TString& trigger, 
			const TString& ppTrigger, 
			Double_t       etaShift=0,
			Bool_t         rebinned=true, 
			Bool_t         empirical=true,
			Bool_t         write=false)
  {
    THStack* pp = GetStack(0, "pp", ppsNN, ppTrigger, rebinned, empirical);
    if (!pp) return;

    THStack* tgt = GetStack(0, system, sNN, trigger, rebinned, empirical);
    if (!tgt) return;

    TH1* ref = static_cast<TH1*>(pp->GetHists()
				 ->FindObject("dndetaForward_all"));
    TGraphErrors* g = H2G(ref, etaShift);
    if (!g) return;

    TGraph* low = 0;
    TGraph* up  = 0;
    // Info("", "Get denominator error bands");
    ErrorGraphs(g, low, up);

    TMultiGraph* outRatio = (write ? new TMultiGraph("fwd","") : 0);
    TMultiGraph* outPP    = (write ? new TMultiGraph("pp","") : 0);
    THStack* ratios = new THStack("ratios", tgt->GetTitle());
    TIter next(tgt->GetHists());
    TH1* h = 0;

    while ((h = static_cast<TH1*>(next()))) {
      TH1* r = HOverG(h, g, low, up); 
      ratios->Add(r);
      if (outRatio) {
	TGraph* rg = H2G(r, 0);
	rg->SetLineColor(h->GetLineColor());
	rg->SetLineStyle(h->GetLineStyle());
	rg->SetFillColor(kBlue-10);
	rg->SetFillStyle(0);
	outRatio->Add(rg, "p");
      }
    }
    if (outPP) outPP->Add(g, "p");
    delete low;
    delete up;

    TString eTgt = SNNString(system, sNN);
    TString epp  = SNNString("pp", ppsNN);
  
    TString tTgt(Form("%s %s %s", system.Data(), eTgt.Data(), trigger.Data()));
    TString tpp(Form("pp %s %s", epp.Data(), ppTrigger.Data()));
    TString t(Form("%s scaled by %s", tTgt.Data(), tpp.Data()));
    ratios->SetTitle(t);

    TCanvas* canvas = MakeCanvas(system, sNN, trigger);
    canvas->SetName(Form("%s_%04d_%s_pp_%04d_%s",
			 system.Data(), sNN, trigger.Data(), 
			 ppsNN, ppTrigger.Data()));
    canvas->SetTitle(t);

    ratios->Draw("nostack");
    ratios->GetXaxis()->SetTitle(Form("#eta%s", 
				   (TMath::Abs(etaShift) < 1e-6 
				    ? "" : "-#eta_{MC}")));
    const char* dNdeta = "1/#it{N} d#it{N}_{ch}/d#it{#eta}";
    ratios->GetYaxis()->SetTitle(Form("#frac{%s|_{%s}}{%s|_{%s}}",
				   dNdeta, tTgt.Data(), dNdeta, tpp.Data()));
    ratios->GetYaxis()->SetTitleOffset(1.5);
    
    TMultiGraph* other   = GetOther(system, sNN, trigger);
    TMultiGraph* ppOther = GetOther("pp", ppsNN, ppTrigger);
    TMultiGraph* outOther = (write ? new TMultiGraph("cen", "") : 0);
    if (other && ppOther && ppOther->GetListOfGraphs()) { 
      // Info("", "Got %d others", other->GetListOfGraphs()->GetEntries());
      TGraph* otherRef = 
	static_cast<TGraph*>(ppOther->GetListOfGraphs()->At(0));
      TGraph* dlow = 0;
      TGraph* dup  = 0;
      // Info("", "Get other denominator error bands");
      if (otherRef->IsA()->InheritsFrom(TGraphAsymmErrors::Class())) 
	ErrorGraphs(static_cast<TGraphAsymmErrors*>(otherRef), dlow, dup);
      else if (otherRef->IsA()->InheritsFrom(TGraphErrors::Class())) 
	ErrorGraphs(static_cast<TGraphErrors*>(otherRef), dlow, dup);
      else 
	ErrorGraphs(otherRef, dlow, dup);

      if (otherRef) {
	outPP->Add(otherRef, "p");
	TIter nextG(other->GetListOfGraphs());
	TGraph* g1 = 0;
	while ((g1 = static_cast<TGraph*>(nextG()))) {
	  TGraph* rg = GOverG(g1, otherRef, dlow, dup, etaShift);
	  if (rg) rg->Draw("p same");
	  outOther->Add(rg, "p");
	}
      }
      delete dlow;
      delete dup;
    }
    else 
      Warning("", "No other data for %s,%d,%s", 
	      system.Data(),sNN,trigger.Data());

    PrintCanvas(canvas,"pdf png");

    if (!write) return;
    TFile* out = TFile::Open(Form("%s/%s.root", PlotPrefix(),
				  canvas->GetName()), "RECREATE");
    outRatio->Write();
    outOther->Write();
    outPP->Write();
    out->Write();
  }
  /** 
   * Symmetrice distribution and scale by pp
   * 
   * @param ppsNN     PP center of mass energy 
   * @param trigger   Numerator trigger 
   * @param ppTrigger PP trigger 
   * @param option    Drawing option 
   * @param rebinned  If true, use rebinned data 
   * @param empirical If true, use data corrected by empirical corr. 
   * @param write     If true, also write to ROOT file 
   */
  static void SymScaleBypp(UShort_t        ppsNN, 
			   const TString&  trigger, 
			   const TString&  ppTrigger, 
			   const Option_t* option="e2",
			   Bool_t          rebinned=true, 
			   Bool_t          empirical=true,
			   Bool_t          write=false)
  {
    TString  system("pPb");
    TString  system2("Pbp");
    UShort_t sNN  = 5023;
    TString  e    = SNNString(system, sNN);

    TString t1(trigger); t1.ReplaceAll("X", "A");
    TString t2(trigger); t2.ReplaceAll("X", "C");
    THStack* s1 = GetStack(0, system,  sNN, t1, rebinned, empirical);
    THStack* s2 = GetStack(0, system2, sNN, t2, rebinned, empirical);
  
    if (!s1 || !s2) return;
    THStack* tgt = Symmetrice(s1, s2); 
    THStack* pp  = Symmetrice(GetStack(0, "pp", ppsNN, ppTrigger, 
				       rebinned, empirical)); 
    if (!pp) return;

    TH1* ref = static_cast<TH1*>(pp->GetHists()
				 ->FindObject("dndetaForward_all"));
    TGraphErrors* g = H2G(ref, 0);
    if (!g) return;

    TGraph* low = 0;
    TGraph* up  = 0;
    // Info("", "Get denominator error bands");
    ErrorGraphs(g, low, up);

    THStack*     ratios    = new THStack("ratios", tgt->GetTitle());
    TMultiGraph* outRatio  = (write ? new TMultiGraph("fwd","") : 0);
    TMultiGraph* outPP     = (write ? new TMultiGraph("pp","") : 0);
    TMultiGraph* outdNdeta = (write ? new TMultiGraph("dNdeta","") : 0);
    TIter next(tgt->GetHists());
    TH1* h = 0;
    while ((h = static_cast<TH1*>(next()))) {
      TH1* r = HOverG(h, g, low, up, true);
      r->SetFillStyle(0);
      r->SetFillColor(kBlue-10);
      Info("", "Adding %s to stack", r->GetName());
      ratios->Add(r, "p");
      if (outRatio) {
	TGraph* rg = H2G(r, 0);
	rg->SetLineColor(h->GetLineColor());
	rg->SetLineStyle(h->GetLineStyle());
	rg->SetFillColor(kBlue-10);
	rg->SetFillStyle(0);
	outRatio->Add(rg, "p");
      }
      if (outdNdeta) {
	TGraph* hg = H2G(h, 0);
	hg->SetLineColor(h->GetMarkerColor());
	outdNdeta->Add(hg, "p");
      }
    }
    if (outPP) {
      outPP->Add(g, "p1");
    }
    delete low;
    delete up;

    TString eTgt = SNNString(system, sNN);
    TString epp  = SNNString("pp", ppsNN);  
    TString tTgt(Form("%s %s %s", system.Data(), eTgt.Data(), trigger.Data()));
    TString tpp(Form("pp %s %s", epp.Data(), ppTrigger.Data()));
    TString t(Form("%s scaled by %s", tTgt.Data(), tpp.Data()));
    ratios->SetTitle(t);

    TCanvas* canvas = MakeCanvas(system, sNN, trigger);
    canvas->SetName(Form("%s_%04d_%s_pp_%04d_%s_sym",
			 system.Data(), sNN, trigger.Data(), 
			 ppsNN, ppTrigger.Data()));
    canvas->SetTitle(t.Data());

    ratios->Draw(Form("nostack %s", option));
    ratios->GetXaxis()->SetTitle("#eta");
    const char* dNdeta = "1/#it{N} d#it{N}_{ch}/d#it{#eta}";
    ratios->GetYaxis()->SetTitle(Form("#frac{%s|_{%s}}{%s|_{%s}}",
				   dNdeta, tTgt.Data(), dNdeta, tpp.Data()));
    ratios->GetYaxis()->SetTitleOffset(1.5);
    if (outRatio && outRatio->GetHistogram()) {
      outRatio->GetHistogram()->SetXTitle(ratios->GetXaxis()->GetTitle());
      outRatio->GetHistogram()->SetYTitle(ratios->GetYaxis()->GetTitle());
    }
    
    TMultiGraph* other    = GetOther(system, sNN, t1);
    TMultiGraph* ppOther  = GetOther("pp", ppsNN, ppTrigger);
    TMultiGraph* outOther = (write ? new TMultiGraph("cen", "") : 0);
    if (other && ppOther && ppOther->GetListOfGraphs()) { 
      // Info("", "Got %d others", other->GetListOfGraphs()->GetEntries());
      TGraph* otherRef = 
	static_cast<TGraph*>(ppOther->GetListOfGraphs()->At(0));
      TGraph* dlow = 0;
      TGraph* dup  = 0;
      if (outPP) {
	outPP->Add(otherRef, "p1");
      }
      // Info("", "Get other denominator error bands");
      if (otherRef->IsA()->InheritsFrom(TGraphAsymmErrors::Class())) 
	ErrorGraphs(static_cast<TGraphAsymmErrors*>(otherRef), dlow, dup);
      else if (otherRef->IsA()->InheritsFrom(TGraphErrors::Class())) 
	ErrorGraphs(static_cast<TGraphErrors*>(otherRef), dlow, dup);
      else 
	ErrorGraphs(otherRef, dlow, dup);

      if (otherRef) {
	TIter nextG(other->GetListOfGraphs());
	TGraph* g1 = 0;
	while ((g1 = static_cast<TGraph*>(nextG()))) {
	  TGraph* rg = GOverG(g1, otherRef, dlow, dup, 0);
	  rg->SetLineColor(g1->GetLineColor());
	  rg->SetLineStyle(g1->GetLineStyle());
	  rg->SetFillColor(kBlue-10);
	  rg->SetFillStyle(0);
	  if (rg) rg->Draw(Form("p same %s", option));
	  if (outOther)  outOther->Add(rg, "p1");
	  if (outdNdeta) {
	    g1->SetLineColor(g1->GetMarkerColor());
	    outdNdeta->Add(g1, "p");
	  }
	}
      }
      delete dlow;
      delete dup;
    }
    else 
      Warning("", "No other data for %s,%d,%s", 
	      system.Data(),sNN,trigger.Data());
    if (outOther && outOther->GetHistogram()) {
      outOther->GetHistogram()->SetXTitle(ratios->GetXaxis()->GetTitle());
      outOther->GetHistogram()->SetYTitle(ratios->GetYaxis()->GetTitle());
    }

    PrintCanvas(canvas, "png pdf");

    if (!write) return;
    TFile* out = TFile::Open(Form("%s/%s.root",PlotPrefix(),
				  canvas->GetName()), "RECREATE");
    out->cd();
    outRatio->Write();
    outOther->Write();
    outdNdeta->Write();
    outPP->Write();
    out->Write();
    
  }
  /** 
   * Append to ouput name 
   * 
   * @param out   Output name 
   * @param what  What to append 
   */
  static void Add2Out(TString& out, const char* what)
  {
    if (!out.IsNull()) out.Append("_");
    out.Append(what);
  }
  /** 
   * Append to ouput name 
   * 
   * @param out   Output name 
   * @param sNN   What to append 
   */
  static void Add2Out(TString& out, UShort_t sNN)
  {
    if (!out.IsNull()) out.Append("_");
    out.Append(Form("%04d", sNN));
  }
  /** 
   * Draw all plots 
   * 
   * @param syss   Systems to plot 
   * @param sNNs   Energies to plot 
   * @param trgS   Triggers to plot 
   * @param flags  Some flags 
   */
  static void DrawAll(const char** syss,
		      UShort_t*    sNNs,
		      const char** trgS,
		      UShort_t     flags=0x0)
  {
    Int_t nSys = 0;
    Int_t nSNN = 0;
    Int_t nTrg = 0;
    const char** pSys   = syss;
    UShort_t*    pSNN   = sNNs;
    const char** pTrg   = trgS;
    const char*  exps[] = { "ALICE", "CMS", "WIP", 0 };
    TString out;
    while (*pSys) { nSys++; Add2Out(out, *pSys); pSys++; }
    while (*pSNN) { nSNN++; Add2Out(out, *pSNN); pSNN++; }
    while (*pTrg) { nTrg++; Add2Out(out, *pTrg); pTrg++; }
    if (nSys <= 0) {
      ::Warning("DrawAll", "No systems specified");
      return;
    }
    if (nSNN <= 0) {
      ::Warning("DrawAll", "No energies specified");
      return;
    }
    if (nTrg <= 0) {
      ::Warning("DrawAll", "No triggers specified");
      return;
    }

    if (nSys != 1 && nSNN != 1 && nTrg != 1) {
      ::Warning("DrawAll", "At least one of sys (%d), sNN (%d), "
		"or trg (%d) must be singular", nSys, nSNN, nTrg);
      return;
    }
    Int_t nHoriz = 0;
    Int_t nVert  = 1;
    if (nSys == 1) { nHoriz = nSNN; }
    if (nSNN == 1) { nHoriz = nSys; }

    const Color_t kAliceBlue   = AliceBlue();
    // const Color_t kAliceRed    = AliceRed();

    Bool_t collapse = (flags & 0x1);
    TCanvas* c = new TCanvas("all", out, nHoriz*500, nVert*500);
    c->SetTopMargin(0.01);
    c->SetRightMargin(0.01);
    if (collapse) 
      c->Divide(nHoriz, nVert, 0, 0);
    else
      c->Divide(nHoriz, nVert);

    // --- Loop over specs -------------------------------------------
    TList          stacks;
    Double_t       max = 0;
    TObjArray      unique;
    const Double_t pp8000e[] = { 0.85, 0.93, 1, 0 };
    Double_t       lw        = (nSNN == 1 ? .6  : .83);
    Double_t       ly1       = (nSNN == 1 ? .7  : .15);
    Double_t       ly2       = (nSNN == 1 ? .95 : .35);
    
    Int_t iPad = 1;
    pSys  = syss;
    while (*pSys) {
      pSNN = sNNs;
      while (*pSNN) {
	TVirtualPad* q = c->cd(iPad++);
	TLegend* leg = 0;
	if (pSys == syss && pSNN == sNNs) {
	  Double_t lx1 = .15;
	  leg =  new TLegend(lx1,
			     ly1 /*0.24*/, lx1+lw, ly2); // , stit);
	  leg->SetNColumns(nTrg == 1 ? 2 : nTrg);
	  leg->SetFillColor(0);
	  leg->SetFillStyle(0);
	  leg->SetBorderSize(0);
	  leg->SetTextColor(kAliceBlue);
	  leg->SetTextFont(42);
	}
	const Double_t* effs = ((((*pSys)[0] == (*pSys)[1]) && *pSNN == 8000)
				? pp8000e : 0);
	  
	TPair* p = GetDataOther(leg, unique, *pSys, *pSNN, trgS, exps,
				"e5", true, true,  effs);
	if (!p) {
	  Warning("", "No data for pad %d", iPad-1);
	  pSNN++;
	  continue;
	}
	if (!p->Key()) {
	  Warning("", "No data for pad %d", iPad-1);
	  pSNN++;
	  continue;
	}
	THStack* stack = static_cast<THStack*>(p->Key());
	stacks.Add(stack);
	if (nSNN == 1) stack->SetTitle(*pSys);
	if (nTrg == 1) stack->SetTitle(Form("%s - %s", *pSys, trgS[0]));
	Double_t lFac = (TString("Pbp").EqualTo(*pSys) ? 1.5 : 1.2);
	Double_t lMax = lFac*stack->GetMaximum("nostack");
	stack->SetMaximum(lMax);
	stack->SetMinimum(.3);
	stack->Draw("nostack");
	stack->GetHistogram()->SetXTitle("#eta");
	stack->GetHistogram()->SetYTitle("1/#it{N} d#it{N}_{ch}/d#it{#eta}");
	max = TMath::Max(lMax, max);
	// Info("", "Maximum: %f (%f)", lMax, max);
	if (p->Value()) {
	  p->Value()->Draw("p");
	}
	if (leg) leg->Draw();

	q->Modified();
	q->Update();
	q->cd();
	
	pSNN++;
      }
      pSys++;
    }

    Int_t nPad = nHoriz * nVert;
    if (collapse) {
      // ::Info("", "Maximum is %f", max);
      TIter iStack(&stacks);
      THStack* stack = 0;
      while ((stack = static_cast<THStack*>(iStack())))
	stack->SetMaximum(max);
      for (iPad = 1; iPad <= nPad; iPad++) {
	TVirtualPad* q = c->cd(iPad);
	q->Modified();
	q->Update();
	q->cd();
      }
    }
    
    // --- Make legend of unique names -------------------------------
    TParameter<int>* u =
      new TParameter<int>("PWG-LF/GEO - work in progress", 20);
    u->SetUniqueID(nSNN == 1 ? kBlack : kRed+2);
    unique.Add(u);

    Double_t lx1  = 1-lw-.02;
    ly2           = TMath::Max(ly2,0.5);
    if (nPad == 1) {
      ly2 = ly1;
      ly1 = 0.5;
    }
      
    TLegend* uleg = MakeUniqueLegend(lx1,ly1,lx1+lw,ly2,unique, nSNN);
    TVirtualPad* q = c->cd(nPad);
    q->SetRightMargin(0.01);
    uleg->Draw();

    PrintCanvas(c, "pdf png");
    c->SaveAs(Form("%s/%s.root", PlotPrefix(), out.Data()));
    
  }
  static TLatex* MakeTitle(Double_t x, Double_t y, 
			   const TString& system, 
			   UShort_t       sNN, 
			   const TString& trigger)
  {
    TString e = SNNString(system, sNN);
    TString t = trigger;
    if (t.Contains("CENT")) { 
      t.ReplaceAll("CENT", "by centrality (");
      t.Append(")");
    }
    TString txt(Form("%s @ %s %s", 
		     system.Data(), e.Data(), t.Data()));
    TLatex* ltx = new TLatex(x, y, txt);
    ltx->SetTextFont(42);
    ltx->SetTextColor(AliceBlue());
    ltx->SetNDC();
    
    ltx->Draw();
    
    return ltx;
  }
  static TLegend* MakeUniqueLegend(Double_t   x1, 
				   Double_t   y1, 
				   Double_t   x2, 
				   Double_t   y2, 
				   TObjArray& unique,
				   Int_t      nSNN)
  {
    const Color_t kAliceBlue   = AliceBlue();

    TLegend* uleg = new TLegend(x1, y1, x2, y2 /*0.23*/);
    uleg->SetNColumns(1 /*2*/);
    uleg->SetFillColor(0);
    uleg->SetFillStyle(0);
    uleg->SetBorderSize(0);
    uleg->SetTextColor(kAliceBlue);
    uleg->SetTextFont(42);

    return MakeUniqueLegend(uleg, unique, nSNN);
  }
  static TLegend* MakeUniqueLegend(TLegend* uleg, TObjArray& unique, Int_t nSNN)
  {
    TIter nextU(&unique);
    TLegendEntry* e = 0;
    TParameter<int>* u = 0;
    while ((u = static_cast<TParameter<int>*>(nextU()))) {
      e = uleg->AddEntry("dummy", u->GetName(), nSNN == 1 ? "p" : "f");
      e->SetMarkerStyle(u->GetVal());
      e->SetMarkerSize(1.7);
      e->SetMarkerColor(u->GetUniqueID());    
      e->SetFillColor(u->GetUniqueID());    
      e->SetFillStyle(1001);    
      e->SetLineColor(kBlack); // u->GetUniqueID());    
    }
    return uleg;
  }			    
  /** 
   * Get the marker size based on marker type 
   * 
   * @param style Style of marker 
   * 
   * @return return relative size 
   */
  static Float_t MarkerSize(Int_t style)
  {
    switch (style) {
    case 20: return 1.6;
    case 21: return 1.5;
    case 22: return 1.7;
    }
    ::Warning("", "Marker style %d maps to default", style);
    return 1.6;
  }
  static const char* CentLimitName(Bool_t isMult, Float_t v)
  {
    if (isMult)                    return Form("%3d", Int_t(v));
    if ((Int_t(v*100) % 100) == 0) return Form("%3d%%", Int_t(v)); 
    if ((Int_t(v*100)  % 10) == 0) return Form("%5.1f%%", v);
    return Form("%6.2f%%", v);
  }
  /** 
   * REturn a pair of a stack (data) and multi-graph (other)
   * 
   * @param leg            Possible legend to fill 
   * @param unique         For re-entrant calling to set unique names
   * @param system         Collision system
   * @param sNN            Collision energy 
   * @param trigs          Triggers to draw 
   * @param exps           Other to get the data for 
   * @param errOpt         Error option 
   * @param rebinned       Rebinned or full
   * @param empirical      Empirical correction or not 
   * @param effs           Efficiencies 
   * 
   * @return Pair or null
   */
  static TPair* GetDataOther(TLegend*        leg,
			     TObjArray&      unique,
			     const TString&  system,
			     UShort_t        sNN,
			     const char**    trigs,
			     const char**    exps,
			     Option_t*       errOpt="e5",
			     Bool_t          rebinned=false,
			     Bool_t          empirical=true,
			     const Double_t* effs=0)
  {

    TPair*          ret    = 0;
    TString         stit   = SNNString(system, sNN);
    THStack*        stack  = new THStack("stack", stit.Data());    
    TList           points;
    Int_t           marker = 20;
    Bool_t          seSeen = false;
    Int_t           fill   = 0;
    Int_t           style  = 3001;
    TLegendEntry*   e      = 0;
    const char**    ptrig  = trigs;
    const Double_t* peff   = effs;
    TString         allT   = "";
    Int_t           verbose= 0;
    // if (leg) leg->SetHeader(stit);
    points.SetOwner(false);
    while ((*ptrig)) {
      if (!allT.IsNull()) allT = ",";
      allT.Append(*ptrig);

      TString trg(*ptrig);
      Bool_t   isCent = trg.BeginsWith("CENT");
      THStack* tmp    = GetStack(0, system, sNN, *ptrig,
				 rebinned, empirical, isCent ? 20 : marker);
      if (!tmp) {
	Warning("", "No data for %s %d %s", system.Data(), sNN, *ptrig);
	ptrig++;
	marker++;
	if (peff) peff++;
	continue;
      }
      Double_t eff = (sNN == 8000 && peff ? *peff : 1);
      TList*   lst = tmp->GetHists();
      TIter    next(lst);
      TObject* o = 0;
      while ((o = next())) { 
	TString n(o->GetName());			
	if (n.Contains("mirror")) continue;
	TH1* h = static_cast<TH1*>(o->Clone());
	if (isCent) {
	  Color_t cOld = h->GetMarkerColor();
	  Color_t cNew = Brighten(cOld);
	  h->SetMarkerColor(Brighten(h->GetMarkerColor()));
	  h->SetFillColor(Brighten(h->GetFillColor()));
	  h->SetLineColor(Brighten(h->GetLineColor()));
	  Printf("%s: old=%d new=%d real=%d",n.Data(),
		 cOld,cNew,h->GetMarkerColor());
	}	
	h->SetDirectory(0);
	h->Scale(eff);
	if (n.Contains("SysError")) {
	  seSeen = true;
	  fill     = h->GetFillColor();
	  style    = h->GetFillStyle();
	  TH1*  hc = static_cast<TH1*>(h->Clone(Form("%s_syserror_C",*ptrig)));
	  TH1*  ha = h;
	  ha->SetName(Form("%s_syserror_A", *ptrig));
	  for (Int_t j = 1; j < h->GetNbinsX(); j++) { 
	    Double_t x = h->GetXaxis()->GetBinCenter(j);
	    TH1* z = (x < 0 ? ha : hc);
	    z->SetBinContent(j, 0);
	    z->SetBinError(j, 0);
	  }
	  
	  stack->Add(ha, errOpt);
	  stack->Add(hc, errOpt);
	  continue;
	} // SysError
	if (!isCent) h->SetName(*ptrig);
	h->SetMarkerSize(MarkerSize(h->GetMarkerStyle()));
	points.Add(h);

	if (leg) {
	  TString nm(*ptrig);
	  if (isCent) {
	    TString ct(h->GetName());
	    // Info("", "Centralty %s", ct.Data());
	    Int_t idx = ct.Index("cent");
	    if (idx != kNPOS) {
	      ct.Remove(0,idx+4);
	      TObjArray* tokens = ct.Tokenize("_");
	      TString    t1     = tokens->At(0)->GetName();
	      TString    t2     = tokens->At(1)->GetName();
	      TString    s1     = t1.Strip(TString::kLeading, '0');
	      TString    s2     = t2.Strip(TString::kLeading, '0');
	      s1.ReplaceAll("d", ".");
	      s2.ReplaceAll("d", ".");
	      TString    c1     = CentLimitName(false, s1.Atof());
	      TString    c2     = CentLimitName(false, s2.Atof());
	      nm = Form("%s-%s", c1.Data(), c2.Data());
	      // Printf("Got name=%s t1=%s t2=%s c1=%s c2=%s -> %s",
	      //        ct.Data(), t1.Data(), t2.Data(), s1.Data(), s2.Data(),
	      //        nm.Data());
	      tokens->Delete();
	      // Printf("%s (%s %s -> %s %s)", nm.Data(),
	      //        t1.Data(), c1.Data(),
	      //        t2.Data(), c2.Data());
		    
	    }
	  }
	  e = leg->AddEntry("dummy", nm, "pl");
	  e->SetMarkerStyle(h->GetMarkerStyle());
	  e->SetMarkerSize(h->GetMarkerSize());
	  if (isCent) {
	    e->SetLineColor(h->GetLineColor());
	    e->SetMarkerColor(h->GetMarkerColor());
	  }
	}
      }
      
      delete tmp;
      ptrig++;
      marker++;
      if (peff) peff++;
    }
    if (points.GetEntries() <= 0) {
      ::Warning("", "No data for %s %d %s", system.Data(), sNN, allT.Data());
      delete stack;
      return ret;
    }

    if (leg) { 
      leg->SetTextFont(42);
      leg->SetTextColor(AliceBlue());
    }
    if (false && seSeen && leg) { 
      e = leg->AddEntry("dummy", "7.6% sys. error", "f");
      e->SetFillColor(fill);
      e->SetFillStyle(style);
      leg->SetNColumns(leg->GetNColumns()+1);
    }
    TIter nextP(&points);
    TH1*  data = 0;
    while ((data = static_cast<TH1*>(nextP()))) 
      stack->Add(data);
    
    stack->SetMaximum(1.2*stack->GetMaximum("nostack"));
    stack->SetMinimum(.3);
    
    ptrig                    = trigs;
    marker                   = 20;
    TParameter<int>*   u     = 0;
    TMultiGraph*       other = new TMultiGraph("others","Others");
    TString            allE  = "";
    while ((*ptrig)) { 
      const char** pexp = exps;
      TString trg(*ptrig);
      Bool_t   isCent = (trg.BeginsWith("CENT"));

      while ((*pexp)) {
	if (ptrig == trigs) {
	  if (!allE.IsNull()) allE.Append(",");
	  allE.Append(*pexp);
	}
	Info("", "System=%s sNN=%d Trigger=%s Exp=%s", system.Data(),
	     sNN, *ptrig, *pexp);
	TMultiGraph* mg = GetOther(system, sNN, *ptrig, *pexp, verbose);
	if (!mg) {
	  pexp++;
	  continue;
	}
	TIter nextG(mg->GetListOfGraphs());
	TGraph* g = 0;
	while ((g = static_cast<TGraph*>(nextG()))) {
	  // g->Draw("p same");
	  other->Add(g, "p");
	  if (!isCent) {
	    g->SetMarkerStyle(marker);
	    g->SetMarkerSize(MarkerSize(marker));
	  }
	  if (isCent) {
	    Color_t cOld = g->GetMarkerColor();
	    Color_t cNew = Brighten(cOld);
	    g->SetMarkerColor(Brighten(g->GetMarkerColor()));
	    g->SetFillColor(Brighten(g->GetFillColor()));
	    g->SetLineColor(Brighten(g->GetLineColor()));
	    Printf("%s: old=%d new=%d real=%d",g->GetName(),
		   cOld,cNew,g->GetMarkerColor());
	  }
	  TString ut(g->GetTitle());
	  if (ut.BeginsWith("PWG-UD/MULT - "))
	    ut = "PWG-UD/MULT - work in progress";

	  if (unique.FindObject(ut)) continue;

	  u = new TParameter<int>(ut,g->GetMarkerStyle());
	  u->SetUniqueID(isCent ? Int_t(kBlack) : g->GetMarkerColor());
	  unique.Add(u);
	  
	}
	pexp++;
      }
      marker++;
      ptrig++;
    }
    if (!other ||
	!other->GetListOfGraphs() ||
	other->GetListOfGraphs()->GetEntries() < 0) {
      ::Warning("", "No data for %s %d %s %s",
		system.Data(), sNN, allT.Data(), allE.Data());
      if (other) delete other;
      other = 0;
    }
    // ::Info("", "Making pair of data (%p) and other (%p)", stack, other);
    ret = new TPair(stack, other);

    return ret;
  }			   
};

#endif
