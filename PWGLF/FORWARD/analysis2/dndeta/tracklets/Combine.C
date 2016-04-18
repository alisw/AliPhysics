#ifndef __CINT__
# include <TFile.h>
# include <THStack.h>
# include <TH1.h>
# include <TH2.h>
# include <TError.h>
# include <TMath.h>
# include <TClass.h>
# include <TCanvas.h>
# include <TSystem.h>
# include <TFitResult.h>
# include <TF1.h>
# include "GraphSysErr.C"
#else
class TFile;
class THStack;
class TH1;
class TH2;
class TCanvas;
class TDirectory;
class TDirectory;
#endif

const Bool_t kCombineLoaded = true;
Int_t cW = 1200;
Int_t cH =  800;

//____________________________________________________________________
void AddPath(const TString& dir, Bool_t prepend=true)
{
  TString d(gSystem->ExpandPathName(dir.Data()));
  gSystem->AddIncludePath("-I%s", d.Data());
  const char* oldPath = gROOT->GetMacroPath();
  gROOT->SetMacroPath(Form(".:%s:%s",
			   prepend ? d.Data() : oldPath,
			   prepend ? oldPath  : d.Data()));
}

//____________________________________________________________________
/** 
 * Define a correlated uncertainty 
 * 
 * @param o    Graph
 * @param name Name
 * @param val  Value (relative)
 * @param c    Color
 */
void MakeCommon(TObject* o, const char* name, Double_t val, Color_t c)
{
  GraphSysErr* gse = static_cast<GraphSysErr*>(o);
  Int_t id = gse->DefineCommon(name,true,val,GraphSysErr::kBox);
  gse->SetSysFillColor(id, c);
  gse->SetSysLineColor(id, c);
}
//____________________________________________________________________
/** 
 * Declare a point-to-point uncertainty 
 * 
 * @param o    Graph
 * @param name Name 
 * @param c    Color
 * 
 * @return Id of Syst.unc.
 */
Int_t MakeP2P(TObject* o, const char* name, Color_t c)
{
  GraphSysErr* gse = static_cast<GraphSysErr*>(o);
  Int_t id = gse->DeclarePoint2Point(name,true,GraphSysErr::kBox);
  gse->SetSysFillColor(id, c);
  gse->SetSysLineColor(id, c);
  return id;
}
  
//====================================================================
Double_t SysEval(Double_t x, Double_t sMin, Double_t sMax, Double_t xMax)
{
  return sMin + TMath::Power(x/xMax, 2)*(sMax-sMin);
}

Double_t CSysEval(Double_t x, Double_t sMin, Double_t sMax)
{
  return SysEval(x, sMin, sMax, 80);
}
Double_t EtaSysEval(Double_t x, Double_t sMin, Double_t sMax)
{
  return SysEval(x, sMin, sMax, 2);
}
  
//____________________________________________________________________
/** 
 * Make a GraphSysErr object
 * 
 * @param g  Graph 
 * @param c1 Least centrality 
 * @param c2 Largest centrality 
 * 
 * @return Newly created GraphSysErr
 */
TObject* MakeGSE(TH1* g, Int_t bin)
{
  if (!g) return 0;

  Double_t c1 = 0;
  Double_t c2 = 100;
  switch (bin) {
  case 0: c1 =  0; c2 =  5; break;
  case 1: c1 =  5; c2 = 10; break;
  case 2: c1 = 10; c2 = 20; break;
  case 3: c1 = 20; c2 = 30; break;
  case 4: c1 = 30; c2 = 40; break;
  case 5: c1 = 40; c2 = 50; break;
  case 6: c1 = 50; c2 = 60; break;
  case 7: c1 = 60; c2 = 70; break;
  case 8: c1 = 70; c2 = 80; break;
  case 9: c1 = 80; c2 = 90; break;
  }    
  if (!gROOT->GetClass("GraphSysErr")) return 0;
  TString  bnn; bnn.Form("%03dd%02d_%03dd%02d",
			 Int_t(c1), Int_t(c1*100)%100, 
			 Int_t(c2), Int_t(c2*100)%100);
  Bool_t  mc      = g->GetMarkerStyle() == 30;
  TString nme(bnn);
  if (mc) nme.Prepend("CENTT_");
  else    nme.Prepend("CENT_");

  Color_t  col     = g->GetMarkerColor();
  // Double_t bg      = (1-c1/100)*(2-0.1)+0.1;
  // Double_t c       = TMath::Power(c1/100,2)*(6.2-0.4)+0.4;
  Double_t bg      = CSysEval(c2, 0.02, 0.001);
  Double_t c       = CSysEval(c2, 0.005, 0.075);
  GraphSysErr* gse = new GraphSysErr(nme,"");  
  gse->SetName(nme);
  gse->SetTitle(Form("%5.1f - %5.1f%%", c1, c2));
  gse->SetKey("author", "PREGHENELLA : 2015");
  gse->SetKey("title", "dNch/deta in PbPb at 5023 GeV");
  gse->SetKey("obskey", "DN/DETARAP");
  gse->SetKey("reackey", "PB PB --> CHARGED X");
  gse->SetKey("laboratory", "CERN");
  gse->SetKey("accelerator", "LHC");
  gse->SetKey("detector", "TRACKLETS");
  gse->SetKey("reference", "ALICE-AN-2830");
  gse->AddQualifier("CENTRALITY IN PCT", Form("%.1f TO %.1f",c1,c2));
  gse->AddQualifier("SQRT(S)/NUCLEON IN GEV", "5023");
  gse->SetXTitle("ETARAP");
  gse->SetYTitle("DN/DETARAP");
  gse->SetMarkerStyle(g->GetMarkerStyle());
  gse->SetMarkerSize(g->GetMarkerSize());
  gse->SetDataOption(GraphSysErr::kNoTick);
  gse->SetMarkerColor(col);
  gse->SetLineColor(col);
  gse->SetFillColor(col);
  gse->SetSumFillColor(col);
  gse->SetSumLineColor(col);
  gse->SetSumOption(GraphSysErr::kBox);
  gse->SetCommonSumFillColor(col);
  gse->SetCommonSumLineColor(col);
  gse->SetCommonSumOption(GraphSysErr::kBox);
  if (!mc) {
    MakeCommon(gse, "Particle composition",	0.01,   col);
    MakeCommon(gse, "Weak decay",		0.01,   col);
    MakeCommon(gse, "pT extrapolation",		0.02,   col);
    MakeCommon(gse, "EG dependence",		0.02,   col);
    MakeCommon(gse, "Background subrtaction",	bg,     col);
    MakeCommon(gse, "Centrality",		c,      col);
  }
  Int_t acc = 0;
  if (!mc) acc = MakeP2P(gse, "Acceptance", col);

  Int_t j =  0;
  for (Int_t i = 1; i <= g->GetNbinsX(); i++) {
    Double_t eta  = g->GetXaxis()->GetBinCenter(i);
    Double_t eEta = g->GetXaxis()->GetBinWidth(i)/2;
    Double_t xo   = TMath::Abs(eta)+eEta;
    if (xo > 2) continue;
    Double_t ea  = 0.02*TMath::Power(xo/2,2);
    gse->SetPoint(j, eta, g->GetBinContent(i));
    gse->SetPointError(j, eEta, eEta);
    gse->SetStatError(j, g->GetBinError(i),g->GetBinError(i));
    if (acc != 0) 
      gse->SetSysError(acc, j, eEta, eEta, ea/100, ea/100);
    j++;
  }
  return gse;
    
  
}

//====================================================================
TObject* GetO(TDirectory* dir, const char* name="result", TClass* cls=0)
{
  if (!dir) {
    Warning("GetHS", "No directory");
    return 0;
  }
  TString par; par = gSystem->DirName(name);
  TString bse; bse = gSystem->BaseName(name);

  if (!par.EqualTo(".")) {
    TDirectory* save = dir;
    dir = save->GetDirectory(par);
    if (!dir) {
      // Warning("", "Failed to get directory %s", par.Data());
      // save->ls();
      return 0;
    }
  }
  TObject* o = dir->Get(bse);
  if (!o) {
    Warning("GetO", "%s not found in %s", name, dir->GetName());
    return 0;
  }
  if (!cls) return o;
  if (!o->IsA()->InheritsFrom(cls)) {
    Warning("GetO", "%s is not a %s!", name, cls->GetName());
    return 0;
  }
  return o;
}

//____________________________________________________________________
THStack* GetHS(TDirectory* dir, const char* name="result")
{
  return static_cast<THStack*>(GetO(dir,name,THStack::Class()));
}
//____________________________________________________________________
TH1* GetH1(TDirectory* dir, const char* name="result")
{
  return static_cast<TH1*>(GetO(dir,name,TH1::Class()));
}
//____________________________________________________________________
TH2* GetH2(TDirectory* dir, const char* name="result")
{
  return static_cast<TH2*>(GetO(dir,name,TH2::Class()));
}

//====================================================================
TH1* Combine(TH1* left, TH1* middle, TH1* right, TDirectory* out,
	     Bool_t cutAtTwo=true)
{
  TH1* hs[] = { left, middle, right, 0 };
  TH1* r    = static_cast<TH1*>(middle->Clone());
  r->SetDirectory(0);
  r->Reset();
  if (cutAtTwo) {
    if (r->GetMarkerStyle() == 20) 
      r->SetMarkerStyle(29);
    else 
      r->SetMarkerStyle(30);
  }

  for (Int_t i = 1; i <= r->GetNbinsX(); i++) {
    Double_t aeta = TMath::Abs(r->GetXaxis()->GetBinCenter(i));
    Double_t weta = r->GetXaxis()->GetBinWidth (i);
    if (cutAtTwo && (aeta + weta/2) > 2) continue; // Kill bins outside +/- 2
    TH1**    ph   = hs;
    Double_t sum  = 0;
    Double_t sumw = 0;
    Int_t    n    = 0;
    while (*ph) {
      TH1* h = *ph;
      ph++;

      Double_t c = h->GetBinContent(i);
      Double_t e = h->GetBinError(i);
      if (c < 1e-6) continue;

      sum   += c;
      sumw  += e * e;
      n     += 1;
    }
    if (n <= 0) continue;
    r->SetBinContent(i, sum/n);
    r->SetBinError  (i, TMath::Sqrt(sumw));
  }
  return r;
}

//____________________________________________________________________
TH2* Combine(TH2* left, TH2* middle, TH2* right, TDirectory* out)
{
  TH2* hs[] = { left, middle, right, 0 };
  TH2* r    = static_cast<TH2*>(middle->Clone());
  r->SetDirectory(0);
  r->Reset();
  if (r->GetMarkerStyle() == 20) 
    r->SetMarkerStyle(29);
  else 
    r->SetMarkerStyle(30);
   

  for (Int_t i = 1; i <= r->GetNbinsX(); i++) {
    Double_t aeta = TMath::Abs(r->GetXaxis()->GetBinCenter(i));
    Double_t weta = r->GetXaxis()->GetBinWidth (i);
    if ((aeta + weta/2) > 2) continue; // Kill bins outside +/- 2
    for (Int_t j = 1; j <= r->GetNbinsY(); j++) {
      TH2**    ph   = hs;
      Double_t sum  = 0;
      Double_t sumw = 0;
      Int_t    n    = 0;
      while (*ph) {
	TH2* h = *ph;
	ph++;
	
	Double_t c = h->GetBinContent(i,j);
	Double_t e = h->GetBinError(i,j);
	if (c < 1e-6) continue;

	sum   += c;
	sumw  += e * e;
	n     += 1;
      }
      if (n <= 0) continue;
      r->SetBinContent(i, sum/n);
      r->SetBinError  (i, TMath::Sqrt(sumw));
    }
  }
  return r;
}

//____________________________________________________________________
void
CombineMap(TH2*        sleft,
	   TH2*        smiddle,
	   TH2*        sright,
	   TDirectory* out)
{
  Int_t nyLeft   = sleft  ->GetYaxis()->GetNbins();
  Int_t nyMiddle = smiddle->GetYaxis()->GetNbins();
  Int_t nyRight  = sright ->GetYaxis()->GetNbins();
  TH2*  maps[]   = { sleft, smiddle, sright };
  TArrayD bins(nyLeft+nyMiddle+nyRight+1);  
  bins[0] = sleft->GetYaxis()->GetXmin();
  Int_t j = 1;
  for (Int_t k = 0; k < 3; k++) {
    TH2* tmp =  maps[k];
    for (Int_t i = 1; i <= tmp->GetYaxis()->GetNbins(); i++, j++)
      bins[j] = tmp->GetYaxis()->GetBinUpEdge(i);
  }

  TH2* ret = new TH2D(smiddle->GetName(), smiddle->GetTitle(),
		      smiddle->GetXaxis()->GetNbins(),
		      smiddle->GetXaxis()->GetXmin(),
		      smiddle->GetXaxis()->GetXmax(),
		      nyLeft+nyMiddle+nyRight, bins.GetArray());
  ret->SetXTitle(smiddle->GetXaxis()->GetTitle());
  ret->SetYTitle(smiddle->GetYaxis()->GetTitle());
  ret->SetZTitle(smiddle->GetZaxis()->GetTitle());  
  ret->SetDirectory(out);
  static_cast<TAttAxis*>(smiddle->GetXaxis())->Copy(*ret->GetXaxis());
  static_cast<TAttAxis*>(smiddle->GetYaxis())->Copy(*ret->GetYaxis());
  static_cast<TAttAxis*>(smiddle->GetZaxis())->Copy(*ret->GetZaxis());

  Int_t base = 0;
  for (Int_t k = 0; k < 3; k++) {
    TH2* tmp =  maps[k];
    for (Int_t i = 1; i <= tmp->GetXaxis()->GetNbins(); i++) {
      for (Int_t j = 1; j <= tmp->GetYaxis()->GetNbins(); j++) {
	ret->SetBinContent(i,base+j,tmp->GetBinContent(i,j));
	ret->SetBinError  (i,base+j,tmp->GetBinError  (i,j));
      }
    }
    base += tmp->GetYaxis()->GetNbins();
  }
}

  
  
  
  
  
  
  
//____________________________________________________________________
void
Combine(THStack*    sleft,
	THStack*    smiddle,
	THStack*    sright,
	THStack*    result,
	TDirectory* out,
	TList*      gses,
	TH1*        cent=0,
	Bool_t      cutAtTwo=true)
{
  TIter    ileft  (sleft  ->GetHists());
  TIter    imiddle(smiddle->GetHists());
  TIter    iright (sright ->GetHists());
  TH1*     hleft   = 0;
  TH1*     hmiddle = 0;
  TH1*     hright  = 0;
  Int_t    cnt     = 0;
  TH1*     mid     = cent ? static_cast<TH1*>(cent->Clone("mid")) : 0;
  if (mid) {
    mid->Reset();
    mid->SetDirectory(out);
    mid->SetXTitle("#eta");
    mid->SetYTitle("d#it{N}_{ch}/d#it{#eta}|_{|#it{#eta}|<0.5}");
    mid->SetTitle("Mid-rapidity");
  }
  Int_t bin = 0;
  while ((hleft    = static_cast<TH1*>(ileft  ())) &&
	 (hmiddle  = static_cast<TH1*>(imiddle())) &&
	 (hright   = static_cast<TH1*>(iright ()))) {
    TH1* h = Combine(hleft, hmiddle, hright, out, cutAtTwo);
    if (mid && !TString(h->GetName()).Contains("truth")) {
      TF1* f = new TF1(Form("f%s", h->GetName()), "pol0", -.5, +.5);
      TFitResultPtr s = h->Fit(f, "QS+", "", -.5, +.5);
      bin++;
      mid->SetBinContent(bin, s->Parameter(0));
      mid->SetBinError  (bin, s->ParError (0));
      Printf("dNch/deta %20s: %6.1f +/- %6.1f (%5.2f)",
	     h->GetName(), s->Parameter(0), s->ParError(0), s->Chi2()/s->Ndf());
    }
    if (result) result->Add(h);
    if (gses) {
      TObject* g = MakeGSE(h, cnt);
      if (h->GetMarkerStyle() == 29) cnt++;
      // Printf("Adding GSE %s -> %s", h->GetName(), g->GetName());
      gses->Add(g);
    }
  }
}  
  
void
Combine(UShort_t flags=0x0, const char* var="none")
{
  const char* fwd = "$ALICE_ROOT/PWGLF/FORWARD/analysis2";
  AddPath("$HOME/GraphSysErr");
  AddPath(TString::Format("%s/gse", fwd), false);
  AddPath(TString::Format("%s/dndeta/tracklets", fwd));
  if (!gROOT->GetClass("GraphSysErr")) 
    gROOT->LoadMacro("$HOME/GraphSysErr/GraphSysErr.C+g");
  
  UShort_t which = flags & 0x3;
  TFile* fleft    = TFile::Open(Form("partial/left_%s_0x%x.root",  var,which),
				"READ");
  TFile* fmiddle  = TFile::Open(Form("partial/middle_%s_0x%x.root",var,which),
				"READ");
  TFile* fright   = TFile::Open(Form("partial/right_%s_0x%x.root", var,which),
				"READ");

  TH1* cent = GetH1(fmiddle, "cent");
  THStack* sleft   = GetHS(fleft);
  THStack* smiddle = GetHS(fmiddle);
  THStack* sright  = GetHS(fright);

  gSystem->mkdir("results",1);
  TFile* out = TFile::Open(Form("results/combine_%s_0x%x.root", var, which),
			   "RECREATE");
  const char* obs  = "\\mathrm{d}N_{\\mathrm{ch}}/\\mathrm{d}\\eta";
  TList*   gses    = new TList;
  THStack* result  = new THStack("result",
				 Form("%s\\hbox{ %s %s}",
				      obs, which == 0x3 ? "combinatorics" :
				      "injection", var));
  Combine(sleft, smiddle, sright, result, out, gses, cent);
  result->Write();
  cent->Write("cent");
  
  if (cent) {
    for (Int_t i = 1; i <= cent->GetNbinsX(); i++) {
      TString name;
      Double_t c1 = cent->GetXaxis()->GetBinLowEdge(i);
      Double_t c2 = cent->GetXaxis()->GetBinUpEdge (i);
      name.Form("cent%03dd%02d_%03dd%02d",
		Int_t(c1), Int_t(c1*100)%100, Int_t(c2), Int_t(c2*100)%100);
      TString sname(name); sname.Append("/summary");
      sleft   = GetHS(fleft,  sname);
      smiddle = GetHS(fmiddle,sname);
      sright  = GetHS(fright, sname);
      if (!sleft || !smiddle || !sright) continue;
      TDirectory* dir = out->mkdir(name);
      dir->cd();
      THStack* summary  = new THStack("summary","");
      Combine(sleft,smiddle,sright,summary, 0, 0);
      summary->Write();
      dir->cd("..");

      TDirectory* det = dir->mkdir("details");
      const char* maps[] = { "realMeas", "realBg", "realSig",
			     "simMeas",  "simBg",  "simSig",
			     "trueGen",  "alpha",  "fiducial", 0 };
      const char** pmap  = maps;
      while (*pmap) {
	sname.Form("%s/details/%s", name.Data(), *pmap);
	TH2* mleft   = GetH2(fleft,  sname);
	TH2* mmiddle = GetH2(fmiddle,sname);
	TH2* mright  = GetH2(fright, sname);
	if (!mleft || !mmiddle || !mright) continue;
	CombineMap(mleft, mmiddle, mright, det);
	pmap++;
      }
      sname.Form("%s/details/deltaInt", name.Data());
      TH1* hleft   = GetH1(fleft,  sname);
      TH1* hmiddle = GetH1(fmiddle,sname);
      TH1* hright  = GetH1(fright, sname);
      if (hleft && hmiddle && hright) {
	TH1* kHist = Combine(hleft, hmiddle, hright, det);
	det->cd();
	kHist->Write();
      }
      

      sname.Form("%s/details/deltas", name.Data());
      sleft   = GetHS(fleft,  sname);
      smiddle = GetHS(fmiddle,sname);
      sright  = GetHS(fright, sname);
      if (!sleft || !smiddle || !sright) continue;
      THStack* deltas  = new THStack("deltas","");
      Combine(sleft,smiddle,sright,deltas, 0, 0, 0, false);
      det->cd();
      deltas->Write();

      TDirectory* avg = dir->mkdir("averages");
      avg->cd();
      const char* avgs[] = { "realAvgMeas", "realAvgBg", "realAvgSig",
			     "simAvgMeas",  "simAvgBg",  "simAvgSig",
			     "truth",       0 };
      const char** pavg  = avgs;
      while (*pavg) {
	sname.Form("%s/averages/%s", name.Data(), *pavg);
	hleft   = GetH1(fleft,  sname);
	hmiddle = GetH1(fmiddle,sname);
	hright  = GetH1(fright, sname);
	if (!hleft || !hmiddle || !hright) continue;
	TH1* avgH = Combine(hleft, hmiddle, hright, avg);
	avgH->Write();
	pavg++;
      }
      
    }
  }
  gSystem->mkdir("plots",1);
  TCanvas* c = new TCanvas("c","c", cW, cH);
  c->SetTopMargin(0.01);
  c->SetRightMargin(0.01); 
  result->Draw("nostack");
  result->GetHistogram()->GetYaxis()->SetTitleOffset(1.5);
  result->GetHistogram()->SetYTitle(obs);
  result->GetHistogram()->SetXTitle("\\eta");
  c->Print(Form("plots/combine_%s_0x%x.png", var, which));

  out->cd();
  if (gses) {
    gses->SetName("container");
    gses->Draw("quad combine stat xbase=2.2");
    gses->Write(gses->GetName(), TObject::kSingleKey);
  }
  Printf("Output stored in %s", out->GetName());
  out->Write();
}


		
//____________________________________________________________________
//
// EOF
// 
