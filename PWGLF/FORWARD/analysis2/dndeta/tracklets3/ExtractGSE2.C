/**
 * @file   PWGLF/FORWARD/analysis2/dndeta/tracklets3/ExtractGSE2.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Wed Apr 27 16:51:39 2016
 * 
 * @brief  Extract GraphSysErr from results 
 * 
 * 
 * @ingroup pwglf_forward_tracklets
 */
#ifndef __CINT__
# include "GraphSysErr.C"
# include <TFile.h>
# include <TList.h>
# include <TMath.h>
# include <TCanvas.h>
# include <TLegend.h>
# include <TSeqCollection.h>
# include <TLegendEntry.h>
# include <THStack.h>
#else
class TSeqCollection;
class TCanvas;
class TGraphErrors;
class THStack;
class TLegend;
class TH1;
#endif

//====================================================================
/** 
 * Find bin number correspondig to centrality range 
 * 
 * @param c1 Least centrality 
 * @param c2 Largest centrality 
 * 
 * @return Bin number 
 */
static Int_t PbPbBin(Double_t c1, Double_t c2)
{
  Double_t c = (c1+c2) / 2;
  if      (c <  5) return 0;
  else if (c < 10) return 1;
  else if (c < 20) return 2;
  else if (c < 30) return 3;
  else if (c < 40) return 4;
  else if (c < 50) return 5;
  else if (c < 60) return 6;
  else if (c < 70) return 7;
  else if (c < 80) return 8;
  else if (c < 90) return 9;
  return                  10;
}
  
//____________________________________________________________________
/** 
 * Get the centrality color for PbPb 
 * 
 * @param c1 Lower edge
 * @param c2 Upper edge
 * 
 * @return Color 
 */
static Color_t PbPbColor(Double_t c1, Double_t c2)
{
  const Color_t cc[] = { kMagenta+2,
			 kBlue+2,
			 kAzure-1, // 10,
			 kCyan+2,
			 kGreen+1,
			 kSpring+5,//+10,
			 kYellow+1,
			 kOrange+5,//+10,
			 kRed+1,
			 kPink+5,//+10,
			 kBlack };
  Int_t bin = PbPbBin(c1,c2);
  return cc[bin];
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
 * @param d   Directory 
 * @param dimen Dimension 
 * @param sNN Collision energy in GeV
 * @param c1  Least centrality 
 * @param c2  Largest centrality 
 * 
 * @return Newly created GraphSysErr
 */
TObject* MakeGSE(TDirectory* d,
		 Int_t dimen,
		 UShort_t sNN,
		 Double_t c1,
		 Double_t c2)
{
  if (!gROOT->GetClass("GraphSysErr")) return 0;
  TString  bin; bin.Form("%03dd%02d_%03dd%02d",
			 Int_t(c1), Int_t(c1*100)%100, 
			 Int_t(c2), Int_t(c2*100)%100);
  Bool_t isAll = (c1+1.e-9 >= c2);
  if (isAll) bin = "all";
  TString sub(bin);
  if (!isAll) sub.Prepend("cent");
  sub.Append(Form("/results%dd/result", dimen));
  TString nme(bin); nme.Prepend("CENT_");
  TH1*    g = GetH1(d, sub);
  if (!g) return 0;

  Double_t eff = g->GetBinContent(0); 
  Printf("Trigger efficiency: %6.4f", eff);
  if (eff < 1e-6) eff = 1;
  g->Scale(eff);
  g->SetBinContent(0,0);
  
  Color_t  col     = PbPbColor(c1,c2); // g->GetMarkerColor();
  // Double_t bg      = (1-c1/100)*(2-0.1)+0.1;
  // Double_t c       = TMath::Power(c1/100,2)*(6.2-0.4)+0.4;
  // Preliminary centrality sys error for Xe-Xe
  Double_t aMin    = sNN == 5440 ? -0.09  : 0.02;
  Double_t pComp   = sNN == 5440 ? 0.002  : 0.01;
  Double_t wDcy    = sNN == 5440 ? 0.004  : 0.01;
  Double_t ptInt   = sNN == 5440 ? 0.004  : 0.02;
  Double_t ptMax   = sNN == 5440 ? 0.010  : 0.00;
  Double_t egDep   = sNN == 5440 ? 0.000  : 0.02;
  Double_t cMin    = sNN == 5440 ? 0.001  : sNN == 5023 ? 0.005 : 0.004;
  Double_t cMax    = sNN == 5440 ? 0.048  : sNN == 5023 ? 0.075 : 0.062;
  Double_t tMax    = sNN == 5440 ? 0.040  : 0;
  Double_t tMin    = sNN == 5440 ? 0.003  : 0;
  Double_t bgMax   = sNN == 5440 ? 0.011  : 0.02;
  Double_t bgMin   = sNN == 5440 ? 0.005  : 0.001;
  Double_t bg      = CSysEval(c2, bgMax, bgMin);
  Double_t c       = CSysEval(c2, cMin, cMax);
  Double_t em1     = sNN == 5440 ? 0.012 : sNN == 5023 ? 0.04 : 0.02;
  Double_t em2     = sNN == 5440 ? 0.040 : 5023 ? 0.04 : 0.02;
  Double_t em      = c1 >= 80 ? em2 : c1 >= 70 ? em1 : 0;
  Double_t tkl     = CSysEval(c2, tMax, tMin);
  GraphSysErr* gse = new GraphSysErr(g->GetNbinsX());
  gse->SetName(nme);
  gse->SetTitle(Form("%5.1f - %5.1f%%", c1, c2));
  //gse->SetKey("author",(sNN == 5023?"PREGHENELLA : 2015":"SHAHOYAN : 2013"));
  gse->SetKey("title", Form("dNch/deta in PbPb at %d GeV", sNN));
  gse->SetKey("obskey", "DN/DETARAP");
  gse->SetKey("reackey", "PB PB --> CHARGED X");
  gse->SetKey("laboratory", "CERN");
  gse->SetKey("accelerator", "LHC");
  gse->SetKey("detector", "TRACKLETS");
  gse->SetKey("reference", sNN==5023 ? "ALICE-AN-2830" : "ALICE-AN-2180");
  if (!isAll)gse->AddQualifier("CENTRALITY IN PCT", Form("%.1f TO %.1f",c1,c2));
  gse->AddQualifier("SQRT(S)/NUCLEON IN GEV", Form("%d", sNN));
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
  MakeCommon(gse, "Particle composition",	pComp,  col);
  MakeCommon(gse, "Weak decay",			wDcy,   col);
  Int_t pte = -1;
  if (ptMax < 1e-3)
    MakeCommon(gse, "pT extrapolation",		ptInt,  col);
  else
    pte = MakeP2P(gse, "pT extrapolation", col);
  MakeCommon(gse, "Background subtraction",	bg,     col);
  MakeCommon(gse, "Centrality",			c,      col);
  MakeCommon(gse, "EM contamination",		em,     col);
  if (egDep > 1e-6) 
    MakeCommon(gse, "EG dependence",		egDep,  col);
  if (sNN==5440)
    MakeCommon(gse, "Material budget",		0.001,  col);
  if (tkl > 1e-6)
    MakeCommon(gse, "Tracklet selection",	tkl,    col);
    
  if (eff) {
    MakeCommon(gse, "TRIGGER", 0.02, col);
  }
  Int_t acc = -1;
  if (aMin > 0)
    acc = MakeP2P(gse, "Acceptance", col);
  else
    MakeCommon(gse, "Acceptance", -aMin, col);
  

  Int_t j =  0;
  for (Int_t i = 1; i <= g->GetNbinsX(); i++) {
    Double_t eta  = g->GetXaxis()->GetBinCenter(i);
    Double_t eEta = g->GetXaxis()->GetBinWidth(i)/2;
    Double_t xo   = TMath::Abs(eta)+eEta;
    if (xo > 2) continue;
    Double_t ea  = aMin*TMath::Power(xo/2,2);
    Double_t ep  = EtaSysEval(xo,ptInt, ptMax);
    gse->SetPoint(j, eta, g->GetBinContent(i));
    gse->SetPointError(j, eEta, eEta);
    gse->SetStatError(j, g->GetBinError(i),g->GetBinError(i));
    if (acc > 0) gse->SetSysError(acc, j, eEta, eEta, ea/100, ea/100);
    if (pte > 0) gse->SetSysError(pte, j, eEta, eEta, ep/100, ep/100);
    j++;
  }
  return gse;
    
  
}

//____________________________________________________________________
/** 
 * Make a GraphSysErr object for simulations 
 * 
 * @param d   Directory 
 * @param dimen Dimension 
 * @param sNN Collision energy in GeV
 * @param c1 Least centrality 
 * @param c2 Largest centrality 
 * 
 * @return Newly created GraphSysErr
 */
TObject* MakeTGSE(TDirectory* d,
		  Int_t dimen,
		  UShort_t sNN,
		  Double_t c1,
		  Double_t c2)
{
  if (!gROOT->GetClass("GraphSysErr")) return 0;
  TString  bin; bin.Form("%03dd%02d_%03dd%02d",
			 Int_t(c1), Int_t(c1*100)%100, 
			 Int_t(c2), Int_t(c2*100)%100);
  Bool_t isAll = (c1+1.e-9 >= c2);
  if (isAll) bin = "all";
  TString sub(bin);
  if (!isAll) sub.Prepend("cent");
  sub.Append(Form("/results%dd/simG",dimen));
  TString nme(bin); nme.Prepend("CENTT_");
  TH1*    g = GetH1(d, sub);
  if (!g) return 0;
  
  Color_t  col     = PbPbColor(c1,c2); // g->GetMarkerColor();
  // Double_t bg      = (1-c1/100)*(2-0.1)+0.1;
  // Double_t c       = TMath::Power(c1/100,2)*(6.2-0.4)+0.4;
  GraphSysErr* gse = new GraphSysErr(g->GetNbinsX());
  gse->SetName(nme);
  gse->SetTitle(Form("%5.1f - %5.1f%%", c1, c2));
  gse->SetKey("author", (sNN == 5023 ? "PREGHENELLA : 2015":"SHAHOYAN : 2013"));
  gse->SetKey("title", Form("dNch/deta in PbPb at %d GeV", sNN));
  gse->SetKey("obskey", "DN/DETARAP");
  gse->SetKey("reackey", "PB PB --> CHARGED X");
  gse->SetKey("laboratory", "CERN");
  gse->SetKey("accelerator", "LHC");
  gse->SetKey("detector", "TRACKLETS");
  gse->SetKey("reference", sNN==5023 ? "ALICE-AN-2830" : "ALICE-AN-2180");
  if (!isAll)gse->AddQualifier("CENTRALITY IN PCT", Form("%.1f TO %.1f",c1,c2));
  gse->AddQualifier("SQRT(S)/NUCLEON IN GEV", Form("%d", sNN));
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
    j++;
  }
  return gse;
    
  
}

//____________________________________________________________________
TObject*
GetO(TDirectory* dir, const char* name, TClass* cls=0)
{
  if (!dir) {
    Warning("GetO", "No directory passed");
    return 0;
  }

  TObject* o = dir->Get(name);
  if (!o) {
    Warning("GetO", "object %s not found in %s",
	    name, dir->GetPath());
    return 0;
  }
  if (!cls) return o;
  if (!o->IsA()->InheritsFrom(cls)) {
    Warning("GetO", "Object %s in %s is not a %s, but a %s",
	    name, dir->GetPath(), cls->GetName(), o->ClassName());
    return 0;
  }
  return o;
}
//____________________________________________________________________
TDirectory* GetD(TDirectory* dir, const char* name)
{
  return static_cast<TDirectory*>(GetO(dir,name,TDirectory::Class()));
}

//____________________________________________________________________
TH1* GetH1(TDirectory* dir, const char* name)
{
  return static_cast<TH1*>(GetO(dir,name,TH1::Class()));
}
//____________________________________________________________________
/**
 * Steering 
 * 
 */
void
ExtractGSE2(const char* input, UShort_t sNN=5023)
{
  if (!gROOT->GetClass("GraphSysErr")) 
    gROOT->LoadMacro("$HOME/GraphSysErr/GraphSysErr.C+g");

  TString base = input; //gSystem->DirName(input); base.ReplaceAll(".root","");
  TFile*  file = TFile::Open(Form("%s/result.root",input), "READ");
  if (!file) {
    Warning("ExtractGSE2", "Failed to open %s/result.root", input);
    return;
  }
  
  Int_t dimen = -1;
  if      (base.EndsWith("unit"))    dimen = 0;
  else if (base.EndsWith("const"))   dimen = 1;
  else if (base.EndsWith("eta"))     dimen = 2;
  else if (base.EndsWith("etaipz"))  dimen = 3;
  if (dimen < 0) {
    Error("ExtractGSE2", "Don't know how to extract dimension from %s",base);
    return;
  }

  TH1* frame = 0;
  TH1* cent  = GetH1(file, "realCent");

  Bool_t   first  = true;
  TList*   stack  = new TList;
  TList*   truths = new TList;
  for (Int_t i = 1; i <= cent->GetNbinsX(); i++) {
    Double_t c1 = cent->GetXaxis()->GetBinLowEdge(i);
    Double_t c2 = cent->GetXaxis()->GetBinUpEdge(i);
    TObject* g  = MakeGSE(file,  dimen, sNN, c1, c2);
    TObject* t  = MakeTGSE(file, dimen, sNN, c1, c2);
    Double_t gmin, gmax, tmin, tmax;
    if (!g) continue;
    stack->Add(g);    
    truths->Add(t);    
    if (first) g->Draw("quad stat combine axis");
    else       g->Draw("quad stat combine");
    if (t)     t->Draw("quad");
    first = false;
    if (!frame) {
      GraphSysErr* gse = static_cast<GraphSysErr*>(g);
      if (gse->GetMulti())
	frame = gse->GetMulti()->GetHistogram();
    }
  }
  if (frame) frame->SetMinimum(1);
	       
  // TString obase(base); obase.Prepend("GSE_");
  
  std::ofstream out(Form("%s/gse.input", base.Data()));
  GraphSysErr::Export(stack, out, "HFC", 2);
  out << "*E" << std::endl;
  out.close();

  TFile* rout = TFile::Open(Form("%s/gse.root", base.Data()), "RECREATE");
  stack->AddAll(truths);
  stack->Write("container", TObject::kSingleKey);
  rout->Write();
  
}

void
ExtractGSE(Int_t flags)
{
  ExtractGSE(Form("MiddNdeta_0x%x.root", flags));
}
//____________________________________________________________________
