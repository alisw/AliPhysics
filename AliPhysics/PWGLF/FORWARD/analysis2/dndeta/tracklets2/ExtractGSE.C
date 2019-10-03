/**
 * @file   PWGLF/FORWARD/analysis2/dndeta/tracklets2/ExtractGSE.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Mon Aug 15 09:38:59 2016
 * 
 * @brief  
 * 
 * 
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
/** 
 * Evaluate a systematic error 
 * 
 * @param x    Where 
 * @param sMin Least error 
 * @param sMax Largest error 
 * @param xMax Largest @f$ x@f$ 
 * 
 * @return Systematic error 
 */
Double_t SysEval(Double_t x, Double_t sMin, Double_t sMax, Double_t xMax)
{
  return sMin + TMath::Power(x/xMax, 2)*(sMax-sMin);
}
/** 
 * Evaluate a centrality dependent systematic error 
 * 
 * @param x    Centrality 
 * @param sMin Least error 
 * @param sMax Largest error 
 * 
 * @return Systematic error 
 */
Double_t CSysEval(Double_t x, Double_t sMin, Double_t sMax)
{
  return SysEval(x, sMin, sMax, 80);
}
/** 
 * Evaluate an @f$ \eta@f$  dependent systematic error 
 * 
 * @param x    @f$ \eta@f$ 
 * @param sMin Least error 
 * @param sMax Largest error 
 * 
 * @return  Systematic error
 */
Double_t EtaSysEval(Double_t x, Double_t sMin, Double_t sMax)
{
  return SysEval(x, sMin, sMax, 2);
}
  
//____________________________________________________________________
/** 
 * Make a GraphSysErr object
 * 
 * @param d  Directory
 * @param c1 Least centrality 
 * @param c2 Largest centrality 
 * 
 * @return Newly created GraphSysErr
 */
TObject* MakeGSE(TDirectory* d, Double_t c1, Double_t c2)
{
  if (!gROOT->GetClass("GraphSysErr")) return 0;
  TString  bin; bin.Form("%03dd%02d_%03dd%02d",
			 Int_t(c1), Int_t(c1*100)%100, 
			 Int_t(c2), Int_t(c2*100)%100);
  TString sub(bin); sub.Prepend("cent"); sub.Append("/dndeta");
  TString nme(bin); nme.Prepend("CENT_");
  TH1*    g = GetH1(d, sub);
  if (!g) return 0;
  
  Color_t  col     = g->GetMarkerColor();
  // Double_t bg      = (1-c1/100)*(2-0.1)+0.1;
  // Double_t c       = TMath::Power(c1/100,2)*(6.2-0.4)+0.4;
  Double_t bg      = CSysEval(c2, 0.02, 0.001);
  Double_t c       = CSysEval(c2, 0.005, 0.075);
  GraphSysErr* gse = new GraphSysErr(g->GetNbinsX());
  gse->SetName(nme);
  gse->SetTitle(Form("%5.1f - %5.1f%%", c1, c2));
  gse->SetKey("author", "PREGHENELLA : 20150");
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
  MakeCommon(gse, "Particle composition",	0.01,   col);
  MakeCommon(gse, "Weak decay",			0.01,   col);
  MakeCommon(gse, "pT extrapolation",		0.02,   col);
  MakeCommon(gse, "EG dependence",		0.02,   col);
  MakeCommon(gse, "Background subrtaction",	bg,     col);
  MakeCommon(gse, "Centrality",			c,      col);
  Int_t acc = MakeP2P(gse, "Acceptance", col);

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
    gse->SetSysError(acc, j, eEta, eEta, ea/100, ea/100);
    j++;
  }
  return gse;
    
  
}

//____________________________________________________________________
/** 
 * Get an object from a directory/collection 
 * 
 * @param dir  Directory/collection 
 * @param name Name of object 
 * @param cls  Optionally, the type 
 * 
 * @return Pointer to read object (if type is specified, only if of
 * the right type), or null
 */
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
/** 
 * Get a directory from a directory/collection 
 * 
 * @param dir  Directory/collection 
 * @param name Name of object 
 * 
 * @return Pointer to read object or null
 */
TDirectory* GetD(TDirectory* dir, const char* name)
{
  return static_cast<TDirectory*>(GetO(dir,name,TDirectory::Class()));
}

//____________________________________________________________________
/** 
 * Get a histogram from a directory/collection 
 * 
 * @param dir  Directory/collection 
 * @param name Name of object 
 * 
 * @return Pointer to read object or null
 */
TH1* GetH1(TDirectory* dir, const char* name)
{
  return static_cast<TH1*>(GetO(dir,name,TH1::Class()));
}
//____________________________________________________________________
/**
 * Steering 
 * 
 * @param input Input file 
 */
void
ExtractGSE(const char* input)
{
  if (!gROOT->GetClass("GraphSysErr")) 
    gROOT->LoadMacro("$HOME/GraphSysErr/GraphSysErr.C+g");

  TString base = gSystem->BaseName(input); base.ReplaceAll(".root","");
  TFile*  file = TFile::Open(input, "READ");
  if (!file) return;

  TH1* cent = GetH1(file, "realCent");

  Bool_t first = true;
  TList* stack = new TList;
  for (Int_t i = 1; i <= cent->GetNbinsX(); i++) {
    Double_t c1 = cent->GetXaxis()->GetBinLowEdge(i);
    Double_t c2 = cent->GetXaxis()->GetBinUpEdge(i);
    TObject* g  = MakeGSE(file, c1, c2);
    if (!g) continue;
    stack->Add(g);    
    if (first) g->Draw("quad stat combine axis");
    else       g->Draw("quad stat combine");
    first = false;
  }

  TString obase(base); obase.Prepend("GSE_");
  
  std::ofstream out(Form("%s.input", obase.Data()));
  GraphSysErr::Export(stack, out, "HFC", 2);
  out << "*E" << std::endl;
  out.close();

  TFile* rout = TFile::Open(Form("%s.root", obase.Data()), "RECREATE");
  stack->Write("container", TObject::kSingleKey);
  rout->Write();
  
}
/** 
 * Create GSE 
 * 
 * @param flags Flags of input file  
 */
void
ExtractGSE(Int_t flags)
{
  ExtractGSE(Form("MiddNdeta_0x%x.root", flags));
}
//____________________________________________________________________
