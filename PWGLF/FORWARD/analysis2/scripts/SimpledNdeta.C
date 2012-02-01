#ifndef __CINT__
# include <TStyle.h>
# include <TFile.h>
# include <TList.h>
# include <TH1.h>
# include <THStack.h>
# include <TLegend.h>
# include <TLegendEntry.h>
# include <TLatex.h>
# include <TLine.h>
# include <TString.h>
# include <TCanvas.h>
# include <TError.h>
# include <TColor.h>
#else
class THStack;
class TAxis;
#endif

/**
 * A simple script to draw results from MakedNdeta.C (or similar)
 * 
 * @ingroup pwg2_forward_scripts
 */
/** 
 * Get a stack from the passed list 
 * 
 * @param list   List to get the stack from 
 * @param name   Name of stack 
 * @param rebin  Optional rebinning - must exists in list 
 * 
 * @return Stack or null
 * 
 * @ingroup pwg2_forward_scripts
 */
THStack*
GetStack(const TList* list, const char* name, Int_t rebin)
{
  if (!list) { 
    Warning("GetStack", "List is null");
    return 0;
  }

  TString n(name);
  if (rebin > 1) n.Append(Form("_rebin%02d", rebin));
  TObject* o = list->FindObject(n);
  if (!o) { 
    Warning("GetStack", "No %s object found in %s", n.Data(), list->GetName());
    return 0;
  }
  return static_cast<THStack*>(o);
}

/** 
 * Get a histogram from a list 
 * 
 * @param list   List 
 * @param name   Name of histogram
 * @param rebin  Rebinning factor
 * 
 * @return Histogram or null
 * 
 * @ingroup pwg2_forward_scripts
 */
TH1*
GetHist(const TList* list, const char* name, Int_t rebin)
{
  if (!list) { 
    Warning("GetStack", "List is null");
    return 0;
  }
  TList* all = static_cast<TList*>(list->FindObject("all"));
  if (!all) { 
    Warning("GetHist", "List all not found in %s", list->GetName());
    return 0;
  }

  TString n(name);
  if (rebin > 1) n.Append(Form("_rebin%02d", rebin));
  TObject* o = all->FindObject(n);
  if (!o) { 
    Warning("GetStack", "No %s object found in %s", n.Data(), all->GetName());
    return 0;
  }
  return static_cast<TH1*>(o);
}

/** 
 * Add histograms from one stack to another 
 * 
 * @param p      Parent stack to add to 
 * @param list   List to look for the stack in 
 * @param name   Name of stack to add 
 * @param rebin  Optional rebinning - must exists in list 
 * 
 * @return Added stack or null
 * 
 * @ingroup pwg2_forward_scripts
 */
THStack*
AddStack(THStack* p, const TList* list, const char* name, Int_t rebin)
{
  THStack* s = GetStack(list, name, rebin);
  if (!s) return 0;
  
  TIter next(s->GetHists());
  TH1*  hist = 0;
  while ((hist = static_cast<TH1*>(next()))) 
    p->Add(hist);
  return s;
}

/** 
 * Build up a centrality legend 
 * 
 * @param c Centrality axis 
 * 
 * @ingroup pwg2_forward_scripts
 */
void
BuildCentLegend(const TAxis* c)
{
  if (!c) return;
  TLegend* l = new TLegend(1.2*gPad->GetLeftMargin(), 
			   .55, .35, 1-gPad->GetTopMargin());
  l->SetFillColor(0);
  l->SetFillStyle(0);
  l->SetBorderSize(0);
  l->SetTextFont(22);
  l->SetHeader("Centralities");
  l->SetTextFont(132);
  
  Int_t    nCol     = gStyle->GetNumberOfColors();
  for (Int_t i = 0; i < c->GetNbins(); i++) {
    UShort_t low    = c->GetBinLowEdge(i+1);
    UShort_t high   = c->GetBinUpEdge(i+1);
    Float_t  fc     = (low+double(high-low)/2) / 100;
    Int_t    icol   = TMath::Min(nCol-1,int(fc * nCol + .5));
    Int_t    col    = gStyle->GetColorPalette(icol);
    TLegendEntry* e = l->AddEntry(Form("dummy%02d", i), 
				  Form("%3d%% - %3d%%", low, high), "p");
    e->SetMarkerColor(col);
				  
  }
  l->Draw();
}

/** 
 * Build a legend.  Histograms a filtered for the same title 
 * 
 * @param stack Stack of histograms 
 * @param c     Centrality axis.  If present, markers are black 
 * 
 * @ingroup pwg2_forward_scripts
 */
void
BuildLegend(const THStack* stack, const TAxis* c)
{
  Double_t x1 = .75, x2 = 1-gPad->GetRightMargin();
  Double_t y1 = .8,  y2 = 1-gPad->GetTopMargin();
  if (!c) { 
    // PP 
    x1 = .4; 
    y1 = .4;
    x2 = .8;
    y2 = .6;
  }
  TLegend* l = new TLegend(x1, y1, x2, y2, "");
  l->SetFillColor(0);
  l->SetFillStyle(0);
  l->SetBorderSize(0);
  l->SetTextFont(132);
  
  TObjArray seen;
  seen.SetOwner();
  TIter next(stack->GetHists());
  TH1* hist = 0;
  while ((hist = static_cast<TH1*>(next()))) { 
    TString n(hist->GetTitle());
    if (n.Contains("mirrored")) continue;
    if (seen.FindObject(n.Data())) continue;
    TObjString* ns = new TObjString(n);
    ns->SetUniqueID(((hist->GetMarkerStyle() & 0xFFFF) << 16) |
		    ((hist->GetMarkerColor() & 0xFFFF) <<  0));
    seen.Add(ns);
  }
  seen.Sort();
  TIter nextu(&seen);
  TObject* s = 0;
  Int_t i = 0;
  while ((s = nextu())) { 
    TLegendEntry* dd = l->AddEntry(Form("data%2d", i++), 
				   s->GetName(), "p");
    Int_t style = (s->GetUniqueID() >> 16) & 0xFFFF;
    Int_t color = (s->GetUniqueID() >>  0) & 0xFFFF;
    dd->SetLineColor(kBlack);
    if (c) dd->SetMarkerColor(kBlack);
    else   dd->SetMarkerColor(color);
    dd->SetMarkerStyle(style);
    Int_t st = dd->GetMarkerStyle();
    if (st == 27 || st == 28 || st == 29 || st == 30 || st == 33 || st == 34) 
      dd->SetMarkerSize(1.4*dd->GetMarkerSize());
  }
  // TLegendEntry* sep = l->AddEntry("s", "_", "h");
  // sep->SetTextSize(0.01);
  TLine* sep = new TLine(0,0,1,1);
  sep->SetLineWidth(1);
  sep->DrawLineNDC(x1+.02, y1-.005, x2-.03, y1-.005);
  // sep->Draw();
  
  TLegend* l2 = new TLegend(x1, y1-.005, x2, y1-.15, "");
  l2->SetFillColor(0);
  l2->SetFillStyle(0);
  l2->SetBorderSize(0);
  l2->SetTextFont(132);
  
  TLegendEntry* d1 = l2->AddEntry("d1", "Data", "p");
  d1->SetLineColor(kBlack);
  d1->SetMarkerColor(kBlack);
  d1->SetMarkerStyle(20);
  TLegendEntry* d2 = l2->AddEntry("d2", "Mirrored data", "p");
  d2->SetLineColor(kBlack);
  d2->SetMarkerColor(kBlack);
  d2->SetMarkerStyle(24);
  
  l->Draw();
  l2->Draw();
}

/** 
 * Add additional information
 *  
 * @param forward  List of info
 * @param prelim   Preliminary mark 
 */
void
AddInformation(TList* forward, bool prelim=true)
{
  Double_t x  = .12;
  Double_t y  = .95;
  Double_t ts = 0.05;
  if (prelim) {
    TLatex* wip = new TLatex(x, y, "Work in progress");
    wip->SetNDC();
    wip->SetTextColor(TColor::GetColor(234,26,46));
    wip->SetTextAlign(13);
    wip->SetTextFont(132);
    wip->SetTextSize(ts);
    wip->Draw();
    y -= ts;
  }

  TObject* sNN = forward->FindObject("sNN");
  TObject* sys = forward->FindObject("sys");
  TObject* trg = forward->FindObject("trigger");
  TObject* vtx = forward->FindObject("vtxAxis");
  TObject* sch = forward->FindObject("scheme");

  TString t(sys->GetTitle());
  Bool_t isPP = t == "pp";

  
  TString s = Form("%s @ #sqrt{s%s}=", 
		   sys->GetTitle(),
		   (isPP ? "" : "_{NN}"));
  Int_t cms = sNN->GetUniqueID();
  if (cms > 1000) s.Append(Form("%5.2fTeV", float(cms)/1000));
  else            s.Append(Form("%3dGeV", cms));
  s.Append(Form(", %s", trg->GetTitle()));

  // if (isPP) { x = .3; y = .3; }
  if (isPP) { 
    x  = 1-gPad->GetRightMargin()-.01;
    y  = 1-gPad->GetTopMargin()-.01;
    ts = .04;
  }
  TLatex* ltx = new TLatex(x, y, s.Data());
  ltx->SetNDC();
  ltx->SetTextColor(TColor::GetColor(41,73,156));
  ltx->SetTextAlign((isPP ? 33 : 13));
  ltx->SetTextFont(132);
  ltx->SetTextSize(ts);
  ltx->Draw();
  y -= ts;
  ltx->DrawLatex(x, y, vtx->GetTitle());
  y -= ts;
  ltx->DrawLatex(x, y, sch->GetTitle());
}  

/** 
 * A function (double Gaussian)
 * 
 * @param xp Independent variables
 * @param pp Parameters 
 * 
 * @return Value of function
 * 
 * @ingroup pwg2_forward_scripts
 */
Double_t myFunc(Double_t* xp, Double_t* pp)
{
  Double_t x  = xp[0];
  Double_t a1 = pp[0];
  Double_t a2 = pp[1];
  Double_t s1 = pp[2];
  Double_t s2 = pp[3];
  return a1*(TMath::Gaus(x, 0, s1) - a2 * TMath::Gaus(x, 0, s2));
}

/** 
 * Make systematic error band 
 * 
 * @param cen     Central result
 * @param fwd     Forward result
 * @param sysErr  Systematic error (fractional)
 * 
 * @return 
 * 
 * @ingroup pwg2_forward_scripts
 */
TH1* 
MakeSysError(const TH1* cen, const TH1* fwd, Double_t sysErr=0.7)
{
  TH1* tmp = static_cast<TH1*>(fwd->Clone("dndetaFitted"));
  tmp->SetMarkerStyle(0);
  tmp->SetFillColor(kGray);
  tmp->SetFillStyle(3001);
  tmp->SetDirectory(0);
  Double_t xlow  = 100;
  Double_t xhigh = 0;
  for (Int_t i = 1; i <= tmp->GetNbinsX(); i++) {
    Double_t cc = cen->GetBinContent(i);
    Double_t cf = fwd->GetBinContent(i);
    Double_t ec = fwd->GetBinError(i);
    Double_t ef = fwd->GetBinError(i);
    Double_t nc = cf;
    Double_t ne = ef;
    if (cc < 0.001 && cf < 0.01) continue;
    xlow  = TMath::Min(tmp->GetXaxis()->GetBinLowEdge(i),xlow);
    xhigh = TMath::Max(tmp->GetXaxis()->GetBinUpEdge(i),xhigh);
    if (cc > 0.001) {
      nc = cc;
      ne = ec;
      if (cf > 0.001) {
	nc  = (cf + cc) / 2;
	ne  = TMath::Sqrt(ec*ec + ef*ef);
      }
    }
    tmp->SetBinContent(i, nc);
    tmp->SetBinError(i, ne);
  }
  TF1* tmpf  = new TF1("tmpf", "gaus", xlow, xhigh);
  tmp->Fit(tmpf, "NQ", "");
  
  TF1* fit = new TF1("f", myFunc, xlow, xhigh, 4);
  fit->SetParNames("a_{1}", "a_{2}", "#sigma_{1}", "#sigma_{2}");
  fit->SetParameters(tmpf->GetParameter(0), 
		     .2, 
		     tmpf->GetParameter(2), 
		     tmpf->GetParameter(2)/4);
  tmp->Fit(fit,"0WQ","");
  for (Int_t i = 1; i <= tmp->GetNbinsX(); i++) {
    Double_t tc = tmp->GetBinContent(i);
    if (tc < 0.01) continue;
    Double_t x = tmp->GetXaxis()->GetBinCenter(i);
    Double_t y = fit->Eval(x);
    tmp->SetBinContent(i, y);
    tmp->SetBinError(i,sysErr*y);
  }
  return tmp;
}

/** 
 * Function to draw the results from forward_dndeta.root file 
 * 
 * @param what     What to draw 
 * @param rebin    Rebinnig.  Note, the data must be present in the file
 * @param filename File to open and draw stuff from >
 * 
 * @ingroup pwg2_forward_scripts
 */
void
SimpledNdeta(Int_t what=0x5, 
	     Int_t rebin=5, const char* filename="forward_dndeta.root")
{
  TFile* file = TFile::Open(filename, "READ");
  if (!file) { 
    Error("SimpledNdeta", "File %s not found", filename);
    return;
  }

  TList* forward = static_cast<TList*>(file->Get("ForwardResults"));
  TList* central = static_cast<TList*>(file->Get("CentralResults"));
  TList* mctruth = static_cast<TList*>(file->Get("MCTruthResults"));

  THStack* all = new THStack("all", "All");
  if (what & 0x1) AddStack(all, forward, "dndeta",      rebin);
  if (what & 0x2) AddStack(all, forward, "dndetaMC",    rebin);
  if (what & 0x1) AddStack(all, central, "dndeta",      rebin);
  if (what & 0x2) AddStack(all, central, "dndetaMC",    rebin);
  if (what & 0x4) AddStack(all, mctruth, "dndeta",      rebin);
  
  TH1* tmp = 0;
  if (what & 0x1) { 
    Double_t sysErr = 0.07;
    TH1* fwd = GetHist(forward, "dndetaForward", rebin);
    TH1* cen = GetHist(central, "dndetaCentral", rebin);
    tmp = MakeSysError(cen, fwd, sysErr);
    // fit = tmpf;
  }
  
  TAxis* centAxis = static_cast<TAxis*>(forward->FindObject("centAxis"));
  Bool_t isPP = centAxis == 0;

  gStyle->SetPalette(1);
  gStyle->SetOptTitle(0);
  TCanvas* c = new TCanvas("dndeta", "dN/deta results", 900, 600);
  c->SetFillColor(0);
  c->SetFillStyle(0);
  c->SetBorderSize(0);
  c->SetBorderMode(0);
  c->SetTopMargin(0.03);
  c->SetRightMargin(0.03);
  
  all->Draw("nostack");
  TAxis* xa = all->GetHistogram()->GetXaxis();
  xa->SetTitleFont(132);
  xa->SetLabelFont(132);
  xa->SetTitle("#eta");
  TAxis* ya = all->GetHistogram()->GetYaxis();
  ya->SetTitleFont(132);
  ya->SetLabelFont(132);
  ya->SetTitle("dN_{ch}/d#eta");

  if (tmp) { 
    tmp->Draw("e5 same");
    all->Draw("same nostack");
  }
  

  // if (fit) fit->Draw("same");
  // if (tmp) tmp->Draw("same");

  BuildCentLegend(centAxis);
  BuildLegend(all, centAxis);

  AddInformation(forward);

  c->SaveAs("dndeta_simple.C");
  c->SaveAs("dndeta_simple.png");
  c->SaveAs("dndeta_simple.root");
}
//
// EOF
//
