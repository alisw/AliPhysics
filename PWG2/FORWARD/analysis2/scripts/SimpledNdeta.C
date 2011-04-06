#ifndef __CINT__
# include <TStyle.h>
# include <TFile.h>
# include <TList.h>
# include <TH1.h>
# include <THStack.h>
# include <TLegend.h>
# include <TLegendEntry.h>
# include <TLine.h>
# include <TString.h>
# include <TCanvas.h>
# include <TError.h>
#else
class THStack;
class TAxis;
#endif

/**
 * A simple script to draw results from MakedNdeta.C (or similar)
 * 
 */
/** 
 * Get a stack from the passed list 
 * 
 * @param list   List to get the stack from 
 * @param name   Name of stack 
 * @param rebin  Optional rebinning - must exists in list 
 * 
 * @return Stack or null
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
 * Add histograms from one stack to another 
 * 
 * @param p      Parent stack to add to 
 * @param list   List to look for the stack in 
 * @param name   Name of stack to add 
 * @param rebin  Optional rebinning - must exists in list 
 * 
 * @return Added stack or null
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
 */
void
BuildLegend(const THStack* stack, const TAxis* c)
{
  TLegend* l = new TLegend(.75, .8, 1-gPad->GetRightMargin(), 
			   1-gPad->GetTopMargin(), "");
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
  }
  // TLegendEntry* sep = l->AddEntry("s", "_", "h");
  // sep->SetTextSize(0.01);
  TLine* sep = new TLine(0,0,1,1);
  sep->SetLineWidth(1);
  sep->DrawLineNDC(.77, .795, 1-gPad->GetRightMargin()-.03, .795);
  // sep->Draw();
  
  TLegend* l2 = new TLegend(.75, .7, 1-gPad->GetRightMargin(), .79, "");
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
 * Function to draw the results from forward_dndeta.root file 
 * 
 * @param rebin    Rebinnig.  Note, the data must be present in the file
 * @param filename File to open and draw stuff from >
 */
void
SimpledNdeta(Int_t rebin=5, const char* filename="forward_dndeta.root")
{
  TFile* file = TFile::Open(filename, "READ");
  if (!file) { 
    Error("SimpledNdeta", "File %s not found", filename);
    return;
  }

  TList* forward = static_cast<TList*>(file->Get("ForwardResults"));
  TList* central = static_cast<TList*>(file->Get("CentralResults"));

  THStack* all = new THStack("all", "All");
  AddStack(all, forward, "dndeta",      rebin);
  AddStack(all, forward, "dndetaMC",    rebin);
  AddStack(all, forward, "dndetaTruth", rebin);
  AddStack(all, central, "dndeta",      rebin);
  
  TAxis* centAxis = static_cast<TAxis*>(forward->FindObject("centAxis"));

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

  BuildCentLegend(centAxis);
  BuildLegend(all, centAxis);

  c->SaveAs("dndeta.C");
}
//
// EOF
//
