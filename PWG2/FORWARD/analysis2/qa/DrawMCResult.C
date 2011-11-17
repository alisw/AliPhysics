/**
 * @file   DrawMCResult.C
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Thu Jul  7 10:57:01 2011
 * 
 * @brief  Script to draw steps (deprecated version - use DrawSteps.C)
 * 
 * 
 * @deprecated Use DrawSteps.C instead
 * @ingroup pwg2_forward_scripts_qa
 */
/** 
 * 
 * 
 * @param forward Forward folder
 * @param sub     Sub-folder to get
 * @param name    Name of stack
 * 
 * @return Stack, or null
 *
 * @deprecated Use DrawSteps.C instead
 * @ingroup pwg2_forward_scripts_qa
 */
THStack*
GetStack(const TList& forward,  const char* sub, const char* name)
{
  TList* lsub = static_cast<TList*>(forward.FindObject(sub));
  if (!lsub) { 
    Warning("GetStack", "Sub list %s not found in %s", sub, forward.GetName());
    return 0;
  }
  THStack* ret = static_cast<THStack*>(lsub->FindObject(name));
  if (!ret) 
    Warning("GetStack" "Stack %s not found in %s", name, sub);
  return ret;
}

/** 
 * Rebin a histogram
 * 
 * @param h      Histogram
 * @param rebin  Rebinning factor
 * 
 * @return Pointer to histogram
 *
 * @deprecated Use DrawSteps.C instead
 * @ingroup pwg2_forward_scripts_qa
 */
TH1* 
Rebin(TH1* h, Int_t rebin)
{
  if (rebin <= 1) return h;
  h->Rebin(rebin);
  h->Scale(1. / rebin);
  return h;
}

/** 
 * Take ratio of two histograms
 * 
 * @param h1  Numerator
 * @param h2  Denominator
 * 
 * @return Newly allocated histogram containg ratio
 *
 * @deprecated Use DrawSteps.C instead
 * @ingroup pwg2_forward_scripts_qa
 */
TH1*
Ratio(const TH1* h1, const TH1* h2)
{
  if (!h1) return;
  if (!h2) return;
  
  TH1* copy = static_cast<TH1*>(h2->Clone("tmp"));
  copy->SetName(Form("%s_%s", h2->GetName(), h1->GetName()));
  copy->SetTitle(Form("%s/%s", h2->GetTitle(), h1->GetTitle()));
  copy->SetDirectory(0);
  copy->Divide(h1);

  return copy;
}
/** 
 * Take ratio of histograms in stacks
 * 
 * @param r   Return stack of ratios
 * @param h1  Numerators
 * @param h2  Denominators
 * 
 * @return How many histograms in the return stack
 *
 * @deprecated Use DrawSteps.C instead
 * @ingroup pwg2_forward_scripts_qa
 */
Int_t 
Ratio(THStack* r, const THStack* h1, const THStack* h2)
{
  if (!h1) return 0;
  if (!h2) return 0;

  int n1 = h1->GetHists()->GetEntries();
  int n2 = h2->GetHists()->GetEntries();
  int nH = 0;
  for (int i = 0; i < n1 && i < n2; i++) { 
    TH1* hh1 = static_cast<TH1*>(h1->GetHists()->At(i));
    TH1* hh2 = static_cast<TH1*>(h2->GetHists()->At(i));
    TH1* h   = Ratio(hh1, hh2);
    if (!h) continue;
    nH++;
    r->Add(h);
  }
  return nH;
}
/** 
 * Draw MC results
 * 
 * @param filename  Input file name
 * @param rebin     Rebinning factor
 * @param ratios    Whether to show ratios
 *
 * @deprecated Use DrawSteps.C instead
 * @ingroup pwg2_forward_scripts_qa
 */
void
DrawMCResult(const char* filename="forward.root", Int_t rebin=1,
	     Bool_t ratios=true)
{
  gStyle->SetPalette(1);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);

  TFile* file = TFile::Open(filename, "READ");
  if (!file) { 
    Error("DrawMCResult", "failed to open %s", filename);
    return;
  }

  TList* forward = static_cast<TList*>(file->Get("Forward"));
  if (!forward) { 
    Error("DrawMCResult", "List Forward not found in %s", filename);
    return;
  }
  THStack* res    = GetStack(*forward, "ringResults", "all");
  THStack* mcRes  = GetStack(*forward, "mcRingResults", "all");
  THStack* deltas = GetStack(*forward, "fmdSharingFilter", "sums");
  THStack* nchs   = GetStack(*forward, "fmdDensityCalculator", "sums");
  THStack* prims  = GetStack(*forward, "fmdCorrector", "sums");
  
  TH1* sumEta = static_cast<TH1*>(forward->FindObject("mcSumEta"));
  if (!sumEta) { 
    Warning("DrawMCResults", "mcSumEta not found in Forward");
  }


  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleStyle(0);
  // gStyle->SetTitleColor(kBlack);


  TCanvas* c = new TCanvas("c", "C", 900, 700);
  c->SetFillColor(0);
  c->SetBorderSize(0);
  c->SetTopMargin(0.05);
  c->SetRightMargin(0.05);

  Double_t y1 = (ratios ? .3 : 0);
  TPad* p1 = new TPad("p1", "p1", 0, y1, 1, 1, 0, 0, 0); 
  p1->SetBottomMargin(ratios ? 0.01 : .10);
  p1->SetFillColor(0);
  p1->SetBorderSize(0);
  p1->SetTopMargin(0.05);
  p1->SetRightMargin(0.05);
  p1->Draw();
  p1->cd();

  THStack* all = new THStack("all", "Analysis steps");
  if (res) {
    res->SetTitle("dN_{ch}/d#eta");
    TIter next(res->GetHists());
    TH1* h = 0;
    while ((h = static_cast<TH1*>(next()))) all->Add(Rebin(h,rebin));
  }
  if (mcRes) {
    mcRes->SetTitle("Track-Refs");
    TIter next(mcRes->GetHists());
    TH1* h = 0;
    while ((h = static_cast<TH1*>(next()))) all->Add(Rebin(h,rebin));
  }
  if (deltas) {
    deltas->SetTitle("#sum_{} #Delta/#Delta_{mip}");
    TIter next(deltas->GetHists());
    TH1* h = 0;
    while ((h = static_cast<TH1*>(next()))) { 
      h->SetMarkerStyle(25);
      all->Add(Rebin(h,rebin));
    }
  }
  if (nchs) {
    nchs->SetTitle("#sum_{} N_{ch,incl}");
    TIter next(nchs->GetHists());
    TH1* h = 0;
    while ((h = static_cast<TH1*>(next()))) { 
      h->SetMarkerStyle(21);
      all->Add(Rebin(h,rebin));
    }
  }
  if (prims) {
    prims->SetTitle("#sum_{} N_{ch,primary}");
    TIter next(prims->GetHists());
    TH1* h = 0;
    while ((h = static_cast<TH1*>(next()))) { 
      h->SetMarkerStyle(22);
      all->Add(Rebin(h,rebin));
    }
  }
  if (sumEta) all->Add(sumEta);
  all->Draw("nostack");
  all->GetHistogram()->SetXTitle("#eta");
  all->GetHistogram()->SetYTitle("signal");
  all->GetHistogram()->GetXaxis()->SetLabelFont(132);
  all->GetHistogram()->GetXaxis()->SetTitleFont(132);
  all->GetHistogram()->GetYaxis()->SetLabelFont(132);
  all->GetHistogram()->GetYaxis()->SetTitleFont(132);
  c->SetGridx();

  TLegend* l = new TLegend(.33, .2, .53, .9);
  TLegendEntry* e = 0;
  l->SetFillColor(0);
  l->SetFillStyle(0);
  l->SetBorderSize(0);
  l->SetNColumns(1);
  l->SetTextFont(132);
  Int_t i = 0;
  if (res) { 
    TIter next(res->GetHists());
    TH1* h = 0;
    while ((h = static_cast<TH1*>(next()))) {
      e = l->AddEntry(Form("dummy%02d", i++),h->GetTitle(),"pl");
      e->SetMarkerStyle(20);
      e->SetMarkerColor(h->GetMarkerColor());
    }
    e = l->AddEntry(Form("dummy%02d", i++),res->GetTitle(),"pl");
    e->SetMarkerStyle(20);
    e->SetMarkerColor(kBlack);
  }
  if (deltas) { 
    e = l->AddEntry(Form("dummy%02d", i++), deltas->GetTitle(),"pl");
    TH1* h = static_cast<TH1*>(deltas->GetHists()->At(0));
    e->SetMarkerStyle(h->GetMarkerStyle());
    e->SetMarkerColor(kBlack);
  }
  if (nchs) { 
    e = l->AddEntry(Form("dummy%02d",i++),nchs->GetTitle(),"pl");
    TH1* h = static_cast<TH1*>(nchs->GetHists()->At(0));
    e->SetMarkerStyle(h->GetMarkerStyle());
    e->SetMarkerColor(kBlack);
  }
  if (prims) { 
    e = l->AddEntry(Form("dummy%02d", i++), prims->GetTitle(),"pl");
    TH1* h = static_cast<TH1*>(prims->GetHists()->At(0));
    e->SetMarkerStyle(h->GetMarkerStyle());
    e->SetMarkerColor(kBlack);
  }

  if (mcRes) { 
    e = l->AddEntry(Form("dummy%02d", i++), mcRes->GetTitle(), "pl");
    TH1* h = static_cast<TH1*>(mcRes->GetHists()->At(0));
    e->SetMarkerStyle(h->GetMarkerStyle());
    e->SetMarkerColor(kBlack);
  }

  if (sumEta) l->AddEntry(sumEta);
  l->Draw();


  if (!ratios) return;

  c->cd();
  TPad* p2 = new TPad("p2", "p2", 0, 0, 1, y1, 0, 0, 0); 
  p2->SetTopMargin(0);
  p2->SetFillColor(0);
  p2->SetBorderSize(0);
  p2->SetRightMargin(0.05);
  p2->Draw();
  p2->cd();

  THStack* rs = new THStack("ratios", "Ratios");
  Int_t nDN = Ratio(rs, deltas, nchs);
  Int_t nNR = Ratio(rs, nchs, res);
  Int_t nRP = Ratio(rs, res, prims);
  rs->Draw("nostack");

  TLegend* ll = new TLegend(.38, .2, .48, .9);
  ll->SetFillColor(0);
  ll->SetFillStyle(0);
  ll->SetBorderSize(0);
  ll->SetNColumns(1);
  ll->SetTextFont(132);
  if (nDN) {
    e = ll->AddEntry("d1", Form("#frac{%s}{%s}", 
					      nchs->GetTitle(), 
					      deltas->GetTitle()), "lp");
    e->SetMarkerStyle(21);
  }
  if (nNR) {
    e = ll->AddEntry("d2", Form("#frac{%s}{%s}", res->GetTitle(), 
				nchs->GetTitle()), "lp");
    e->SetMarkerStyle(20);
  }
  if (nRP) {
    e = ll->AddEntry("d3", Form("#frac{%s}{%s}", prims->GetTitle(), 
			    res->GetTitle()), "lp");
    e->SetMarkerStyle(22);
  }
  ll->Draw();
  

}
 
//
// EOF
//
