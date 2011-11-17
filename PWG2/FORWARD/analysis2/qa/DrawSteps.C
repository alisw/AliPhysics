/**
 * @file   DrawSteps.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Thu Nov 17 11:34:01 2011
 * 
 * @brief  
 * 
 * 
 * @defgroup pwg2_forward_scripts_qa Quality Assurance scripts
 * @ingroup pwg2_forward_scripts
 */
/** 
 * Get a stack 
 * 
 * @param forward   Input list
 * @param sub       Sub-list
 * @param name      Name of stack
 * 
 * @return A stack or null
 * 
 * @ingroup pwg2_forward_scripts_qa
 */
THStack*
GetStack(const TList& forward,  const char* sub, const char* name)
{
  TList* lsub = &forward;

  if (sub && sub[0] != '\0') 
    lsub = static_cast<TList*>(forward.FindObject(sub));

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
 * @return Histogram
 * 
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
 * Ratio of two histograms 
 * 
 * @param h1 numerator
 * @param h2 denominator
 * 
 * @return Ratio
 * 
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
 * Ratio all histograms in stacks 
 * 
 * @param r  Result
 * @param h1 Numerators
 * @param h2 Denominators 
 * 
 * @return Number of histograms 
  * 
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
 * Add a histogram to the all stack
 * 
 * @param all         Stack
 * @param h           Histogram
 * @param singleStep  Showing individual steps?
 * 
 * @ingroup pwg2_forward_scripts_qa
*/
void
AddToAll(THStack* all, const TH1* h, Bool_t singleStep)
{
  TH1* copy = static_cast<TH1*>(h->Clone(Form("%s_copy", h->GetName())));
  copy->SetDirectory(0);
  if (singleStep) { 
    copy->SetMarkerColor(kGray);
    copy->SetLineColor(kGray);
  }
  all->Add(copy);
}

/** 
 * Dim an entry
 * 
 * @param thisId  This step
 * @param step    Current step
 * @param e       Entry in legend 
 * 
 * @ingroup pwg2_forward_scripts_qa
 */
void
DimEntry(Int_t thisId, Int_t step, TLegendEntry* e)
{
  
  Int_t col = (thisId == step || step <= 0) ? kBlack : kGray;
  e->SetMarkerColor(col);
  e->SetLineColor(col);
  e->SetTextColor(col);
}

/** 
 * Draw a step
 * 
 * @param deltas   From energy loss
 * @param nchs     After 2nd correction
 * @param prims    Primaries
 * @param dndeta   Result 
 * @param step     Step number 
 * 
 * @ingroup pwg2_forward_scripts_qa
 */
void
DrawStep(THStack* deltas, THStack* nchs, THStack* prims, 
	 TH1*     dndeta, Int_t step)
{
  THStack* all = new THStack("all", "Analysis steps");
  if (step > 0) all->SetTitle(Form("Step %d", step));

  if (deltas) {
    deltas->SetTitle("#sum_{} #Delta/#Delta_{mip}");
    TIter next(deltas->GetHists());
    TH1* h = 0;
    while ((h = static_cast<TH1*>(next()))) { 
      h->SetMarkerStyle(25);
      // Info("DrawStep", "Adding %s", h->GetName());
      AddToAll(all, h, step>0);
    }
  }
  if (nchs) {
    nchs->SetTitle("#sum_{} N_{ch,incl}");
    TIter next(nchs->GetHists());
    TH1* h = 0;
    while ((h = static_cast<TH1*>(next()))) { 
      h->SetMarkerStyle(21);
      // Info("DrawStep", "Adding %s", h->GetName());
      AddToAll(all, h, step>0);
    }
  }
  if (prims) {
    prims->SetTitle("#sum_{} N_{ch,primary}");
    TIter next(prims->GetHists());
    TH1* h = 0;
    while ((h = static_cast<TH1*>(next()))) { 
      h->SetMarkerStyle(22);
      // Info("DrawStep", "Adding %s", h->GetName());
      AddToAll(all, h, step>0);
    }
  }
  if (dndeta) {
    dndeta->SetTitle("1/N dN_{ch}/d#eta");
    dndeta->SetMarkerStyle(20);
    dndeta->SetMarkerColor(kBlack);
    // Info("DrawStep", "Adding %s", dndeta->GetName());
    AddToAll(all, dndeta, step>0);
  }

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
  if (deltas) { 
    TIter next(deltas->GetHists());		
    TH1*  h = 0;
    while ((h = static_cast<TH1*>(next()))) {
      e = l->AddEntry(Form("dummy%02d", i++),h->GetTitle(),"pl");
      e->SetMarkerStyle(20);
      e->SetMarkerColor(h->GetMarkerColor());
    }
    e = l->AddEntry(Form("dummy%02d", i++), deltas->GetTitle(),"pl");
    TH1* h = static_cast<TH1*>(deltas->GetHists()->At(0));
    e->SetMarkerStyle(h->GetMarkerStyle());
    DimEntry(1, step, e);
  }
  if (nchs) { 
    e = l->AddEntry(Form("dummy%02d",i++),nchs->GetTitle(),"pl");
    TH1* h = static_cast<TH1*>(nchs->GetHists()->At(0));
    e->SetMarkerStyle(h->GetMarkerStyle());
    DimEntry(2, step, e);
  }
  if (prims) { 
    e = l->AddEntry(Form("dummy%02d", i++), prims->GetTitle(),"pl");
    TH1* h = static_cast<TH1*>(prims->GetHists()->At(0));
    e->SetMarkerStyle(h->GetMarkerStyle());
    DimEntry(3, step, e);
  }
  if (dndeta) { 
    e = l->AddEntry(Form("dummy%02d", i++), dndeta->GetTitle(),"pl");
    e->SetMarkerStyle(dndeta->GetMarkerStyle());
    DimEntry(4, step, e);
  }
  l->Draw();

  TString what;
  if (step > 0) {
    switch (step) { 
    case 1: 
      deltas->Draw("same nostack"); 
      what = "After merging";
      break;
    case 2: 
      nchs->Draw("same nostack"); 
      what = "After particle counting";
      break;
    case 3: 
      prims->Draw("same nostack"); 
      what = "After corrections";
      break;
    case 4: 
      dndeta->Draw("same"); 
      what = "After normalisation";
      break;
    default: 
      Error("DrawSteps", "Unknown step: %d (must be in 1-4)");
      break;
    }
  }
  TLatex* ltx = new TLatex(.95, .85, what);
  ltx->SetNDC();
  ltx->SetTextSize(.07);
  ltx->SetTextAlign(33);
  ltx->SetTextFont(132);
  ltx->Draw();
}

/** 
 * Draw steps
 * 
 * @param filename Input file 
 * @param single   Whether to show individial steps 
 * 
 * @ingroup pwg2_forward_scripts_qa
 */
void DrawSteps(const char* filename="forward.root", Bool_t single=true)
{
  gStyle->SetPalette(1);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);

  TFile* file = TFile::Open(filename, "READ");
  if (!file) { 
    Error("DrawMCResult", "failed to open %s", filename);
    return;
  }
  const char* fname2 = "forward_dndeta.root";
  TFile* file2 = TFile::Open(fname2, "READ");
  if (!file2) { 
    Error("DrawSteps", "File %s not found", fname2);
  }

  TList* forward = static_cast<TList*>(file->Get("Forward"));
  if (!forward) { 
    Error("DrawMCResult", "List Forward not found in %s", filename);
    return;
  }
  TList* forwardRes = (file2 ? 
		       static_cast<TList*>(file2->Get("ForwardResults")) :
		       0);
  TList* forwardAll = (forwardRes ? 
		       static_cast<TList*>(forwardRes->FindObject("all")) :
		       0);
		       
  
  // THStack* res    = GetStack(*forward, "ringResults", "all");
  // THStack* mcRes  = GetStack(*forward, "mcRingResults", "all");
  THStack* deltas = GetStack(*forward, "fmdSharingFilter", "sums");
  THStack* nchs   = GetStack(*forward, "fmdDensityCalculator", "sums");
  THStack* prims  = GetStack(*forward, "fmdCorrector", "sums");
  TH1*     dndeta = (forwardAll ? 
		     static_cast<TH1*>(forwardAll->FindObject("dndetaForward")):
		     0);

  Info("DrawSteps", "Got steps deltas=%p, nchs=%p, prims=%p, dndeta=%p",
       deltas, nchs, prims, dndeta);


  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleX(.7);
  gStyle->SetTitleY(.95);
  gStyle->SetTitleH(.1);
  gStyle->SetTitleW(.25);
  gStyle->SetOptTitle(1);
  // gStyle->SetTitleColor(kBlack);


  
  if (!single) { 
    TCanvas* c = new TCanvas("c", "C", 900, 700);
    c->SetFillColor(0);
    c->SetBorderSize(0);
    c->SetTopMargin(0.05);
    c->SetRightMargin(0.05);

    DrawStep(deltas, nchs, prims, dndeta, 0);
    c->SaveAs("steps_all.png");
    return;
  }
  Int_t nSteps = 0;
  if (deltas) nSteps++;
  if (nchs)   nSteps++;
  if (prims)  nSteps++;
  if (dndeta) nSteps++;

  Int_t w = (nSteps >= 4 ? 1100 :  700);
  Int_t h = (nSteps >= 4 ? 800  : 1100);

  TCanvas* c = new TCanvas("c", "C", w, h);
  c->SetFillColor(0);
  c->SetBorderSize(0);
  c->SetTopMargin(0.05);
  c->SetRightMargin(0.05);

  if (nSteps >= 4) 
    c->Divide(2,(nSteps+1)/2,0,0);
  else 
    c->Divide(1,nSteps,0,0);
  
  for (Int_t i=1; i<=nSteps; i++) { 
    TVirtualPad* p = c->cd(i);
    p->SetFillColor(0);
    p->SetFillStyle(0);
    p->SetBorderSize(0);
    p->SetGridx();
    p->SetGridy();

    DrawStep(deltas, nchs, prims, dndeta, i);
  }
  c->SaveAs("steps_comic.png");
}
//
// EOF
//

    
