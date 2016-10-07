/**
 * @file   YieldExpectations.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Fri Oct  7 12:01:42 2016
 * 
 * @brief  Do back of the envelope calculation of effect of reweighing
 * 
 * @ingroup pwglf_forward_tracklets
 * 
 */
#include "AliTrackletAODUtils.C"
#ifndef __CINT__
# include <vector>
# include <TFile.h>
# include <TError.h>
# include <TCollection.h>
# include <TString.h>
# include <TH1.h>
# include <TH2.h>
# include <THStack.h>
# include <TClass.h>
# include <TBrowser.h>
# include <TCanvas.h>
# include <TLegend.h>
# include <TLegendEntry.h>
# include <TStyle.h>
#else
class TFile;
class TClass;
class TH1;
class TH2;
class TDirectory;
class THStack;
class TCollection;
class TBrowser;
class TCanvas;
class TVirtualPad;
class TLegend;
#endif
/**
 * 
 * 
 * @ingroup pwglf_forward_tracklets
 */
struct Calculations : public AliTrackletAODUtils 
{
  struct One
  {
    THStack* fStack;
    TH1*     fK0S;
    TH1*     fKpm;
    TH1*     fLambda;
    TH1*     fXi;
    TH1*     fSigma;
    TH1*     fOther;
    TH1*     fYields;
    
    One(const char* filename,
	const char* title,
	Color_t     color,
	Double_t    offset,
	Double_t    width=.15,
	Double_t    c1=0,
	Double_t    c2=0)
      : fStack(0),
	fK0S(0),
	fKpm(0),
	fLambda(0),
	fXi(0),
	fSigma(0),
	fOther(0),
	fYields(0)
    {
      TString    cntN = CentName(c1,c2);
      TFile*     file = OpenFile(filename);
      Container* top  = GetC(file, "MidRapidityMCResults");
      Container* cent = GetC(top,  cntN);
      Container* gen  = GetC(cent, "generated");
      Container* mix  = GetC(gen,  "mix");

      fStack  = new THStack("dists", title);
      fYields = new TH1D("yields", title, 6, 0, 6);
      fYields->SetMarkerSize(1.5);
      ModHist(fYields, 0, color, offset, width);
      
      TIter    next(mix);
      TObject* o = 0;
      while ((o = next())) {
	if (!o->IsA()->InheritsFrom(TH1::Class())) continue;
	TH1*    t = static_cast<TH1*>(o);

	TString nme(t->GetName());
	if (!nme.BeginsWith("eta_")) continue;

	nme.ReplaceAll("eta_","");
	Int_t pdg = nme.Atoi();

	TH1** h = 0;
	switch (pdg) {
	case 310:   h = &fK0S;    break;
	case 321:   h = &fKpm;    break;
	case 3122:  h = &fLambda; break;
	case 3322:  h = &fXi;     break; // Old, wrong code
	  // case 3312:  h = &fXi;     break; // New, right code
	case 3112:  h = &fSigma;  break;
	case 3212:  h = &fSigma;  break; // Old, wrong code
	  // case 3222:  h = &fSigma;  break; // New, right code
	default:    h = &fOther;  break;
	}
	if (!h) {
	  Printf("Unknown PDG=%d", pdg);
	  continue;
	}
	if (!*h) {
	  *h = static_cast<TH1*>(t->Clone());
	  (*h)->Reset();
	}
	(*h)->Add(t);
      }
      fOther ->SetTitle("Other");
      fOther ->SetMarkerStyle(20); fOther ->SetMarkerSize(1.2);
      fK0S   ->SetMarkerStyle(21); fK0S   ->SetMarkerSize(1);
      fKpm   ->SetMarkerStyle(25); fKpm   ->SetMarkerSize(1);
      fLambda->SetMarkerStyle(23); fLambda->SetMarkerSize(1.2);
      fXi    ->SetMarkerStyle(22); fXi    ->SetMarkerSize(1.2);
      fSigma ->SetMarkerStyle(26); fSigma ->SetMarkerSize(1.2);
      fStack->Add(fK0S   );
      fStack->Add(fKpm   );
      fStack->Add(fLambda);
      fStack->Add(fXi    );
      fStack->Add(fSigma );
      fStack->Add(fOther );

      if (color == kMagenta+1) {	
	fK0S   ->Scale(1.52233);
	fKpm   ->Scale(1.41178); 
	fLambda->Scale(2.75002); 
	fXi    ->Scale(3.24110); 
	fSigma ->Scale(2.75002); 
      }	fOther ->Scale(1);       
      
      TIter nextH(fStack->GetHists());
      TH1*  hist = 0;
      Int_t bin  = 1;
      while ((hist = static_cast<TH1*>(nextH()))) {
	CalcYield(hist, bin++);
	ModHist(hist, 0, color, offset, width);
      }
    }
    void ModHist(TH1*        h,
		 const char* title,
		 Color_t     color,
		 Double_t    offset,
		 Double_t    width)
    {
      // Printf("%10s offset=%f  width=%f", h->GetTitle(), offset, width);
      h->SetDirectory(0);
      h->SetMarkerColor(color);
      h->SetLineColor(color);
      h->SetFillColor(color);
      h->SetBarOffset(offset);
      h->SetBarWidth(width);

      TString nme(Form("%s%s",h->GetName(),title)); nme.ReplaceAll("-","");
      h->SetName(nme);

      if (h == fYields) return;
      TAxis*   axis  = h->GetXaxis();
      Double_t shift = axis->GetBinWidth(1)*(offset-0.5);
      axis->SetLimits(axis->GetXmin()+shift, axis->GetXmax()+shift);
      if (title) h->SetTitle(Form("%s #minus %s", h->GetTitle(), title));
    }
    void CalcYield(TH1* h, Int_t bin)
    {
      Double_t err;
      Double_t yie = h->IntegralAndError(1,h->GetNbinsX(),err);
      fYields->SetBinContent(bin, yie);
      fYields->SetBinError  (bin, err);
      fYields->GetXaxis()->SetBinLabel(bin, h->GetTitle());
    }
    void AddToStacks(THStack* dists, THStack* yields,
		     TLegend* ld,    TLegend* ly)
    {
      TIter next(fStack->GetHists());
      TH1*  hist = 0;
      while ((hist = static_cast<TH1*>(next()))) {
	dists->Add(hist);
	if (!ld) continue;
	TLegendEntry* e = ld->AddEntry("dummy", hist->GetTitle(), "p");
	e->SetMarkerStyle(hist->GetMarkerStyle());
	e->SetMarkerSize(1.5*hist->GetMarkerSize());
      }
      
      yields->Add(fYields, "hist bar text90");
      TLegendEntry* e = ly->AddEntry("dummy", fYields->GetTitle(), "f");
      e->SetFillColor(fYields->GetFillColor());
      e->SetFillStyle(1001);      
    }      
  };
  void Run()
  {
    One*  stk    = new One("stk.root",   "Reduced",   kRed+1,     0.15);
    One*  expe   = new One("stk.root",   "Expected",  kMagenta+1, 0.35);
    One*  wstk   = new One("wstk.root",  "Reweighed", kBlue+1,    0.55);
    One*  hijing = new One("hijing.root","As-is",     kGreen+1,   0.75);
    One*  ones[] = { stk, expe, wstk, hijing, 0 };
    One** ptr    = ones;

    gStyle->SetPaintTextFormat("5.0f");
    gStyle->SetErrorX(0.2);
    TCanvas* c = new TCanvas("c","c",1200,900);
    c->Divide(2,1);
    TVirtualPad* p = c->GetPad(1);
    TVirtualPad* q = c->GetPad(2);
    p->SetTopMargin(0.01);
    p->SetRightMargin(0.01);
    p->SetLogy();
    p->SetTicks();
    q->SetTopMargin(0.01);
    q->SetRightMargin(0.01);
    q->SetLeftMargin(0.13);
    q->SetLogy();
    q->SetTicks();

    TLegend* le = new TLegend(p->GetLeftMargin()+.01, .75,
			      1-p->GetRightMargin(), .85);
    TLegend* ly = new TLegend(q->GetLeftMargin()+.01, .75,
			      1-q->GetRightMargin(),
			      1-q->GetTopMargin());
    le->SetBorderSize(0);
    le->SetFillStyle(0);
    le->SetNColumns(6);
    ly->SetBorderSize(0);
    ly->SetFillStyle(0);
    
    THStack* dists  = new THStack("dists",  "");
    THStack* yields = new THStack("yields", "");

    while (*ptr) {
      (*ptr)->AddToStacks(dists, yields, (*ptr == stk ? le : 0), ly);
      ptr++;
    }
    p->cd();
    dists->Draw("nostack");
    dists->GetHistogram()->SetXTitle("#it{#eta}");
    dists->GetHistogram()->SetYTitle("d#it{N}_{X}/d#it{#eta}");    
    le->Draw();

    q->cd();
    yields->SetMaximum(2*yields->GetMaximum("nostack"));
    yields->Draw("nostack");
    yields->GetHistogram()->GetXaxis()->SetLabelSize(0.07);
    yields->GetHistogram()->GetYaxis()->SetLabelSize(0);
    yields->GetHistogram()->GetYaxis()->SetTitleOffset(1.2);    
    yields->GetHistogram()->SetYTitle(Form("#int_{-2}^{+2}d#it{#eta} %s",
					   dists->GetHistogram()->GetYaxis()
					   ->GetTitle()));
    ly->Draw();
  }
    
    
};



/** 
 * 
 * 
 * @ingroup pwglf_forward_tracklets
 */
void YieldExpectations()
{
  Calculations c;
  c.Run();
}
