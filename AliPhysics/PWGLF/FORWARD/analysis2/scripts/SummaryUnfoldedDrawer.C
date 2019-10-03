#include "SummaryDrawer.C"
#include <TMultiGraph.h>
#include <TKey.h>
#include <TList.h>

/**
 * Class to draw summary of unfolding
 * 
 * @ingroup pwglf_forward_multdist
 */
struct SummaryUnfoldedDrawer : public SummaryDrawer
{

  SummaryUnfoldedDrawer() : SummaryDrawer() {}
  //____________________________________________________________________
  void Run(const char* fname)
  {
    TString filename(fname);
    TFile* file = TFile::Open(filename,"READ");
    if (!file) { 
      Error("Run", "Failed to open \"%s\"", filename.Data());
      return;
    }
    
    TString pdfName(filename);
    pdfName.ReplaceAll(".root", ".pdf");
    CreateCanvas(pdfName, 0); // what & kLandscape);

    // --- Title page ------------------------------------------------
    DrawTitlePage(file);
    
    // --- Loop over all keys in the file ----------------------------
    TIter next(file->GetListOfKeys());
    TKey* k  = 0;
    while ((k = static_cast<TKey*>(next()))) {
      file->cd();

      TString clName = k->GetClassName();
      TClass* cl     = gROOT->GetClass(clName);
      if (!cl) { 
	// Warning("Run", "Ignoring object %s of unknown type %s", 
	//         k->GetName(), clName.Data());
	continue;
      }
      if (!cl->InheritsFrom(TDirectory::Class())) {
	// Warning("Run", "Ignoring object %s of type %s", 
	//         k->GetName(), clName.Data());
	continue;
      }
      file->cd(k->GetName());
      ProcessType(gDirectory);
    }
    file->cd();
    CloseCanvas();
  }
  //____________________________________________________________________
  void DrawTitlePage(const TDirectory* file)
  {
    fBody->cd();

    TLatex* ltx = new TLatex(.5, .7, "Raw P(#it{N}_{ch}) #rightarrow "
			     "corrected P(#it{N}_{ch})");
    ltx->SetNDC();
    ltx->SetTextSize(0.07);
    ltx->SetTextAlign(22);
    ltx->Draw();

    Double_t y = .6;
    
    Double_t save = fParName->GetTextSize();
    fParName->SetTextSize(0.03);
    fParVal->SetTextSize(0.03);

    TObject* method = GetObject(file, "method");
    TObject* sys    = GetObject(file, "sys");
    TObject* sNN    = GetObject(file, "sNN");
    TObject* trig   = GetObject(file, "trigger");
    Double_t regP;
    Double_t minIpZ;
    Double_t maxIpZ;
    Bool_t   self;
    GetParameter(file, "regParam", regP);
    GetParameter(file, "minIpZ", minIpZ);
    GetParameter(file, "maxIpZ", maxIpZ);
    GetParameter(file, "self", self);

    DrawParameter(y, "Consistency check", self ? "yes" : "no");
    DrawParameter(y, "System", sys->GetTitle());
    DrawParameter(y, "#sqrt{s_{NN}}", sNN->GetTitle());
    DrawParameter(y, "Trigger", trig->GetTitle());
    DrawParameter(y, "Method", method->GetTitle());
    DrawParameter(y, "Reg. param.", Form("%g", regP));
    DrawParameter(y, "IP_{z} range", Form("%+5.2fcm - %+5.2fcm", 
					  minIpZ, maxIpZ));


    PrintCanvas("Title page");
    fParName->SetTextSize(save);
    fParVal->SetTextSize(save);
  }
  //____________________________________________________________________
  void ProcessType(TDirectory* d) 
  { 
    Printf(" ProcessType: %s", d->GetPath());
    // d->ls();
    MakeChapter(d->GetName()); 
    static TRegexp regex("[pm][0-9]d[0-9]*_[pm][0-9]d[0-9]*");

    // --- Loop over bins in this directory --------------------------
    TIter next(d->GetListOfKeys());
    TKey* k  = 0;
    while ((k = static_cast<TKey*>(next()))) {
      d->cd();
      TString clName = k->GetClassName();
      TClass* cl     = gROOT->GetClass(clName);
      if (!cl) { 
	// Warning("ProcessType", "Ignoring object %s of unknown type %s", 
	// k->GetName(), clName.Data());
	continue;
      }
      if (!cl->InheritsFrom(TDirectory::Class())) {
	// Warning("ProcessType", "Ignoring object %s of type %s", 
	//  k->GetName(), clName.Data());
	continue;
      }
      TString n(k->GetName());
      if (n.Index(regex) == kNPOS) { 
	Warning("ProcessType", "Ignoring non-bin directory %s", n.Data());
	continue;
      }
      d->cd(n);
      ProcessBin(gDirectory);
    }
    d->cd();
    DrawResults(d);
  }
  //____________________________________________________________________
  void DrawResults(TDirectory* d) 
  { 
    THStack* c = GetStack(d, "corrected");
    if (!c) {
      Warning("DrawResults", "Stack of corrected results not found!");
      return;
    }
    DrawInPad(fBody, 0, c, "nostack", kLogy);
    c->GetXaxis()->SetTitle("#it{N}_{ch}");
    c->GetYaxis()->SetTitle("P(#it{N}_{ch})");

    TObject* alice = d->Get("alice");
    TObject* cms   = d->Get("cms");
    if (cms)   cms->Draw();
    if (alice) alice->Draw();

    if (alice || cms) { 
      TObject* dummy = 0;
      TLegend* l = new TLegend(0.15, 0.12, .4, .3, "", "NDC");
      l->SetFillColor(0);
      l->SetFillStyle(0);
      l->SetBorderSize(0);
      TLegendEntry* e = l->AddEntry(dummy, "This work", "F");
      e->SetFillColor(kRed+2);
      e->SetFillStyle(1001);
      if (alice) { 
	e = l->AddEntry(dummy, "ALICE", "F");
	e->SetFillColor(kPink+1);
	e->SetFillStyle(1001);
      }
      if (cms) { 
	e = l->AddEntry(dummy, "CMS", "F");
	e->SetFillColor(kGreen+2);
	e->SetFillStyle(1001);
      }
      l->Draw();
    }
    PrintCanvas(" Results");

    THStack* s = GetStack(d, "ratios");
    if (s) {
      fBody->SetLogy(false);
      DrawInPad(fBody, 0, s, "nostack", kGridy);
      s->GetYaxis()->SetTitle("(this - other)/other");
      s->GetXaxis()->SetTitle("#it{N}_{ch}");

      TObject* dummy = 0;
      TLegend* l = new TLegend(0.65, 0.12, .95, .3, "", "NDC");
      TLegendEntry* e = 0;
      l->SetFillColor(0);
      l->SetFillStyle(0);
      l->SetBorderSize(0);
      if (alice) { 
	e = l->AddEntry(dummy, "to ALICE", "F");
	e->SetFillColor(kPink+1);
	e->SetFillStyle(1001);
      }
      if (cms) { 
	e = l->AddEntry(dummy, "to CMS", "F");
	e->SetFillColor(kGreen+2);
	e->SetFillStyle(1001);
      }
      l->Draw();
      if (s->GetMinimum() > -1) s->SetMinimum(-1);
      
      PrintCanvas("  Ratios");
    }
  }
  //____________________________________________________________________
  void ProcessBin(TDirectory* d) 
  {
    Printf("  ProcessBin: %s", d->GetPath());

    TString tmp(d->GetName());
    tmp.ReplaceAll("p", "+");
    tmp.ReplaceAll("m", "-");
    tmp.ReplaceAll("d", ".");
    tmp.ReplaceAll("_", " ");
    TObjArray* tokens = tmp.Tokenize(" ");
    if (!tokens || tokens->GetEntriesFast() < 2) { 
      Error("Other2Stack", "Failed to decode eta range from %s", 
	    d->GetName());
      if (tokens) tokens->Delete();
      return;
    }
    Double_t eta1 = static_cast<TObjString*>(tokens->At(0))->String().Atof();
    Double_t eta2 = static_cast<TObjString*>(tokens->At(1))->String().Atof();
    tokens->Delete();

    fBody->Divide(1,5,0,0);
    TVirtualPad* p = fBody->cd(5);
    p->SetRightMargin(0.15);
    p->SetTopMargin(0.05);
    p = fBody->cd(4);
    p->SetBottomMargin(0.15);
    THStack* all = GetStack(d,"all");
    if (!all) {
      Warning("ProcessBin", "Argh! All stack not found!");
      return;
    }
    DrawInPad(fBody,1,all, "nostack", kLogy);
    DrawInPad(fBody,2,GetH1(d,"ratioCorrTruth"), "", 0);
    DrawInPad(fBody,3,GetH1(d,"ratioUnfAcc"), "", 0);
    DrawInPad(fBody,4,GetH1(d,"triggerVertex"), "", 0);
    DrawInPad(fBody,5,GetH2(d,"response"), "colz", kLogz);
    all->GetXaxis()->SetTitle("#it{N}_{ch}");
    all->GetYaxis()->SetTitle("P(#it{N}_{ch})");

    PrintCanvas(Form(" %+5.2f<eta<%+5.2f", eta1, eta2));
    
    DrawSteps(all, eta1, eta2);
  }
  //____________________________________________________________________
  void DrawSteps(THStack* stack, Double_t e1, Double_t e2) 
  {
    fBody->Divide(2, 3, 0, 0);
    TList* hists = stack->GetHists();
    Int_t  nHist = hists->GetEntries(); 
    
    // Make a background stack 
    THStack* bg = static_cast<THStack*>(stack->Clone("bg"));
    bg->SetTitle();
    
    for (Int_t i = 0; i < nHist; i++) { 
      // Loop over histograms and set the saved color
      TH1* h        = static_cast<TH1*>(bg->GetHists()->At(i));
      h->SetMarkerColor(kGray+1);
      h->SetFillColor(kGray);
      h->SetLineColor(kGray);

      TList* lf = h->GetListOfFunctions();
      if (lf) { 
	TObject* ll = lf->FindObject("legend");
	if (ll) lf->Remove(ll);
      }
    }
    const char* txt[] = { "MC 'truth'", 
			  "Selected MC 'truth'", 
			  "Measured", 
			  "Unfolded", 
			  "Corrected" };
    // Now loop again, this time drawing the stack
    for (Int_t i = 0; i < nHist; i++) { 
      
      DrawInPad(fBody, i+1, bg, "nostack", kLogy);
      DrawInPad(fBody, i+1, hists->At(i), "same hist p e", kLogy);
      gPad->SetGridx();
      gPad->SetGridy();
      bg->GetXaxis()->SetTitle("#it{N}_{ch}");
      bg->GetYaxis()->SetTitle("P(#it{N}_{ch})");

      TLatex* l = new TLatex(.95, .95, Form("Step %d", i));
      l->SetNDC();
      l->SetTextAlign(33);
      l->SetTextSize(0.06);
      l->Draw();
      l->DrawLatex(.95, .88, txt[i]);
      
    }
    PrintCanvas(Form(" %+5.2f<eta<%+5.2f - Steps", e1, e2));
  }
};
//
// EOF
//

